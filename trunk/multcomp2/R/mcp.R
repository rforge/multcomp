
### multiple comparison procedures for generalized linear models
### `object' requires methods for `model.matrix', `model.frame', `terms', 
### `coef' and `vcov'
### `hypotheses' is a named list of characters (`type' arguments in `linhypo')

mcp <- function(object, hypotheses = NULL, 
                alternative = c("two.sided", "less", "greater")) {

    ### extract model matrix, frame and terms
    mm <- try(model.matrix(object))
    if (inherits(mm, "try-error"))
        stop("no ", sQuote("model.matrix"), " method for ", 
             sQuote("object"), " found!")

    mf <- try(model.frame(object))
    if (inherits(mf, "try-error"))
        stop("no ", sQuote("model.frame"), " method for ", 
             sQuote("object"), " found!")

    tm <- try(terms(object))
    if (inherits(tm, "try-error"))
        stop("no ", sQuote("terms"), " method for ", 
             sQuote("object"), " found!")

    beta <- try(coef(object))
    if (inherits(beta, "try-error"))
        stop("no ", sQuote("coef"), " method for ",
             sQuote("object"), " found!")
    sigma <- try(vcov(object))
    if (inherits(sigma, "try-error"))
        stop("no ", sQuote("vcov"), " method for ",
             sQuote("object"), " found!")       

    alternative <- match.arg(alternative)

    df <- 0
    if (class(object)[1] %in% c("aov", "lm")) {
        class(object) <- "lm"
        df <- summary(object)$df[2]
    }


    ### OK! You know what you want!
    if (is.matrix(hypotheses)) { 
        if(ncol(hypotheses) != length(beta))
            stop("dimensions of ", sQuote("hypotheses"), " and ", sQuote("coef(object)"),
                 "don't match")

         RET <- list(object = object, 
                hypotheses = hypotheses, beta = beta, sigma = sigma,
                type = "user-defined",
                alternative = alternative, df = df)
         class(RET) <- "mcp"
         return(RET)
     }

    ### factors and contrasts
    contrasts <- attr(mm, "contrasts")
    factors <- attr(tm, "factors")
   
    ### linear hypotheses
    if (!is.list(hypotheses) || is.null(names(hypotheses)))
        stop(sQuote("hypotheses"), "is not a named list")
    nhypo <- names(hypotheses)
    checknm <- nhypo %in% rownames(factors) & sapply(mf[nhypo], is.factor)
    if (!all(checknm)) 
        stop("Factor", nhypo[!checknm], "not found!")
    for (nm in nhypo) {
        if (is.character(hypotheses[[nm]]))
            hypotheses[[nm]] <- linhypo(mf[[nm]], hypotheses[[nm]])
    }

    ### transform linear hypotheses using model contrasts
    hypo <- vector(mode = "list", length = length(nhypo))
    names(hypo) <- nhypo

    for (nm in nhypo) {
        ### extract contrast matrix for each factor
        if (is.character(contrasts[[nm]])) {
            C <- do.call(contrasts[[nm]], 
                         list(n = nlevels(mf[[nm]])))
        } else {
            C <- contrasts[[nm]]
        }
        ### and transform the original linear hypothesis
        Kstar <- hypotheses[[nm]] %*% C
        pos <- factors[nm,] == 1
        ### average over interaction terms (if any)
        if (sum(pos) > 1) {
            Kinter <- c()
            for (i in which(pos)[-1]) {
                k <- sum(attr(mm, "assign") == i) / ncol(Kstar)
                if (sum(factors[1:which(rownames(factors) == nm), i]) == 1) {
                    Kinter <- cbind(Kinter, Kstar[,rep(1:ncol(Kstar), k), drop = FALSE] / (k + 1))
                } else {
                    Kinter <- cbind(Kinter, 
                        Kstar[,rep(1:ncol(Kstar), rep(k, ncol(Kstar))), drop = FALSE] / (k + 1))
                }
            }
            Kstar <- cbind(Kstar, Kinter)
        }
        hypo[[nm]] <- list(K = Kstar,
                           where = attr(mm, "assign") %in% which(factors[nm,] == 1),
                           type = paste(attr(hypotheses[[nm]], "type"), "(", nm, ")", sep = ""))
    }

    ### create matrix of all transformed linear hypotheses
    M <- matrix(0, nrow = sum(sapply(hypo, function(x) nrow(x$K))),
                    ncol = ncol(mm))
    colnames(M) <- colnames(mm)

    count <- 1
    for (h in hypo) {
        M[count:(count + nrow(h$K) - 1), h$where] <- h$K
        count <- count + nrow(h$K)
    }
    rownames(M) <- unlist(lapply(hypo, function(x) rownames(x$K)))


    ### create `mcp' object
    RET <- list(object = object, 
                hypotheses = M, beta = beta, sigma = sigma,
                type = paste(sapply(hypo, function(x) x$type), collapse = ";"),
                alternative = alternative,
                df = df)
    class(RET) <- "mcp"
    return(RET)
}

summary.mcp <- function(object, distribution = adjusted(), ...) {

    pq <- distribution(object)
    object$mcp <- cbind(pq$coefficients, pq$sigma, pq$tstat, 
                        pq$pvalues)
    colnames(object$mcp) <- c("Estimate", "Std. Error",
        ifelse(object$df == 0, "z value", "t value"), "p value")
    attr(object$mcp, "type") <- pq$type
    class(object) <- c("summary.mcp", "mcp")
    return(object)
}

confint.mcp <- function(object, parm, level = 0.95, ...) {

    pq <- pqmcp(object)
    object$confint <- cbind(pq$coefficients, 
                            pq$qfunction(conf.level = level, ...))
    colnames(object$confint) <- c("Estimate", "lwr", "upr")
    attr(object$confint, "conf.level") <- level
    class(object) <- c("confint.mcp", "mcp")
    return(object)
}
