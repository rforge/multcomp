
### multiple comparison procedures for generalized linear models
### `object' requires methods for `model.matrix', `model.frame', `terms', 
### `coef' and `vcov'
### `hypotheses' is a named list of characters (`type' arguments in `contrMat')

mcp <- function(object, hypotheses = NULL) {

    ### extract model matrix, frame and terms
    mm <- try(model.matrix(object))
    if (inherits(mm, "try-error"))
        stop("no", sQuote("model.matrix"), "method for", 
             sQuote("object"), "found!")

    mf <- try(model.frame(object))
    if (inherits(mf, "try-error"))
        stop("no", sQuote("model.frame"), "method for", 
             sQuote("object"), "found!")

    tm <- try(terms(object))
    if (inherits(tm, "try-error"))
        stop("no", sQuote("terms"), "method for", 
             sQuote("object"), "found!")

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
                    Kinter <- cbind(Kinter, Kstar[,rep(1:ncol(Kstar), k)] / (k + 1))
                } else {
                    Kinter <- cbind(Kinter, 
                        Kstar[,rep(1:ncol(Kstar), rep(k, ncol(Kstar)))] / (k + 1))
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
    beta <- try(coef(object))
    if (inherits(beta, "try-error"))
        stop("no", sQuote("coef"), "method for",
             sQuote("object"), "found!")
    sigma <- try(vcov(object))
    if (inherits(sigma, "try-error"))
        stop("no", sQuote("vcov"), "method for",
             sQuote("object"), "found!")       

    RET <- list(object = object, 
                hypotheses = M, beta = beta, sigma = sigma,
                type = paste(sapply(hypo, function(x) x$type), collapse = ";"))
    class(RET) <- "mcp"
    return(RET)
}

print.mcp <- function(x, ...) {
    cat("\n\t", "Multiple Comparison Procedures for Linear Hypotheses\n")
    cat("\t", "Estimates: ")
    print(x$hypotheses %*% x$beta)
}

summary.mcp <- function(object, logical = FALSE, ...) {

    ### use multivariate t distribution for linear models,
    ### normal distribution otherwise
    df <- 0
    asympt <- TRUE
    if (inherits(object$object, "lm")) {
        class(object$object) <- "lm"
        df <- summary(object$object)$df[2]
        asympt <- FALSE
    }

    ### OK, we are done, call the work horse ...
    csimtest(estpar = object$beta,
            df = df,
            covm = object$sigma,
            cmatrix = object$hypotheses,
            ctype = object$type,
            ttype = ifelse(logical, "logical", "free"),
            asympt = asympt,
            ...)
}
    

confint.mcp <- function(object, parm, level = 0.95, ...) {

    ### use multivariate t distribution for linear models,
    ### normal distribution otherwise
    df <- 0
    asympt <- TRUE
    if (inherits(object$object, "lm")) {
        class(object$object) <- "lm"
        df <- summary(object$object)$df[2]
        asympt <- FALSE
    }

    ### OK, we are done, call the work horse ...
    csimint(estpar = object$beta, 
            df = df, 
            covm = object$sigma, 
            cmatrix = object$hypotheses, 
            ctype = object$type, 
            conf.level = level, 
            asympt = asympt, 
            ...) 
}
