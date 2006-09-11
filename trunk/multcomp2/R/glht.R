
### simultaneous inference procedures for general linear hypothesis
### in a fitted generalized linear model `model'
###
### `model' requires methods for `model.matrix', `model.frame', `terms', 
### `coef' and `vcov'
### `K' is either a named list of characters or matrices _or_ a matrix.

glht <- function(model, K, m = 0,
                 alternative = c("two.sided", "less", "greater")) {

    ### extract model matrix, frame and terms
    mm <- try(model.matrix(model))
    if (inherits(mm, "try-error"))
        stop("no ", sQuote("model.matrix"), " method for ", 
             sQuote("model"), " found!")

    mf <- try(model.frame(model))
    if (inherits(mf, "try-error"))
        stop("no ", sQuote("model.frame"), " method for ", 
             sQuote("model"), " found!")

    tm <- try(terms(model))
    if (inherits(tm, "try-error"))
        stop("no ", sQuote("terms"), " method for ", 
             sQuote("model"), " found!")

    beta <- try(coef(model))
    if (inherits(beta, "try-error"))
        stop("no ", sQuote("coef"), " method for ",
             sQuote("model"), " found!")

    sigma <- try(vcov(model))
    if (inherits(sigma, "try-error"))
        stop("no ", sQuote("vcov"), " method for ",
             sQuote("model"), " found!")       

    alternative <- match.arg(alternative)

    df <- 0
    if (class(model)[1] %in% c("aov", "lm")) {
        class(model) <- "lm"
        df <- summary(model)$df[2]
    }

    if (!any(c(is.character(K),
               is.matrix(K),
               is.list(K))))
        stop(sQuote("K"), " is neither character, matrix nor list")

    ### interpret K in terms of coef(model)
    if (is.character(K)) {
        tmp <-  chr2K(K, names(beta))
        K <- tmp$K
        if (m != 0)
            warning(sQuote("m"), " is given in both ", sQuote("K"),
                    " and ", sQuote("m"), " -- the latter is ignored")
        m <- tmp$m
    }
    if (is.matrix(K)) { 
        if (ncol(K) != length(beta))
            stop("dimensions of ", sQuote("K"), " and ", 
                 sQuote("coef(model)"), "don't match")

        if (length(m) == 1) m <- rep(m, nrow(K))
        if (length(m) != nrow(K))
            stop("dimensions of ", sQuote("K"), " and ", 
                 sQuote("m"), "don't match")

        RET <- list(model = model, K = K, m = m,
                    beta = beta, sigma = sigma, df = df,
                    alternative = alternative,
                    type = "user-defined")
        class(RET) <- "glht"
        return(RET)
    }

    ### interpret K in terms of factor levels
    ### extract factors and contrasts
    contrasts <- attr(mm, "contrasts")
    factors <- attr(tm, "factors")
    intercept <- attr(tm, "intercept") != 0
   
    ### linear hypotheses
    if (!is.list(K) || is.null(names(K)))
        stop(sQuote("K"), "is not a named list")
    nhypo <- names(K)
    checknm <- nhypo %in% rownames(factors) & sapply(mf[nhypo], is.factor)
    if (!all(checknm)) 
        stop("Factor(s) ", nhypo[!checknm], " not found!")
    mnull <- m
    for (nm in nhypo) {
        if (is.character(K[[nm]])) {

            Kchr <- function(kch) {
                ### check if kch is suitable as `type' argument to `contrMat'
                types <- eval(formals(contrMat)$type)
                pm <- pmatch(kch, types)
                ### if yes, compute K from `contrMat'
                if (!is.na(pm)) {
                    tmpK <- contrMat(table(mf[[nm]]), type = types[pm])
                    tmpm <- rep(mnull, nrow(tmpK))
                } else {
                    ### if not, interpret kch as an expression
                    tmp <-  chr2K(kch, levels(mf[[nm]]))
                    tmpK <- tmp$K
                    tmpm <- tmp$m
                }
                if (is.null(rownames(tmpK)))
                    rownames(tmpK) <- 1:nrow(tmpK)
                list(K = tmpK, m = tmpm)
            }
            
            tmp <- lapply(K[[nm]], Kchr)
            K[[nm]] <- do.call("rbind", lapply(tmp, function(x) x$K))
            if (length(m) == 1) {
                m <- drop(unlist(sapply(tmp, function(x) x$m)))
            } else {
                m <- c(m, drop(unlist(sapply(tmp, function(x) x$m))))
            }
        }
    }

    ### transform linear hypotheses using model contrasts
    hypo <- vector(mode = "list", length = length(nhypo))
    names(hypo) <- nhypo

    for (nm in nhypo) {
        ### extract contrast matrix for each factor from model fit
        if (is.character(contrasts[[nm]])) {
            C <- do.call(contrasts[[nm]], 
                         list(n = nlevels(mf[[nm]])))
        } else {
            C <- contrasts[[nm]]
        }
        ### and transform the original linear hypotheses 
        ### K beta to K C beta^* 
        if (intercept) {
            Kstar <- K[[nm]] %*% C
        } else {
            ### model.matrix has `contrasts' argument even if no intercept
            ### was fitted and the contrast actually hasn't been applied
            Kstar <- K[[nm]]
        }
        pos <- factors[nm,] == 1
        ### average over interaction terms (if any)
        if (sum(pos) > 1) {
            Kinter <- c()
            for (i in which(pos)[-1]) {
                k <- sum(attr(mm, "assign") == i) / ncol(Kstar)
                if (sum(factors[1:which(rownames(factors) == nm), i]) == 1) {
                    Kinter <- cbind(Kinter, 
                        Kstar[,rep(1:ncol(Kstar), k), drop = FALSE] / (k + 1))
                } else {
                    Kinter <- cbind(Kinter, 
                        Kstar[,rep(1:ncol(Kstar), rep(k, ncol(Kstar))), 
                              drop = FALSE] / (k + 1))
                }
            }
            Kstar <- cbind(Kstar, Kinter)
        }
        hypo[[nm]] <- list(K = Kstar,
            where = attr(mm, "assign") %in% which(factors[nm,] == 1),
            type = paste(attr(K[[nm]], "type"), "(", nm, ")", sep = ""))
    }

    ### create matrix of all linear hypoheses
    Ktotal <- matrix(0, nrow = sum(sapply(hypo, function(x) nrow(x$K))),
                     ncol = ncol(mm))
    colnames(Ktotal) <- colnames(mm)

    count <- 1
    for (h in hypo) {
        Ktotal[count:(count + nrow(h$K) - 1), h$where] <- h$K
        count <- count + nrow(h$K)
    }
    rownames(Ktotal) <- unlist(lapply(hypo, function(x) rownames(x$K)))

    if (length(m) == 1) m <- rep(m, nrow(Ktotal))
    if (length(m) != nrow(Ktotal))
        stop("dimensions of ", sQuote("K"), " and ", 
             sQuote("m"), "don't match")

    ### create `glht' model
    RET <- list(model = model, K = Ktotal, m = m, 
                beta = beta, sigma = sigma, df = df,
                alternative = alternative,
                type = paste(sapply(hypo, function(x) x$type), 
                             collapse = ";"))
    class(RET) <- "glht"
    return(RET)
}

coef.glht <- function(object, null = FALSE, ...) 
{
    if (null) return(object$m)
    betahat <- coef(object$model)
    if (inherits(object$model, "lmer"))
        betahat <- coeflmer(object$model)
    drop(object$K %*% betahat)
}

vcov.glht <- function(object, ...) 
    object$K %*% tcrossprod(as.matrix(vcov(object$model)), object$K)

summary.glht <- function(object, test = adjusted(), ...) 
{
    pq <- test(object)
    object$mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, 
                           pq$pvalues)
    attr(object$mtests, "error") <- attr(pq$pvalues, "error")
    colnames(object$mtests) <- c("Estimate", "Std. Error",
        ifelse(object$df == 0, "z value", "t value"), "p value")
    attr(object$mtests, "type") <- pq$type
    class(object) <- c("summary.glht", "glht")
    return(object)
}

confint.glht <- function(object, parm, level = 0.95, ...) 
{
    pq <- pqglht(object)
    ci <- pq$qfunction(conf.level = level, ...)
    object$confint <- cbind(pq$coefficients, ci)
    colnames(object$confint) <- c("Estimate", "lwr", "upr")
    attr(object$confint, "conf.level") <- level
    attr(object$confint, "calpha") <- attr(ci, "calpha")
    attr(object$confint, "error") <- attr(ci, "error")
    class(object) <- c("confint.glht", "glht")
    return(object)
}
