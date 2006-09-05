
contrMat <- function(n, ...) {
    x <- factor(rep(1:length(n), n))
    if (!null(names(n)))
        levels(x) <- names(n)
    .Deprecated("linhypo", package = "multcomp")
    linhypo(x)
}

simint <- function(object, ...) UseMethod("simint")
simtest <- function(object, ...) UseMethod("simtest")

csimint <- function(estpar, df, covm, cmatrix=NULL, ctype="user-defined",
                    conf.level=0.95,
                    alternative=c("two.sided","less","greater"), asympt=FALSE,
                    eps=0.001, maxpts=1000000)
{
    if (is.null(cmatrix)) cmatrix <- diag(length(estpar))
    object <- list(object = NULL, hypotheses = cmatrix, beta = estpar, sigma = covm,
                   type = ctype, alternative = match.arg(alternative),
                   df = ifelse(is.null(df) || asympt, 0, df))
    class(object) <- "mcp"
    .Deprecated("mcp", package = "multcomp")
    confint(object, level = conf.level, abseps = eps, maxpts = maxpts)
}

simint.default <- function(object, 
    type = c("Dunnett", "Tukey", "Sequen", "AVE", "Changepoint", "Williams", "Marcus",
             "McDermott"), 
    cmatrix = NULL, conf.level = 0.95,
    alternative = c("two.sided","less","greater"), 
    eps = 0.001, maxpts = 1e+06, whichf) {

    if (!is.null(cmatrix)) {
        linhypo <- cmatrix
    } else {
        linhypo <- list(match.arg(type))
        names(linhypo) <- whichf
    }
    tmcp <- mcp(object, hypotheses = linhypo, alternative = match.arg(alternative))
    .Deprecated("mcp", package = "multcomp")
    confint(tmcp, level = conf.level, abseps = eps, maxpts = maxpts)
}

simint.formula <- function(formula, data=list(), subset, na.action, ...)
{
    cl <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(cl), 0)
    cl <- cl[c(1, m)]
    cl[[1]] <- as.name("lm")
    object <- eval(cl, parent.frame())
    addargs <- list(...)
    whichf <- addargs$whichf
    if (is.null(whichf)) {
        mm <- model.frame(object)
        whichf <- names(mm)[sapply(mm, class) == "factor"]
    }
    addargs$whichf <- whichf
    addargs$object <- object
    do.call("simint.default", addargs)
}

simint.lm <- function(object, psubset = NULL, ...) {

    beta <- coef(object)
    if (is.null(psubset)) 
        simin.default(object, cmatrix = diag(length(beta)), ...)
    
    psubset <- which(beta %in% beta[psubset])
    simtest.default(object, cmatrix = diag(length(beta))[psubset,])
}


csimtest <- function(estpar, df, covm, cmatrix=NULL, ctype="user-defined",
                     ttype=c("free","logical"),
                     alternative=c("two.sided","less","greater"), asympt=FALSE,
                     eps=0.001, maxpts=1000000)
{
    if (is.null(cmatrix)) cmatrix <- diag(length(estpar))
    object <- list(object = NULL, hypotheses = cmatrix, beta = estpar, sigma = covm,
                   type = ctype, alternative = match.arg(alternative),
                   df = ifelse(is.null(df) || asympt, 0, df))
    class(object) <- "mcp"
    ttype <- match.arg(ttype)
    if (ttype == "free")
        distr <- adjusted("free")
    if (ttype == "logical")
        distr <- adjusted("Westfall")
    .Deprecated("mcp", package = "multcomp")
    summary(object, distribution = distr, abseps = eps, maxpts = maxpts)
}

simtest.default <- function(object, 
    type = c("Dunnett", "Tukey", "Sequen", "AVE", "Changepoint", "Williams", "Marcus",
             "McDermott"), 
    ttype = c("free", "logical"),
    cmatrix = NULL, conf.level = 0.95,
    alternative = c("two.sided","less","greater"), 
    eps = 0.001, maxpts = 1e+06, whichf) {

    if (!is.null(cmatrix)) {
        linhypo <- cmatrix
    } else {
        linhypo <- list(match.arg(type))
        names(linhypo) <- whichf
    }
    tmcp <- mcp(object, hypotheses = linhypo, alternative = match.arg(alternative))
    .Deprecated("mcp", package = "multcomp")
    ttype <- match.arg(ttype)
    if (ttype == "free")
        distr <- adjusted("free")
    if (ttype == "logical")
        distr <- adjusted("Westfall")
    summary(tmcp, distribution = distr, abseps = eps, maxpts = maxpts)
}

simtest.formula <- function(formula, data=list(), subset, na.action, ...)
{
    cl <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(cl), 0)
    cl <- cl[c(1, m)]
    cl[[1]] <- as.name("lm")
    object <- eval(cl, parent.frame())
    addargs <- list(...)
    whichf <- addargs$whichf
    if (is.null(whichf)) {
        mm <- model.frame(object)
        whichf <- names(mm)[sapply(mm, class) == "factor"]
    }
    addargs$whichf <- whichf
    addargs$object <- object
    do.call("simtest.default", addargs)
}

simtest.lm <- function(object, psubset = NULL, ...) {

    beta <- coef(object)
    if (is.null(psubset)) 
        simin.default(object, cmatrix = diag(length(beta)), ...)
    
    psubset <- which(beta %in% beta[psubset])
    simtest.default(object, cmatrix = diag(length(beta))[psubset,])
}

summary.summary.mcp <- function(object, ...) invisible(x)
summary.confint.mcp <- function(object, ...) invisible(x)
