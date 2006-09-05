
simint <- function(object, ...) UseMethod("simint")
simtest <- function(object, ...) UseMethod("simtest")

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
    .Deprecated("mcp", package = "multcomp2")
    confint(tmcp, level = conf.level, abseps = eps, maxpts = maxpts)
}

simint.formula <- function(formula, data=list(), subset, na.action, ...)
{
    cl <- match.call(expand.dots = FALSE)
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
    .Deprecated("mcp", package = "multcomp2")
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
