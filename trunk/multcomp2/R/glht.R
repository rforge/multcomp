
### extract coefficients, covariance matrix and 
### degrees of freedom (if available) from `model'
coefvcovdf <- function(model) {

    ### extract coefficients and their covariance matrix
    beta <- try(coef(model))
    if (inherits(beta, "try-error"))
        stop("no ", sQuote("coef"), " method for ",
             sQuote("model"), " found!")
    if (inherits(model, "lmer"))
        beta <- coeflmer(model)

    sigma <- try(vcov(model))
    if (inherits(sigma, "try-error"))
        stop("no ", sQuote("vcov"), " method for ",
             sQuote("model"), " found!")       

    ### check if a linear model was supplied
    df <- 0
    if (class(model)[1] %in% c("aov", "lm")) {
        class(model) <- "lm"
        df <- summary(model)$df[2]
    }

    list(beta = beta, sigma = sigma, df = df)
}

### general linear hypotheses
glht <- function(model, linfct, ...) UseMethod("glht", linfct)

### K coef(model) _!alternative_ m
glht.matrix <- function(model, linfct, m = 0, 
    alternative = c("two.sided", "less", "greater"), ...) {

    ### extract coefficients and their covariance matrix, df
    tmp <- coefvcovdf(model)
    beta <- tmp$beta

    alternative <- match.arg(alternative)
    if (!is.numeric(m))
        stop(sQuote("m"), " is not a numeric vector")

    if (ncol(linfct) != length(beta))
        stop(sQuote("ncol(linfct)"), " is not equal to ", 
             sQuote("length(coef(model))"))

    if (is.null(colnames(linfct)))
        colnames(linfct) <- names(beta)

    if (is.null(rownames(linfct)))
        rownames(linfct) <- 1:nrow(linfct)

    if (!is.numeric(m)) stop(sQuote("m"), " is not a numeric vector")
    if (length(m) == 1) m <- rep(m, nrow(linfct))
    if (length(m) != nrow(linfct))
        stop(sQuote("nrow(linfct)"), " is not equal to ",
             sQuote("length(m)"))

    RET <- list(model = model, K = linfct, m = m,
                beta = beta, sigma = tmp$sigma, df = tmp$df,
                alternative = alternative,
                type = "user-defined")
    class(RET) <- "glht"
    RET
}

glht.character <- function(model, linfct, ...) {
    ### extract coefficients and their covariance matrix
    beta <- try(coef(model))
    if (inherits(beta, "try-error"))
        stop("no ", sQuote("coef"), " method for ",
             sQuote("model"), " found!")

    tmp <- chrlinfct2matrix(linfct, names(beta))
    return(glht(model, linfct = tmp$K, m = tmp$m, 
                alternative = tmp$alternative))
}

glht.expression <- function(model, linfct, ...) 
    glht(model, as.character(linfct), ...)

### multiple comparison procedures
glht.mcp <- function(model, linfct, ...) {

    ### extract factors and contrast matrices from `model'
    tmp <- mcp2matrix(model, linfct = linfct)

    args <- list(model = model, linfct = tmp$K)
    if (!is.null(tmp$alternative))
        args$alternative <- tmp$alternative
    if (any(tmp$m != 0))
        args$m <- tmp$m
    args <- c(args, list(...))

    return(do.call("glht", args))
}
