
### general linear hypotheses
glht <- function(model, linfct, ...) UseMethod("glht", linfct)

### K coef(model) _!alternative_ rhs
glht.matrix <- function(model, linfct, 
    alternative = c("two.sided", "less", "greater"), rhs = 0, 
    df = NULL, ...) {

    ### extract coefficients and their covariance matrix, df
    tmp <- coefvcovdf(model)
    beta <- tmp$beta

    alternative <- match.arg(alternative)
    if (!is.numeric(rhs))
        stop(sQuote("rhs"), " is not a numeric vector")

    if (ncol(linfct) != length(beta))
        stop(sQuote("ncol(linfct)"), " is not equal to ", 
             sQuote("length(coef(model))"))

    if (is.null(colnames(linfct)))
        colnames(linfct) <- names(beta)

    if (is.null(rownames(linfct)))
        rownames(linfct) <- 1:nrow(linfct)

    if (length(rhs) == 1) rhs <- rep(rhs, nrow(linfct))
    if (length(rhs) != nrow(linfct))
        stop(sQuote("nrow(linfct)"), " is not equal to ",
             sQuote("length(rhs)"))

    RET <- list(model = model, K = linfct, m = rhs,
                beta = beta, sigma = tmp$sigma, 
                df = ifelse(is.null(df), tmp$df, df),
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
    return(glht(model, linfct = tmp$K, rhs = tmp$m, 
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
        args$rhs <- tmp$m
    args <- c(args, list(...))

    return(do.call("glht", args))
}
