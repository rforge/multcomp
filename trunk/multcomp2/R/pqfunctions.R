
# $Id: simint.R,v 1.52 2005/07/25 15:25:05 hothorn Exp $

pqmcp <- function(object) 
{

    beta <- object$beta
    sigma <- object$sigma
    linhypo <- object$hypotheses
    df <- object$df
    alternative <- object$alternative
    p <- length(beta)
    covm  <- linhypo %*% sigma %*% t(linhypo)
    d     <- 1/sqrt(diag(covm))
    if (length(d) > 1)
      d <- diag(d)              
    cr    <- d %*% covm %*% d

    ests  <- linhypo %*% beta
    ses   <- sqrt(diag(covm))
    tvals <- ests/ses
    dim   <- ncol(cr)

    pfunction <- function(type = c("univariate", "Bonferroni", "adjusted"), ...) {

        type <- match.arg(type) 

        pfct <- function(q) {
            switch(alternative, "two.sided" = {
                      low <- rep(-abs(q), dim)
                      upp <- rep( abs(q), dim)
               }, "less" = {
                      low <- rep(      q, dim)
                      upp <- rep(    Inf, dim)
               }, "greater" = {
                      low <- rep(   -Inf, dim)
                      upp <- rep(      q, dim)
               })
               pmvt(lower = low, upper = upp, df = df, corr = cr, ...)
        }

        switch(alternative, "two.sided" = {
            if (df > 0) pvals <- 2*(1-pt(abs(tvals),df))     
            else        pvals <- 2*(1-pnorm(abs(tvals)))
        }, "less" = {
            if (df > 0) pvals <- pt(tvals,df) 
            else        pvals <- pnorm(tvals)
        }, "greater" = {
            if (df > 0) pvals <- 1-pt(tvals,df)
            else        pvals <- 1-pnorm(tvals)
        })

        if (type == "univariate") {
            return(pvals)
        }

        if (type == "Bonferroni")
            return(pmin(1, dim * pvals))

        if (type == "adjusted")
            return(1 - apply(tvals, 1, pfct))
    }

    qfunction <- function(conf.level, ...) {

        tail <- switch(alternative, "two.sided" = "both.tails",
                                    "less"      = "upper.tail",
                                    "greater"   = "lower.tail")
        calpha <- qmvt(conf.level, df = df, corr = cr, tail = tail, 
                       ...)$quantile

        switch(alternative, "two.sided" = {  
            LowerCL <- ests - calpha*ses
            UpperCL <- ests + calpha*ses
        }, "less" = {
            LowerCL <- rep(-Inf, dim)
            UpperCL <- ests - calpha*ses
        }, "greater" = {
            LowerCL <- ests - calpha*ses
            UpperCL <- rep( Inf, dim)
        })

        cint <- cbind(LowerCL, UpperCL)
        colnames(cint) <- c("lower", "upper")
        attr(cint, "conf.level") <- conf.level
        return(cint)
    }
    RET <- list(pfunction = pfunction, qfunction = qfunction,
                coefficients = ests, sigma = ses, tstat = tvals)
    class(RET) <- "pqmcp"
    RET
}

univariate <- function() 
{
    function(object) {
        RET <- pqmcp(object)
        RET$pvalues <- RET$pfunction("univariate")
        RET$type <- "univariate"
        RET
    }
}

Bonferroni <- function() {
    function(object) {
        RET <- pqmcp(object)
        RET$pvalues <- RET$pfunction("Bonferroni")
        RET$type <- "Bonferroni"
        RET
    }
}

adjusted <- function(...) {
    function(object) {
        RET <- pqmcp(object)
        RET$pvalues <- RET$pfunction("adjusted", ...)
        RET$type <- "adjusted"
        RET
    }
}
