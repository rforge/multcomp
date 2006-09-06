
# $Id: simint.R,v 1.52 2005/07/25 15:25:05 hothorn Exp $

pqglht <- function(object) 
{

    beta <- object$beta
    sigma <- object$sigma
    K <- object$K
    m <- object$m
    df <- object$df
    alternative <- object$alternative

    p <- length(beta)
    covm <- K %*% sigma %*% t(K)
    d <- 1/sqrt(diag(covm))
    if (length(d) > 1) d <- diag(d)              
    cr    <- d %*% covm %*% t(d)

    betahat  <- K %*% beta
    ses <- sqrt(diag(covm))
    tstat <- (betahat - m) / ses
    dim <- ncol(cr)

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
            if (df > 0) pvals <- 2*(1-pt(abs(tstat),df))     
            else        pvals <- 2*(1-pnorm(abs(tstat)))
        }, "less" = {
            if (df > 0) pvals <- pt(tstat,df) 
            else        pvals <- pnorm(tstat)
        }, "greater" = {
            if (df > 0) pvals <- 1-pt(tstat,df)
            else        pvals <- 1-pnorm(tstat)
        })

        if (type == "univariate") {
            return(pvals)
        }

        if (type == "Bonferroni")
            return(pmin(1, dim * pvals))

        if (type == "adjusted")
            return(1 - apply(tstat, 1, pfct))
    }

    qfunction <- function(conf.level, ...) {

        tail <- switch(alternative, "two.sided" = "both.tails",
                                    "less"      = "upper.tail",
                                    "greater"   = "lower.tail")
        calpha <- qmvt(conf.level, df = df, corr = cr, tail = tail, 
                       ...)$quantile

        switch(alternative, "two.sided" = {  
            LowerCL <- betahat - calpha*ses
            UpperCL <- betahat + calpha*ses
        }, "less" = {
            LowerCL <- rep(-Inf, dim)
            UpperCL <- betahat - calpha*ses
        }, "greater" = {
            LowerCL <- betahat - calpha*ses
            UpperCL <- rep( Inf, dim)
        })

        cint <- cbind(LowerCL, UpperCL)
        colnames(cint) <- c("lower", "upper")
        attr(cint, "conf.level") <- conf.level
        return(cint)
    }
    RET <- list(pfunction = pfunction, qfunction = qfunction,
                coefficients = betahat, sigma = ses, tstat = tstat)
    class(RET) <- "pqglht"
    RET
}

univariate <- function() 
{
    function(object) {
        RET <- pqglht(object)
        RET$pvalues <- RET$pfunction("univariate")
        RET$type <- "univariate"
        RET
    }
}

adjusted <- function(type = c("free", "Bonferroni", "Shaffer", "Westfall"), ...) 
{
    type <- match.arg(type)

    ### simple Bonferroni-adjustment
    if (type == "Bonferroni") {
        return(function(object) {
            RET <- pqghlt(object)
            RET$pvalues <- RET$pfunction("Bonferroni")
            RET$type <- "Bonferroni"
            RET
        })
    }

    ### usual max-type adjustment over all linear hypotheses
    if (type == "free") {
        return(function(object) {
            RET <- pqglht(object)
            RET$pvalues <- RET$pfunction("adjusted", ...)
            RET$type <- type
            RET
        })
    }

    ### Westfall (1997, JASA): constraints and correlations 
    ### or
    ### Shaffer (1886, JASA): constraints
    return(function(object) {
        RET <- pqglht(object)
        tstat <- switch(object$alternative, 
                        "less" = RET$tstat,
                        "greater" = -RET$tstat,
                        "two.sided" = -abs(RET$tstat))
        C <- object$K
        Corder <- C[order(tstat), ]
        ms <- maxsets(Corder)
        p <- sapply(ms, function(x) {
           max(sapply(x, function(s) {
                object$K <- Corder[s, , drop = FALSE]
                min(pqglht(object)$pfunction(ifelse(type == "Westfall", 
                                                   "adjusted", "Bonferroni"), 
                                            ...))
            }))
        })
        for (i in 2:length(p))
            p[i] <- max(p[i-1], p[i])
        ### <FIXME> what happens in case of ties??? </FIXME> ###
        RET$pvalues <- p[rank(tstat)]
        RET$type <- type
        RET
    })
}
