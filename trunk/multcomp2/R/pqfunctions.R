
# $Id: simint.R,v 1.52 2005/07/25 15:25:05 hothorn Exp $

pqglht <- function(object) 
{
    betahat <- coef(object)
    covm <- vcov(object)
    m <- coef(object, null = TRUE)
    df <- object$df

    ses <- sqrt(diag(covm))
    tstat <- (betahat - m) / ses
    cr <- cov2cor(covm)
    dim <- ncol(cr)

    pfunction <- function(type = c("univariate", "Bonferroni", "adjusted"), 
                          ...) {

        type <- match.arg(type) 

        pfct <- function(q) {
            switch(object$alternative, "two.sided" = {
                      low <- rep(-abs(q), dim)
                      upp <- rep( abs(q), dim)
               }, "less" = {
                      low <- rep(q, dim)
                      upp <- rep(Inf, dim)
               }, "greater" = {
                      low <- rep(-Inf, dim)
                      upp <- rep(q, dim)
               })
               pmvt(lower = low, upper = upp, df = df, corr = cr, ...)
        }

        switch(object$alternative, "two.sided" = {
            if (df > 0) pvals <- 2*(1 - pt(abs(tstat),df))     
            else        pvals <- 2*(1 - pnorm(abs(tstat)))
        }, "less" = {
            if (df > 0) pvals <- pt(tstat,df) 
            else        pvals <- pnorm(tstat)
        }, "greater" = {
            if (df > 0) pvals <- 1 - pt(tstat,df)
            else        pvals <- 1 - pnorm(tstat)
        })

        if (type == "univariate")
            return(pvals)

        if (type == "Bonferroni")
            return(pmin(1, dim * pvals))

        if (type == "adjusted") {
            ret <- numeric(length(tstat))
            error <- 0
            for (i in 1:length(tstat)) {
                tmp <- pfct(tstat[i])
                if (error < attr(tmp, "error")) 
                    error <- attr(tmp, "error")
                ret[i] <- tmp
            }
            ret <- 1 - ret
            attr(ret, "error") <- error
            return(ret)
        }
    }

    qfunction <- function(conf.level, ...) {

        tail <- switch(object$alternative, "two.sided" = "both.tails",
                                    "less"      = "upper.tail",
                                    "greater"   = "lower.tail")
        calpha <- qmvt(conf.level, df = df, corr = cr, tail = tail, 
                       ...)
        error <- calpha$estim.prec
        calpha <- calpha$quantile

        switch(object$alternative, "two.sided" = {  
            LowerCL <- betahat - calpha * ses
            UpperCL <- betahat + calpha * ses
        }, "less" = {
            LowerCL <- rep(-Inf, dim)
            UpperCL <- betahat - calpha * ses
        }, "greater" = {
            LowerCL <- betahat - calpha * ses
            UpperCL <- rep( Inf, dim)
        })

        cint <- cbind(LowerCL, UpperCL)
        colnames(cint) <- c("lower", "upper")
        attr(cint, "conf.level") <- conf.level
        attr(cint, "calpha") <- calpha
        attr(cint, "error") <- error
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

adjusted <- function(type = c("free", "Bonferroni", "Shaffer", "Westfall"), 
                     ...) 
{
    type <- match.arg(type)

    ### simple Bonferroni-adjustment
    if (type == "Bonferroni") {
        return(function(object) {
            RET <- pqglht(object)
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
        m <- coef(object, null = TRUE)
        tstat <- switch(object$alternative, 
                        "less" = RET$tstat,
                        "greater" = -RET$tstat,
                        "two.sided" = -abs(RET$tstat))
        C <- object$K
        Corder <- C[order(tstat), , drop = FALSE]
        Cm <- m[order(tstat)]
        ms <- maxsets(Corder)
        error <- 0
        p <- sapply(ms, function(x) {
           max(sapply(x, function(s) {
                object$K <- Corder[s, , drop = FALSE]
                object$m <- Cm[s]
                tmp <- pqglht(object)$pfunction(ifelse(type == "Westfall", 
                                                   "adjusted", "Bonferroni"), 
                                            ...)
                tmperr <- attr(tmp, "error")
                if (!is.null(tmperr) && tmperr > error)
                    error <<- tmperr
                min(tmp)
            }))
        })
        for (i in 2:length(p))
            p[i] <- max(p[i-1], p[i])
        ### <FIXME> what happens in case of ties??? </FIXME> ###
        RET$pvalues <- p[rank(tstat)]
        attr(RET$pvalues, "error") <- error
        RET$type <- type
        RET
    })
}

### copied from package MASS  
MPinv <- function (X, tol = sqrt(.Machine$double.eps))
{
    if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X)))
        stop("X must be a numeric or complex matrix")
    if (!is.matrix(X))
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X))
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    if (all(Positive)) 
        RET <- Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))   
    else if (!any(Positive))
        RET <- array(0, dim(X)[2:1])
    else RET <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
        t(Xsvd$u[, Positive, drop = FALSE]))
    return(list(MPinv = RET, rank = sum(Positive)))
}

global <- function(type = c("Chisq", "F")) {
    type <- match.arg(type)
    
    fct <- function(object) {

        RET <- pqglht(object)
        betahat <- RET$coefficients
        m <- coef(object, null = TRUE)
        covm <- vcov(object)

        tmp <- betahat - m
        MP <- MPinv(covm)
        SSH <- t(tmp) %*% MP$MPinv %*% tmp

        q <- MP$rank
        df <- df.residual(object$model)
        if (is.null(df)) {
            type <- "Chisq"
            warning(sQuote("df.residual"), " is not available for ",
                    sQuote("model"), " performing F test")
        }

        if (type == "Chisq") {
            pval <- pchisq(SSH, q, lower.tail = FALSE)
        } else {
            pval <- pf(SSH/q, q, df, lower.tail = FALSE)
        }
        RET$pvalue  <- pval
        RET$type <- type
        RET$SSH <- SSH
        RET$fstat <- SSH/q
        RET$df <- c(q, df)
        class(RET) <- "summary.glht.global"
        return(RET)
    }
    return(fct)
}
