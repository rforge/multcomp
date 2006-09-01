
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

adjusted <- function(type = c("free", "Bonferroni", "Shaffer", "Westfall"), ...) 
{
    type <- match.arg(type)

    ### simple Bonferroni-adjustment
    if (type == "Bonferroni") {
        return(function(object) {
            RET <- pqmcp(object)
            RET$pvalues <- RET$pfunction("Bonferroni")
            RET$type <- "Bonferroni"
            RET
        })
    }

    ### usual max-type adjustment over all linear hypotheses
    if (type == "free") {
        return(function(object) {
            RET <- pqmcp(object)
            RET$pvalues <- RET$pfunction("adjusted", ...)
            RET$type <- "adjusted (free)"
            RET
        })
    }

    ### Westfall (1997, JASA): constraints and correlations 
    ### or
    ### Shaffer (1886, JASA): constraints
    return(function(object) {
        RET <- pqmcp(object)
        tstat <- switch(object$alternative, 
                        "less" = RET$tstat,
                        "greater" = -RET$tstat,
                        "two.sided" = -abs(RET$tstat))
        C <- object$hypotheses
        Corder <- C[order(tstat), ]
        ms <- maxsets(Corder)
        p <- sapply(ms, function(x) {
           max(sapply(x, function(s) {
                object$hypotheses <- Corder[s, , drop = FALSE]
                min(pqmcp(object)$pfunction(ifelse(type == "Westfall", 
                                                   "adjusted", "Bonferroni"), 
                                            ...))
            }))
        })
        for (i in 2:length(p))
            p[i] <- max(p[i-1], p[i])
        ### <FIXME> what happens in case of ties??? </FIXME> ###
        RET$pvalues <- p[rank(tstat)]
        RET$type <- ifelse(type == "Westfall", 
            "adjusted (Westfall -- constraints and correlations)",
            "adjusted (Shaffer -- constraints with Bonferroni)")
        RET
    })
}

### compute all possible (ordered) subsets of the index set K
### cf. Westfall (1997, Section 3.2)
allsubsets <- function(K) 
{
    if (length(K) == 0) return(list(NULL))
    if (length(K) == 1) return(list(K))
    ret <- as.list(K)
    for (i in 1:(length(K)-1)) {
        tmp <- allsubsets(K[-(1:i)])
        for (j in 1:length(tmp))
            tmp[[j]] <- c(K[i], tmp[[j]])
        ret <- c(ret, tmp)
    }
    ret
}

### check if any of C[,1:(min(K)-1)] is in column space of C[,K]
### cf. Westfall (1997, Section 3.2)
checkCS <- function(K, C) 
{
    if (length(K) == ncol(C)) return(TRUE)
    CK <- C[,K,drop = FALSE]
    Cj <- C[,1:(min(K)-1), drop = FALSE]
    tmp <- Cj - (CK %*% ginv(CK) %*% Cj)
    all(colSums(tmp^2) > .Machine$double.eps)
}

### remove redundant index sets
rmsets <- function(sets) 
{
    if (length(sets) == 1) return(sets)
    rm <- logical(length(sets))
    for (j in 1:(length(sets) - 1)) {
        set <- sets[[j]]
        rm[j] <- any(sapply(sets[(j + 1):length(sets)], function(x) 
                            all(set %in% x)))
    }
    sets[!rm]
}

### compute maximal sets of linear hypotheses
### cf. Westfall (1997, Section 3.2)
maxsets <- function(hypotheses) 
{
    C <- t(hypotheses)
    k <- ncol(C)
    p <- nrow(C)
    S <- 1:k
    ret <- vector(mode = "list", length = k)

    for (j in S) {
        tmp <- allsubsets(S[-(1:j)])
        for (i in 1:length(tmp))
            tmp[[i]] <- c(j, tmp[[i]])
        if (length(tmp) > 1 || length(tmp[[1]]) > 1)
            tmp <- c(j, tmp)
        ret[[j]] <- tmp[sapply(tmp, checkCS, C = C)]
        
    }
    lapply(ret, rmsets)
}
