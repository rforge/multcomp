# $Id: simint.R,v 1.52 2005/07/25 15:25:05 hothorn Exp $

csimint <- function(estpar, df, covm, cmatrix=NULL, ctype="user-defined",
                    conf.level=0.95,
                    alternative=c("two.sided","less","greater"), asympt=FALSE,
                    eps=0.001, maxpts=1000000, ...)
{
    if (!is.vector(estpar) & !is.matrix(estpar)) stop("estpar not a vector")
    p <- length(estpar)
    if (missing(df) & !asympt) {
      stop("df is missing")
    } else {
      if (missing(df)) df <- 0
      if (!all.equal(df - floor(df), 0)) stop("df not an integer")
    }
    if (!is.matrix(covm)) {
      if (length(covm) == 1) 
        covm <- as.matrix(covm)
      else 
        stop("covm is not a matrix")
    }
    cm <- cmatrix
    if (ctype !="user-defined") cmatrix <- NULL

    alternative <- match.arg(alternative)

    if (asympt) df <- 0                          
 
    covm  <- cm %*% covm %*% t(cm)
    d     <- 1/sqrt(diag(covm))                            
    if (length(d) > 1)
      d <- diag(d)              
    cr    <- d %*% covm %*% d

    ests  <- cm %*% estpar
    ses   <- sqrt(diag(covm))
    tvals <- ests/ses
    dim   <- ncol(cr)

    # compute the p-values

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
           },)
           pmvt(lower=low, upper=upp, df=df, corr=cr,
                abseps=eps, maxpts=maxpts)
    }

    switch(alternative, "two.sided" = {
        if (df>0) rawp <- 2*(1-pt(abs(tvals),df))     
        else      rawp <- 2*(1-pnorm(abs(tvals)))
    }, "less" = {
        if (df>0) rawp <- pt(tvals,df) 
        else      rawp <- pnorm(tvals)
    }, "greater" = {
       if (df>0) rawp <- 1-pt(tvals,df)
       else      rawp <- 1-pnorm(tvals)
    },)

    adjp <- 1-apply(tvals, 1, pfct)
    bonp <- pmin(1,dim*rawp)

    # and the simultaneous confidence intervals

    tail <- switch(alternative, "two.sided" = "both.tails",
                                "less"      = "upper.tail",
                                "greater"   = "lower.tail")
    calpha <- qmvt(conf.level, df = df, corr = cr, abseps = eps / 10, 
                   maxpts = maxpts, tail = tail)$quantile

    switch(alternative, "two.sided" = {  
        LowerCL <- ests - calpha*ses
        UpperCL <- ests + calpha*ses
    }, "less" = {
        LowerCL <- rep(-Inf, dim)
        UpperCL <- ests - calpha*ses
    }, "greater" = {
        LowerCL <- ests - calpha*ses
        UpperCL <- rep( Inf, dim)
    },)

    cint <- cbind(LowerCL, UpperCL)
    colnames(cint) <- c("lower", "upper")
    attr(cint, "conf.level") <- conf.level

    RET <- list(cmatrix = cm, ctype = ifelse(is.null(cmatrix), ctype, NA), 
                estimate = ests, sd = ses, statistics = tvals,
                p.value.raw = rawp, p.value.bon = bonp,
                p.value.adj = adjp, conf.int = cint, eps=eps, calpha=calpha,
                asympt = asympt, alternative = alternative)
    class(RET) <- "hmtest"
    RET
}
