
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
    tmp <- Cj - (CK %*% MPinv(CK)$MPinv %*% Cj)
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
maxsets <- function(K) 
{
    C <- t(K)
    k <- ncol(C)
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

### oh dear!
### Cox models don't have any intercept ...
model.matrix.coxph <- function(object, ...) 
{
    mm <- model.matrix.default(object)
    at <- attributes(mm)
    mm <- mm[,-1]
    at$dim[2] <- at$dim[2] - 1
    at$dimnames[[2]] <- at$dimnames[[2]][-1]
    at$assign <- at$assign[-1]
    attributes(mm) <- at
    mm
}

### some methods of lmer objects
model.matrix.lmer <- function(object, ...) {
    x <- object@X
    if (is.null(x))
        stop("models of class ", sQuote("lmer"), 
             " need to be fitted with argument ", 
             sQuote("model = TRUE"))
    x
}

model.frame.lmer <- function(object, ...) {
    x <- object@frame
    if (is.null(x))
        stop("models of class ", sQuote("lmer"), 
             " need to be fitted with argument ", 
             sQuote("model = TRUE"))
    x
}

terms.lmer <- function(object, ...) {
    x <- object@terms
    if (is.null(x))
        stop("models of class ", sQuote("lmer"), 
             " need to be fitted with argument ", 
             sQuote("model = TRUE"))
    x
}

coeflmer <- function(object, ...) {
    x <- object@fixef
    names(x) <- rownames(vcov(object))
    x
}

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
    sigma <- as.matrix(sigma)

    ### check if a linear model was supplied
    df <- 0
    if (class(model)[1] %in% c("aov", "lm")) {
        class(model) <- "lm"
        df <- summary(model)$df[2]
    }

    list(beta = beta, sigma = sigma, df = df)
}

### modified from package MASS  
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
