
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

### handle contrasts specified via (character) expressions
### extract left hand side of an expression
lhs <- function(ex) {

    if (length(ex) != 1)
        stop("expression is not of length 1")

    if (length(ex[[1]]) != 3 || ex[[1]][[1]] != "=")
        stop("expression has not three element or does not contain ", sQuote("="))

    return(ex[[1]][[2]])
}

### extract right hand side of an expression
rhs <- function(ex) {

    if (length(ex) != 1)
        stop("expression is not of length 1")

    if (length(ex[[1]][[3]]) == 2)
        return(-ex[[1]][[3]][[2]])

    rhs <- ex[[1]][[3]]
    if (!is.numeric(rhs) || length(rhs) > 1)
        stop("right hand side of expression is not a scalar numeric")

    rhs
}

### extract coefficients for factor levels
coefs <- function(ex) {

    if (length(ex) == 1)
        return(list(coef = 1, level = ex))

    if (length(ex) == 2)
        return(list(coef = -1, level = ex[[2]]))

    if (length(ex) == 3) {

        if (ex[[1]] == "*")
            return(list(coef = ex[[2]], level = ex[[3]]))

        cf <- coefs(ex[[3]])
        if (ex[[1]] == "-")
            cf$coef <- cf$coef * (-1)

        return(cf)
    }
}
 
### convert an expression of factor levels to a linear hypothesis
expression2K <- function(ex, y) {

    ex <- parse(text = ex)
    m <- rhs(ex)
    x <- lhs(ex)
    K <- rep.int(0, length(y))
    while(TRUE) {
        cf <- coefs(x)
        if (!any(cf$level == y))
            stop("level \"", cf$level, "\" not found")
        K[y == cf$level] <- cf$coef
        if (length(x) > 1 && !is.numeric(x[[2]])) {
            x <- x[[2]]   
        } else {
            break
        }
    }
    return(list(K = K, m = m))
}

### convert a vector of character expressions of 
### factor levels to a linear hypothesis
chr2K <- function(ex, y) {

    if (!is.character(y))
        stop(sQuote("y"), " is not a a character vector")
    if (class(ex) != "character")
        stop(sQuote(ex), " is not a character vector")

    K <- c()
    m <- c()
    for (x in ex) {
        tmp <- expression2K(parse(text = x), y)
        K <- rbind(K, tmp$K)
        m <- c(m, tmp$m)
    }
    colnames(K) <- levels(y)
    rownames(K) <- ex
    list(K = K, m = m)
}
