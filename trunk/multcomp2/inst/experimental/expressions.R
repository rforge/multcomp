
### extract coefficients and variable names
coefs <- function(ex) {

    if (length(ex) == 1)
        return(list(coef = 1, var = as.character(ex)))

    if (length(ex) == 2)
        return(list(coef = -1, var = as.character(ex[[2]])))

    if (length(ex) == 3) {

        if (ex[[1]] == "*")
            return(list(coef = eval(ex[[2]]), var = as.character(ex[[3]])))

        cf <- coefs(ex[[3]])
        if (ex[[1]] == "-") 
            cf$coef <- cf$coef * (-1)

        return(cf)
    }
}

### extract coefficients and variable names
expression2coef <- function(ex) {

    cf <- c()
    nm <- c()

    m <- rhs(ex)
    x <- lhs(ex)
    alternative <- as.character(ex[[1]][[1]])
    alternative <- switch(alternative, "<=" = "greater",
                          ">=" = "less",
                          "=" = "two.sided")

    while(TRUE) {

        tmp <- coefs(x)
        cf <- c(cf, tmp$coef)
        nm <- c(nm, tmp$var)

        ### x == "A"
        if (is.name(x)) break
        ### x == "-1"
        if (is_num(x)) break
        ### x == "3 * A"
        if (length(x) == 3 && is_num(x[[2]])) break
        ### x == "-3 * A"
        if (length(x) == 3 && is_num(x[[2]])) break
        x <- x[[2]]   
    }
    return(list(K = cf, nm = nm, m = m, alternative = alternative))
}

chr2matrix <- function(ex, var) {

    alternative <- c()
    K <- c()
    m <- c()
    for (x in ex) {
        tmp <- expression2coef(parse(text = x))
        ct <- rep(0, length(var))
        for (n in tmp$nm)
            ct[var == n] <- tmp$K[tmp$nm == n]
        K <- rbind(K, ct)
        m <- c(m, tmp$m)
        alternative <- c(alternative, tmp$alternative)
    }
    colnames(K) <- var
    rownames(K) <- ex
    list(K = K, m = m, alternative = unique(alternative))
}

is_num <- function(x) {
    opt <- options("show.error.messages")
    options(show.error.messages = FALSE)
    tr <- try(eval(x))
    options(opt)
    if (!inherits(tr, "try-error")) return(is.numeric(tr))
    return(FALSE)
}

### handle contrasts specified via (character) expressions
### extract left hand side of an expression
lhs <- function(ex) {

    if (length(ex) != 1)
        stop("expression is not of length 1")

    if (length(ex[[1]]) != 3 || !(as.character(ex[[1]][[1]]) %in% c("<=", ">=", "=")))
        stop("expression has not three element or does not contain ", sQuote("<=, >=, ="))


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
