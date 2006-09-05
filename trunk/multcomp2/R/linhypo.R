# $Id$

linhypo <- function(x, type = c("Dunnett", "Tukey", "Sequen", "AVE",
                                "Changepoint", "Williams", "Marcus",
                                "McDermott"), base = 1) {

    if (!is.factor(x)) stop(sQuote("x"), "is not a factor!")
    if (nlevels(x) < 2) stop("less than 2 groups!")
    if (any(tabulate(x) < 2)) 
        stop("less than 2 observations in at least one group!")
    k <- nlevels(x)
    n <- tabulate(x)
    if (base < 1 || base > k) stop("base is not between 1 and ", k)

    RET <- c()
    rnames <- c()
    if (!is.null(levels(x)))
        varnames <- levels(x)
    else 
        varnames <- 1:length(n)

    kindx <- 1:k
    type <- match.arg(type)

    switch(type, "Dunnett" = {
        for(i in kindx[-base])
            RET <- rbind(RET, as.numeric(kindx == i) - as.numeric( kindx == base))
        rnames <- paste(varnames[kindx[-base]], "-", varnames[base], sep="")
    }, "Tukey" = {
        for (i in 1:(k-1)) {
            for(j in (i+1):k) {
                RET  <- rbind(RET, as.numeric(kindx==j)-as.numeric(kindx==i))
                rnames <- c(rnames, paste(varnames[j], "-", varnames[i],
                                          sep=""))
            }
        }
    }, "Sequen" =  {
        for (i in 2:k) {
            RET  <- rbind(RET, as.numeric(kindx==i)-as.numeric(kindx==i-1))
            rnames <- c(rnames, paste(varnames[i], "-", varnames[i-1],
                                      sep=""))
        }
    }, "AVE" = {
        help <- c(1,  -n[2:k]/sum(n[2:k]))
        RET <- rbind(RET, help)
        for (i in 2:(k-1)) {
            x <- sum(n[1:(i-1)])+sum(n[(i+1):k])
            help <- c(-n[1:(i-1)]/x, 1, -n[(i+1):k]/x)
            RET <- rbind(RET, help)
        }
        help <- c(-n[1:(k-1)]/sum(n[1:(k-1)]), 1)
        RET  <- rbind(RET, help)
        rnames <- paste("C", 1:nrow(RET))
    }, "Changepoint" = {
        for (i in 1:(k-1)) {
            help <- c(-n[1:i]/sum(n[1:i]), n[(i+1):k]/sum(n[(i+1):k]))
            RET <- rbind(RET, help)
        }
        rnames <- c(rnames, paste("C", 1:nrow(RET), sep=""))
    }, "Williams" = {
        for (i in 1:(k-2)) {
            help <-  c(-1, rep(0, k-i-1), n[(k-i+1):k]/sum(n[(k-i+1):k]))
            RET <- rbind(RET, help)
        }
        help <- c(-1, n[2:k]/sum(n[2:k]))
        RET <- rbind(RET, help)
        rnames <- c(rnames, paste("C", 1:nrow(RET), sep=""))
    }, "Marcus" = {
        cm1 <- matrix(0, nrow=k-1, ncol=k)
        cm2 <- cm1
        for (i in 1:(k-1)) {
            cm1[i,(i+1):k] <- n[(i+1):k]/sum(n[(i+1):k])
            cm2[i,1:i] <- n[1:i]/sum(n[1:i])
        }
        row <- k*(k-1)/2
        index <- 1
        for (i in 1:(k-1)) {
            for (j in 1:i) {
                help <- cm1[i,]-cm2[j,]
                RET <- rbind(RET, help)
                index <- index+1
            }
        }
        rnames <- c(rnames, paste("C", 1:nrow(RET), sep=""))
     }, "McDermott" = {
         for(i in 1:(k-2)) {
             help  <- c(-n[1:i]/sum(n[1:i]), 1, rep(0, k-i-1))
             RET <- rbind(RET, help)
         }
         help <- c(-n[1:(k-1)]/sum(n[1:(k-1)]), 1)
         RET  <- rbind(RET, help)
         rnames <- c(rnames, paste("C", 1:nrow(RET), sep=""))
    })
    rownames(RET) <- rnames
    colnames(RET) <- varnames
    attr(RET, "type") <- type
    class(RET) <- "linhypo"
    RET
}

print.linhypo <- function(x, ...) {

    cat("\n\t", paste("Matrix of Linear Hypotheses (", 
                      attr(x, "type"), ")\n\n", sep = ""))
    attr(x, "type") <- NULL
    class(x) <- "matrix"
    print(x)
}
