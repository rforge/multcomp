
print.glht <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n\t", "General Linear Hypotheses\n\n")
    x <- cbind(coef(x), coef(x, null = TRUE))
    colnames(x) <- c("Estimate", "Hypothesis")
    cat("Linear Hypotheses:\n")
    print(x, digits = digits)
    cat("\n")
    invisible(x)
}

print.summary.glht <- function(x, digits = max(3, getOption("digits") - 3), 
                              ...) 
{
    cat("\n\t", "Multiple Tests for General Linear Hypotheses\n\n")
    cat("Fit: ")
    if (inherits(x$model, "lmer")) {
        print(x$model@call)
    } else {
        print(x$model$call)
    }
    cat("\n")
    type <- attr(x$mtests, "type")
    attr(x$mtests, "type") <- NULL
    ### print p values according to simulation precision
    error <- attr(x$mtests, "error")
    if (!is.null(error) && error > .Machine$double.eps) {
        sig <- which.min(abs(1 / error - (10^(1:10))))
        sig <- 1 / (10^sig)
    } else {
        sig <- .Machine$double.eps
    }
    cat("Linear Hypotheses:\n")
    printCoefmat(x$mtests, digits = digits, 
                 has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)
    switch(type, 
        "univariate" = cat("(Univariate p values reported)"),
        "Bonferroni" = cat("(Bonferroni-adjusted p values reported)"),
        "free" = cat("(Adjusted p values reported)"),
        "Shaffer" = cat("(Adjusted p values reported -- Shaffer method)"),
        "Westfall" = cat("(Adjusted p values reported -- Westfall method)")
    )
    cat("\n\n")
    invisible(x)                    
}

print.confint.glht <- function(x, digits = max(3, getOption("digits") - 3), 
                              ...) 
{
    cat("\n\t", "Simultaneous Confidence Intervals for General Linear Hypotheses\n\n")
    level <- attr(x$confint, "conf.level")
    attr(x$confint, "conf.level") <- NULL
    cat("Fit: ")
    if (inherits(x$model, "lmer")) {
        print(x$model@call)
    } else {
        print(x$model$call)
    }
    cat("\n")
    error <- attr(x$confint, "error")
    if (!is.null(error) && error > .Machine$double.eps)
        digits <- min(digits, which.min(abs(1 / error - (10^(1:10)))))
    cat("Estimated Quantile =", round(attr(x$confint, "calpha"), digits))
    cat("\n\n")
    cat("Linear Hypotheses:\n")
    print(format(x$confint, nsmall = digits, digits = digits), quote = FALSE)
    cat("\n")
    cat(paste(level * 100, 
              "% family-wise confidence level\n", sep = ""), "\n\n")
    invisible(x)
}

print.contrMat <- function(x, digits = max(3, getOption("digits") - 3), ...) {

    cat("\n\t", "General Linear Hypotheses Matrix\n")
    cat("\t\t", attr(x, "type"), "Contrasts\n\n")
    attr(x, "type") <- NULL
    class(x) <- "matrix"  
    print(x, digits = digits)
    invisible(x)
}

print.summary.glht.global <- function(x, 
    digits = max(3, getOption("digits") - 3), ...) {

    print.glht(x, digits = digits)
    cat("Global Test:\n")
    if (x$gtest$type == "Chisq") {
        pr <- data.frame(x$gtest$SSH, x$gtest$df[1], x$gtest$pval)
        names(pr) <- c("Chisq", "DF", "Pr(>Chisq)")
    }
    if (x$gtest$type == "F") {
        pr <- data.frame(x$gtest$fstat, x$gtest$df[1], x$gtest$df[2], 
                         x$gtest$pval)
        names(pr) <- c("F", "DF1", "DF2", "Pr(>F)")
    }
    print(pr, digits = digits)
    invisible(x)
}
