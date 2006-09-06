
print.glht <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\n\t", "Multiple Comparison Procedures for Linear Hypotheses\n\n")
    x <- x$K %*% x$beta
    colnames(x) <- "Linear Hypotheses:"
    print(x, digits = digits)
    cat("\n")
    invisible(x)
}

print.summary.glht <- function(x, digits = max(3, getOption("digits") - 3), 
                              ...) 
{
    cat("\n\t", "Multiple Tests for General Linear Hypotheses\n\n")
    cat("Fit: ")
    print(x$model$call)
    cat("\n")
    type <- attr(x$mtests, "type")
    attr(x$mtests, "type") <- NULL
    cat("Linear Hypotheses:\n")
    printCoefmat(x$mtests, digits = digits, 
                 has.Pvalue = TRUE, P.values = TRUE)
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
    print(x$object$call)
    cat("\n")
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
