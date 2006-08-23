
print.mcp <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\n\t", "Multiple Comparison Procedures for Linear Hypotheses\n\n")
    x <- x$hypotheses %*% x$beta
    colnames(x) <- "Linear Hypotheses:"
    print(x, digits = digits)
    cat("\n")
    invisible(x)
}

print.summary.mcp <- function(x, digits = max(3, getOption("digits") - 3), 
                              ...) 
{
    cat("\n\t", "Multiple Comparison Procedures for Linear Hypotheses\n\n")
    cat("Fit: ")
    print(x$object$call)
    cat("\n")
    type <- attr(x$mcp, "type")
    attr(x$mcp, "type") <- NULL
    cat("Linear Hypotheses:\n")
    printCoefmat(x$mcp, digits = digits, 
                 has.Pvalue = TRUE, P.values = TRUE)
    switch(type, 
        "univariate" = cat("(Univariate p values reported)"),
        "Bonferroni" = cat("(Bonferroni-adjusted p values reported)"),
        "adjusted" = cat("(Adjusted p values reported)"))
    cat("\n\n")
    invisible(x)                    
}

print.confint.mcp <- function(x, digits = max(3, getOption("digits") - 3), 
                              ...) 
{
    cat("\n\t", "Multiple Comparison Procedures for Linear Hypotheses\n\n")
    level <- attr(x$confint, "conf.level")
    attr(x$confint, "conf.level") <- NULL
    cat("Fit: ")
    print(x$object$call)
    cat("\n")
    cat("Linear Hypotheses:\n")
    print(x$confint, digits = digits)    
    cat("\n")
    cat(paste(level * 100, 
              "% family-wise confidence level\n", sep = ""), "\n\n")
    invisible(x)
}
