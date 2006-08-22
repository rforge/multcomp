
print.mcp <- function(x, ...) {
    cat("\n\t", "Multiple Comparison Procedures for Linear Hypotheses\n\n")
    x <- x$hypotheses %*% x$beta
    colnames(x) <- "Estimate"
    print(x)
}

print.summary.mcp <- function(x, ...) {

    cat("\n\t", "Multiple Comparison Procedures for Linear Hypotheses\n\n")
    cat("Fit: ")
    print(x$object$call)
    cat("\n")
    type <- attr(x$mcp, "type")
    attr(x$mcp, "type") <- NULL
    print(x$mcp)
    cat("\n")
    switch(type, 
        "raw" = cat("\nRaw p values reported"),
        "Bonferroni" = cat("\nBonferroni-adjusted p values reported"),
        "adjusted" = cat("Adjusted p values reported"))
    cat("\n\n")
    invisible(x)                    
}

print.confint.mcp <- function(x, ...) {

    cat("\n\t", "Multiple Comparison Procedures for Linear Hypotheses\n\n")
    level <- attr(x$confint, "conf.level")
    attr(x$confint, "conf.level") <- NULL
    cat("Fit: ")
    print(x$object$call)
    cat("\n")
    print(x$confint)    
    cat("\n")
    cat(paste(level * 100, 
              "% family-wise confidence level\n", sep = ""), "\n\n")
    invisible(x)

}
