
### uhhh -- mainly copy and paste from plot.TukeyHSD
plot.confint.mcp <- function(x, ...) {

    xi <- x$confint
    yvals <- nrow(xi):1
    plot(c(xi[, "lwr"], xi[, "upr"]), rep.int(yvals, 2), 
         type = "n", axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, ...)
    axis(2, at = nrow(xi):1, labels = dimnames(xi)[[1]], 
         srt = 0, ...)
    abline(h = yvals, lty = 1, lwd = 0, col = "lightgray")
    abline(v = 0, lty = 2, lwd = 0, ...)
    segments(xi[, "lwr"], yvals, xi[, "upr"], yvals, ...)
    segments(as.vector(xi), rep.int(yvals - 0.1, 3), as.vector(xi), 
        rep.int(yvals + 0.1, 3), ...)
    title(main = paste(format(100 * attr(x$confint, "conf.level"), 
          2), "% family-wise confidence level\n", sep = ""), 
          xlab = "Linear Hypotheses")
    box()
}
