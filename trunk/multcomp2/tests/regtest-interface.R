
library("multcomp2")

testdata <- data.frame(y = rnorm(21), 
                       f1 <- factor(c(rep(c("A", "B", "C"), 7))),
                       f2 <- factor(c(rep("D", 10), rep("E", 11))),
                       x <- rnorm(21))

# one-way ANOVA
coef(amod <- aov(y ~ f1, data = testdata))
glht(amod, K = list(f1 = "Dunnett"))

# and a continuous covariable: ANCOVA
coef(lmod <- lm(y ~ f1 + x, data = testdata))
glht(lmod, K = list(f1 = "Dunnett"))

# ANCOVA with an additional factor as covariable
coef(lmod <- lm(y ~ f1 + f2 + x, data = testdata))
glht(lmod, K = list(f1 = "Dunnett"))

# and with interaction terms
coef(lmod <- lm(y ~ f1 + f2 + f2:f1 + x, data = testdata))
glht(lmod, K = list(f1 = "Dunnett"))

# with contrasts as expressions
glht(lmod, K = list(f1 = c("B - A = 0", "C - A = 0")))

tmp <- multcomp2:::chr2K(c(l1 = "x1 - x2 = 2", 
                          l2 = "x2 + 3 * x3 = 1"), 
                          paste("x", 1:3, sep = ""))

stopifnot(max(abs(tmp$K - rbind(c(1, -1, 0), c(0, 1, 3)))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(tmp$m - c(2, 1))) < sqrt(.Machine$double.eps))
