
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

y <- gl(3, 3)
levels(y) <- paste("x", 1:3, sep = "")
max(abs(multcomp2:::chr2K(c(l1 = "x1 - x2 = 2", l2 = "x2 + 3 * x3 = 1"), y)$K - 
        rbind(c(1, -1, 0), c(0, 1, 3))))

