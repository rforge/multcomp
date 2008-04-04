

data("alpha", package = "coin")
library("multcomp")
library("sandwich")

amod <- aov(elevel ~ alength, data = alpha)
confint(glht(amod, linfct = mcp(alength = "Tukey")))

class(amod) <- "lm"
confint(glht(amod, linfct = mcp(alength = "Tukey"), vcov = sandwich))

