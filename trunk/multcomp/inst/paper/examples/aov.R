
data("alpha", package = "coin")
library("multcomp")

amod <- aov(elevel ~ alength, data = alpha)
confint(glht(amod, linfct = mcp(alength = "Tukey")))
