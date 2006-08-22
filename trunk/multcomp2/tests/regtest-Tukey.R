
library("multcomp2")

data("warpbreaks")
fm1 <- aov(breaks ~ wool * tension, data = warpbreaks)

TukeyHSD(fm1, "tension", ordered = FALSE)
confint(mcp(fm1, hypotheses = list(tension = "Tukey")))
summary(mcp(fm1, hypotheses = list(tension = "Tukey")))

TukeyHSD(fm1, "wool", ordered = FALSE)
confint(mcp(fm1, hypotheses = list(wool = "Tukey")))
summary(mcp(fm1, hypotheses = list(wool = "Tukey")))

