
library("multcomp2")

set.seed(290875)

data("warpbreaks")
fm1 <- aov(breaks ~ wool * tension, data = warpbreaks)

TukeyHSD(fm1, "tension", ordered = FALSE)
confint(glht(fm1, K = list(tension = "Tukey")))
summary(glht(fm1, K = list(tension = "Tukey")))

TukeyHSD(fm1, "wool", ordered = FALSE)
confint(glht(fm1, K = list(wool = "Tukey")))
summary(glht(fm1, K = list(wool = "Tukey")))
