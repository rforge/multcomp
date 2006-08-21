
source("linhypo.R")
source("mcp.R")
source("simint.R")
source("simtest.R")

library(multcomp)
data(waste)

mod2way <- aov(waste ~ envir * temp, data = waste)
a <- mcp(mod2way, list(temp = "Tukey"))#, envir = "Dunnett"))
confint(a)

TukeyHSD(mod2way, "temp")

