
source("../R/linhypo.R")
source("../R/mcp.R")
source("../R/pqfunctions.R")
source("../R/print.R")

library(multcomp)
data(waste)

mod2way <- aov(waste ~ envir * temp, data = waste)
a <- mcp(mod2way, list(temp = "Dunnett", envir = "Dunnett"))
summary(a)
confint(a)

TukeyHSD(mod2way, "temp")


mod2way <- aov(waste ~ temp, data = waste)
a <- mcp(mod2way, list(temp = "Dunnett"y)

x <- 1:100
df <- data.frame(y = x + rnorm(length(x)), x = x)

mod <- lm(y ~ x)

confint(mcp(mod, hypo = cbind(1, seq(from = 10, to = 100, by = 10))))
