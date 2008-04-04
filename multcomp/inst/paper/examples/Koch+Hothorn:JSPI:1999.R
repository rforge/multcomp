
library("multcomp")

### Table 5
dose <- c("control", "low10", "med50", "high100")
mice <- data.frame(ndead = c(4, 1, 6, 8), 
                   nalive = c(36, 19, 14, 12), 
                   dose = factor(dose, levels = dose))

gmod <- glm(cbind(ndead, nalive) ~ dose, data = mice, family = binomial())

lh <- glht(gmod, linfct = mcp(dose = "Dunnett"), alternative = "greater")
summary(lh)
confint(lh)

