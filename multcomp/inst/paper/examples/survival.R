
load(url("http://www.imbe.med.uni-erlangen.de/~hothorn/AML_Bullinger.rda", open = "r"))
risk <- rep(0, nrow(clinical))
rlev <- levels(clinical[, "Cytogenetic.group"])
risk[clinical[, "Cytogenetic.group"] %in% rlev[c(7,8,4)]] <- "low"
risk[clinical[, "Cytogenetic.group"] %in% rlev[c(5, 9)]] <- "intermediate"
risk[clinical[, "Cytogenetic.group"] %in% rlev[-c(4,5, 7,8,9)]] <- "high"
risk <- as.factor(risk)
library("survival")
library("multcomp")
plot(survfit(Surv(time, event) ~ risk, data = clinical))
cmod <- coxph(Surv(time, event) ~ risk, data = clinical)
gh <- glht(cmod, linfct = mcp(risk = "Tukey"))
summary(gh)
confint(gh)
plot(confint(gh))


