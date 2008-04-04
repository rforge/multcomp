
data("alzheimer", package = "coin")
library("multcomp")

y <- alzheimer$disease == "Alzheimer's"

gmod <- glm(y ~ smoking * gender, data = alzheimer, family = binomial())
plot(confint(glht(gmod, linfct = mcp(smoking = "Sequen"))))

a <- cbind(levels(alzheimer$smoking), "Female")
b <- cbind(levels(alzheimer$smoking), "Male")
d <- rbind(a, b)
colnames(d) <- c("smoking", "gender")
rownames(d) <- paste(d[,1], d[,2], sep = ":")
d <- as.data.frame(d)

K <- model.matrix(~ smoking * gender, data = d)

ci <- confint(glht(gmod, linfct = K))
ci$confint <- apply(ci$confint, 2, binomial()$linkinv)
plot(ci)
