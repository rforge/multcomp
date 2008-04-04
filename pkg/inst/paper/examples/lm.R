
bmod <- lm(DEXfat ~ ., data = bodyfat)
p <- length(coef(bmod))
K <- diag(p)[-1,]
rownames(K) <- names(coef(bmod))[-1]
summary(glht(bmod, lin = K))
summary(bmod)

