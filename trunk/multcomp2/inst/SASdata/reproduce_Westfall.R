
library("multcomp2")

rda <- list.files(pattern = "\\.rda")
sapply(rda, function(x) load(file = x, env = .GlobalEnv))

### weights loss data, page 47, program 3.2

  amod <- aov(wloss ~ diet, data = wloss)
  amod

  gh <- glht(amod, list(diet = "Tukey"))

  # page 49 / 50
  confint(gh)

  amod <- aov(wloss ~ diet - 1, data = wloss)
  K <- diag(nlevels(wloss$diet))
  rownames(K) <- levels(wloss$diet)
  gh <- glht(amod, K)

  # page 61
  confint(gh)


### tox data, page 56, program 3.7

  tox$g <- as.factor(tox$g)
  amod <- aov(gain ~ g, data = tox)
  amod

  # page 56
  gh <- glht(amod, list(g = "Dunnett"))
  confint(gh)

  # page 59
  gh <- glht(amod, list(g = "Dunnett"), alternative = "less")
  confint(gh)

### coupon data, page 62, program 3.11

  coupon$discount <- as.factor(coupon$discount)
  amod <- aov(purchase ~ discount - 1, data = coupon)

  gh <- glht(amod, K = rbind("linear" = c(-3, -1,  1,  3),
                             "quad" =  c(-2,  2,  2, -2),
                             "cubic" = c(-1,  3, -3,  1)))

  # page 63
  summary(gh)

### recover data, page 66, program 4.1

  amod <- aov(minutes ~ blanket, data = recover)

  gh <- glht(amod, list(blanket = "Tukey"))

  # page 68
  confint(gh)

  # page 75
  summary(gh)

  gh <- glht(amod, list(blanket = "Dunnett"))

  # page 78
  confint(gh)

  # page 79
  summary(gh)

  gh <- glht(amod, list(blanket = "Dunnett"), alternative = "less")

  # page 80
  confint(gh, level = 0.9)

  # page 80
  summary(gh)

  # page 80
  amod <- aov(minutes ~ blanket - 1, data = recover)
  confint(glht(amod, K = diag(4)), level = 0.9)

### house prices, page 84, program 5.1

  amod <- aov(price ~ location + sqfeet + age, data = house)
  gh <- glht(amod, list(location = "Tukey"))

  # page 85
  confint(gh)

  # page 96
  summary(gh)
  summary(gh, test = univariate())

### rat growth data, page 99, program 5.9

  amod <- aov(W4 ~ trt + I(W0 - W3), data = ratgrwth)

  gh <- glht(amod, list(trt = "Dunnett"), alternative = "less")

  summary(gh)
  confint(gh)

### Alzheimer data, page 103, program 5.12

  alz$therapy <- as.factor(alz$therapy)
  amod <- aov(score ~ therapy * since + age, data = alz)

  gh <- glht(amod, K = list(therapy = "Tukey"))

  confint(gh)

### litter data, page 109, program 6.2.2

  litter$dose <- as.factor(litter$dose)
  amod <- aov(weight ~ dose + gesttime + number, data = litter)

  K <- rbind("cont-low"  = c(1, -1,  0,  0),
             "cont-mid"  = c(1,  0, -1,  0),
             "cont-high" = c(1,  0,  0, -1))

### hourse regression line

  houseA <- subset(house, location == "A")

  lmod <- lm(price ~ sqfeet, data = houseA)
  K <- cbind(1, grid <- seq(from = 1000, to = 3000, by = 200))
  rownames(K) <- grid

  gh <- glht(lmod, K = K)

  confint(gh)

