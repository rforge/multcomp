
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

### house regression line

  houseA <- subset(house, location == "A")

  lmod <- lm(price ~ sqfeet, data = houseA)
  K <- cbind(1, grid <- seq(from = 1000, to = 3000, by = 200))
  rownames(K) <- grid

  gh <- glht(lmod, K = K)

  confint(gh)

### patient satisfaction, page 125, program 6.10

  pat_sat <- pat_sat[order(pat_sat$severe),]
  lmod <- lm(satisf ~ age + severe + anxiety, data = pat_sat)
  K <- cbind(1, mean(pat_sat$age), pat_sat$severe, mean(pat_sat$anxiety))

  gh <- glht(lmod, K = K)

  ci <- confint(gh)

  # page 127
  plot(pat_sat$severe, ci$confint[,"Estimate"], 
       xlab = "Severity", ylab = "Satisfaction", type = "b", 
       ylim = c(30, 80), xlim = c(45, 60))
  lines(pat_sat$severe, ci$confint[,"lwr"], lty = 2)
  lines(pat_sat$severe, ci$confint[,"upr"], lty = 2)

### tire data, page 127, program 6.12

  amod <- aov(cost ~ make + make:mph, data = tire)

  x <- seq(from = 10, to = 75, by = 5)
  K <- cbind(0, 1, x, -x)
  rownames(K) <- x

  # page 129
  gh <- glht(amod, K = K)

  confint(gh)

### cholesterol data, page 153, program 8.1

  amod <- aov(response ~ trt, data = cholesterol)
  
  gh <- glht(amod, K = list(trt = c("B - A = 0",
                                    "C - A = 0",
                                    "D - A = 0",
                                    "E - A = 0",
                                    "C - B = 0",
                                    "D - B = 0",
                                    "E - B = 0",
                                    "D - C = 0",
                                    "E - C = 0",
                                    "E - D = 0")))

  # page 171
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall"))

  gh <- glht(amod, K = list(trt = c("D - E = 0",
                                    "C - E = 0",
                                    "C - D = 0",
                                    "3 * B - C - D - E = 0",
                                    "3 * A - C - D - E = 0")))

  # page 172
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall"))

### waste data, page 177, program 9.1

  waste$temp <- as.factor(waste$temp)
  waste$envir <- as.factor(waste$envir)
  amod <- aov(waste ~ temp * envir, data = waste)

  # page 179
  confint(glht(amod, K = list(temp = "Tukey")))
  confint(glht(amod, K = list(envir = "Tukey")))

  gh <- glht(amod, K = list(temp = "Tukey", envir = "Tukey"))
 
  # page 181
  confint(gh)
  
  # page 186
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall"))

### drug data, page 187, program 9.7

  drug$drug <- as.factor(drug$drug)
  amod <- aov(response ~ drug * disease, data = drug)

  confint(glht(amod, K = list(drug = "Tukey")))

### detergents data, page 189, program 9.8

  detergent$block <- as.factor(detergent$block)
  detergent$detergent <- as.factor(detergent$detergent)
  amod <- aov(plates ~ block + detergent, data = detergent)

  gh <- glht(amod, K = list(detergent = "Tukey"))

  # page 190
  confint(gh)

  # page 192
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall"))

###   