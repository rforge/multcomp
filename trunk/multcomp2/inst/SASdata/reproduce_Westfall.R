
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

  # page 100
  summary(gh)
  confint(gh)

### Alzheimer data, page 103, program 5.12

  alz$therapy <- as.factor(alz$therapy)
  amod <- aov(score ~ therapy * since + age, data = alz)

  gh <- glht(amod, K = list(therapy = "Tukey"))
  gh$K[,8:11] <- gh$K[,8:11] * 10
  gh$K[,"age"] <- mean(alz$age)

  confint(gh)

### litter data, page 109, program 6.2.2

  litter$dose <- as.factor(litter$dose)
  amod <- aov(weight ~ dose + gesttime + number, data = litter)

  K <- rbind("cont-low"  = c(1, -1,  0,  0),
             "cont-mid"  = c(1,  0, -1,  0),
             "cont-high" = c(1,  0,  0, -1),
              otrend = c(1.5, 0.5, -0.5, -1.5) / 2,
              atrend = c(0, 5, 50, 500) - mean(c(0, 5, 50, 500)),
              ltrend = -(log(1:4) - mean(log(1:4))))
  K["atrend",] <- K["atrend",] / -max(K["atrend",])

  gh <- glht(amod, K = list(dose = K))

  # page 110
  summary(gh, test = univariate())

  # page 111
  gh$alternative <- "greater"
  summary(gh, test = univariate())
  summary(gh)
  confint(gh)

  # page 174
  gh$alternative <- "greater"
  summary(gh, test = adjusted("Westfall"))

### house regression line

  houseA <- subset(house, location == "A")

  lmod <- lm(price ~ sqfeet, data = houseA)
  K <- cbind(1, grid <- seq(from = 1000, to = 3000, by = 200))
  rownames(K) <- grid

  gh <- glht(lmod, K = K)

  # page 123
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

  amod <- aov(cost ~ make + make:mph - 1, data = tire)

  x <- seq(from = 10, to = 70, by = 5)
  K <- cbind(1, -1, x, -x)
  rownames(K) <- x

  gh <- glht(amod, K = K)

  # page 129
  confint(gh)

### cholesterol data, page 153, program 8.1

  cholesterol$trt <- as.factor(rev(as.integer(cholesterol$trt)))
  levels(cholesterol$trt) <- LETTERS[1:5]
  amod <- aov(response ~ trt - 1, data = cholesterol)
  
  gh <- glht(amod, K = list(trt = "Tukey"))

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
  
  # page 186 CHECK!
  amod <- aov(waste ~ temp + envir, data = waste)

  gh <- glht(amod, K = list(temp = "Tukey", envir = "Tukey"),
             alternative = "greater")
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

### pigs data, page 195, program 9.13

  pigs$pen <- as.factor(pigs$pen)
  pigs$feed <- as.factor(pigs$feed)
  amod <- aov(gain ~ pen + feed * sex + initial, data = pigs)

  S <- matrix(c(1, -1), ncol = 2, dimnames = list("F-M", c("F", "M")))
  gh <- glht(amod, K = list(feed = "Tukey", sex = S))
  gh$K <- rbind(gh$K, "initial" = c(rep(0, 10), 1))
  gh$m <- c(gh$m, 0)

  # page 194
  confint(gh)

  # page 195
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall"))

### respiratory, page 196, program 9.14

         ### three-way ANOVA without intercept
       amod <- aov(Score ~ Treatment:AgeGroup:InitHealth - 1, data = respiratory)
       
       ### compute weighted linear hypotheses in several steps 
       ### overall active vs. placebo
       CA  <- c(13,  0, 11,  0, 13,  0, 17,  0)
       CP  <- c( 0, 14,  0, 12,  0, 19,  0, 12)
       CA  <- CA/sum(CA)
       CP  <- CP/sum(CP)
       C1  <- CP-CA

       ### for older subgroup only
       CAO <- c(13,  0,  0,  0, 13,  0,  0,  0) 
       CPO <- c( 0, 14,  0,  0,  0, 19,  0,  0) 
       CAO <- CAO/sum(CAO)
       CPO <- CPO/sum(CPO)
       C2  <- CPO - CAO

       ### for younger subgroup only 
       CAY <- c(0,  0, 11,  0,  0,  0, 17,  0) 
       CPY <- c(0,  0,  0, 12,  0,  0,  0, 12) 
       CAY <- CAY/sum(CAY)
       CPY <- CPY/sum(CPY)
       C3  <- CPY - CAY

       ### subgroup with inital good health
       CAG <- c(13,  0, 11,  0,  0,  0,  0,  0) 
       CPG <- c( 0, 14,  0, 12,  0,  0,  0,  0) 
       CAG <- CAG/sum(CAG)
       CPG <- CPG/sum(CPG)
       C4  <- CPG - CAG

       ### subgroup with inital poor health
       CAP <- c(0,  0,  0,  0, 13,  0, 17,  0 ) 
       CPP <- c(0,  0,  0,  0,  0, 19,  0, 12 ) 
       CAP <- CAP/sum(CAP)
       CPP <- CPP/sum(CPP)
       C5  <- CPP - CAP

       ### all 4 subgroup combinations of age and initial health condition 
       C6  <- c(-1,  1,  0,  0,  0,  0,  0,  0)
       C7  <- c( 0,  0,  0,  0, -1,  1,  0,  0)
       C8  <- c( 0,  0, -1,  1,  0,  0,  0,  0)
       C9  <- c( 0,  0,  0,  0,  0,  0, -1,  1)

       ### matrix of linear hypotheses
       C <- rbind(C1, C2, C3, C4, C5, C6, C7, C8, C9)   
       rownames(C) <- c("Overall", "Older", "Younger", "Good Init", "Poor Init",
                        "Old x Good", "Old x Poor", "Young x Good", "Young x Poor") 

       ### set up one-sided multiple comparisons
       rht <- glht(amod, K = C, alternative="less")

       ### see Westfall et al. (1999, page 198)
       summary(rht, test = univariate())
       summary(rht, test = adjusted("Shaffer"))

### wine data, page 199, program 9.16

  wine$Purchase <- as.factor(wine$Purchase)
  wine$CustomerType <- as.factor(wine$CustomerType)
  wine$light <- as.factor(wine$light)
  wine$music <- as.factor(wine$music)

  amod <- glm(Purchase ~ CustomerType + light + music + CustomerType:light + CustomerType:music
              + light:music + CustomerType:light:music + handle + examine, data = wine,
              family = binomial())

### wloss data, page 205, program 10.1

  library("lme4")
  lmod <- lmer(wloss ~ diet + (1 | i), data = wloss, model = TRUE)
  
  gh <- glht(lmod, list(diet = "Tukey"))

  # page 205
  confint(gh)

  # page 207 / 208
  summary(gh)

