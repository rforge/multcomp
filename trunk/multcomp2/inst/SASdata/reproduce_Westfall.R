
library("multcomp2")

set.seed(290875)

rda <- list.files(pattern = "\\.rda")
sapply(rda, function(x) load(file = x, env = .GlobalEnv))

### weights loss data, page 47

  amod <- aov(wloss ~ diet, data = wloss)
  amod

  gh <- glht(amod, mcp(diet = "Tukey"))

  # page 49 / 50
  confint(gh)

  amod <- aov(wloss ~ diet - 1, data = wloss)
  K <- diag(nlevels(wloss$diet))
  rownames(K) <- levels(wloss$diet)
  gh <- glht(amod, K)

  # page 61
  confint(gh)


### tox data, page 56

  tox$g <- as.factor(tox$g)
  amod <- aov(gain ~ g, data = tox)
  amod

  # page 56
  gh <- glht(amod, mcp(g = "Dunnett"))
  confint(gh)

  # page 59
  gh <- glht(amod, mcp(g = "Dunnett"), alternative = "less")
  confint(gh)


### coupon data, page 62

  coupon$discount <- as.factor(coupon$discount)
  amod <- aov(purchase ~ discount - 1, data = coupon)

  gh <- glht(amod, linfct = mcp(discount = rbind(
                                    linear = c(-3, -1,  1,  3),
                                    quad =  c(-2,  2,  2, -2),
                                    cubic = c(-1,  3, -3,  1))))

  # page 63
  summary(gh)


### recover data, page 66

  amod <- aov(minutes ~ blanket, data = recover)

  gh <- glht(amod, mcp(blanket = "Tukey"))

  # page 68
  confint(gh)

  # page 75
  summary(gh)

  gh <- glht(amod, mcp(blanket = "Dunnett"))

  # page 78
  confint(gh)

  # page 79
  summary(gh)

  gh <- glht(amod, mcp(blanket = "Dunnett"), alternative = "less")

  # page 80
  confint(gh, level = 0.9)

  # page 80
  summary(gh)

  # page 80
  amod <- aov(minutes ~ blanket - 1, data = recover)
  confint(glht(amod, linfct = diag(4)), level = 0.9)


### house prices, page 84

  amod <- aov(price ~ location + sqfeet + age, data = house)
  gh <- glht(amod, mcp(location = "Tukey"))

  # page 85
  confint(gh)

  # page 96
  summary(gh)
  summary(gh, test = univariate())


### rat growth data, page 99

  amod <- aov(w4 ~ trt + I(w0 - w3), data = ratgrwth)

  gh <- glht(amod, mcp(trt = "Dunnett"), alternative = "less")

  # page 100
  summary(gh)
  confint(gh)


### Alzheimer data, page 103

  alz$therapy <- as.factor(alz$therapy)
  amod <- aov(score ~ therapy * since + age, data = alz)

  gh <- glht(amod, linfct = mcp(therapy = "Tukey"))
  gh$K[,8:11] <- gh$K[,8:11] * 10
  gh$K[,"age"] <- mean(alz$age)

  confint(gh)


### litter data, page 109

  litter$dose <- as.factor(litter$dose)
  amod <- aov(weight ~ dose + gesttime + number, data = litter)

  K <- rbind("cont-low"  = c(1, -1,  0,  0),
             "cont-mid"  = c(1,  0, -1,  0),
             "cont-high" = c(1,  0,  0, -1),
              otrend = c(1.5, 0.5, -0.5, -1.5) / 2,
              atrend = c(0, 5, 50, 500) - mean(c(0, 5, 50, 500)),
              ltrend = -(log(1:4) - mean(log(1:4))))
  K["atrend",] <- K["atrend",] / -max(K["atrend",])

  gh <- glht(amod, linfct = mcp(dose = K))

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


### house data -- regression line

  houseA <- subset(house, location == "A")

  lmod <- lm(price ~ sqfeet, data = houseA)
  K <- cbind(1, grid <- seq(from = 1000, to = 3000, by = 200))
  rownames(K) <- grid

  gh <- glht(lmod, linfct = K)

  # page 123
  confint(gh)


### patient satisfaction, page 125

  pat_sat <- pat_sat[order(pat_sat$severe),]
  lmod <- lm(satisf ~ age + severe + anxiety, data = pat_sat)
  K <- cbind(1, mean(pat_sat$age), pat_sat$severe, mean(pat_sat$anxiety))

  gh <- glht(lmod, linfct = K)

  ci <- confint(gh)

  # page 127
  plot(pat_sat$severe, ci$confint[,"Estimate"], 
       xlab = "Severity", ylab = "Satisfaction", type = "b", 
       ylim = c(30, 80), xlim = c(45, 60))
  lines(pat_sat$severe, ci$confint[,"lwr"], lty = 2)
  lines(pat_sat$severe, ci$confint[,"upr"], lty = 2)


### tire data, page 127

  amod <- aov(cost ~ make + make:mph - 1, data = tire)

  x <- seq(from = 10, to = 70, by = 5)
  K <- cbind(1, -1, x, -x)
  rownames(K) <- x

  gh <- glht(amod, linfct = K)

  # page 129
  confint(gh)


### cholesterol data, page 153

  cholesterol$trt <- as.factor(rev(as.integer(cholesterol$trt)))
  levels(cholesterol$trt) <- LETTERS[1:5]
  amod <- aov(response ~ trt - 1, data = cholesterol)
  
  gh <- glht(amod, linfct = mcp(trt = "Tukey"))

  # page 171
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall"))

  gh <- glht(amod, linfct = mcp(trt = c("D - E = 0",
                                    "C - E = 0",
                                    "C - D = 0",
                                    "3 * B - C - D - E = 0",
                                    "3 * A - C - D - E = 0")))

  # page 172
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall"))


### waste data, page 177

  waste$temp <- as.factor(waste$temp)
  waste$envir <- as.factor(waste$envir)
  amod <- aov(waste ~ temp * envir, data = waste)

  # page 179
  confint(glht(amod, linfct = mcp(temp = "Tukey")))
  confint(glht(amod, linfct = mcp(envir = "Tukey")))

  gh <- glht(amod, linfct = mcp(temp = "Tukey", envir = "Tukey"))
 
  # page 181
  confint(gh)
  
  # page 186 CHECK!
  amod <- aov(waste ~ temp + envir, data = waste)

  gh <- glht(amod, linfct = mcp(temp = "Tukey", envir = "Tukey"),
             alternative = "greater")
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall", abseps = 0.1))


### drug data, page 187

  drug$drug <- as.factor(drug$drug)
  amod <- aov(response ~ drug * disease, data = drug)

  confint(glht(amod, linfct = mcp(drug = "Tukey")))


### detergents data, page 189

  detergent$block <- as.factor(detergent$block)
  detergent$detergent <- as.factor(detergent$detergent)
  amod <- aov(plates ~ block + detergent, data = detergent)

  gh <- glht(amod, linfct = mcp(detergent = "Tukey"))

  # page 190
  confint(gh)

  # page 192
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall"))


### pigs data, page 195

  pigs$pen <- as.factor(pigs$pen)
  pigs$feed <- as.factor(pigs$feed)
  amod <- aov(gain ~ pen + feed * sex + initial, data = pigs)

  S <- matrix(c(1, -1), ncol = 2, dimnames = list("F-M", c("F", "M")))
  gh <- glht(amod, linfct = mcp(feed = "Tukey", sex = S))
  gh$K <- rbind(gh$K, "initial" = c(rep(0, 10), 1))
  gh$m <- c(gh$m, 0)

  # page 194
  confint(gh)

  # page 195
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall"))


### respiratory, page 196, program 9.14

  amod <- aov(score ~ treatment:agegroup:inithealth - 1, data = respiratory)
       
  CA  <- c(13,  0, 11,  0, 13,  0, 17,  0)
  CP  <- c( 0, 14,  0, 12,  0, 19,  0, 12)
  CA  <- CA/sum(CA)
  CP  <- CP/sum(CP)
  C1  <- CP-CA

  CAO <- c(13,  0,  0,  0, 13,  0,  0,  0) 
  CPO <- c( 0, 14,  0,  0,  0, 19,  0,  0) 
  CAO <- CAO/sum(CAO)
  CPO <- CPO/sum(CPO)
  C2  <- CPO - CAO

  CAY <- c(0,  0, 11,  0,  0,  0, 17,  0) 
  CPY <- c(0,  0,  0, 12,  0,  0,  0, 12) 
  CAY <- CAY/sum(CAY)
  CPY <- CPY/sum(CPY)
  C3  <- CPY - CAY

  CAG <- c(13,  0, 11,  0,  0,  0,  0,  0) 
  CPG <- c( 0, 14,  0, 12,  0,  0,  0,  0) 
  CAG <- CAG/sum(CAG)
  CPG <- CPG/sum(CPG)
  C4  <- CPG - CAG

  CAP <- c(0,  0,  0,  0, 13,  0, 17,  0 ) 
  CPP <- c(0,  0,  0,  0,  0, 19,  0, 12 ) 
  CAP <- CAP/sum(CAP)
  CPP <- CPP/sum(CPP)
  C5  <- CPP - CAP

  C6  <- c(-1,  1,  0,  0,  0,  0,  0,  0)
  C7  <- c( 0,  0,  0,  0, -1,  1,  0,  0)
  C8  <- c( 0,  0, -1,  1,  0,  0,  0,  0)
  C9  <- c( 0,  0,  0,  0,  0,  0, -1,  1)

  C <- rbind(C1, C2, C3, C4, C5, C6, C7, C8, C9)   
  rownames(C) <- c("Overall", "Older", "Younger", "Good Init", "Poor Init",
                   "Old x Good", "Old x Poor", "Young x Good", "Young x Poor") 

  gh <- glht(amod, linfct = C, alternative="less")

  # page 198
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))

### wine data, page 199

  wine$purchase <- as.factor(wine$purchase)
  wine$customertype <- as.factor(wine$customertype)
  wine$light <- as.factor(wine$light)
  wine$music <- as.factor(wine$music)

  amod <- glm(purchase ~ customertype + light + music + customertype:light + customertype:music
              + light:music + customertype:light:music + handle + examine, data = wine,
              family = binomial())


### wloss data, page 205

  library("lme4")
  lmod <- lmer(wloss ~ diet + (1 | i), data = wloss, model = TRUE)
  
  gh <- glht(lmod, mcp(diet = "Tukey"))

  # page 205
  confint(gh)

  # page 207 / 208
  summary(gh)


### detergent data, page 211

  detergent$block <- as.factor(detergent$block)
  detergent$detergent <- as.factor(detergent$detergent)
  lmod <- lmer(plates ~ detergent + (1 | block), data = detergent, 
               model = TRUE)

  gh <- glht(lmod, mcp(detergent = "Tukey"))

  # page 211
  confint(gh)


### waste data

  waste$temp <- as.factor(waste$temp)
  waste$envir <- as.factor(waste$envir)
  lmod <- lmer(waste ~ temp + (1 | envir) + (1 | envir : temp),
               data = waste)

  gh <- glht(lmod, mcp(temp = "Tukey"))

  # page 213
  confint(gh)


### halothane data, page 214

  lmod <- lmer(rate ~ treatment + (1 | dog), data = halothane, 
               model = TRUE)

  gh <- glht(lmod, linfct = mcp(treatment = "Tukey"))

  # page 215
  confint(gh)

  gh <- glht(lmod, 
      linfct = mcp(treatment = c("Tukey", 
            Halo = "-0.5 * HA - 0.5 * HP + 0.5 * LA + 0.5 * LP = 0",
            CO2 = "0.5 * HA -0.5 * HP + 0.5 * LA -0.5 * LP = 0",
            Interaction = "HA - HP - LA + LP = 0")))

  # page 217
  summary(gh, test = univariate())
  summary(gh, test = adjusted("Shaffer"))
  summary(gh, test = adjusted("Westfall"))
  

### multipleendpoints data, page 218

  ##multipleendpoints$endpoint <- as.factor(multipleendpoints$endpoint)
  ##lmod <- lmer(y ~ treatment:endpoint + (1 | subject), 
  ##             data = multipleendpoints, model = TRUE)

### obesity, page 220

### heart, 222
