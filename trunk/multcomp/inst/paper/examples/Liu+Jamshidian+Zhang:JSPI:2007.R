
library("multcomp")

### systolic blood pressure data (source: Kleinbaum et al.)
bp <- structure(list(gender = structure(as.integer(c(2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1)), .Label = c("male", "female"), class = "factor"), sbp =
as.integer(c(144,
138, 145, 162, 142, 170, 124, 158, 154, 162, 150, 140, 110, 128,
130, 135, 114, 116, 124, 136, 142, 120, 120, 160, 158, 144, 130,
125, 175, 158, 185, 152, 159, 176, 156, 184, 138, 172, 168, 176,
164, 154, 124, 142, 144, 149, 128, 130, 138, 150, 156, 134, 134,
174, 174, 158, 144, 139, 180, 165, 172, 160, 157, 170, 153, 148,
140, 132, 169)), age = as.integer(c(39, 45, 47, 65, 46, 67, 42,
67, 56, 64, 56, 59, 34, 42, 48, 45, 17, 20, 19, 36, 50, 39, 21,
44, 53, 63, 29, 25, 69, 41, 60, 41, 47, 66, 47, 68, 43, 68, 57,
65, 57, 61, 36, 44, 50, 47, 19, 22, 21, 38, 52, 41, 18, 51, 55,
65, 33, 23, 70, 56, 62, 51, 48, 59, 40, 35, 33, 26, 61))), .Names =
c("gender",
"sbp", "age"), row.names = c("1", "2", "3", "4", "5", "6", "7",
"8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
"19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
"30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",
"41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51",
"52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62",
"63", "64", "65", "66", "67", "68", "69"), class = "data.frame")

### fit linear models for both gender simultaneously
lmod <- lm(sbp ~ gender * age, data = bp)
coef(lmod)

### define grid / linear functions
age <- seq(from = 17, to = 47, by = 2)
K <- cbind(0, 1, 0, age)
rownames(K) <- paste("age", age, sep = "")

### two-sided linear hypothesis
lh <- glht(lmod, linfct = K)
summary(lh, test = Ftest())

### simultaneous confidence bands
sci <- confint(lh, level = 0.99)
sci ### quantile: 2.969 +/- 0.003, page 61

### reproduce Figure 1
plot(age, coef(lh), ylim = c(-30, 2), type = "b", 
     ylab = expression(X^T * (b[F] - b[M])))
lines(age, sci$confint[,"upr"])
lines(age, sci$confint[,"lwr"])
abline(h = 0, lty = 2)

### add unadjusted point-wise confidence sets
uci <- confint(lh, level = 0.99, adjusted = FALSE)

lines(age, uci$confint[,"upr"], lty = 2)
lines(age, uci$confint[,"lwr"], lty = 2)

### one-sided: upper band only
lh <- glht(lmod, linfct = K, alternative = "less")

sci <- confint(lh, level = 0.99)
sci ### quantile: -2.711, page 64

### reproduce Figure 4
plot(age, coef(lh), ylim = c(-8, 2), type = "b",
     ylab = expression(X^T * (b[F] - b[M])))
lines(age, sci$confint[,"upr"])
abline(h = 0, lty = 2)

uci <- confint(lh, level = 0.99, adjusted = FALSE)
lines(age, uci$confint[,"upr"], lty = 2)

### fitness data (source: SAS User's Guide 1990)

fitness <- structure(list(age = as.integer(c(44, 40, 44, 42, 38, 47, 40,
43, 44, 38, 44, 45, 45, 47, 54, 49, 51, 51, 48, 49, 57, 54, 52,
50, 51, 54, 51, 57, 49, 48, 52)), weight = c(89.47, 75.07, 85.84,
68.15, 89.02, 77.45, 75.98, 81.19, 81.42, 81.87, 73.03, 87.66,
66.45, 79.15, 83.12, 81.42, 69.63, 77.91, 91.63, 73.37, 73.37,
79.38, 76.32, 70.87, 67.25, 91.63, 73.71, 59.08, 76.32, 61.24,
82.78), oxygen = c(44.609, 45.313, 54.297, 59.571, 49.874, 44.811,
45.681, 49.091, 39.442, 60.055, 50.541, 37.388, 44.754, 47.273,
51.855, 49.156, 40.836, 46.672, 46.774, 50.388, 39.407, 46.08,
45.441, 54.625, 45.118, 39.203, 45.79, 50.545, 48.673, 47.92,
47.467), runtime = c(11.37, 10.07, 8.65, 8.17, 9.22, 11.63, 11.95,
10.85, 13.08, 8.63, 10.13, 14.03, 11.12, 10.6, 10.33, 8.95, 10.95,
10, 10.25, 10.08, 12.63, 11.17, 9.63, 8.92, 11.08, 12.88, 10.47,
9.93, 9.4, 11.5, 10.5), restpulse = as.integer(c(62, 62, 45,
40, 55, 58, 70, 64, 63, 48, 45, 56, 51, 47, 50, 44, 57, 48, 48,
67, 58, 62, 48, 48, 48, 44, 59, 49, 56, 52, 53)), runpulse =
as.integer(c(178,
185, 156, 166, 178, 176, 176, 162, 174, 170, 168, 186, 176, 162,
166, 180, 168, 162, 162, 168, 174, 156, 164, 146, 172, 168, 186,
148, 186, 170, 170)), maxpulse = as.integer(c(182, 185, 168,
172, 180, 176, 180, 170, 176, 186, 168, 192, 176, 164, 170, 185,
172, 168, 164, 168, 176, 165, 166, 155, 172, 172, 188, 155, 188,
176, 172))), .Names = c("age", "weight", "oxygen", "runtime",
"restpulse", "runpulse", "maxpulse"), class = "data.frame", row.names =
c("1",
"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
"14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
"25", "26", "27", "28", "29", "30", "31"))

### define cutpoint
fitness$runcut <- fitness$runtime > 12
table(fitness$runcut) ### is this OK?

### fit linear model
lmod <- lm(oxygen ~ runcut + runcut:(age + runpulse) - 1, data = fitness)

### estimate b2 - b1
K <- rbind(c(1, -1, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 1, -1))
b2mb1 <- glht(lmod, linfct = K)
coef(b2mb1) ### is this OK?

### set up X^top (b2 - b1)
age <- seq(from = 38, to = 45, by = 1)
pulse <- seq(from = 160, to = 186, length = length(age))
X <- cbind(1, rep(age, rep(length(age), length(age))),
              rep(pulse, length(age)))
rownames(X) <- paste("age", X[,2], "pulse", round(X[,3], 1), sep = "")

### does not match with Figure 3
lh <- glht(b2mb1, linfct = X)
z <- matrix(coef(lh), nrow = length(age))
z
persp(x = age, y = pulse, z = z, tick = "detailed")
