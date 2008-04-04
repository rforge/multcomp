
### read in data
rda <- list.files(pattern = "\\.rda")
sapply(rda, function(x) load(file = x, env = .GlobalEnv))

### modify
ffoo <- function(x, n) factor(x, labels = paste(n, sort(unique(x)), sep = ""))
tox$g <- ffoo(tox$g, "g")

coupon$discount <- ffoo(coupon$discount, "d")

alz$therapy <- ffoo(alz$therapy, "t")

litter$dose <- ffoo(litter$dose, "d")

cholesterol$trt <- as.factor(rev(as.integer(cholesterol$trt)))
levels(cholesterol$trt) <- LETTERS[1:5]

waste$temp <- ffoo(waste$temp, "t")
waste$envir <- ffoo(waste$envir, "e")

drug$drug <- ffoo(drug$drug, "d")

detergent$block <- ffoo(detergent$block, "b")
detergent$detergent <- ffoo(detergent$detergent, "d")

pigs$pen <- ffoo(pigs$pen, "p")
pigs$feed <- ffoo(pigs$feed, "f")

wine$purchase <- ffoo(wine$purchase, "p")
wine$customertype <- ffoo(wine$customertype, "c")
wine$light <- ffoo(wine$light, "l")
wine$music <- ffoo(wine$music, "m")

rm("rda")
rm("ffoo")
x <- ls()
save(list = x, file = "Westfall_et_al.rda")
