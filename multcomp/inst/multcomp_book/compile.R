
library("utils")
library("tools")
Rnw <- list.files(pattern = "Rnw")
sapply(Rnw, function(f) Sweave(f))
texi2dvi("Multiple_Comparisons_in_R.tex", pdf = TRUE, clean = FALSE)
