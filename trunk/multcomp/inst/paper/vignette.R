
tex <- readLines("gpsi.tex")
rnw1 <- readLines("illustrations.Rnw")
rnw2 <- readLines("trees.Rnw")

tex <- c("%%\\VignetteIndexEntry{Simultaneous Inference in General Parametric Models}", 
         "%%\VignetteDepends{multcomp,mboost,survival,robustbase,lme4}", 
         "%%\\usepackage{Sweave}", tex)

w1 <- grep("input*\\{illustrations\\}", tex)         
tex <- c(tex[1:(w1-1)], rnw1, tex[(w1 + 1):length(tex)])

w1 <- grep("input*\\{trees\\}", tex)         
tex <- c(tex[1:(w1-1)], rnw2, tex[(w1 + 1):length(tex)])

writeLines(tex, "generalsiminf.Rnw")
