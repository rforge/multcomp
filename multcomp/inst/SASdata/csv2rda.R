
csv <- list.files(pattern = ".csv")

sapply(csv, function(x) {
    print(x)
    tmp <- read.table(x, header = TRUE, sep = ",")
    print(colnames(tmp))
    ### tmp <- tmp[, !(colnames(tmp) %in% c("i", "j"))]
    names(tmp) <- tolower(names(tmp))
    name <- tolower(gsub(".csv", "", x))
    eval(parse(text = paste(name, " <- tmp")))
    prompt(name = name)
    save(list = name, file = tolower(gsub("csv", "rda", x)))
})

