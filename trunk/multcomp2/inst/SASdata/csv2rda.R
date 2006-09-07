
csv <- list.files(pattern = ".csv")

sapply(csv, function(x) {
    print(x)
    tmp <- read.table(x, header = TRUE, sep = ",")
    name <- tolower(gsub(".csv", "", x))
    eval(parse(text = paste(name, " <- tmp")))
    prompt(name = name)
    save(list = name, file = tolower(gsub("csv", "rda", x)))
})

