library(methylKit)

args <- commandArgs(TRUE)

processBismarkAln(
    location = args[1],
    sample.id = args[2],
    assembly = "hg38",
    save.folder = getwd(),
    save.context = c("CpG"),
    read.context = "CpG",
    nolap = FALSE,
    mincov = 5,
    minqual = 20,
    phred64 = FALSE,
    treatment = NULL,
    save.db = FALSE
)
