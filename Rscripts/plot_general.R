library(ggplot2)
argv <- commandArgs(trailingOnly = T)

file_trimmed <- argv[1]
output <- argv[2]

a <- read.table(file_trimmed, header=TRUE)
tiff(output)
qplot(type, data=a, weight=count, geom="bar")
dev.off()
