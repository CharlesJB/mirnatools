library(ggplot2)
argv <- commandArgs(trailingOnly = T)

file_trimmed <- argv[1]
output <- argv[2]

a <- read.table(file_trimmed, header=FALSE)
colnames(a) <- c("position", "count")

tiff(output)
qplot(position, count, data = a, geom = "line")
dev.off()

