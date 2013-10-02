library(ggplot2)
argv <- commandArgs(trailingOnly = T)

file_trimmed <- argv[1]
output <- argv[2]

a <- read.table(file_trimmed, header=TRUE)
tiff(output, height=nrow(a)*60)
qplot(type, data=a, weight=percent, geom="bar") + facet_wrap(~ name, ncol=1)
dev.off()
