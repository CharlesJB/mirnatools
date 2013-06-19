library(ggplot2)

argv <- commandArgs(trailingOnly = T)
file_summary <- argv[1]
file_output <- argv[2]

data <- read.table(file_summary, header=TRUE)
data <- data[with(data, order(-Count)),]

tiff(file_output)
g <- ggplot(data, aes(x=RNA_specie, y=Count))
g + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust=1, size = 12))
dev.off()
