drawPlot <- function(type, inputFilename, outputBasename) {
    if (type == "jpeg") {
	filename = paste(outputBasename, ".jpeg", sep="")
        jpeg(filename)
    }
    if (type == "tiff") {
	filename = paste(outputBasename, ".tiff", sep="")
        tiff(filename)
    }
    data <- read.table(inputFilename)
    names(data)[1] <- "miRNA name"
    names(data)[2] <- "Percentage of total miRNA"
    pie(data[,2], labels = data[,1], main = "Distribution of the most abundant miRNAs", clockwise=TRUE)
    dev.off()
}

argv <- commandArgs(trailingOnly = T)

if (length(argv) == 3) {
	type <- argv[1]
	inputFilename <- argv[2]
	outputBasename <- argv[3]

	drawPlot(type, inputFilename, outputBasename)

} else {

	cat("\nUsage:\n")
	cat("percent_pie_chart.R <type> <inputFilename> <outputBasename>\n")
	cat("    type: type of file (jpeg, tiff, pdf, etc...)\n")
	cat("    inputFilename: File containing the name and proportion of miRNA to draw.\n")
	cat("    outputBasename: basename for the output file.\n")

}
