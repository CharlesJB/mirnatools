library(sfsmisc)

drawPlot <- function(type, inputFilename, outputBasename) {
    data <- read.table(inputFilename)
    x <- seq(1:nrow(data))
    y <- data[,2]

    if (type == "jpeg") {
filename = paste(outputBasename, ".jpeg", sep="")
        jpeg(filename)
    }
    if (type == "tiff") {
filename = paste(outputBasename, ".tiff", sep="")
        tiff(filename)
    }
    # Draw empty plot
    plot(x, y, type = "n", xaxt = "n", yaxt="n", xlab="", ylab="", log = "y", col="blue")

    # Add axis infos
    eaxis(1, padj=-0.5, cex.axis=0.8)
    eaxis(2, padj=-0.5, cex.axis=0.8)

    # Add text
    mtext(side=1, text="MicroRNAs detected in human platelet", line=2.0)
    mtext(side=2, text="MicroRNA expression level (per million reads)", line=2.0)
    mtext(side=3, text="Profiling of microRNAs", line=1.2, cex=1.3)

    # Add grid
    grid(lwd=1.5)

    # Add lines
    lines(x, y, col="red", lwd=2.2)
    lines(x, rep(0.4, length(x)), col="black")

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
cat("profiling.R <type> <inputFilename> <outputBasename>\n")
cat(" type: type of file (jpeg, tiff, pdf, etc...)\n")
cat(" inputFilename: File containing the relative count of mature miRNA\n")
cat(" outputBasename: basename for the output file.\n")

}
