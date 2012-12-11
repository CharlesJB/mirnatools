library(sfsmisc)

data <- read.table("mature_combinedMatches_relativeCount.txt")
x <- seq(1:nrow(data))
y <- data[,2]

drawPlot <- function(type) {
    if (type == "jpeg") {
        jpeg("Profiling_miRNA.jpeg")
    }
    if (type == "tiff") {
        tiff("Profiling_miRNA.tiff")
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

drawPlot("jpeg")
drawPlot("tiff")
