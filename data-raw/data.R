produce_mirbase <- function(version, dest_dir = ".") {
    # Download database
    filename <- "mature.fa.gz"
    url <- "ftp://mirbase.org/pub/mirbase"
    url <- paste(url, version, filename, sep = "/")
    filename <- paste(dest_dir, filename, sep = "/")
    if (! file.exists(filename)) {
        download.file(url, filename)
    }

    # Import in R
    Biostrings::readRNAStringSet(filename)
}

precalculate_occurences <- function(seqlength, txbd, bsgenome) {
    results <- list()
    # 1. Generate all possible subsequences
    all_sequences <- generate_all_sequences(seqlength)

    # 2. Extract the UTR regions
    utr3 <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
    utr5 <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)

    # 3. Extract the UTR sequences
    utr3 <- unlist(Biostrings::getSeq(bsgenome, utr3))
    utr3_revC <- Biostrings::reverseComplement(utr3)
    utr5 <- unlist(Biostrings::getSeq(bsgenome, utr5))
    utr5_revC <- Biostrings::reverseComplement(utr5)

    # 4. Calculate the number of occurences
    count_occurences <- function(string, sequences) {
        sum(Biostrings::vcountPattern(string, sequences))
    }
    results[["3UTR"]] <- lapply(all_sequences, count_occurences, utr3)
    results[["3UTR_RevC"]] <- lapply(all_sequences, count_occurences, utr3_revC)
    results[["5UTR"]] <- lapply(all_sequences, count_occurences, utr5)
    results[["5UTR_RevC"]] <- lapply(all_sequences, count_occurences, utr5_revC)
    results <- lapply(results, unlist)
    results <- lapply(results, function(x) { names(x) <- all_sequences; x })
    results
}


