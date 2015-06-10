generate_all_sequences <- function(len) {
    dna <- c("A", "C", "G", "T")
    l <- rep(list(dna), len)
    apply(expand.grid(l), 1, paste, collapse = "")
}

get_occurence_count <- function(sequence, genome, utr) {
    occurence_count[[genome]][[utr]][[sequence]]
}

get_occurence_percent <- function(sequence, genome, utr) {
   oc <- get_occurence_count(sequence, genome, utr)
   all_oc <- occurence_count[[genome]][[utr]]
   sum(oc > all_oc) / length(all_oc) * 100
}

get_mirbase_count <- function(position, substring_length, genome, utr,
                              reverse_complement = FALSE) {
    current_mirbase <- subset_mirbase(genome, reverse_complement)
    start <- position
    end <- position + substring_length - 1
    substrings <- substring(current_mirbase, start, end)
    counts <- vapply(substrings, get_occurence_count,
           genome = genome, utr = utr,
           FUN.VALUE = numeric(1))
    mean(counts)
}

get_mirbase_percent <- function(position, sequence, substring_length, genome,
                                utr, reverse_complement = FALSE) {
    current_mirbase <- subset_mirbase(genome, reverse_complement)
    start <- position
    end <- position + substring_length - 1
    substrings <- substring(current_mirbase, start, end)
    counts <- vapply(substrings, get_occurence_count,
           genome = genome, utr = utr,
           FUN.VALUE = numeric(1))
    sequence_count <- get_occurence_count(sequence, genome, utr)
    sum(sequence_count > counts) / length(counts) * 100
}

subset_mirbase <- function(genome, reverse_complement) {
    name <- ""
    if (genome == "mm10") {
        name <- "mmu"
    }
    i <- grepl(name, names(mirbase))
    if (reverse_complement) {
        Biostrings::reverseComplement(as(mirbase[i], "DNAStringSet"))
    } else {
        as(mirbase[i], "DNAStringSet")
    }
}
