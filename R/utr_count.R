#' Count number of occurence of DNA substring in UTR regions
#'
#' This function will look for the number of occurence of every substrings of
#' length \code{substring_length} in the UTR regions of a given genome.
#'
#' @param sequence_name A \code{character} vector of length 1 corresponding to
#'                      name of the current sequence being analysed.
#' @param sequence A \code{character} vector or a \code{DNAString} object
#'                 corresponding to the miRNA sequence to analyze.
#' @param substring_length The length of each substring to use to find perfect
#'                         matches in the UTR regions.
#' @param utr Which UTR region to use? "5UTR", "3UTR" or "both".
#'            Default: "5UTR".
#' @param reverse_complement Use the reverse complement of the sequence?
#'                           Default: FALSE.
#'
#' @return A \code{data.frame} with the following columns:
#'     \itemize{
#'         \item name The sequence_name argument value.
#'         \item sequence The sequence argument value.
#'         \item substring The current_substring.
#'         \item substring The position of the substring in the sequence.
#'         \item count The number of perfect matches for the substring in the
#'                     UTR regions.
#'         \item percent The percentage of all the possible substrings of
#'                       \code{substring_length} that have a smaller number of
#'                       perfect matches in the UTR regions.
#'         \item mirbase_count The number of perfect matches for the all the
#'                             miRBase substrings at the same position of its
#'                             respective sequence in the specified genome.
#'         \item mirbase_percent The percentage of all miRBase substrings at
#'                               the same position in its respective sequence
#'                               that have a smaller number of perfect matches
#'                               in the UTR region of the specified genome.
#'         \item reverse_complement The value of the reverse complement
#'                                  argument.
#'     }
#'
#' @export
utr_count <- function(sequence_name, sequence, substring_length, genome, utr,
                      reverse_complement = FALSE) {
    if (!is.character(sequence_name) | nchar(sequence_name) == 0) {
        msg <- "sequence_name should be a character string with more than one "
        msg <- paste0(msg, "character.")
        stop(msg)
    }
    if (!is.character(sequence) & !is(sequence, "DNAString") &
        ! is(sequence, "RNAString")) {
        msg <- "sequence argument should be a character vector or a DNAString"
        msg <- paste(msg, "object.")
        stop(msg)
    }
    if (is.character(sequence)) {
        valids <- c("A", "a", "C", "c", "G", "g", "T", "t", "U", "u", "N", "n")
        split_sequence <- unlist(strsplit(sequence, ""))
        check <- vapply(split_sequence, function(x) x %in% valids, logical(1))
        if (!all(check)) {
            msg <- "sequence argument contains invalid values (not A, C, G, T"
            msg <- paste(msg, "or N).")
            stop(msg)
        }
    }
    if (!is.numeric(substring_length) |
        round(substring_length) != substring_length) {
        msg <- "substring_length should be a numeric vector of length 1 with no"
        msg <- paste(msg, "decimals.")
    }
    if (nchar(sequence) < substring_length) {
        msg <- "sequence argument length should be greater of equal to"
        msg <- paste(msg, "substring_length.")
        stop(msg)
    }
    if (!is.character(utr) | ! utr %in% c("5UTR", "3UTR", "both")) {
        msg <- "utr argument must be \"5UTR\", \"3UTR\" or \"both\"."
        stop(msg)
    }
    if (! is.logical(reverse_complement)) {
        msg <- "reverse_complement argument must be a logical"
        stop(msg)
    }
    # 1. Generate substrings
    if (reverse_complement) {
        sequence <- as(sequence, "DNAString")
        sequence <- Biostrings::reverseComplement(sequence)
    }
    sequence <- as.character(sequence)
    start <- seq(1, nchar(sequence) - substring_length + 1, 1)
    stop <- seq(6, nchar(sequence), 1)
    substrings <- substring(sequence, start, stop)

    # 2. Produce data.frame
    position <- seq_along(substrings)
    count <- vapply(substrings, get_occurence_count,
                    genome = genome, utr = utr,
                    FUN.VALUE = numeric(1))
    percent <- vapply(substrings, get_occurence_percent,
                    genome = genome, utr = utr,
                    FUN.VALUE = numeric(1))
    mirbase_count <- function(i) {
        get_mirbase_count(position = i, substring_length = substring_length,
                          genome = genome, utr = utr,
                          reverse_complement = reverse_complement)
    }
    mirbase_count <- vapply(position, mirbase_count, FUN.VALUE = numeric(1))
    mirbase_percent <- function(i) {
        get_mirbase_percent(position = i, sequence = substrings[i],
                            substring_length = substring_length,
                            reverse_complement = reverse_complement,
                            genome = genome, utr = utr)
    }
    mirbase_percent <- vapply(position, mirbase_percent,
                              FUN.VALUE = numeric(1))

    data.frame(name = sequence_name,
               sequence = sequence,
               substring = substrings,
               position = position,
               count = count,
               percent = percent,
               mirbase_count = mirbase_count,
               mirbase_percent = mirbase_percent,
               reverse_complement = reverse_complement,
               row.names = NULL)
}
