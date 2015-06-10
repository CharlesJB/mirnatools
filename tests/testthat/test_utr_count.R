library(mirnatools)

context("Test for function utr_count")

test_that("sequence_name argument is in a valid format", {
    msg <- "sequence_name should be a character string with more than one "
    msg <- paste0(msg, "character.")
    expect_error(utr_count(1), msg)
    expect_error(utr_count(""), msg)
})

test_that("sequence argument is in a valid format", {
    msg <- "sequence argument should be a character vector or a DNAString"
    msg <- paste(msg, "object.")
    expect_error(utr_count("abc", 1), msg)
    expect_error(utr_count("abc", factor("AAAAAA")), msg)
})

test_that("sequence argument contains valid values", {
    msg <- "sequence argument contains invalid values (not A, C, G, T"
    msg <- paste("or N).")
    expect_error(utr_count("abc", "X"), msg)
})

test_that("substring_length argument is in a valid format", {
    msg <- "substring_length should be a numeric vector of length 1 with no"
    msg <- paste(msg, "decimals.")
    expect_error(utr_count("abc", "AAAAAA", "A"))
    expect_error(utr_count("abc", "AAAAAA", 1.2))
})

test_that("sequence argument is long enough", {
    msg <- "sequence argument length should be greater of equal to"
    msg <- paste(msg, "substring_length.")
    expect_error(utr_count("abc", "", 6), msg)
    expect_error(utr_count("abc", "A", 6), msg)
    expect_error(utr_count("abc", "AAAAA", 6), msg)
})

test_that("utr argument is valid", {
    msg <- "utr argument must be \"5UTR\", \"3UTR\" or \"both\"."
    expect_error(utr_count("abc", "AAAAAA", 6, "mm10", 1), msg)
    expect_error(utr_count("abc", "AAAAAA", 6, "mm10", "a"), msg)
})

test_that("reverse_complement argument is valid", {
    msg <- "reverse_complement argument must be a logical."
    expect_error(utr_count("abc", "AAAAAA", 6, "mm10", "5UTR", 1), msg)
    expect_error(utr_count("abc", "AAAAAA", 6, "mm10", "5UTR", "a"), msg)
})
