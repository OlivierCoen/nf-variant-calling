#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
library(optparse)
library(dplyr)


options(error = traceback)


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

get_args <- function() {
    option_list <- list(
        make_option("--RO", dest="RO_file"),
        make_option("--AO", dest="AO_file"),
        make_option("--design", dest="design_file"),
        make_option("--out", dest="output_file")
    )

    args <- parse_args(OptionParser(
        option_list = option_list
        ))
    return(args)
}


get_sample_lists <- function(design) {
    grouped_design <- design %>%
        unique() %>%
      arrange(population) %>%
      group_by(phenotype) %>%
      summarise(samples = list(sample), .groups = "drop")

  sample_lists <- grouped_design$samples
  return(sample_lists)
}


get_cmh_pval <- function(R0, R1, A0, A1, correct = TRUE) {
    # matrices [n_snp × n_pop]
    # R for reference allele, A for alternative allele (but this can be swapped, as long as the same allele is the reference)
    # 0 for phenotype 0, 1 for phenotype 1 (whichever these are)
    # correct = applies Yates contiguity correction

    n <- R0 + R1 + A0 + A1

    row_sum1 <- R0 + R1
    row_sum2 <- A0 + A1
    col_sum1 <- R0 + A0
    col_sum2 <- R1 + A1

    expected <- (row_sum1 * col_sum1) / n
    V <- (row_sum1 * row_sum2 * col_sum1 * col_sum2) / (n^2 * (n - 1))

    delta <- abs(rowSums(R0 - expected, na.rm = TRUE))
    yates <- 0
    if(correct) {
        yates <- pmin(delta, 0.5)
    }

    numerator <- (delta-yates)^2
    denominator <- rowSums(V, na.rm = TRUE)

    stat <- numerator / denominator
    p <- pchisq(stat, 1, lower.tail = FALSE)
    p[!is.finite(p)] <- NA_real_

    return(p)
}


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


main <- function() {

    args <- get_args()

    RO <- read.csv(args$RO_file)
    AO <- read.csv(args$AO_file)
    design <- read.csv(args$design_file)

    sample_lists <- get_sample_lists(design)

    if (length(sample_lists) != 2) {
        message("Exactly two phenotypes needed here!")
        quit(save = "no", status = 1)
    }

    samples_pheno_1 <- sample_lists[[1]]
    samples_pheno_2 <- sample_lists[[2]]

    message("Samples for phenotype 1: ", paste(samples_pheno_1, collapse = ", "))
    message("Samples for phenotype 2: ", paste(samples_pheno_2, collapse = ", "))

    p_values <- get_cmh_pval(
      R0 = RO[,cbind(samples_pheno_1)],
      R1 = RO[,cbind(samples_pheno_2)],
      A0 = AO[,cbind(samples_pheno_1)],
      A1 = AO[,cbind(samples_pheno_2)],
      correct = TRUE
    )

    p_values <- p.adjust(p_values, method = "fdr")

    write.table(p_values, file = args$output_file, row.names = FALSE, col.names = FALSE)
}


#####################################################
# ENTRYPOINT
#####################################################

main()
