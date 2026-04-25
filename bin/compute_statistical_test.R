#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("arrow"))
library(argparse)
library(dplyr)
library(arrow)


options(error = traceback)

NB_MONTE_CARLO_SIMULATIONS <- 10000


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

ALLOWED_METHODS <- c("cmh", "fet", "chisq")

get_args <- function() {
    parser <- ArgumentParser(description='Compute stat test')
    parser$add_argument("--method", choices = ALLOWED_METHODS, help="Method to use")
    parser$add_argument("--RO", dest="RO_file", help="Ref counts")
    parser$add_argument("--AO", dest="AO_file", help="Alt counts")
    parser$add_argument("--design", dest="design_file", help="Design")
    parser$add_argument("--out", dest="output_file", help="Output file")
    
    args <- parser$parse_args()
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

cast_to_numeric <- function(df){
  return(df %>% mutate(across(everything(), as.numeric)))
}

compute_fet_from_contingency <- function(mat) {
  # Compute the fischer exact test from a contingency table
  # if they are more than 2 columns (phenotypes), use Monte Carlo simulation
  if ( ncol(mat) > 2 ) {
    res <- fisher.test(
      mat, 
      alternative = "two.sided",
      simulate.p.value = TRUE, 
      B = NB_MONTE_CARLO_SIMULATIONS
    )
  } else {
      res <- fisher.test(mat, alternative = "two.sided")
  }
  return(res$p.value)
}

compute_chisq_from_contingency <- function(mat) {
  # Compute the chi-square test from a contingency table
  res <- chisq.test(mat)
  return(res$p.value)
}


compute_test <- function(mathod, R0, A0, sample_lists) {
    # matrices [n_snp × n_pop]
    # Apply stat test row-wise (per SNP), aggregating across populations

    n_snp <- nrow(R0)
    
    # grouping reference and alternative allele by phenotype
    R <- list()
    A <- list()
    for ( pheno in seq_along(sample_lists) ) {
      R[[pheno]] <- R0[, sample_lists[[pheno]]]
      A[[pheno]] <- A0[, sample_lists[[pheno]]]
    }

    p_values <- vapply(seq_len(n_snp), function(i) {
      
        # Aggregate counts across populations for this SNP
        mat_elements <- list()
        for ( pheno in seq_along(R) ) {
          ref_pheno <- sum(R[[pheno]][i, ], na.rm = TRUE)
          alt_pheno <- sum(A[[pheno]][i, ], na.rm = TRUE)
          mat_elements[[pheno]] <- c(ref_pheno, alt_pheno)
        }
        mat <- matrix(unlist(mat_elements), nrow = 2)
        
        # Example for a 2x2 contingency table:
        #               pheno 1            pheno 2
        # ref allele    r[[pheno1]]        r[[pheno2]]
        # alt allele    a[[pheno1]]        a[[pheno2]]

        if (any(is.na(mat)) || sum(mat) == 0) return(NA_real_)
        
        # add pseudocount to avoid zeroes (TODO: see how to improve that)
        mat  = mat + 1
        
        if (mathod == "fet") {
          p_value <- compute_fet_from_contingency(mat)
        } else if (mathod == "chisq") {
          p_value <- compute_chisq_from_contingency(mat)
        }
        
        p_value
        
    }, numeric(1))

    return(p_values)
}


compute_cochran_mantel_haenszel_test <- function(R0, R1, A0, A1, correct = TRUE) {
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
  
    RO <- arrow::read_parquet(args$RO_file, as_tibble=FALSE)
    AO <- arrow::read_parquet(args$AO_file, as_tibble=FALSE)
    design <- read.csv(args$design_file)

    sample_lists <- get_sample_lists(design)
    
    for (i in seq_along(sample_lists)) {
      message(paste("Samples for phenotype", i, ":", sample_lists[[i]]))
    }

    RO <- cast_to_numeric(RO)
    AO <- cast_to_numeric(AO)
    
    if ( args$method == "cmh" ) {
      
        if (length(sample_lists) != 2) {
            message("Exactly two phenotypes needed here!")
            quit(save = "no", status = 1)
        }
        samples_pheno_1 <- sample_lists[[1]]
        samples_pheno_2 <- sample_lists[[2]]
        
        p_values <- compute_cochran_mantel_haenszel_test(
          R0 = RO[,samples_pheno_1],
          R1 = RO[,samples_pheno_2],
          A0 = AO[,samples_pheno_1],
          A1 = AO[,samples_pheno_2],
          correct = TRUE
        )
        
    } else if ( args$method %in% c("fet", "chisq") ) {
      
      p_values <- compute_test(args$method, RO, AO, sample_lists)
      
    } else { 
      error(paste("Method not recognised:", args$method))
    }

    p_values <- p.adjust(p_values, method = "fdr")
    
    write.table(p_values, file = args$output_file, row.names = FALSE, col.names = FALSE)
}


#####################################################
# ENTRYPOINT
#####################################################

main()
