cmh_pval <- function(R0, R1, A0, A1, correct = TRUE) {
  # matrices [n_snp × n_pop]
  # R for reference allele, A for alternative allele (but this can be swapped, as long as the same allele is the reference)
  # 0 for phenotype 0, 1 for phenotype 1 (whichever these are)
  # correct = applies Yates contiguity correction

  n <- R0 + R1 + A0 + A1

  row_sum1 <- R0 + R1
  row_sum2 <- A0 + A1
  col_sum1 <- R0 + A0
  col_sum2 <- R1 + A1

  Expected <- (row_sum1 * col_sum1) / n
  V <- (row_sum1 * row_sum2 * col_sum1 * col_sum2) / (n^2 * (n - 1))

  delta <- abs(rowSums(R0 - Expected, na.rm = TRUE))
  yates = 0
  if(correct) {
    yates = pmin(delta, 0.5)
  }

  numerator <- (delta-yates)^2
  denominator <- rowSums(V, na.rm = TRUE)

  stat <- numerator / denominator
  p <- pchisq(stat, 1, lower.tail = FALSE)
  p[!is.finite(p)] <- NA_real_
  return(p)
}
