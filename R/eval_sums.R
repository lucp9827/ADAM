#' Evaluate Performance Metrics for Statistical Tests
#'
#' This function calculates Type I error rate (FPP), sensitivity (TPP), 
#' and false discovery proportion (FDP) based on the provided p-values and 
#' group indices for I0 (true negatives) and I1 (true positives). 
#'
#' @param pval A numeric vector of p-values.
#' @param I0_id Indices of true negatives (I0).
#' @param I1_id Indices of true positives (I1).
#' @param alpha Significance level for determining significant p-values. Default is 0.05.
#'
#' @return A numeric vector with three elements:
#'   \describe{
#'     \item{\code{FPP}}{False positive proportion (Type I error rate).}
#'     \item{\code{TPP}}{True positive proportion (sensitivity).}
#'     \item{\code{FDP}}{False discovery proportion.}
#'   }
#'
#' @examples
#' \dontrun{
#'   pvals <- c(0.01, 0.02, 0.5, 0.04, 0.3)
#'   I0 <- c(3, 5)
#'   I1 <- c(1, 2, 4)
#'   metrics <- eval_sums(pvals, I0, I1, alpha = 0.05)
#' }
#'
#' @export
eval_sums <- function(pval, I0_id, I1_id, alpha = 0.05) {
  TN <- sum(pval[I0_id] >= alpha, na.rm = TRUE)
  FP <- sum(pval[I0_id] < alpha, na.rm = TRUE)
  TP <- sum(pval[I1_id] < alpha, na.rm = TRUE)
  FN <- sum(pval[I1_id] >= alpha, na.rm = TRUE)
  
  FPP <- FP / (FP + TN)
  TPP <- TP / (TP + FN)
  
  if ((FP + TP) == 0) {
    FDP <- 0
  } else {
    FDP <- FP / (FP + TP)
  }
  
  return(c(FPP, TPP, FDP))
}
