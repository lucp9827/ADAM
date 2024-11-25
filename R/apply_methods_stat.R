#' apply_methods_stat: Statistical Method Application Function
#'
#' This function applies multiple statistical methods to a given dataset and returns
#' the results in a structured list. It handles errors gracefully by using `tryCatch`.
#'
#' @param data A dataset on which to apply the statistical methods.
#'
#' @return A list containing results from various statistical methods.
#' @examples
#' \dontrun{
#'   results <- apply_methods_stat(data)
#' }
#' @export
apply_methods_stat <- function(data) {

  # ALDEx2 (no test statistic)
  res_aldex2 <- tryCatch({aldex_method(data)}, error = function(e) { cat("Error") })

  res_aldex2_t <- res_aldex2$ALDEx2_t
  res_aldex2_w <- res_aldex2$ALDEx2_w
  res_aldex2 <- res_aldex2$total

  # CornCob
  res_corncob <- tryCatch({corncob_method(data)}, error = function(e) { cat("Error") })

  # ANCOMBC
  res_ancombc <- tryCatch({ANCOMBC_method(data)}, error = function(e) { cat("Error") })

  # LinDA
  res_LinDA <- tryCatch({LinDA_test(data)}, error = function(e) { cat("Error") })

  # Wilcoxon test (TSS normalization)
  res_wilcox <- tryCatch({WilcoxTest(data)}, error = function(e) { cat("Error") })

  # Wilcoxon test (CLR normalization)
  res_wilcox_clr <- tryCatch({WilcoxTest_clr(data)}, error = function(e) { cat("Error") })

  # R-sign methods
  res_Rsignmethods <- tryCatch({RsignRI_test(data, model_R = "All")}, error = function(e) { cat("Error") })

  if (is.null(res_Rsignmethods)) {
    res_Rsignmethods <- tryCatch({RsignRI_test_sparse(data, model_R = "All")}, error = function(e) { cat("Error") })
  }

  res_marg <- res_Rsignmethods$Marginal
  res_cond <- res_Rsignmethods$Conditional
  res_ri <- res_Rsignmethods$RI

  # DESeq2
  res_deseq <- tryCatch({deseq2_method(data)}, error = function(e) { cat("Error") })

  # Output
  res <- list(
    res_aldex2_t = res_aldex2_t,
    res_aldex2_w = res_aldex2_w,
    res_corncob = res_corncob,
    res_ancombc = res_ancombc,
    res_wilcox = res_wilcox,
    res_marginal = res_marg,
    res_conditional = res_cond,
    res_RI = res_ri,
    res_linda = res_LinDA,
    res_deseq = res_deseq,
    res_wilcox_clr = res_wilcox_clr
  )

  return(res)
}
