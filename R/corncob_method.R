#' Corncob Method
#'
#' This function applies the Corncob method for differential abundance analysis using likelihood ratio tests.
#'
#' @param data A phyloseq object containing count data and metadata.
#'
#' @return A data frame with test statistics, p-values, and adjusted p-values.
#'
#' @examples
#' \dontrun{
#'   results <- corncob_method(data)
#' }
#'
#' @importFrom corncob differentialTest
#' @export
corncob_method <- function(data) {
  res <- corncob::differentialTest(
    formula = ~ group,
    phi.formula = ~ group,
    formula_null = ~ 1,
    phi.formula_null = ~ group,
    test = "LRT", boot = FALSE,
    data = data,
    fdr_cutoff = 0.05
  )
  
  res_stat <- sapply(1:ntaxa(data), function(i) {
    if (length(res[["all_models"]][[i]]) != 1) {
      res[["all_models"]][[i]][["coefficients"]][2, 3]
    } else {
      NA
    }
  })
  
  res_pval <- res[["p"]]
  res_adjP <- p.adjust(as.numeric(res_pval), method = "BH")
  names(res_pval) <- names(res_adjP) <- names(res_stat) <- rownames(data)
  
  out <- data.frame(stat = res_stat, pval = res_pval, adjP = res_adjP)
  return(out)
}
