#' ANCOM-BC Method
#'
#' This function applies the ANCOM-BC method for differential abundance analysis.
#'
#' @param data A phyloseq object containing count data and metadata.
#'
#' @return A data frame with test statistics, p-values, and adjusted p-values.
#'
#' @examples
#' \dontrun{
#'   results <- ANCOMBC_method(data)
#' }
#'
#' @importFrom ANCOMBC ancombc
#' @export
ANCOMBC_method <- function(data) {
  res <- ANCOMBC::ancombc(data, formula = "group", prv_cut = 0.05)
  
  res_pval <- res[["res"]][["p_val"]][["group"]]
  names(res_pval) <- res[["res"]][["p_val"]][["taxon"]]
  res_adjP <- p.adjust(res_pval, method = "BH")
  
  res_stat <- res[["res"]][["W"]][["group"]]
  names(res_stat) <- res[["res"]][["p_val"]][["taxon"]]
  
  out <- data.frame(stat = res_stat, pval = res_pval, adjP = res_adjP)
  rownames(out) <- res[["res"]][["taxon"]]
  return(out)
}
