#' ALDEx2 Method
#'
#' This function applies the ALDEx2 method for differential abundance analysis using compositional data.
#'
#' @param data A phyloseq object containing count data and metadata.
#'
#' @return A list containing t-test and Wilcoxon test results, as well as the full ALDEx2 output.
#'
#' @examples
#' \dontrun{
#'   results <- aldex_method(data)
#' }
#'
#' @importFrom ALDEx2 aldex
#' @export
aldex_method <- function(data) {
  counts <- data.frame(otu_table(data))
  group <- sample_data(data)$group
  
  res <- ALDEx2::aldex(counts, group, test = 't')
  
  res_stat <- res$effect
  rownames(res) -> names(res_stat)
  
  res_pvalue_t <- res$we.ep
  res_pvalue_w <- res$wi.ep
  res_adjP_t <- p.adjust(res_pvalue_t, method = "BH")
  res_adjP_w <- p.adjust(res_pvalue_w, method = "BH")
  
  out_t <- data.frame(stat = res_stat, pval = res_pvalue_t, adjP = res_adjP_t)
  out_w <- data.frame(stat = res_stat, pval = res_pvalue_w, adjP = res_adjP_w)
  
  out <- list(ALDEx2_t = out_t, ALDEx2_w = out_w, total = res)
  return(out)
}
