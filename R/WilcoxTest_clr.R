#' Wilcoxon Rank-Sum Test with CLR Normalization
#'
#' This function applies the Wilcoxon rank-sum test to a dataset with CLR normalization.
#'
#' @param data A phyloseq object containing count data and metadata.
#'
#' @return A data frame with test statistics, p-values, and adjusted p-values.
#'
#' @examples
#' \dontrun{
#'   results <- WilcoxTest_clr(data)
#' }
#'
#' @importFrom microbiome transform
#' @importFrom phyloseq get_variable
#' @importFrom stats wilcox.test
#' @export
WilcoxTest_clr <- function(data) {
  data_norm <- microbiome::transform(data, transform = "clr")
  counts <- as(otu_table(data_norm), "matrix")
  groupBinary <- get_variable(data_norm, "group") == get_variable(data_norm, "group")[1L]
  
  stat <- apply(counts, 1, function(x) wilcox.test(x[groupBinary], x[!groupBinary], exact = FALSE)$statistic)
  pVals <- apply(counts, 1, function(x) wilcox.test(x[groupBinary], x[!groupBinary], exact = FALSE)$p.value)
  adjPVals <- p.adjust(pVals, method = "BH")
  
  out <- data.frame(stat = stat, pval = pVals, adjP = adjPVals)
  return(out)
}
