#' Wilcoxon Rank-Sum Test
#'
#' This function applies the Wilcoxon rank-sum test to a dataset without normalization.
#'
#' @param data A phyloseq object containing count data and metadata.
#'
#' @return A data frame with test statistics, p-values, and adjusted p-values.
#'
#' @examples
#' \dontrun{
#'   results <- WilcoxTest(data)
#' }
#'
#' @importFrom phyloseq transform_sample_counts get_variable
#' @importFrom stats wilcox.test
#' @export
WilcoxTest <- function(data) {
  data <- transform_sample_counts(data, function(x) x / sum(x))
  counts <- as(otu_table(data), "matrix")
  groupBinary <- get_variable(data, "group") == get_variable(data, "group")[1L]
  
  stat <- apply(counts, 1, function(x) wilcox.test(x[groupBinary], x[!groupBinary], exact = FALSE)$statistic)
  pVals <- apply(counts, 1, function(x) wilcox.test(x[groupBinary], x[!groupBinary], exact = FALSE)$p.value)
  adjPVals <- p.adjust(pVals, method = "BH")
  
  out <- data.frame(stat = stat, pval = pVals, adjP = adjPVals)
  return(out)
}
