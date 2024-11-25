#' DESeq2 Method
#'
#' This function applies the DESeq2 method for differential abundance analysis on a given dataset.
#'
#' @param data A phyloseq object containing count data and metadata.
#'
#' @return A data frame with DESeq2 results, including test statistics, p-values, and adjusted p-values.
#'
#' @examples
#' \dontrun{
#'   results <- deseq2_method(data)
#' }
#'
#' @importFrom phyloseq phyloseq_to_deseq2
#' @importFrom DESeq2 estimateSizeFactors DESeq results
#' @export
deseq2_method <- function(data) {
  if (!taxa_are_rows(data)) {
    data <- t(data)
  }
  
  data@sam_data$group <- factor(data@sam_data$group)
  
  dds <- phyloseq::phyloseq_to_deseq2(data, design = as.formula("~ group"))
  dds <- DESeq2::estimateSizeFactors(dds, type = 'poscounts')
  
  test_deseq2 <- DESeq2::DESeq(dds)
  
  res <- DESeq2::results(test_deseq2, independentFiltering = FALSE)
  res_pval <- res$pvalue
  res_adjP <- p.adjust(as.numeric(res$pvalue), method = "BH")
  res_stat <- res$stat
  
  out <- data.frame(stat = res_stat, pval = res_pval, adjP = res_adjP)
  return(out)
}
