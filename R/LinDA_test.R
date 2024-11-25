#' LinDA Test
#'
#' This function applies the LinDA method from the MicrobiomeStat package for differential abundance analysis.
#'
#' @param data A phyloseq object containing count data and metadata.
#'
#' @return A data frame with LinDA test statistics, p-values, and adjusted p-values.
#'
#' @examples
#' \dontrun{
#'   results <- LinDA_test(data)
#' }
#'
#' @importFrom MicrobiomeStat linda
#' @export
LinDA_test <- function(data) {
  otu_matrix <- as.data.frame(otu_table(data))
  group <- data.frame(group = sample_data(data)$group)
  
  linda.obj <- MicrobiomeStat::linda(otu_matrix, group, formula = '~group', alpha = 0.05)
  pVals <- linda.obj[["output"]][["group"]][["pvalue"]]
  adjPVals <- p.adjust(pVals, method = "BH")
  
  stat <- linda.obj[["output"]][["group"]][["stat"]]
  out <- data.frame(pval = pVals, adjP = adjPVals, stat = stat)
  return(out)
}
