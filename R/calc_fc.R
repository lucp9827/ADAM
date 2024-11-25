#' Calculate Log Fold Changes between Groups
#' This function calculates log fold changes (LFC) for OTUs between two groups within a phyloseq object.
#' It supports two normalization methods: TSS (Total Sum Scaling) and CLR (Centered Log-Ratio).
#'
#' @param data_1 A phyloseq object representing the first group for analysis.
#' @param data_2 A phyloseq object representing the second group for analysis.
#' @param norm_method Character string specifying the normalization method. Choose either "TSS" or "CLR".
#' @return A numeric vector with log fold changes for each OTU between the two groups.
#' @examples
#' # Assuming you have two phyloseq objects for two groups:
#' # log_fc <- calc_fc(group_1, group_2, "TSS")
#' @export
calc_fc = function(data_1, data_2, norm_method){
  
  # Calculate row means for each OTU in both groups
  mean_1 = rowMeans(otu_table(data_1))
  mean_2 = rowMeans(otu_table(data_2))
  
  # Calculate log fold changes based on normalization method
  if (norm_method == "TSS") {
    logfoldchanges = abs(log2(mean_1 / mean_2))
  } else if (norm_method == "CLR") {
    logfoldchanges = abs((mean_1 - mean_2) / log(2))
  } else {
    stop("Invalid normalization method. Choose either 'TSS' or 'CLR'.")
  }
  
  return(logfoldchanges)
}