#' perm_boot: Permutation and Bootstrap Sampling Function
#'
#' This function performs permutation and bootstrap sampling on training data
#' and applies statistical methods to the modified data. The results are saved 
#' in a list.
#'
#' @param data_training Training data for the procedure.
#' @param B Number of permutation/bootstrap datasets to generate.
#' @param nSamples Number of samples of the original dataset.
#' @param path Directory path to save the results (if saving functionality is enabled).
#' @param dataset Name of the dataset for identification purposes.
#' 
#' @return A list of intermediate results obtained after applying statistical methods.
#' @examples
#' \dontrun{
#'   res <- perm_boot(data_training, B = 100, nSamples = 50, path = "results/", dataset = "example")
#' }
#' 
#' @importFrom signtrans Trim
#' @importFrom phyloseq tax_table otu_table sample_names
#' @export
perm_boot <- function(data_training, B, nSamples, path, dataset) {
  
  # To save results for each (B-times) permutation/bootstrap dataset
  res_intermediate <- list()
  
  # Trim training data
  data_trim <- signtrans::Trim(data_training, minPrev = 0.05)
  data_modified <- data_trim
  
  # Identify index of I0/I1 OTUs
  I0_id <- which(tax_table(data_trim)[, 3] == "I_0")
  I1_id <- which(tax_table(data_trim)[, 3] == "I_1")
  
  set.seed(9827)
  for (i in seq_len(B)) {
    
    counts_trim <- data.frame(otu_table(data_trim))
    
    # Permutation method
    for (j in I0_id) {
      counts_trim[j, ] <- sample(counts_trim[j, ], replace = FALSE)
    }
    
    # Bootstrap method
    for (otu_index in I1_id) {
      otu_values_1 <- counts_trim[otu_index, seq_len(nSamples / 2)]
      otu_values_2 <- counts_trim[otu_index, seq((nSamples / 2) + 1, nSamples)]
      
      # Perform bootstrapping within the OTU
      bootstrap_ind_1 <- sample(as.numeric(otu_values_1), replace = TRUE)
      bootstrap_ind_2 <- sample(as.numeric(otu_values_2), replace = TRUE)
      vector_values <- c(bootstrap_ind_1, bootstrap_ind_2)
      
      cnt <- 1
      while (!mean(vector_values >= 1) >= 0.05 && cnt < 100) {
        bootstrap_ind_1 <- sample(as.numeric(otu_values_1), replace = TRUE)
        bootstrap_ind_2 <- sample(as.numeric(otu_values_2), replace = TRUE)
        cnt <- cnt + 1
      }
      counts_trim[otu_index, ] <- c(bootstrap_ind_1, bootstrap_ind_2)
    }
    
    colnames(counts_trim) <- sample_names(data_trim)
    data_modified@otu_table <- otu_table(counts_trim, taxa_are_rows = TRUE)
    save(data_modified, file = paste0(path, "/Variations/Modified_Test_data", i, ".RData"))
    
    res_intermediate[[i]] <- apply_methods_stat(data_modified)
  }
  
  save(res_intermediate, file = paste0(path, "/res_perm_boot", dataset, ".RData"))
  
  return(res_intermediate)
}
