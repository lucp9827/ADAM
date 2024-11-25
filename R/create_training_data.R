#' Create Test Dataset for Microbiome Analysis
#'
#' This function creates test data based on sampling OTUs from different strata.
#'
#' @param data A phyloseq object containing the microbiome data.
#' @param TPR True Positive Rate for the I1 group (numeric).
#' @param TPR_B True Positive Rate for the B group (numeric).
#' @param perc_strata Sampling fraction (numeric).
#' @param treshold_I1 Threshold for I1 group (numeric).
#' @param treshold_I01 Threshold for I0 group (numeric).
#' @param treshold_I02 Threshold for between groups (numeric).
#' @param OTU_replace Boolean indicating whether to replace sampled OTUs.
#' @param combinations Data frame of combinations (data frame).
#' @param combinations_I0 Data frame of I0 combinations (data frame).
#' @param combinations_I1 Data frame of I1 combinations (data frame).
#' @param combinations_B Data frame of B combinations (data frame).
#' @param path Directory path to save training data (character string).
#' @param k Index or identifier for the training data (numeric).
#' @return A phyloseq object containing the training data.
#' @import dplyr
#' @import phyloseq
#' @import microbiome
#' @export
create_training_data <- function(data, TPR, TPR_B, perc_strata,norm_method='TSS',
                                 threshold_I1, threshold_I01, threshold_I02,
                                 OTU_replace,
                                 combinations,
                                 combinations_I0,
                                 combinations_I1,
                                 combinations_B, path, k) {
  
  nOtu <- ntaxa(data)
  keep_otu <- c()
  i <- 1
  
  # Number of Pseudo OTUs to sample (per strata)
  size_I0 <- (1 - TPR - TPR_B) * nOtu
  size_I1 <- TPR * nOtu
  size_B <- TPR_B * nOtu
  
  size_I0_strata <- size_I0 * perc_strata
  size_I1_strata <- size_I1 * perc_strata
  size_B_strata <- size_B * perc_strata
  
  # Sampling Pseudo OTUs to create training data
  set.seed(9827) 
  while(length(keep_otu) <= nOtu) {
    
    if (i == 1) {
      # Sample rows for each strata 
      sampled_data_I0 <- combinations_I0 %>%
        group_by(strata) %>%
        sample_n(size = ceiling(size_I0_strata), replace = OTU_replace)
      
      sampled_data_I1 <- combinations_I1 %>%
        group_by(strata) %>%
        sample_n(size = ceiling(size_I1_strata), replace = OTU_replace)
      
      sampled_data_B <- combinations_B %>%
        group_by(strata) %>%
        sample_n(size = ceiling(size_B_strata), replace = OTU_replace)
      
      sampled_data <- rbind(sampled_data_I0, sampled_data_I1, sampled_data_B)
      
      # Remaining pairs
      remaining_df_I0 <- anti_join(combinations_I0, sampled_data_I0, by = "name")
      remaining_df_I1 <- anti_join(combinations_I1, sampled_data_I1, by = "name")
      remaining_df_B <- anti_join(combinations_B, sampled_data_B, by = "name")
      
      pair_id <- which(combinations$name %in% sampled_data$name)
      
    } else {
      sampled_data_I0_extra <- c()
      sampled_data_B_extra <- c()
      sampled_data_I1_extra <- c()
      
      if(length(keep_id_t01) < size_I0) {
        n_I0 <- (size_I0 - length(keep_id_t01)) * perc_strata
        sampled_data_I0_extra <- combinations_I0 %>%
          group_by(strata) %>%
          sample_n(size = ceiling(n_I0), replace = OTU_replace)
      }
      
      if (length(keep_id_I1) < size_I1) {
        n_I1 <- (size_I1 - length(keep_id_I1)) * perc_strata
        sampled_data_I1_extra <- combinations_I1 %>%
          group_by(strata) %>%
          sample_n(size = ceiling(n_I1), replace = OTU_replace)
      }
      
      if (length(keep_id_B) < size_B) {
        n_IB <- (size_B - length(keep_id_B)) * perc_strata
        sampled_data_B_extra <- combinations_B %>%
          group_by(strata) %>%
          sample_n(size = ceiling(n_IB), replace = OTU_replace)
      }
      
      sampled_data <- rbind(sampled_data_I0_extra, sampled_data_B_extra, sampled_data_I1_extra)
      pair_id_extra <- which(combinations$name %in% sampled_data$name)
      pair_id <- c(which(combinations$name %in% keep_otu), pair_id_extra)
    }
    
    # Construct count table 
    group_0 <- subset_samples(data, group == 0)
    counts <- data.frame(otu_table(group_0))  
    
    combinations_subset <- combinations[pair_id, ] 
    
    counts_training_tmp <- cbind(counts[combinations_subset$otu1, ],
                                 counts[combinations_subset$otu2, ])
    
    rownames(counts_training_tmp) <- make.unique(combinations_subset$name)
    colnames(counts_training_tmp) <- c(paste0(colnames(counts), "-grp1"),
                                       paste0(colnames(counts), "-grp2"))
    
    # Create phyloseq object
    otuTab <- otu_table(counts_training_tmp, taxa_are_rows = TRUE)
    tax_data <- tax_table(data.frame(make.unique(combinations_subset$name)))
    rownames(tax_data) <- make.unique(combinations_subset$name)
    sample_data <- sample_data(data.frame(group = c(rep(0, ncol(counts)), rep(1, ncol(counts)))))
    rownames(sample_data) <- colnames(counts_training_tmp)
    data_tmp <- merge_phyloseq(otuTab, sample_data, tax_data)
    
    # Normalize & calculate lfc of temporary training data 
    if (norm_method == "TSS") {
      # TSS normalization
      data_tmp_norm <- transform_sample_counts(data_tmp, function(x) x / sum(x))
    } else if (norm_method == "CLR") {
      # CLR normalization
      data_tmp_norm <- microbiome::transform(data_tmp, transform = "clr")
    } else {
      stop("Invalid normalization method. Choose either 'TSS' or 'CLR'.")
    }
    
    lfc_tmp <- calc_fc(subset_samples(data_tmp_norm, group == 1),
                       subset_samples(data_tmp_norm, group == 0),
                       norm_method)
    
    # Define I1/B OTUs of temporary training data
    keep_id_I1 <- which(abs(lfc_tmp) > threshold_I1)
    keep_id_B <- which((abs(lfc_tmp) > threshold_I02) & (abs(lfc_tmp) < threshold_I1))
    
    # Define I0 based on two thresholds
    keep_id_t02 <- which(abs(lfc_tmp) < threshold_I02)
    select <- names(lfc_tmp[keep_id_t02])
    I0_t02 <- prune_taxa(select, data_tmp_norm)  
    
    # Permute I0
    group_permuted <- sample(sample_data(I0_t02)$group)
    sample_data(I0_t02)$group <- group_permuted
    
    lfc_tmp <- calc_fc(subset_samples(I0_t02, group == 1),
                       subset_samples(I0_t02, group == 0), norm_method)
    
    keep_id_t01 <- which(abs(lfc_tmp) < threshold_I01)
    
    # OTU to keep 
    keep_otu <- c(names(keep_id_I1), names(keep_id_t01), names(keep_id_B))
    
    i <- i + 1
  }
  
  # Create phyloseq object of training data 
  pair_id <- which(combinations$name %in% keep_otu)
  combinations_final <- combinations[pair_id, ] 
  
  counts_training <- cbind(counts[combinations_final$otu1, ],
                           counts[combinations_final$otu2, ])
  
  rownames(counts_training) <- combinations_final$name
  colnames(counts_training) <- c(paste0(colnames(counts), "-grp1"),
                                 paste0(colnames(counts), "-grp2"))
  
  # Phyloseq object
  otuTab <- otu_table(counts_training, taxa_are_rows = TRUE)
  tax_data <- tax_table(data.frame(combinations_final$name))
  rownames(tax_data) <- combinations_final$name
  sample_data <- sample_data(data.frame(group = c(rep(0, ncol(counts)), rep(1, ncol(counts)))))
  rownames(sample_data) <- c(paste0(colnames(counts), "-grp1"),
                             paste0(colnames(counts), "-grp2"))
  
  data_training <- merge_phyloseq(otuTab, sample_data, tax_data)
  
  if (norm_method == "TSS") {
    # TSS normalization
    data_training_norm <- transform_sample_counts(data_training, function(x) x / sum(x))
  } else if (norm_method == "CLR") {
    # CLR normalization
    data_training_norm <- microbiome::transform(data_training, transform = "clr")
  } else {
    stop("Invalid normalization method. Choose either 'TSS' or 'CLR'.")
  }
  
  
  lfc_final <- calc_fc(subset_samples(data_training_norm, group == 0),
                       subset_samples(data_training_norm, group == 1),norm_method)
  
  lfc_final <- data.frame(lfc_final)
  lfc_final[abs(lfc_final$lfc_final) <= threshold_I01, 'Threshold'] <- "I_0"
  lfc_final[lfc_final$lfc_final > threshold_I1, 'Threshold'] <- "I_1"
  lfc_final[lfc_final$lfc_final < -threshold_I1, 'Threshold'] <- "B"
  
  data_training@tax_table@.Data=cbind(data_training@tax_table@.Data[,1],lfc_final$lfc_final,lfc_final$Threshold)
  
  # Save the training data
  save(data_training, file = paste0(path, "/Training_data", k, ".RData"))
  
  return(data_training)
}
