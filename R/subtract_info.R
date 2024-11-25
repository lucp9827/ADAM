#' Subtract information from dataset based on thresholds
#'
#' This function processes a microbiome dataset by normalizing the data using CLR (centered log-ratio) or TSS (Total Sum Scaling) transformation. 
#' It calculates log fold changes (lfc) between two groups and applies thresholds to determine true positive rates (TPR). 
#' The function then creates strata for the OTUs based on relative abundances and thresholds.
#'
#' @param data A phyloseq object containing OTU data for analysis.
#' @param threshold_I1 Numeric threshold for I1. If not provided, calculated from log fold changes.
#' @param threshold_I01 Numeric threshold for I01. If not provided, calculated from log fold changes.
#' @param threshold_I02 Numeric threshold for I02. If not provided, calculated from log fold changes.
#' @param threshold_value Logical indicating whether to compute the thresholds from the data (default is FALSE).
#' @param norm_method Character, either "TSS" or "CLR" to specify the normalization method (default: "TSS").
#' @param log_file Path to a log file where data-driven thresholds are logged (default: "log.txt").
#' @param perc_strata Numeric, indicating the fraction of data to be used per strata (default: 0.1).
#' @return A list with breakpoints for the different strata (I0, I1, between), TPR values, sample counts, and thresholds.
#' @examples
#' data <- example_data()  # Your data loading code here
#' result <- subtract_info(data)
#' @export
subtract_info <- function(data, 
                          threshold_I1=NULL, 
                          threshold_I01=NULL, 
                          threshold_I02=NULL, 
                          threshold_value=FALSE, 
                          norm_method = "TSS",
                          log_file = "test_log.txt",
                          perc_strata = 0.1) {
  
  # Get number of OTUs and samples
  nOtu = ntaxa(data)
  nSamples = nsamples(data)
  
  # Normalize data based on the specified normalization method
  if (norm_method == "TSS") {
    # TSS normalization
    data_normalized <- transform_sample_counts(data, function(x) x / sum(x))
  } else if (norm_method == "CLR") {
    # CLR normalization
    data_normalized <- microbiome::transform(data, transform = "clr")
  } else {
    stop("Invalid normalization method. Choose either 'TSS' or 'CLR'.")
  }
  
  # Subset data by groups
  group_0_norm = subset_samples(data_normalized, group == 0)
  group_1_norm = subset_samples(data_normalized, group == 1)
  
  # Calculate log fold changes (lfc) based on normalized data
  logfoldchanges_norm = calc_fc(group_0_norm, group_1_norm,norm_method)
  
  # Initialize log message
  log_message = ""
  
  # If thresholds are not provided, calculate based on the data
  if (threshold_value == TRUE) {
    
    threshold_I1 = median(logfoldchanges_norm[is.finite(logfoldchanges_norm)], na.rm = TRUE) +
      (3 * mad(logfoldchanges_norm[is.finite(logfoldchanges_norm)], constant = 1, na.rm = TRUE))
    threshold_I01 = median(logfoldchanges_norm[is.finite(logfoldchanges_norm)], na.rm = TRUE)
    threshold_I02 = median(logfoldchanges_norm[is.finite(logfoldchanges_norm)], na.rm = TRUE) +
      1.5 * mad(logfoldchanges_norm[is.finite(logfoldchanges_norm)], constant = 1, na.rm = TRUE)
    
    log_message <- paste("Data-driven thresholds were applied based on log fold changes:",
                         "\nI1 threshold (high impact): Median log fold change + 3 * MAD =", threshold_I1,
                         "\nI01 threshold (baseline): Median log fold change =", threshold_I01,
                         "\nI02 threshold (moderate impact): Median log fold change + 1.5 * MAD =", threshold_I02,
                         sep = "\n")
    
    # Print the message to console
    cat(log_message, "\n")
    
    # Ensure the log file is created if it doesn't exist
    if (!file.exists(log_file)) {
      file.create(log_file)
    }
    
    # Write the message to log file
    writeLines(log_message, log_file)
  }
  
  # Calculate true positive rates (TPR)
  TPR = sum(logfoldchanges_norm > threshold_I1, na.rm = TRUE) / nOtu
  TPR_B = sum((logfoldchanges_norm < threshold_I1) & (logfoldchanges_norm > threshold_I01), na.rm = TRUE) / nOtu
  
  # Sum relative abundance per OTU
  otu_sum = taxa_sums(data_normalized)
  
  # Combine results: log fold changes and sum of relative abundance
  df = data.frame('lfc' = logfoldchanges_norm, 'Sum_RA' = otu_sum)
  
  # Add threshold group information to df
  df <- df %>%
    mutate(I = ifelse(is.na(abs(lfc)), "I_1",
                      ifelse(abs(lfc) < threshold_I01, "I_0",
                             ifelse(abs(lfc) > threshold_I1, "I_1", "I_between"))))
  
  # Order data by sum relative abundance (Sum_RA)
  df_I0 = df[df$I == 'I_0',]
  df_I1 = df[df$I == 'I_1',]
  df_B = df[df$I == 'I_between',]
  
  ordered_I0 <- df_I0[order(df_I0$Sum_RA, decreasing = TRUE),]
  ordered_I1 <- df_I1[order(df_I1$Sum_RA, decreasing = TRUE),]
  ordered_B <- df_B[order(df_B$Sum_RA, decreasing = TRUE),]
  
  # Determine breakpoints based on ordered sumRA
  breakpoints_I0 <- seq(0, nrow(ordered_I0), by = nrow(ordered_I0) * perc_strata)
  breakpoints_I1 <- seq(0, nrow(ordered_I1), by = nrow(ordered_I1) * perc_strata)
  breakpoints_B <- seq(0, nrow(ordered_B), by = nrow(ordered_B) * perc_strata)
  
  ordered_I0$strata <- cut(1:nrow(ordered_I0), breaks = breakpoints_I0, labels = FALSE, include.lowest = TRUE)
  ordered_I1$strata <- cut(1:nrow(ordered_I1), breaks = breakpoints_I1, labels = FALSE, include.lowest = TRUE)
  ordered_B$strata <- cut(1:nrow(ordered_B), breaks = breakpoints_B, labels = FALSE, include.lowest = TRUE)
  
  # Summarize breakpoints for different strata
  break_I0 <- ordered_I0 %>%
    group_by(strata) %>%
    summarize(min_Sum_RA = min(Sum_RA, na.rm = TRUE),
              max_Sum_RA = max(Sum_RA, na.rm = TRUE))
  
  break_I1 <- ordered_I1 %>%
    group_by(strata) %>%
    summarize(min_Sum_RA = min(Sum_RA, na.rm = TRUE),
              max_Sum_RA = max(Sum_RA, na.rm = TRUE))
  
  break_B <- ordered_B %>%
    group_by(strata) %>%
    summarize(min_Sum_RA = min(Sum_RA, na.rm = TRUE),
              max_Sum_RA = max(Sum_RA, na.rm = TRUE))
  
  # Output
  out = list(break_I0 = break_I0,
             break_I1 = break_I1,
             break_B = break_B,
             TPR = TPR,
             TPR_B = TPR_B,
             nSamples = nSamples,
             threshold_I1 = threshold_I1,
             threshold_I01 = threshold_I01,
             threshold_I02 = threshold_I02)
  
  return(out)
}
