#' Perform Pairwise Comparisons of OTUs
#'
#' This function performs pairwise comparisons of OTUs in a microbiome dataset. It calculates various characteristics
#' of pseudo OTU pairs based on normalized data and stratifies them according to defined thresholds.
#'
#' @param data A phyloseq object containing the OTU data for analysis.
#' @param break_I0 A data frame containing breakpoints for the I0 strata.
#' @param break_I1 A data frame containing breakpoints for the I1 strata.
#' @param break_B A data frame containing breakpoints for the between strata.
#' @param threshold_I1 Numeric threshold for high impact OTUs.
#' @param threshold_I01 Numeric threshold for baseline impact OTUs.
#' @return A list containing:
#' \item{combinations}{A data frame of pairwise OTU comparisons and their characteristics.}
#' \item{combinations_I0}{Subset of combinations that fall into the I0 strata.}
#' \item{combinations_I1}{Subset of combinations that fall into the I1 strata.}
#' \item{combinations_B}{Subset of combinations that fall into the between strata.}
#' @examples
#' # result <- pairwise_comp(data, break_I0, break_I1, break_B, threshold_I1, threshold_I01)
#' @export
pairwise_comp = function(data, 
                         break_I0, 
                         break_I1, 
                         break_B, 
                         threshold_I1, 
                         threshold_I01,
                         norm_method='TSS') {
  
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
  
  # Subset data for group 0
  group_0_norm = subset_samples(data_normalized, group == 0)
  
  # Count table as data frame
  counts_normalized = data.frame(otu_table(group_0_norm))
  
  # Make all possible combinations of OTUs
  nOtu = ntaxa(data)
  combinations <- data.frame(combn(1:nOtu, 2, simplify = TRUE))
  
  # Calculate various characteristics of these pseudo OTU pairs
  for (i in 1:ncol(combinations)) {
    
    mean_1 = mean(as.numeric(counts_normalized[combinations[1, i], ]))
    mean_2 = mean(as.numeric(counts_normalized[combinations[2, i], ]))
    
    if (norm_method == "TSS") {
      logfoldchange = abs(log2(mean_1 / mean_2))
    } else if (norm_method == "CLR") {
      logfoldchange = abs((mean_1 - mean_2) / log(2))
    } else {
      stop("Invalid normalization method. Choose either 'TSS' or 'CLR'.")
    }
    
    # Store results in combinations data frame
    combinations[3, i] = logfoldchange
    combinations[4, i] = mean(c(as.numeric(counts_normalized[combinations[1, i], ]),
                                as.numeric(counts_normalized[combinations[2, i], ])))
    combinations[5, i] = var(c(as.numeric(counts_normalized[combinations[1, i], ]),
                               as.numeric(counts_normalized[combinations[2, i], ])))
    combinations[6, i] = mean(c(as.numeric(counts_normalized[combinations[1, i], ]),
                                as.numeric(counts_normalized[combinations[2, i], ])) == 0)
    combinations[7, i] = sum(c(as.numeric(counts_normalized[combinations[1, i], ]),
                               as.numeric(counts_normalized[combinations[2, i], ])))
  }
  
  # Adjust to long format and assign row names
  rownames(combinations) = c("otu1", "otu2", "lfc", "mean", "var", "zerofraction", "Sum_RA")
  combinations = data.frame(t(combinations))
  
  # Name pseudo OTUs
  combinations$name = paste(combinations$otu1, combinations$otu2, sep = '-')
  rownames(combinations) = combinations$name
  
  # Divide pairwise comparisons based on thresholds
  combinations <- combinations %>%
    mutate(I = ifelse(is.na(abs(lfc)), "I_1",
                      ifelse(abs(lfc) < threshold_I01, "I_0",
                             ifelse(abs(lfc) > threshold_I1, "I_1", "I_between"))))
  
  # Stratify pairwise comparisons
  combinations_I0 = combinations[combinations$I == 'I_0', ]
  combinations_I1 = combinations[combinations$I == 'I_1', ]
  combinations_B = combinations[combinations$I == 'I_between', ]
  
  combinations_I0 = combinations_I0[order(combinations_I0$Sum_RA, decreasing = TRUE), ]
  combinations_I1 = combinations_I1[order(combinations_I1$Sum_RA, decreasing = TRUE), ]
  combinations_B = combinations_B[order(combinations_B$Sum_RA, decreasing = TRUE), ]
  
  # Stratification based on breakpoints for combinations_B
  for (i in 1:nobs(break_B)) {
    if (i == 1) {
      combinations_B[(combinations_B$Sum_RA > break_B$max_Sum_RA[i + 1]), 'strata'] = i
    } else {
      if (i == nobs(break_B)) {
        combinations_B[(combinations_B$Sum_RA <= break_B$max_Sum_RA[i]), 'strata'] = i
      } else {
        combinations_B[(combinations_B$Sum_RA <= break_B$max_Sum_RA[i] &
                          combinations_B$Sum_RA > break_B$max_Sum_RA[i + 1]), 'strata'] = i
      }
    }
  }
  
  # Stratification for combinations_I1
  for (i in 1:nobs(break_I1)) {
    if (i == 1) {
      combinations_I1[(combinations_I1$Sum_RA > break_I1$max_Sum_RA[i + 1]), 'strata'] = i
    } else {
      if (i == nobs(break_I1)) {
        combinations_I1[(combinations_I1$Sum_RA <= break_I1$max_Sum_RA[i]), 'strata'] = i
      } else {
        combinations_I1[(combinations_I1$Sum_RA <= break_I1$max_Sum_RA[i] &
                           combinations_I1$Sum_RA > break_I1$max_Sum_RA[i + 1]), 'strata'] = i
      }
    }
  }
  
  # Stratification for combinations_I0
  for (i in 1:nobs(break_I0)) {
    if (i == 1) {
      combinations_I0[(combinations_I0$Sum_RA > break_I0$max_Sum_RA[i + 1]), 'strata'] = i
    } else {
      if (i == nobs(break_I0)) {
        combinations_I0[(combinations_I0$Sum_RA <= break_I0$max_Sum_RA[i]), 'strata'] = i
      } else {
        combinations_I0[(combinations_I0$Sum_RA <= break_I0$max_Sum_RA[i] &
                           combinations_I0$Sum_RA > break_I0$max_Sum_RA[i + 1]), 'strata'] = i
      }
    }
  }
  
  out = list(combinations = combinations,
             combinations_I0 = combinations_I0,
             combinations_I1 = combinations_I1,
             combinations_B = combinations_B)
  
  return(out)
}
