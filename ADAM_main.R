#' ADAM - Adaptive Thresholds Main Function
#'
#' This is the main function of the ADAM (Adaptive Thresholds) workflow. It performs preprocessing,
#' Test data creation, bootstrapping, method selection, and final analysis to identify the
#' best statistical method for differential abundance analysis.
#'
#' @param data A dataset containing bulk simulation data or OTU table to analyze.
#' @param norm_method Normalization method to use. Default is 'TSS'.
#' @param perc_strata Proportion of strata for breaking data. Default is 0.2.
#' @param OTU_replace Logical indicating whether to replace OTUs during training data generation. Default is TRUE.
#' @param threshold_value Logical indicating whether to use thresholding in the analysis. Default is TRUE.
#' @param alpha Significance level for FDR control. Default is 0.05.
#' @param B Number of bootstrap/permutation datasets. Default is 2.
#' @param path File path to save results. Default is the current working directory.
#' @param dataset A numeric identifier for the dataset. Default is 1.
#'
#' @return A list containing the final results of the ADAM workflow.
#'
#' @examples
#' \dontrun{
#'   library(dplyr)
#'   path <- getwd()
#'   dataset <- 1
#'   perc_strata <- 0.2
#'   OTU_replace <- TRUE
#'   norm_method <- 'TSS'
#'   B <- 2
#'   threshold_value <- TRUE
#'   alpha <- 0.05
#'   data <- sim.data.bulk.p_B[[1]]
#'
#'   final_results <- ADAM_main(
#'     data = data,
#'     norm_method = norm_method,
#'     perc_strata = perc_strata,
#'     OTU_replace = OTU_replace,
#'     threshold_value = threshold_value,
#'     alpha = alpha,
#'     B = B,
#'     path = path,
#'     dataset = dataset
#'   )
#' }
#'
#' @importFrom dplyr %>%
#' @export
ADAM_main <- function(
    data,
    norm_method = 'TSS',
    perc_strata = 0.2,
    OTU_replace = TRUE,
    threshold_value = TRUE,
    alpha = 0.05,
    B = 1,
    path = getwd(),
    dataset = 1
) {
  # Step 1: Obtain information from original data
  break_info <- subtract_info(
    data = data,
    threshold_value = threshold_value,
    norm_method = norm_method,
    perc_strata = perc_strata
  )

  # Step 2: Pairwise comparisons of taxa
  pairwise_info <- pairwise_comp(
    data = data,
    break_I0 = break_info$break_I0,
    break_I1 = break_info$break_I1,
    break_B = break_info$break_B,
    threshold_I1 = break_info$threshold_I1,
    threshold_I01 = break_info$threshold_I01,
    norm_method = norm_method
  )

  # Step 3: Create test data
  training_data <- create_training_data(
    data = data,
    TPR = break_info$TPR,
    TPR_B = break_info$TPR_B,
    perc_strata = perc_strata,
    norm_method = norm_method,
    threshold_I1 = break_info$threshold_I1,
    threshold_I01 = break_info$threshold_I01,
    threshold_I02 = break_info$threshold_I02,
    OTU_replace = OTU_replace,
    combinations = pairwise_info$combinations,
    combinations_I0 = pairwise_info$combinations_I0,
    combinations_I1 = pairwise_info$combinations_I1,
    combinations_B = pairwise_info$combinations_B,
    path = path,
    k = dataset
  )

  # Step 4: Perform permutation/bootstrap to create variations of test data
  res_intermediate <- perm_boot(
    data_training = training_data,
    B = B,
    nSamples = break_info$nSamples,
    path = path,
    dataset = dataset
  )

  # Step 5: Choose the best-performing method
  res_best_method <- choose_best_method(
    data_training = training_data,
    res_intermediate = res_intermediate,
    alpha = alpha,
    B = B,
    k = dataset,
    path = path
  )

  # Step 6: Analyze the original data using the best method
  final <- analyze_original_data(
    data = data,
    method = res_best_method$res_best['Method'],
    res_best_method = res_best_method,
    alpha = alpha,
    eval_type = "Sim"
  )

  return(final)
}

