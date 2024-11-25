#' Choose Best Method
#'
#' This function evaluates multiple statistical methods based on their performance on permutation/bootstrap datasets.
#' It calculates metrics such as Type I error rate, sensitivity, and FDR, and selects the best-performing method.
#'
#' @param data_training A phyloseq object containing training data.
#' @param res_intermediate A list of results from intermediate analysis using different methods.
#' @param B Number of bootstrap/permutation datasets.
#' @param k An identifier for the training dataset.
#' @param path Directory path to save the results.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{res_all}}{Summary metrics for all methods.}
#'     \item{\code{res_best}}{The best-performing method based on the evaluation criteria.}
#'     \item{\code{res_all_tmp}}{Detailed results for all methods.}
#'   }
#'
#' @examples
#' \dontrun{
#'   best_method <- choose_best_method(data_training, res_intermediate, B = 100, k = 1, path = "results/")
#' }
#'
#' @importFrom signtrans Trim
#' @importFrom phyloseq otu_table tax_table
#' @export
choose_best_method <- function(data_training, res_intermediate, alpha, B, k, path) {

  # Trim training data
  data_trim <- signtrans::Trim(data_training, minPrev = 0.05)
  counts_trim <- data.frame(otu_table(data_trim))

  # Identify indices of I0/I1 OTUs
  I0_id <- which(tax_table(data_trim)[, 3] == "I_0")
  I1_id <- which(tax_table(data_trim)[, 3] == "I_1")

  # Initialize result data frames for each method
  result_aldex <- data.frame()
  result_corncob <- data.frame()
  result_ancombc <- data.frame()
  result_wilcox <- data.frame()
  result_wilcox_clr <- data.frame()
  result_marg <- data.frame()
  result_RI <- data.frame()
  result_LinDA <- data.frame()
  result_deseq2 <- data.frame()

  # Evaluate each method across B datasets
  for (i in seq_len(B)) {
    res_aldex <- res_intermediate[[i]][["res_aldex2_w"]]
    if (!is.null(res_aldex)) {
      result_aldex[i, 1:3] <- eval_sums(res_aldex$adjP, I0_id, I1_id, alpha)
    }

    res_corncob <- res_intermediate[[i]][["res_corncob"]]
    if (!is.null(res_corncob)) {
      result_corncob[i, 1:3] <- eval_sums(res_corncob$adjP, I0_id, I1_id, alpha)
    }

    res_ancombc <- res_intermediate[[i]][["res_ancombc"]]
    if (!is.null(res_ancombc)) {
      result_ancombc[i, 1:3] <- eval_sums(res_ancombc$adjP, I0_id, I1_id, alpha)
    }

    res_wilcox <- res_intermediate[[i]][["res_wilcox"]]
    if (!is.null(res_wilcox)) {
      result_wilcox[i, 1:3] <- eval_sums(res_wilcox$adjP, I0_id, I1_id, alpha)
    }

    res_wilcox_clr <- res_intermediate[[i]][["res_wilcox_clr"]]
    if (!is.null(res_wilcox_clr)) {
      result_wilcox_clr[i, 1:3] <- eval_sums(res_wilcox_clr$adjP, I0_id, I1_id, alpha)
    }

    res_linda <- res_intermediate[[i]][["res_linda"]]
    if (!is.null(res_linda)) {
      result_LinDA[i, 1:3] <- eval_sums(res_linda$adjP, I0_id, I1_id, alpha)
    }

    res_deseq <- res_intermediate[[i]][["res_deseq"]]
    if (!is.null(res_deseq)) {
      result_deseq2[i, 1:3] <- eval_sums(res_deseq$adjP, I0_id, I1_id, alpha)
    }

    # Reference taxa evaluation adjustments
    res_marg <- res_intermediate[[i]][["res_marginal"]]
    names_I0 <- rownames(counts_trim[I0_id, ])
    names_I1 <- rownames(counts_trim[I1_id, ])
    I0_id_adj <- which(rownames(res_marg) %in% names_I0)
    I1_id_adj <- which(rownames(res_marg) %in% names_I1)

    if (!is.null(res_marg)) {
      result_marg[i, 1:3] <- eval_sums(res_marg$adjP, I0_id_adj, I1_id_adj, alpha)
    }

    res_RI <- res_intermediate[[i]][["res_RI"]]
    if (!is.null(res_RI)) {
      result_RI[i, 1:3] <- eval_sums(res_RI$adjP, I0_id_adj, I1_id_adj, alpha)
    }
  }

  # Combine results into a single data frame
  Methods <- c("Aldex2", "Corncob", "ANCOM-bc", "Wilcoxon", "Wilcoxon_clr",
               "R-sign_marg", "R-sign_RI", "LinDA", "DESeq2")

  res_all <- rbind(
    colMeans(result_aldex, na.rm = TRUE),
    colMeans(result_corncob, na.rm = TRUE),
    colMeans(result_ancombc, na.rm = TRUE),
    colMeans(result_wilcox, na.rm = TRUE),
    colMeans(result_wilcox_clr, na.rm = TRUE),
    colMeans(result_marg, na.rm = TRUE),
    colMeans(result_RI, na.rm = TRUE),
    colMeans(result_LinDA, na.rm = TRUE),
    colMeans(result_deseq2, na.rm = TRUE)
  )

  res_all_tmp <- list(
    result_aldex,
    result_corncob,
    result_ancombc,
    result_wilcox,
    result_wilcox_clr,
    result_marg,
    result_RI,
    result_LinDA,
    result_deseq2
  )

  res_all <- cbind(res_all, Methods)
  colnames(res_all) <- c("Type I error rate", "Sensitivity", "FDR", "Method")

  # Apply FDR control
  uplim <- 0.07
  res_all_control <- res_all[res_all[, "FDR"] <= uplim, ]

  if (!is.null(dim(res_all_control)[1])) {
    res_best <- res_all_control[res_all_control[, "Sensitivity"] == max(res_all_control[, "Sensitivity"], na.rm = TRUE), ]
  } else {
    res_best <- res_all_control
  }

  if (!is.null(dim(res_best)[1])) {
    id <- sample(1:nrow(res_best))[1]
    res_best <- res_best[id, ]
  }

  out <- list(
    res_all = res_all,
    res_best = res_best,
    res_all_tmp = res_all_tmp
  )

  save(out, file = paste0(path, "res_trainingdata_", k, ".RData"))

  return(out)
}
