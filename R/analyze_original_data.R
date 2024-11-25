#' Analyze Original Data
#'
#' This function analyzes the original data using the best-performing method, evaluates the results, 
#' and provides performance metrics such as Type I error rate, sensitivity, and FDR.
#'
#' @param data A phyloseq object containing the original dataset.
#' @param method The statistical method to apply, as a string.
#' @param res_best_method A data frame containing the results of the best-performing method, 
#'   including metrics such as Type I error rate, sensitivity, and FDR.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{final_res}}{Performance metrics (Type I error rate, sensitivity, FDR).}
#'     \item{\code{method}}{The statistical method used.}
#'   }
#'
#' @examples
#' \dontrun{
#'   results <- analyze_original_data(data, method = "DESeq2", res_best_method = res_best_method)
#' }
#'
#' @importFrom signtrans Trim
#' @importFrom phyloseq tax_table
#' @importFrom stats optimise
#' @export
best_method = function(data, method){
  if (method=='Aldex2'){
    
    res = tryCatch({aldex_method(data)$ALDEx2_w}, error=function(e){cat("Error")})
    
  }
  if (method=='Corncob'){
    
    res = tryCatch({corncob_method(data)}, error=function(e){cat("Error")})
  }
  if (method=='ANCOM-bc'){
    
    res = tryCatch({ANCOMBC_method(data)}, error=function(e){cat("Error")})
  }
  if (method=='Wilcoxon'){
    
    res = tryCatch({WilcoxTest(data)}, error=function(e){cat("Error")})
  }
  
  if (method=='Wilcoxon_clr'){
    
    res = tryCatch({WilcoxTest_clr(data)}, error=function(e){cat("Error")})
  }
  
  if (method=='R-sign_marg'){
    
    res = tryCatch({RsignRI_test(data,model_R='Marginal')}, error=function(e){cat("Error")})
    
    if (is.null(res)){
      res = tryCatch({RsignRI_test_sparse(data,model_R="Marginal")}, error=function(e){cat("Error")})
    }
    
  }
  if (method=='R-sign_RI'){
    
    res = tryCatch({RsignRI_test(data,model_R="RI")}, error=function(e){cat("Error")})
    
    
    if (is.null(res)){
      res = tryCatch({RsignRI_test_sparse(data,model_R="RI")}, error=function(e){cat("Error")})
    }
    
    
  }
  
  if (method=='LinDA'){
    
    res = tryCatch({LinDA_test(data)}, error=function(e){cat("Error")})
    
  }
  
  if (method=='DESeq2'){
    
    res = tryCatch({deseq2_method(data)}, error=function(e){cat("Error")})
    
  }
  return(res)
}
analyze_original_data <- function(data, method, res_best_method,alpha, eval_type="Sim") {
  
  # Trim the original dataset
  data_org_trim <- signtrans::Trim(data, minPrev = 0.05)
  
  # Apply the best method
  test <- best_method(data_org_trim, method)
  
  if (is.null(test)) {
    res_best_method <- res_best_method$res_all[method != res_best_method$res_all[, "Method"], ]
    
    # FDR control threshold
    uplim <- 0.07
    res_all_control <- res_best_method[as.numeric(res_best_method[, "FDR"]) <= uplim, ]
    
    if (!is.null(dim(res_all_control)[1])) {
      res_best <- res_all_control[res_all_control[, "Sensitivity"] == max(res_all_control[, "Sensitivity"], na.rm = TRUE), ]
    } else {
      res_best <- res_all_control
    }
    
    if (!is.null(dim(res_best)[1])) {
      id <- sample(1:nrow(res_best), 1)
      res_best <- res_best[id, ]
    }
    
    test <- best_method(data_org_trim, res_best[4])
    method <- res_best[4]
  }
  
  if (eval_type=="Sim"){
  # Evaluate results
  DE.ind <- tax_table(data_org_trim)[, 1]
  DA_OTU <- which(DE.ind == TRUE)
  NDA_OTU <- which(DE.ind == FALSE)
  
  eval_otu <- rownames(test)
 
  
  if (!is.null(test$adjP)) {
    if (!grepl("sign", method)) {
      final_res <- eval_sums(test$adjP, NDA_OTU, DA_OTU,alpha)
    } else {
      DE.ind_adj <- DE.ind[rownames(DE.ind) %in% rownames(test)]
      DA_OTU <- which(DE.ind_adj == TRUE)
      NDA_OTU <- which(DE.ind_adj == FALSE)
      final_res <- eval_sums(test$adjP, NDA_OTU, DA_OTU,alpha)
    }
  } else {
    final_res <- "no result"
  }
  }else{
    final_res <- "no result"
  }
  out <- list(
    test_res= test,
    final_res = final_res,
    method = method
  )
  return(out)
}
