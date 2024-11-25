#' R-sign Methods for Sparse Data
#'
#' This function applies R-sign methods for sparse data using the RioNorm2-adjusted hk_find_A function.
#'
#' @param data A phyloseq object containing count data and metadata.
#' @param model_R The R-sign model to apply ("All", "Marginal", or "RI").
#'
#' @return A list or data frame with test statistics, p-values, and adjusted p-values for the selected model.
#'
#' @examples
#' \dontrun{
#'   results <- RsignRI_test_sparse(data, model_R = "All")
#' }
#'
#' @importFrom phyloseq prune_taxa merge_phyloseq sample_sums
#' @importFrom stats p.adjust
#' @importFrom signtrans R_sign
#' @export
hk_find_A <- function (OTU_table, min_avg_counts = 5) 
{
  nsamples = dim(OTU_table)[2]
  samobs = apply(OTU_table, 1, function(x) sum(x != 0))
  #hk_pool = OTU_table[samobs >= nsamples * 0.8, ] # Origineel
  hk_pool = OTU_table[samobs >= nsamples * 0.6, ]
  avg_count = apply(hk_pool, 1, mean)
  hk_pool = hk_pool[avg_count >= min_avg_counts, ]
  nOTUs = dim(hk_pool)[1]
  OTU_ID = rownames(hk_pool)
  ratio_var = matrix(0, nrow = nOTUs, ncol = nOTUs)
  for (i in 1:(nOTUs - 1)) {
    for (j in (i + 1):nOTUs) {
      mul = (hk_pool[i, ] * hk_pool[j, ])
      ind = unname(unlist(mul)) != 0
      ratio_var[i, j] = var(unlist(log(hk_pool[i, ][ind]/hk_pool[j, ][ind])))
      ratio_var[j, 1] <- ratio_var[i, j]
    }
  }
  ratio_var[lower.tri(ratio_var, diag = TRUE)] <- 0
  dist = ratio_var[ratio_var > 0]
  First = TRUE
  library(igraph)
  
  #quans = seq(0.01, 0.05, 0.005) # Origineel
  quans = seq(0.01, 0.3, 0.005)
  
  for (q in 1:length(quans)) {
    print(q)
    FindIndex <- length(quans)    #the largest one
    nodes = data.frame(seq(1, nOTUs, 1), OTU_ID)
    colnames(nodes) = c("order", "OTU_ID")
    links = matrix(0, nrow = 1, ncol = 3)
    dist = ratio_var[ratio_var > 0]
    h = quantile(dist, probs = quans[q])
    order = seq(1, nOTUs, 1)
    for (i in 1:dim(ratio_var)[1]) {
      check = ratio_var[i, ]
      ind = order[check < h & check > 0]
      if (length(ind) > 0) {
        dist = check[ind]
        rep = replicate(length(ind), i)
        record = cbind(rep, ind, dist)
        links = rbind(links, record)
      }
    }
    
    if (dim(links)[1] > 3)
    {
      FindIndex = q
      break
    }
  }
  
  riOTUs <- NA
  
  for (q in FindIndex:length(quans)) {
    print(q)
    nodes = data.frame(seq(1, nOTUs, 1), OTU_ID)
    colnames(nodes) = c("order", "OTU_ID")
    links = matrix(0, nrow = 1, ncol = 3)
    dist = ratio_var[ratio_var > 0]
    h = quantile(dist, probs = quans[q])
    order = seq(1, nOTUs, 1)
    for (i in 1:dim(ratio_var)[1]) {
      check = ratio_var[i, ]
      ind = order[check < h & check > 0]
      if (length(ind) > 0) {
        dist = check[ind]
        rep = replicate(length(ind), i)
        record = cbind(rep, ind, dist)
        links = rbind(links, record)
      }
    }
    
    links = links[-1, ]
    net <- graph_from_data_frame(d = links, vertices = nodes, directed = F)
    list_cliques = cliques(net)
    clique_size = sapply(list_cliques, length)
    curr_largest_cliques = largest_cliques(net)
    curr_length = length(curr_largest_cliques[[1]])
    curr_largest_cliques = matrix(unlist(curr_largest_cliques), ncol = curr_length, byrow = T)
    if (curr_length < 3) {
      next
    }
    if (First == TRUE) {
      prev_largest_cliques = curr_largest_cliques
      prev_length = dim(prev_largest_cliques)[2]
      First = FALSE
      next
    }
    if (curr_length == prev_length) {
      riOTUs = prev_largest_cliques
      break
    }
    print(prev_largest_cliques)
    print(curr_largest_cliques)
    updated_cliques = matrix(0, nrow = 1, ncol = curr_length)
    for (i in 1:dim(curr_largest_cliques)[1]) {
      for (j in 1:dim(prev_largest_cliques)[1]) {
        if (sum(prev_largest_cliques[j, ] %in% curr_largest_cliques[i, ]) == prev_length) {
          updated_cliques = rbind(updated_cliques, curr_largest_cliques[i, ])
          break
        }
      }
    }
    if (dim(updated_cliques)[1] == 1) {
      riOTUs = prev_largest_cliques
      break
    }
    else if (dim(updated_cliques)[1] == 2) {
      prev_largest_cliques = matrix(updated_cliques[-1, ], nrow = 1)
      prev_length = dim(prev_largest_cliques)[2]
    }
    else {
      prev_largest_cliques = updated_cliques[-1, ]
      prev_length = dim(prev_largest_cliques)[2]
    }
  }
  
  if(length (riOTUs) == 1)
  {print ("No riOTUs are identified at 10th percentile")}
  
  if(length(riOTUs) > 1)
  {
    riOTUs = riOTUs[order(unique(riOTUs))]
    #riOTUs = riOTUs[order(riOTUs)]
    riOTUs_ID = rownames(hk_pool)[riOTUs]
    riOTUs_count = hk_pool[riOTUs, ]
    size_factor = colSums(riOTUs_count)
    return(list(riOTUs_ID = riOTUs_ID, size_factor = size_factor))
  }
}
RsignRI_test_sparse <- function(data, model_R) {
  OTU_table = data.frame(otu_table(data))
  
  test = hk_find_A(OTU_table)
  
  ref = test$riOTUs_ID
  
  ref_phy = prune_taxa(ref,data)
  
  ref_median = apply(data.frame(otu_table(ref_phy),row.names = NULL),2,median)
  
  `%notin%` <- Negate(`%in%`)
  data_final = subset(otu_table(data),rownames(otu_table(data)) %notin% ref)
  
  data_final <- merge_phyloseq(data_final, tax_table(data), sample_data(data))
  
  
  R_sign_results_all<-list()
  stat_Marginal_R <- data.frame()
  stat_Conditional_R <- data.frame()
  stat_RI_R <- data.frame()
  
  pval_Marginal_R <- data.frame()
  pval_Conditional_R <- data.frame()
  pval_RI_R <- data.frame()
  
  if (length(ref)!=0){
    for (j in (1:ntaxa(data_final))){
      # Define a dataframe with raw count of 1 taxon (column 1) and the average reference frame (column 2)
      count=as.vector(data_final@otu_table[j,])
      
      med_count = c()
      if (median(count)==0){
        med_count = min(count[count!=0])
      } else{
        med_count = median(count)
      }
      count_med = median(count)
      db <- data.frame(counts=count,ref=ref_median*(med_count/median(ref_median)),row.names = NULL) # db[,1]= count taxa, db[,2] = mean reference frame
      
      
      # Compute Library size
      
      libsize <- sample_sums(data_final)
      
      # Define starting sample data
      
      sample_data(data_final)$libsize<-libsize
      sample_data(data_final)$group<-as.numeric(sample_data(data_final)$group)
      
      
      # Define model to be tested 
      formula_R=RATIO~as.factor(group)+libsize
      
      # APPly R-sign methods
      R_sign_results_all[[j]]<-signtrans::R_sign(formula_R,db,startdata=data_final,Method=model_R)
      
      if (model_R =='All'){
        pval_Marginal_R <- rbind(pval_Marginal_R ,R_sign_results_all[[j]]$Marginal$Pval)
        pval_Conditional_R <- rbind(pval_Conditional_R,R_sign_results_all[[j]]$Conditional$Pval)
        pval_RI_R <- rbind(pval_RI_R,R_sign_results_all[[j]]$RI$Pval)
        
        stat_Marginal_R <- rbind(stat_Marginal_R ,R_sign_results_all[[j]]$Marginal$Teststatistic)
        stat_Conditional_R <- rbind(stat_Conditional_R,R_sign_results_all[[j]]$Conditional$Teststatistic)
        stat_RI_R <- rbind(stat_RI_R,R_sign_results_all[[j]]$RI$Teststatistic)
      }
      if (model_R =='Marginal'){
        pval_Marginal_R <- rbind(pval_Marginal_R ,R_sign_results_all[[j]]$Pval)
        stat_Marginal_R <- rbind(stat_Marginal_R ,R_sign_results_all[[j]]$Teststatistic)
      }
      if (model_R =='RI'){
        pval_RI_R <- rbind(pval_RI_R,R_sign_results_all[[j]]$Pval)
        stat_RI_R <- rbind(stat_RI_R,R_sign_results_all[[j]]$Teststatistic)
      }
      
      
      
      
    }}
  
  if (model_R =='All'){
    adjPVals_marg <- p.adjust(unlist(pval_Marginal_R),method='BH') 
    adjPVals_cond <- p.adjust(unlist(pval_Conditional_R),method='BH')
    adjPVals_ri<- p.adjust(unlist(pval_RI_R),method='BH')
    
    Marginal = data.frame("pval" = unlist(pval_Marginal_R), "adjP" = adjPVals_marg, "stat"=unlist(stat_Marginal_R))
    Conditional = data.frame("pval" = unlist(pval_Conditional_R), "adjP" = adjPVals_cond, "stat"=unlist(stat_Conditional_R))
    RI = data.frame("pval" = unlist(pval_RI_R), "adjP" = adjPVals_ri, "stat"=unlist(stat_RI_R))
    
    rownames(otu_table(data_final)) -> rownames(Marginal)
    rownames(otu_table(data_final)) -> rownames(Conditional)
    rownames(otu_table(data_final)) -> rownames(RI)
    
    
    out = list(Marginal = Marginal,
               Conditional = Conditional,
               RI = RI)
  }
  if (model_R =='Marginal'){
    adjPVals_marg <- p.adjust(unlist(pval_Marginal_R),method='BH') 
    
    Marginal = data.frame("pval" = unlist(pval_Marginal_R), "adjP" = adjPVals_marg, "stat"=unlist(stat_Marginal_R))
    
    rownames(otu_table(data_final)) -> rownames(Marginal)
    
    
    
    out =  Marginal
  }
  if (model_R =='RI'){
    
    adjPVals_ri<- p.adjust(unlist(pval_RI_R),method='BH')
    
    RI = data.frame("pval" = unlist(pval_RI_R), "adjP" = adjPVals_ri, "stat"=unlist(stat_RI_R))
    
    rownames(otu_table(data_final)) -> rownames(RI)
    
    
    out = RI
  }
  
  
  
  return(out)
}