#1 ) load pRes 

#2) Add new results to the pRes

helper_add_newMethod_to_pRes <- function(pRes, path_newRes){
  
  newMethod <- data.frame(fread(path_newRes))
  
  # if('sum_nMut' %in% colnames(newMethod)){
  #   newMethod <- newMethod[newMethod$sum_nMut > 1,]
  # }
  
  idx_pVal <- which(colnames(newMethod) %in% c("pValues", "p_value",
                                               "raw_p", "pvals"))
  idx_ID <- which(colnames(newMethod) %in% c("binID", "PCAWG_ID"))
  
  newMethod <- data.frame("PCAWG_IDs" = newMethod[,idx_ID], 
                          "p_value" =  newMethod[,idx_pVal],
                          "q_value" = rep(NA, nrow(newMethod)),
                          "fdr" =  rep(NA, nrow(newMethod)))
  
  newMethod <- newMethod[which(!is.na(newMethod$p_value)),]
  
  pRes$newMethod <- newMethod
  class(pRes) = "pRes"
  pRes
}

add_newMethod_to_pRes <- function(pRes, PATHs_newResults, newRESULTS){
  
  for (i in 1:length(PATHs_newResults)) {
    pRes <- helper_add_newMethod_to_pRes(pRes, PATHs_newResults[i])
    names(pRes) <- c(names(pRes)[1:(length(pRes) -1)], newRESULTS[i])
  }
  
  pRes
  
}



#3) and convert pRes to two diferent pRes
split_cod_nonCod <- function(dat){
  
  dat <- dat[order(dat$p_value, decreasing = F),]
  
  dat_cds <- dat[grep("gc19_pc.cds", dat$PCAWG_IDs),]
  
  pattern <- paste(c("enhancers", "gc19_pc.3utr", "gc19_pc.5utr", "gc19_pc.promCore"), collapse = "|") #, "gc19_pc.ss"
  dat_nc <- dat[grep(pattern, dat$PCAWG_IDs), ]
  
  list("CDS" = dat_cds, "non_coding" = dat_nc)
  
}

split_pRes <- function(pRes){
  cds_nc_pRes <- lapply(pRes, split_cod_nonCod)
  cds_pRes <- lapply(cds_nc_pRes, function(s){
    s$CDS
  })
  
  nc_pRes <- lapply(cds_nc_pRes, function(s){
    s$non_coding
  })
  
  list('CDS_pRes' = cds_pRes, 'nc_pRes' = nc_pRes)
}


# 4) subset splitted pRes to common genes

check_n_elems_method <- function(df, mask_threshold, n_valid_IDs){
  df <- df[!grepl("lncrna.ncrna|lncrna.promCore|gc19_pc.ss|smallrna.ncrna|lncrna.ss|mirna_mat|mirna_pre|mirna.prom", 
                  df$PCAWG_IDs), ]
  
  if (!is.data.frame(df)) {
    stop("Input is not a dataframe")
  }
  
  if (nrow(df) >= mask_threshold * n_valid_IDs) {
    return(df)
  }
}




subset_2_commonIDs <- function(pRes, base_method, mask_threshold = .9){
  
  if (base_method %in% names(pRes)) {
    dat <- pRes[[base_method]]
    
  } else {
    print(paste0(base_method, ' does not exist, we are using DriverPower for selecting common IDs'))
    dat <- pRes[['DriverPower']]
  }
  
  
  valid_IDs <- dat$PCAWG_IDs
  
  pRes <- lapply(pRes, function(s){
    s[which(s$PCAWG_IDs %in% valid_IDs),]
  })
  reference_length <- nrow(dat)
  n_valid_IDs = ceiling(reference_length*mask_threshold)
  
  n_checked_pRes <- lapply(pRes, check_n_elems_method, mask_threshold, n_valid_IDs)
  
  null_idx <- unlist(lapply(n_checked_pRes, is.null))
  n_checked_pRes <- n_checked_pRes[which(!null_idx)]
  
  common_ids <- Reduce(intersect, lapply(n_checked_pRes, function(s){
    s$PCAWG_IDs
  }))
  
  pRes_subset <- lapply(n_checked_pRes, function(s){
    s[which(s$PCAWG_IDs %in% common_ids),]
  })
  
  pRes_subset
}



# 5) recalculate fdr for all methods
recalculate_fdr <- function(pRes){
  
  lapply(pRes, function(s){
    s$fdr <- p.adjust(s$p_value, method = 'fdr')
    s[order(s$p_value, decreasing = F),]
  })
}

# 6) determine significant elements, gold standards, and define true positives
getSigninficantElemens <- function(dat, sigCriteria) {
  
  if (sigCriteria$method == "fdr") {
    if (nrow(dat) != 0) {
      dat$positive <- dat$fdr < sigCriteria$thr
    }
  } else if (sigCriteria$method == "fixedNumberOfElems") {
    
    if (nrow(dat) != 0) {
      dat$positive <- c(rep(TRUE, sigCriteria$num),
                        rep(FALSE, nrow(dat) - sigCriteria$num))
    }
  }
  return(dat)
}


retrieve_goldstandards <- function(ann_PCAWG_ID, based_on, elemType, tissue){ 
  
  ann_PCAWG_ID <- as.data.frame(ann_PCAWG_ID[which(ann_PCAWG_ID$type_of_element %in% c("enhancers", "gc19_pc.3utr", "gc19_pc.5utr", "gc19_pc.cds", "gc19_pc.promCore")),])
  
  if(based_on == "in_CGC_new"){
    goldStandard <- ann_PCAWG_ID[which(ann_PCAWG_ID$in_CGC_new) ,]$PCAWG_IDs
    
  } else if(based_on == "in_CGC_literature"){
    goldStandard <- ann_PCAWG_ID[which(ann_PCAWG_ID$in_CGC_literature),]$PCAWG_IDs
    
  } else if(based_on == "in_CGC_afterPCAWG"){
    goldStandard <- ann_PCAWG_ID[which(ann_PCAWG_ID$in_CGC_afterPCAWG),]$PCAWG_IDs
    
  } else if(based_on == "in_oncoKB"){
    goldStandard <- ann_PCAWG_ID[which(ann_PCAWG_ID$in_oncoKB),]$PCAWG_IDs
    
  } else if(based_on == "in_pcawg"){
    goldStandard <- ann_PCAWG_ID[which(ann_PCAWG_ID$in_pcawg),]$PCAWG_IDs
    
  } else if(based_on == "in_pcawg_tissue"){
    goldStandard <- ann_PCAWG_ID[which(ann_PCAWG_ID$in_pcawg & grepl(tissue, ann_PCAWG_ID$tissues)),]$PCAWG_IDs
    
  } else if(based_on == "any"){
    goldStandard <- ann_PCAWG_ID[which((ann_PCAWG_ID$in_pcawg) | (ann_PCAWG_ID$in_CGC_literature) | (ann_PCAWG_ID$in_oncoKB) | (ann_PCAWG_ID$in_CGC_new)),]$PCAWG_IDs
    
  }
  
  if (elemType == 'CDS') {
    GS = goldStandard[grep('gc19.pc.cds', goldStandard)]
  } else {
    GS = goldStandard[which(!grepl('gc19.pc.cds', goldStandard))]
  }
  GS
}

define_significant_criteria <- function(method, n){
  
  if (method == "fixedNumberOfElems") {
    sigCriteria = list(num=n, method = method)
    
  } else if (method == "fdr") {
    
    sigCriteria = list(thr = n, method = method)
  }
  
  sigCriteria
}

define_element_type <- function(binID_vector){
  
  s <- strsplit(binID_vector, "[::]")
  GenomicElement <- unlist(lapply(s, function(x){x[1]}))
  GenomicElement
  
}


annotate_pRes <- function(pRes, sigCriteria, ann_PCAWG_ID, based_on, tissue) {
  
  sigHits <- lapply(pRes, getSigninficantElemens, sigCriteria)
  
  elemType <- unique(ifelse(unique(unlist(lapply(pRes, function(s){
    unique(define_element_type(s$PCAWG_IDs))
  }))) == 'gc19_pc.cds', 'CDS', 
  'NC'))
  
  gold_standard <- retrieve_goldstandards(ann_PCAWG_ID, based_on, elemType, tissue)
  print(length(gold_standard))
  
  TPs <- lapply(sigHits, function(s){
    s$is_goldStandard <- s$PCAWG_IDs %in% gold_standard
    s$is_TP <- (s$PCAWG_IDs %in% gold_standard) & s$positive
    s
  })
  
  TPs
}

# 8) calculate AUC, AUPR, and other measures
calculate_F1_score <- function(precision, recall){
  
  (2*precision*recall)/(precision+recall)
}

computeStat_1 <- function(pRes, sigCriteria, ann_PCAWG_ID, based_on, tissue) {
  
  comp_pRes <- annotate_pRes(pRes, sigCriteria, ann_PCAWG_ID, based_on, tissue)
  
  TP_CDS <- lapply(comp_pRes, function(s){
    s[which(s$is_TP),"PCAWG_IDs"]
  })
  
  elemType <- unique(ifelse(unique(unlist(lapply(pRes, function(s){
    unique(define_element_type(s$PCAWG_IDs))
  }))) == 'gc19_pc.cds', 'CDS', 
  'NC'))
  goldStandards <- retrieve_goldstandards(ann_PCAWG_ID, based_on, elemType, tissue)
  
  nTPs <- unlist(lapply(comp_pRes,
                        function(s) {sum(s$is_TP)}))
  
  nHits <- unlist(lapply(comp_pRes, 
                         function(s) {sum(s$positive)}))
  
  Recalls <- nTPs/length(unique(goldStandards))
  
  PRECISIONs <- nTPs/nHits
  
  
  F1 <- calculate_F1_score(PRECISIONs, Recalls)
  
  n_elemnts <- unlist(lapply(comp_pRes,
                             function(s) {nrow(s)}))
  
  res_table <- data.frame(n_elemnts, nTPs, nHits, PRECISIONs,
                          Recalls, F1)
  # rownames(res_table) <- names(comp_pRes)
  
  res_table
}


prepare_AUPR <- function(pRes, sigCriteria, ann_PCAWG_ID, based_on, tissue){
  
  comp_pRes <- annotate_pRes(pRes, sigCriteria, ann_PCAWG_ID, based_on, tissue)
  
  predicted <- lapply(comp_pRes, function(s){ 1- s$p_value})
  predicted <- predicted[which(unlist(lapply(predicted,
                                             function(s) length(s))) != 0)]
  
  label <- lapply(comp_pRes, function(s){ s$is_goldStandard})
  label <- label[which(unlist(lapply(label,
                                     function(s) length(s))) != 0)]
  
  list("predicted" = predicted, "label" = label)
  
}  

computeStat_2 <- function(pRes, sigCriteria, ann_PCAWG_ID, based_on, tissue) {
  
  AU_required_list <- prepare_AUPR(pRes, sigCriteria, ann_PCAWG_ID, based_on, tissue)
  
  AUC = AUPR = c()
  PR_curve = list()
  for (i in 1:length(AU_required_list$predicted)) {
    
    tmpAUC <- auc(factor(AU_required_list$label[[i]], levels = c(TRUE, FALSE)), 
                  AU_required_list$predicted[[i]])
    tmpAUPR <- pr.curve(scores.class0 = AU_required_list$predicted[[i]], 
                        weights.class0 = AU_required_list$label[[i]], curve=TRUE)
    
    AUC <- c(AUC, tmpAUC)
    AUPR <- c(AUPR, tmpAUPR$auc.integral)
    PR_curve <- append(PR_curve, list(tmpAUPR$curve))
  }
  
  
  AU <- data.frame(AUC, AUPR)
  AU$method <- names(AU_required_list$label)
  
  list("table" = AU, "curve" = PR_curve)
}

table_measures <- function(preprocessed_pRes, sigCriteria, ann_PCAWG_ID, based_on, tissue){
  
  stat2 = computeStat_2(preprocessed_pRes, sigCriteria, ann_PCAWG_ID, based_on, tissue)
  stat1 = computeStat_1(preprocessed_pRes, sigCriteria, ann_PCAWG_ID, based_on, tissue)
  stat1$method <- rownames(stat1)
  
  result <- left_join(stat1, stat2$table, by = 'method')
  # result <- result[order(result$AUPR, decreasing = T),]
  result <- result[,c("method", "nTPs", "nHits", "PRECISIONs", "Recalls", "F1", "AUC", "AUPR", "n_elemnts")]
  
  list('table' = result, 'curve' = stat2$curve)
}

compare_methods <- function(pRes, df, method_name, ann_PCAWG_ID, based_on, tissue){
  
  nHits_row <- c()
  for(k in 1:nrow(df)){
    sigCr <- define_significant_criteria("fixedNumberOfElems", df$nHits[k])
    st1_nTP_desired <- computeStat_1(pRes[method_name], sigCr, ann_PCAWG_ID, based_on, tissue)
    nHits_row <- c(nHits_row, st1_nTP_desired$nTPs)
    
  }
  sigCr2 <- define_significant_criteria("fixedNumberOfElems",
                                        df[which(df$method==method_name), which(colnames(df)=="nHits")]
  )
  st_nTP_nHitsDesired <- computeStat_1(pRes, sigCr2, ann_PCAWG_ID, based_on, tissue)
  nHits_desired <- st_nTP_nHitsDesired$nTPs
  
  new_cols <- cbind(nHits_row, nHits_desired)
  
  df <- cbind(df, new_cols)
  col_order <- c("method","nTPs","nHits", "nHits_row", "nHits_desired",
                 "PRECISIONs", "Recalls", "F1", "AUC", "AUPR", "n_elemnts" )
  df[,col_order]
}

preprocess_pRes <- function(path_pRes, PATHs_newResults, newRESULTS,
                            element, base_method, selected_methods = NULL){
  load(path_pRes)
  pRes <- pRes[which(!names(pRes) %in% c('compositeDriver'))]
  
  print('pRes loaded')
  
  
  # if (!is.na(selected_methods)){
  #   if(!is.null(selected_methods)){
  #     pRes <- pRes[which(names(pRes) %in% selected_methods)]
  #   } else {
      idx <- names(pRes)
  #   }
  # 
  # }
  
  print('idx...')
  
  
  if (length(PATHs_newResults) != 0){
    pRes <- add_newMethod_to_pRes(pRes, PATHs_newResults, newRESULTS) 
  }
  
  print('add new method')
  
  # if(is.null(selected_methods)){
  #   pRes <- pRes[which(!names(pRes) %in% idx)]
  # }
  
  
  preprocessed_pRes <- prepare_pRes_helper(pRes, element, base_method)
  preprocessed_pRes
}

prepare_pRes_helper <- function(pRes, element, base_method = 'iDriver', mask_threshold = .9){
  x = split_pRes(pRes)
  print('**** 1 ******')
  if (element == 'CDS') {
    preprocessed_pRes <- subset_2_commonIDs(x$CDS_pRes, base_method, mask_threshold)
  } else if (element == 'NC') {
    preprocessed_pRes <- subset_2_commonIDs(x$nc_pRes, base_method, mask_threshold)
  }
  print('**** 2 ******')
  preprocessed_pRes <- recalculate_fdr(preprocessed_pRes)
  preprocessed_pRes
}

save_Measures <- function(path_save, newRESULTS, PATHs_newResults, path_pRes, 
                          sig_definition_method,
                          sig_definition_threshold,
                          path_ann_PCAWG_ID, tissue, based_on,
                          element, base_method, selected_methods = NULL, compareRank_new = list(FALSE, NA), tissue_goldStd_pcawg = NULL){
  ann_PCAWG_ID <- fread(path_ann_PCAWG_ID)
  sigCriteria <- define_significant_criteria(sig_definition_method, sig_definition_threshold)
  preprocessed_pRes <- preprocess_pRes(path_pRes, PATHs_newResults, newRESULTS,
                                       element, base_method, selected_methods)
    
  Measures <- table_measures(preprocessed_pRes, sigCriteria, ann_PCAWG_ID, based_on,  tissue = tissue_goldStd_pcawg)
  
  tables <- Measures$table
  # curve <- Measures$curve
  
  if (compareRank_new[[1]]) {
    tables <- compare_methods(preprocessed_pRes, tables, compareRank_new[[2]], ann_PCAWG_ID, based_on, tissue = tissue_goldStd_pcawg)
  }
  
  tables <- tables[order(tables$AUPR, decreasing = T),]
  
  dir.create(paste0(path_save, "/tables/"), 
             showWarnings = F,
             recursive = T)
  write.csv(tables, file = paste0(path_save, 
                                  "/tables/", tissue, "_table_GoldStd_basedon_",
                                  based_on, "_", element, "_",
                                  sigCriteria$method, ".csv"))
  
  # save(curve, file = paste0(path_save, tissue, 
  #                           "/tables/Measures_GoldStd_basedon_", based_on,
  #                           "_", element, "_",sigCriteria$method, ".RData"))
  
}

create_pRes_MyRes <- function(PATH_newResults, newRES){
  
  res <- lapply(PATH_newResults, fread)
  names(res) <- newRES
  class(res) = "pRes"
  
  res <- lapply(res, function(s){
    s <- data.frame(s)
    s <- s[,which(colnames(s) %in% c("PCAWG_ID", "pvals"))]
    colnames(s) = c("PCAWG_IDs",  "p_value")
    s
  })
  
  
  res
}

save_Measures2 <- function(path_save, newRESULTS, PATHs_newResults, 
                          sig_definition_method,
                          sig_definition_threshold,
                          path_ann_PCAWG_ID, tissue, based_on,
                          element, base_method, selected_methods = NULL, 
                          compareRank_new = list(FALSE, NA), tissue_goldStd_pcawg = NULL){
  ann_PCAWG_ID <- fread(path_ann_PCAWG_ID)
  sigCriteria <- define_significant_criteria(sig_definition_method, sig_definition_threshold)
  
  pRes <- create_pRes_MyRes(PATHs_newResults, newRESULTS)
  
  preprocessed_pRes <- prepare_pRes_helper(pRes, element)
  
  
  Measures <- table_measures(preprocessed_pRes, sigCriteria, ann_PCAWG_ID, based_on,  tissue = tissue_goldStd_pcawg)
  
  tables <- Measures$table
  # curve <- Measures$curve
  
  if (compareRank_new[[1]]) {
    tables <- compare_methods(preprocessed_pRes, tables, compareRank_new[[2]], ann_PCAWG_ID, based_on, tissue = tissue_goldStd_pcawg)
  }
  
  tables <- tables[order(tables$AUPR, decreasing = T),]
  
  dir.create(paste0(path_save, "/tables/"), 
             showWarnings = F,
             recursive = T)
  write.csv(tables, file = paste0(path_save, 
                                  "/tables/", tissue, "_table_GoldStd_basedon_",
                                  based_on, "_", element, "_",
                                  sigCriteria$method, ".csv"))
  
  # save(curve, file = paste0(path_save, tissue, 
  #                           "/tables/Measures_GoldStd_basedon_", based_on,
  #                           "_", element, "_",sigCriteria$method, ".RData"))
  
}
prepare_data_kTop <- function(path_save_tables_rmHypers, k = 5){
  nonHyperCancers <- c(
    "Biliary-AdenoCA",  "Bladder-TCC",# "Bone-Leiomyo",
    "Bone-Osteosarc", "Pancan-no-skin-melanoma-lymph",
    "Breast-AdenoCa",   "CNS-GBM",   "CNS-Medullo",
    #            "CNS-Oligo",   "Cervix-SCC", 
    "CNS-PiloAstro",   "ColoRect-AdenoCA",
    "Eso-AdenoCa",  "Head-SCC" ,  "Kidney-RCC",   "Pan_Cancer",
    "Liver-HCC", "Lung-AdenoCA", "Lung-SCC", "Lymph-BNHL",
    "Lymph-CLL",  "Myeloid-MPN", "Ovary-AdenoCA",
    "Panc-AdenoCA",     "Panc-Endocrine",   "Prost-AdenoCA",
    "Skin-Melanoma",
    "Stomach-AdenoCA",  "Thy-AdenoCA", "Uterus-AdenoCA"
  ) #  for "Cervix-SCC" there is just one method to compare with 
  
  # nonHyperCancers <- c("Pan_Cancer", "Pancan-no-skin-melanoma-lymph" , "ColoRect-AdenoCA","CNS-GBM", 
  #                      "Head-SCC", "Skin-Melanoma", "Lung-AdenoCA", "Lung-SCC",
  #                      "Stomach-AdenoCA", "Uterus-AdenoCA" ,  "Biliary-AdenoCA",
  #                      "Eso-AdenoCa", "Lymph-BNHL")
  
  final_df <- c()
  for (k in 1:5) {
    
    k_df_elem <- c()
    for(element_type in c('CDS','NC')){
      
      LIST = create_methods_tables_allCohorts(hyperCancers = c(), nonHyperCancers, 
                                              path_save_tables_withHypers = NA, 
                                              path_save_tables_rmHypers, measure = 'AUPR', k, 
                                              element_type)
      k_df <- LIST$Table
      k_df$K_top <- k
      k_df$elem <- element_type
      k_df_elem <- rbind(k_df_elem, k_df)
    }
    final_df <- rbind(final_df, k_df_elem)
  }
  final_df
}

# Define the split function
split_cod_nonCod_saveSigs <- function(df){
  
  df <- df[order(df$pvals, decreasing = F), ]
  
  # Ensure PCAWG_ID is a character column before applying grep
  df$PCAWG_ID <- as.character(df$PCAWG_ID)
  
  # Split into coding and non-coding
  df_cds <- df[grep("gc19_pc.cds", df$PCAWG_ID), ]
  
  # Define pattern for non-coding elements
  pattern <- paste(c("enhancers", "gc19_pc.3utr", "gc19_pc.5utr", "gc19_pc.promCore"), collapse = "|")
  df_nc <- df[grep(pattern, df$PCAWG_ID), ]
  
  # Return as a list
  return(list("CDS" = df_cds, "non_coding" = df_nc))
}




get_sigHits_cohort <- function(path_file, sig_threshold = .1){
  
  df <- fread(path_file)  # fread directly returns a data.table
  
  # Filter rows where sum_nMut > 1
  df <- df[df$sum_nMut > 1, ]
  
  # Split into CDS and non-coding
  x <- split_cod_nonCod_saveSigs(df)
  
  # Check results
  head(x$CDS)  # Check CDS subset
  head(x$non_coding)  # Check non-coding subset
  
  x <- lapply(x, function(s){
    s$qvals <- p.adjust(s$pvals, method = 'fdr')
    s
  })
  
  
  
  CDS <- x$CDS
  hits_cds <- CDS[which(CDS$qvals < sig_threshold), ]
  
  NC <- x$non_coding
  hits_nc <- NC[which(NC$qvals < sig_threshold), ]
  
  
  rbind(hits_nc, hits_cds)
}
