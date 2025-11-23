rm(list = ls())

library(data.table)
library(openxlsx)
library(dplyr)
library(pROC)
library(PRROC)

source("benchmark/functions_annotate_elements.R")

path_procPCAWG_res <- "../extdata/procInput/PCAWG_results/"
path_proccessedGE <- "../extdata/procInput/"
tissue <-  "Pancan-no-skin-melanoma-lymph"
path_save_benchRes <- "../extdata/output/benchmark/iDriver_LmntSp/perLmnt/"

########### 1) load the pRes object and annotated PCAWG_IDs ###########
ann_PCAWG_ID_complement <- fread(paste0(path_proccessedGE, "ann_PCAWG_ID_complement.csv"))
load(paste0(path_procPCAWG_res, tissue, ".RData"))
# pRes <- pRes["DriverPower"]
#################################################
add_newMethod_to_pRes <- function(pRes, path_newRes){
  
  if(is.data.frame(path_newRes)){
    newMethod <- path_newRes
  } else if (is.character(path_newRes)) {
    newMethod <- data.frame(fread(path_newRes))
  }
  
  idx_ID <- which(colnames(newMethod) %in% c("binID", "PCAWG_ID", 
                                             "element_ID", "id"))
  idx_pVal <- which(colnames(newMethod) %in% c("pValues", "p_value",
                                               "raw_p", "pvals",
                                               "p.value"))
  
  newMethod <- data.frame("PCAWG_IDs" = newMethod[,idx_ID], 
                          "p_value" =  newMethod[,idx_pVal],
                          "q_value" = rep(NA, nrow(newMethod)))
  
  # all_testGE_IDs <- fread("../extdata/rawInput/PCAWG_test_genomic_elements.bed12.gz")$V4
  # removed_IDs <- all_testGE_IDs[which(!all_testGE_IDs %in% newMethod$PCAWG_IDs)]
  # removed_df <- data.frame("PCAWG_IDs" = removed_IDs, 
  #                          "p_value" =  rep(1, length(removed_IDs)),
  #                          "q_value" = rep(NA, length(removed_IDs)))
  # newMethod <- rbind(newMethod, removed_df)
  newMethod$p_value <- ifelse(is.na(newMethod$p_value), 1, newMethod$p_value)
  newMethod$fdr =  p.adjust(newMethod$p_value, "fdr")
  
  pRes$newMethod <- newMethod
  class(pRes) = "pRes"
  pRes
}


newRESULTS <- c("iDriver", "sampling", "sampling_raw", "sampling_impo", "LmnSp_raw_impo")
PATHs_newResults <- c("../extdata/output/individual_level/withCADD/Pan_Cancer/all_indvLvl_GBM_nBinom.csv",
                      "../extdata/output/individual_level/tmp_sampling/Pan_Cancer/CLT_corrected_all_results.csv",
                      "../extdata/output/individual_level/tmp_sampling/Pan_Cancer/raw/gc19_pc.cdstmp_test_iDriverSampling.csv",
                      "../extdata/output/individual_level/tmp_sampling/Pan_Cancer/gc19_pc.cdstmp_test_iDriverSampling.csv",
                      "../extdata/output/individual_level/iDriverImpoSamp_rawScores_elemSp/raw_elemSp_all_results.csv")


for (i in 1:length(PATHs_newResults)) {
  pRes <- add_newMethod_to_pRes(pRes, PATHs_newResults[i])
  names(pRes) <- c(names(pRes)[1:(length(pRes) -1)], newRESULTS[i])
} 

pRes_annotated <- lapply(pRes, function(s){
  left_join(s, ann_PCAWG_ID_complement)
})
##############################
split_df_LmntTypes <- function(df) {
  
  Lmnts <- c("enhancers", "gc19_pc.3utr", "gc19_pc.5utr", "gc19_pc.cds",
             "gc19_pc.promCore", "gc19_pc.ss")
  # , "lncrna.ncrna", "lncrna.promCore", "lncrna.ss")
  
  df_list <- list()
  for(lmnt in Lmnts){
    df_sub <- df[grep(lmnt, df$PCAWG_IDs),]
    df_list <- append(df_list, list(df_sub))
  }
  names(df_list) <- Lmnts
  df_list
}


pRes2 <- lapply(pRes_annotated, split_df_LmntTypes) 

##############
define_significant_criteria <- function(method, n){
  
  if (method == "fixedNumberOfElems") {
    sigCriteria = list(num=n, method = method)
  } else if (method == "fdr") {
    sigCriteria = list(thr = n, method = method)
  }
  sigCriteria
}
sig_definition_methods <- c("fdr", "fixedNumberOfElems")
sig_definition_thresholds <- c(.1, 200)

j=1
sigCriteria <- define_significant_criteria(sig_definition_methods[j],
                                           sig_definition_thresholds[j])


getSigninficantElemens <- function(df, sigCriteria) {
  
  df <- df[order(df$p_value, decreasing = F),]
  
  if (sigCriteria$method == "fdr") {
    if (nrow(df) != 0) {
      df$positive <- df$fdr < sigCriteria$thr
    }
  } else if (sigCriteria$method == "fixedNumberOfElems") {
    if (nrow(df) != 0) {
      df$positive <- c(rep(TRUE, sigCriteria$num),
                       rep(FALSE, nrow(df) - sigCriteria$num))
    }
  }
  df
}


#############
retrieve_goldstandards <- function(ann_PCAWG_ID, based_on){ #based_on:
  # can be: in_CGC_literature, in_CGC_new, in_oncoKB or in_pcawg or one_of_all
  
  lmnts_except_lncRNAs <- c("enhancers", "gc19_pc.3utr", "gc19_pc.5utr", 
                            "gc19_pc.promCore", "gc19_pc.ss", "gc19_pc.cds")
  if(based_on == "in_CGC_new"){
    goldStandard1 <- ann_PCAWG_ID[which((ann_PCAWG_ID$in_CGC_new) &
                                          (ann_PCAWG_ID$type_of_element %in%
                                             lmnts_except_lncRNAs)),]$PCAWG_IDs
    
  } else if(based_on == "in_CGC_literature"){
    goldStandard1 <- ann_PCAWG_ID[which((ann_PCAWG_ID$in_CGC_literature) &
                                          (ann_PCAWG_ID$type_of_element %in%
                                             lmnts_except_lncRNAs)),]$PCAWG_IDs
  } else if(based_on == "in_pcawg"){
    goldStandard1 <- ann_PCAWG_ID[which((ann_PCAWG_ID$in_pcawg) &
                                          (ann_PCAWG_ID$type_of_element %in%
                                             lmnts_except_lncRNAs)),]$PCAWG_IDs
  } else if(based_on == "in_oncoKB"){
    goldStandard1 <- ann_PCAWG_ID[which((ann_PCAWG_ID$in_oncoKB) &
                                          (ann_PCAWG_ID$type_of_element %in%
                                             lmnts_except_lncRNAs)),]$PCAWG_IDs
  } else if(based_on == "all_Databases"){
    goldStandard1 <- ann_PCAWG_ID[which((ann_PCAWG_ID$in_pcawg | 
                                           ann_PCAWG_ID$in_CGC_literature |
                                           ann_PCAWG_ID$in_CGC_new |
                                           ann_PCAWG_ID$in_oncoKB) &
                                          (ann_PCAWG_ID$type_of_element %in%
                                             lmnts_except_lncRNAs)),]$PCAWG_IDs
  }
  # goldStandard2 <- ann_PCAWG_ID[which(((ann_PCAWG_ID$in_CLC) &
  #                                          (ann_PCAWG_ID$type_of_element %in% 
  #                                             c("lncrna.ncrna","lncrna.promCore")))),]$PCAWG_IDs
  goldStandard2 <- ""
  c(goldStandard1, goldStandard2)
}

ann_PCAWG_ID <- ann_PCAWG_ID_complement
based_on = "all_Databases"
annotate_pRes <- function(ann_PCAWG_ID, based_on, pRes2){
  GoldStandard_set <- retrieve_goldstandards(ann_PCAWG_ID, based_on)
  sigHits <- lapply(pRes2, function(s){
    lapply(s, function(x){
      x1 = getSigninficantElemens(x, sigCriteria)
      x1$is_goldStandard <- x1$PCAWG_IDs %in% GoldStandard_set
      x1$is_TP <- (x1$PCAWG_IDs %in% GoldStandard_set) & (x1$positive)
      x1
    })
  })
}

ann_pRes <- annotate_pRes(ann_PCAWG_ID, based_on, pRes2)

calculate_F1_score <- function(precision, recall){
  
  (2*precision*recall)/(precision+recall)
}

retrieve_nTP_nHits <- function(df, GoldStandard_set){
  nTPs <- sum(df$is_TP)
  nHits <- sum(df$positive)
  Lmnt_type <- na.omit(unique(df$type_of_element))
  RECALL <- nTPs/length(grep(Lmnt_type, GoldStandard_set))
  PRECISION <- nTPs/nHits
  n_elemnts <- nrow(df)
  F1 <- calculate_F1_score(PRECISION, RECALL)
  Stat_1 <- c(nTPs, nHits, Lmnt_type, RECALL, PRECISION, F1, n_elemnts)
  names(Stat_1) <- c("nTPs", "nHits", "Lmnt_type", "RECALL", "PRECISION",
                     "F1", "n_elemnts")
  Stat_1
}
computeStat_1 <- function(ann_pRes, ann_PCAWG_ID, based_on){
  GoldStandard_set <- retrieve_goldstandards(ann_PCAWG_ID, based_on)
  Stat_1s <- lapply(ann_pRes, function(s){
    lapply(s, function(x){
      if(nrow(x) != 0){
        retrieve_nTP_nHits(x, GoldStandard_set)
      }
    })
  })
  
  lapply(Stat_1s, function(s){
    data.frame(do.call(rbind, s))
  })
  
}

Stat_1 <- computeStat_1(ann_pRes, ann_PCAWG_ID, "in_CGC_new")
#############
prepare_AUPR_df <- function(df){
  predicted = 1-df$p_value
  label = df$is_goldStandard
  
  list("label" = label, "predicted" = predicted)
}
prepare_AUPR <- function(ann_pRes){
  lapply(ann_pRes, function(s){
    lapply(s, prepare_AUPR_df)
  })
}
AU_required_list=prepare_AUPR(ann_pRes)

stat2_curve_table <- lapply(AU_required_list, function(s){
  lapply(s, function(x){
    if(length(x$predicted) != 0){
      prd=x$predicted
      lbl=x$label
      AUC = auc(factor(lbl), prd)
      AUPRs = pr.curve(scores.class0 = x$predicted, weights.class0 = x$label,
                       curve=TRUE)
      AUPR = AUPRs$auc.integral
      PR_curve <- AUPRs$curve
      
      list("table" = data.frame("AUC" = AUC, "AUPR" = AUPR),
           "curve" = PR_curve)
    }
    
  })
})

extract_Stat_2_tables <- function(stat2_curve_table){
  tables <- lapply(stat2_curve_table, function(s){
    s2 = lapply(s, function(x){
      x$table
    })
    y = data.frame(do.call(rbind, s2))
    y$Lmnt_type <- rownames(y)
    y
  })
}

Stat_2 <- extract_Stat_2_tables(stat2_curve_table)

pres_methods <- names(Stat_1)

measure_list <- list()
for (m in pres_methods) {
  measures <- left_join(Stat_1[[m]], Stat_2[[m]], by = "Lmnt_type")
  measures$method <- rep(m, nrow(measures))
  measure_list <- append(measure_list, list(measures))
}

all_measures <- do.call(rbind, measure_list)
all_measures <- all_measures %>% group_by(Lmnt_type) %>% 
  arrange(desc(AUPR),.by_group = TRUE)

save_name <- "LmntSp_iDriver_allDataBases"
DIR <- paste0(path_save_benchRes, tissue, "/",save_name, "/tables/")
dir.create(DIR, showWarnings = F, recursive = T)
fwrite(all_measures, file = paste0(DIR, sigCriteria$method,"_", 
                                   sigCriteria$thr, ".csv"), sep = ",")



curves <- lapply(stat2_curve_table, function(s){
  s2 = lapply(s, function(x){
    x$curve
  })
})
########################