rm(list = ls())

library(data.table)
library(dplyr)


prepare_predRate_df <- function(path_predRate, path_test_Y){
  predRate <- fread(path_predRate)
  test_y <- fread(path_test_Y)
  df <- left_join(test_y, predRate, by = 'binID')
  predCol <- which(colnames(df) %in% c('predRate', 'pred_rates'))
  colnames(df)[predCol] <- 'predRate'
  N = df$N
  len = df$length
  
  predRate = df$predRate
  df$raw_nPred = predRate * len * N
  # df$obsRate = df$nMut/(df$length * df$N)
  # df$nPred <- pmin(df$nMut, df$raw_nPred)
  df$nPred <- df$raw_nPred
  df
}

# select_cohort <- function(path_donorInfo, 
#                           cohort, exclude_lymph_melanoma = TRUE,
#                           exclude_hyper_mutated = TRUE,
#                           only_hyperMut_pancancer = FALSE){
#   
#   donorInfo <- fread(path_donorInfo)
#   
#   if (only_hyperMut_pancancer) {
#     donorInfo <- donorInfo[which(donorInfo$HyperMut_donor),]
#   } else {
#     if (exclude_lymph_melanoma) {
#       exceptions <- c("Skin-Melanoma", "SKCM-US",
#                       "Lymph-NOS", 
#                       "Lymph-CLL", "CLLE-ES",
#                       "Lymph-BNHL", "MALY-DE", "DLBC-US")
#       donorInfo <- donorInfo[-which(donorInfo$cohort1 %in% exceptions),] # 2279 donors
#     }
#     
#     if (exclude_hyper_mutated) {
#       donorInfo <- donorInfo[which(!donorInfo$HyperMut_donor),] 
#     }
#     
#     if (!cohort %in% c('Pan_Cancer', 'Pancan-no-skin-melanoma-lymph')) {
#       donorInfo <- donorInfo[which(donorInfo$cohort1 == cohort),]
#     } 
#   }
#   
#   
#   
#   if ('sampleID' %in% colnames(donorInfo)) {
#     donorInfo <- donorInfo[,c('sampleID', "D_id","freq" )]
#     colnames(donorInfo) <- c('sampleID', 'donor_id', 'totalMutNrs')
#   } else {
#     donorInfo <- donorInfo[,c("D_id","freq" )]
#     colnames(donorInfo) <- c('donor_id', 'totalMutNrs')
#   }
#   
#   donorInfo
#   
#   
#   
# }
select_cohort <- function(path_donorInfo, 
                          cohort, exclude_lymph_melanoma = TRUE,
                          exclude_hyper_mutated = TRUE,
                          only_hyperMut_pancancer = FALSE){
  
  donorInfo <- fread(path_donorInfo)
  
  if (only_hyperMut_pancancer) {
    donorInfo <- donorInfo[which(donorInfo$HyperMut_donor),]
  } else {
    if (exclude_lymph_melanoma) {
      exceptions <- c("Skin-Melanoma", "SKCM-US",
                      "Lymph-NOS", 
                      "Lymph-CLL", "CLLE-ES",
                      "Lymph-BNHL", "MALY-DE", "DLBC-US")
      donorInfo <- donorInfo[-which(donorInfo$cohort1 %in% exceptions),] # 2279 donors
    }
    
    if (exclude_hyper_mutated) {
      donorInfo <- donorInfo[which(!donorInfo$HyperMut_donor),] 
    }
    
    if (!cohort %in% c('Pan_Cancer', 'Pancan-no-skin-melanoma-lymph')) {
      if (grepl('subsample', cohort)) {
        donorInfo <- donorInfo[which(donorInfo$subsample_id == cohort),]
      } else {
        donorInfo <- donorInfo[which(donorInfo$cohort1 == cohort),]
      }
      
    } 
  }
  
  
  
  if ('sampleID' %in% colnames(donorInfo)) {
    donorInfo <- donorInfo[,c('sampleID', "D_id","freq" )]
    colnames(donorInfo) <- c('sampleID', 'donor_id', 'totalMutNrs')
  } else {
    donorInfo <- donorInfo[,c("D_id","freq" )]
    colnames(donorInfo) <- c('donor_id', 'totalMutNrs')
  }
  
  donorInfo
  
  
  
}





create_pElemDF <- function(path_donorInfo, cohort, path_predRate, path_test_Y,
                           exclude_lymph_melanoma,  exclude_hyper_mutated, save_name){
  
  if (cohort == 'onlyHyperMuts') {
    donorInfo <- select_cohort(path_donorInfo = path_donorInfo,
                               only_hyperMut_pancancer = TRUE)
  } else {
    donorInfo <- select_cohort(path_donorInfo,
                               cohort,
                               exclude_lymph_melanoma,
                               exclude_hyper_mutated)
  }
  
  
  # create the mutData long format for the cohort of interest and scores of interest
  mutData <- fread(path_MutData)
  mutData <- left_join(donorInfo, mutData)
  
  # Add pElement to mutData
  predRates <- prepare_predRate_df(path_predRate, path_test_Y)
  mutData <- data.frame(mutData)
  df <- mutData[,c('donor_id', 'totalMutNrs')]
  df <- df[!duplicated(df),]
  
  all_cohort_mutsNr <- sum(df$totalMutNrs)
  predRates$pElement <- predRates$nPred/all_cohort_mutsNr
  
  mutData <- left_join(mutData, predRates, by = c('PCAWG_ID' = 'binID'))
  
  pElem_df <- mutData[,c("PCAWG_ID", "nMut", "length", "N", "predRate", "nPred", "pElement")]
  pElem_df <- pElem_df[!duplicated(pElem_df),]
  
  dir_path <- dirname(save_name)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  fwrite(pElem_df, file = save_name, sep = '\t')
}


# ##################### RUN ##########################################
# cohorts <- c('Pan_Cancer',"Liver-HCC",  'Pancan-no-skin-melanoma-lymph', 
#              "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
#              "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
#              "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
#              "Breast-AdenoCa","Biliary-AdenoCA", "Eso-AdenoCa", "CNS-GBM",         
#              "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
#              "Breast-LobularCa", "Thy-AdenoCA", "Myeloid-MPN", "Bone-Leiomyo",    
#              "Lymph-NOS", "CNS-Medullo", "Myeloid-AML", "CNS-Oligo", "Cervix-SCC",
#              "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Epith", "Bone-Osteosarc",
#              "Cervix-AdenoCA", "Breast-DCIS", "Bone-Cart", "Myeloid-MDS" )
# 
# path_MutData <- '../extdata/procInput/iDriverInputs/mutData.tsv'
# path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
# exclude_hyper_mutated = TRUE
# 
# 
# for (cohort in cohorts) {
#   print(cohort)
#   
#   path_test_Y <- paste0('../extdata/procInput/BMRs_2024/observed/', cohort, '/test_y.tsv')
#   
#   exclude_lymph_melanoma = ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
#   
#   
#   for (model in c('eMET')) { # 'GBM', 
#     if (model == 'GBM') {
#       path_predRate <- paste0('../extdata/procInput/BMRs_2024/observed/', cohort, '/GBM/GBM_predTest.tsv')
#       save_name <- paste0('../extdata/procInput/BMRs_2024/observed/', cohort, '/pElems/GBM_orig_pElmens.tsv')
#     } else if (model == 'eMET') {
#       path_predRate <- paste0('../extdata/procInput/BMRs_2024/observed/', cohort, '/eMET/eMET_100_predTest.tsv')
#       save_name = paste0('../extdata/procInput/BMRs_2024/observed/', cohort, '/pElems/eMET_orig_pElmens.tsv')
#     }
#     
#     create_pElemDF(path_donorInfo, cohort, path_predRate, path_test_Y,
#                    exclude_lymph_melanoma,  exclude_hyper_mutated, save_name)
#   }
#   print('==========================================================')
# }

# ##################### with hyperMutated samples ##################### 
# cohorts=c('Pan_Cancer', 'Pancan-no-skin-melanoma-lymph'
#   # 'Skin-Melanoma',  'Stomach-AdenoCA'  , 'Uterus-AdenoCA', 'Lymph-BNHL' ,
#   # 'ColoRect-AdenoCA',  'Eso-AdenoCa', 'Head-SCC',     'Lung-AdenoCA', 'Lung-SCC',
#   # 'Biliary-AdenoCA',  'CNS-GBM'
# )
# 
# 
# path_MutData <- '../extdata/procInput/iDriverInputs/mutData.tsv'
# path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
# exclude_hyper_mutated = FALSE
# 
# 
# for (cohort in cohorts) {
#   print(cohort)
#   
#   path_test_Y <- paste0('../extdata/procInput/BMRs_2024/observed_withHyperMutated/', cohort, '/test_y.tsv')
#   
#   exclude_lymph_melanoma = ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
#   
#   
#   for (model in c('eMET', 'GBM')) { # 'GBM', 
#     if (model == 'GBM') {
#       path_predRate <- paste0('../extdata/procInput/BMRs_2024/observed_withHyperMutated/', cohort, '/GBM/GBM_predTest.tsv')
#       save_name <- paste0('../extdata/procInput/BMRs_2024/observed_withHyperMutated/', cohort, '/pElems/GBM_orig_pElmens.tsv')
#     } else if (model == 'eMET') {
#       path_predRate <- paste0('../extdata/procInput/BMRs_2024/observed_withHyperMutated/', cohort, '/eMET/eMET_100_predTest.tsv')
#       save_name = paste0('../extdata/procInput/BMRs_2024/observed_withHyperMutated/', cohort, '/pElems/eMET_orig_pElmens.tsv')
#     }
#     
#     create_pElemDF(path_donorInfo, cohort, path_predRate, path_test_Y,
#                    exclude_lymph_melanoma,  exclude_hyper_mutated, save_name)
#   }
#   print('==========================================================')
# }
# 
# 
# 
# #####################(sim) pan-cancer  ##########################################
# # path_predRate <- '../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/Pancan-no-skin-melanoma-lymph/eMET/eMET_100_predTest.tsv'
# # path_test_Y <- '../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/Pancan-no-skin-melanoma-lymph/test_y.tsv'
# path_MutData <- '../extdata/procInput/iDriverInputs/simulated/Sanger/mutData.tsv'
# path_donorInfo <- '../extdata/procInput/iDriverInputs/simulated/Sanger/donorInfo.tsv'
# exclude_hyper_mutated = TRUE
# 
# # exclude_lymph_melanoma = TRUE
# 
# # cohort <- 'Pancan-no-skin-melanoma-lymph'
# # path_save <- paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/', cohort, '/pElems/eMET_orig_pElmens.tsv')
# 
# # create_pElemDF(path_donorInfo, cohort, path_predRate, path_test_Y,
# #                exclude_lymph_melanoma,  exclude_hyper_mutated, path_save)
# 
# cohorts <- c(
#   'Pan_Cancer',"Liver-HCC",  'Pancan-no-skin-melanoma-lymph',
#              "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
#              "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
#              "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
#              "Breast-AdenoCa","Biliary-AdenoCA", "Eso-AdenoCa", "CNS-GBM",
#              "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
#              "Breast-LobularCa", "Thy-AdenoCA", "Myeloid-MPN",
#   # "Bone-Leiomyo",  "Lymph-NOS","Bone-Cart",   "Myeloid-MDS"
#   "CNS-Medullo", "Myeloid-AML", "CNS-Oligo", "Cervix-SCC",
#              "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Epith", "Bone-Osteosarc",
#              "Cervix-AdenoCA", "Breast-DCIS" )
# 
# 
# for (cohort in cohorts) {
#   print(cohort)
#   
#   path_test_Y <- paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/', cohort, '/test_y.tsv')
#   
#   exclude_lymph_melanoma = ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
#   
#   
#   for (model in c('eMET', 'GBM')) { 
#     if (model == 'GBM') {
#       
#       path_predRate <- paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/', cohort, '/GBM/GBM_predTest.tsv')
#       save_name <- paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/', cohort, '/pElems/GBM_orig_pElmens.tsv')
#     } else if (model == 'eMET') {
#       path_predRate <- paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/', cohort, '/eMET/eMET_100_predTest.tsv')
#       save_name = paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/', cohort, '/pElems/eMET_orig_pElmens.tsv')
#     }
#     
#     create_pElemDF(path_donorInfo, cohort, path_predRate, path_test_Y,
#                    exclude_lymph_melanoma,  exclude_hyper_mutated, save_name)
#     
#   }
#   print('==========================================================')
# }
# 
# #####################(sim) with hyperMuts pan-cancer  ##########################################
# # path_predRate <- '../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/Pancan-no-skin-melanoma-lymph/eMET/eMET_100_predTest.tsv'
# # path_test_Y <- '../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/Pancan-no-skin-melanoma-lymph/test_y.tsv'
# path_MutData <- '../extdata/procInput/iDriverInputs/simulated/Sanger/mutData.tsv'
# path_donorInfo <- '../extdata/procInput/iDriverInputs/simulated/Sanger/donorInfo.tsv'
# exclude_hyper_mutated = FALSE
# 
# 
# cohorts <- c(
#   # 'Pan_Cancer', 'Pancan-no-skin-melanoma-lymph',
#   "ColoRect-AdenoCA", "CNS-GBM",
#   "Head-SCC", "Skin-Melanoma", "Lung-AdenoCA",  "Lung-SCC",
#   "Stomach-AdenoCA",  "Uterus-AdenoCA", "Biliary-AdenoCA",
#     "Lymph-BNHL" # "Eso-AdenoCa",
#   )
# 
# 
# for (cohort in cohorts) {
#   print(cohort)
#   
#   path_test_Y <- paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_withHyperMuts/', cohort, '/test_y.tsv')
#   
#   exclude_lymph_melanoma = ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
#   
#   
#   for (model in c('eMET', 'GBM')) { 
#     if (model == 'GBM') {
#       
#       path_predRate <- paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_withHyperMuts/', cohort, '/GBM/GBM_predTest.tsv')
#       save_name <- paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_withHyperMuts/', cohort, '/pElems/GBM_orig_pElmens.tsv')
#     } else if (model == 'eMET') {
#       path_predRate <- paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_withHyperMuts/', cohort, '/eMET/eMET_100_predTest.tsv')
#       save_name = paste0('../extdata/procInput/BMRs_2024/simulated/Sanger_withHyperMuts/', cohort, '/pElems/eMET_orig_pElmens.tsv')
#     }
#     
#     create_pElemDF(path_donorInfo, cohort, path_predRate, path_test_Y,
#                    exclude_lymph_melanoma,  exclude_hyper_mutated, save_name)
#     
#   }
#   print('==========================================================')
# }
# 
# 
