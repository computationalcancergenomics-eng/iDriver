rm(list = ls())
library(data.table)
library(rtracklayer)
library(dplyr)
source('tidy/functions_iDriver.R')


map_muts <- function(gr, path_bed){
  
  testGE_gr <- import.bed(path_bed)
  
  names(testGE_gr) = testGE_gr$name
  
  GenomicElement <- unlist(lapply(strsplit(names(testGE_gr), "[::]"), function(x){x[1]}))
  
  lenElement <- width(testGE_gr)
  
  mcols(testGE_gr) <- DataFrame(mcols(testGE_gr), GenomicElement, lenElement)
  
  callable_GE <- split(testGE_gr, names(testGE_gr))
  lenElement_new <- sum(width(callable_GE))
  
  ov <- findOverlaps(gr, callable_GE, ignore.strand = TRUE)
  
  idx_gr <- queryHits(ov)
  idx_callable_GE <- subjectHits(ov)
  
  gr_mappedGE <- gr[idx_gr]
  
  GE_col=unique(mcols(callable_GE, level = "within")[,"GenomicElement"])[idx_callable_GE]
  EL_col = sum(mcols(callable_GE, level = "within")[,"lenElement"])[idx_callable_GE]
  name_col = names(callable_GE)[idx_callable_GE]
  
  mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE),
                                  GenomicElement=GE_col,
                                  elemenntLength = EL_col,
                                  name = name_col)
  
  idx_completely_black <- which(lenElement_new == 0)
  
  
  
  if(length(idx_completely_black) != 0){
    lenElement_new <- lenElement_new[-idx_completely_black]
    callable_GE <- callable_GE[-idx_completely_black]
  }
  
  
  list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE,
       lenAllTests = lenElement_new)
  
}


create_response <- function(all_elements, n_donors){
  
  obs_Muts <- all_elements$MAF_GE_mapped
  df <- as.data.frame(DataFrame(obs_Muts))
  df <- df[,c('name', 'D_id')]
  
  
  resTab <- df %>% 
    count(name, D_id) %>%
    filter(n != 0) %>%
    group_by(binID = name) %>%
    summarise(nMut = sum(n), nSample = n(), .groups = "drop")
  
  idx_lenElement <- resTab$binID
  all_lmnt_length <- all_elements$lenAllTests
  ElementLength <- all_lmnt_length[idx_lenElement]
  resTab$length <- ElementLength
  
  resTab$N <- n_donors
  
  all_binIds <- names(all_elements[["lenAllTests"]])
  nonMut_ids <- all_binIds[which(!all_binIds %in% resTab$binID)]
  res_nonMutated <- data.frame(cbind('binID' = nonMut_ids, 
                                     'nMut' = rep(0, length(nonMut_ids)),
                                     'nSample' = rep(0, length(nonMut_ids)),
                                     'length' = all_lmnt_length[!names(all_lmnt_length) %in% resTab$binID], 
                                     'N' = rep(n_donors, length(nonMut_ids))))
  
  resTab <- rbind(resTab, res_nonMutated)
  
  return(resTab)
}



save_responseTable <- function(gr, path_bed, path_out, save_name){
  
  all_elements <- map_muts(gr, path_bed)
  
  # print("************************ mutations mapped to the genomic intervals************************")
  n_donors <- length(unique(gr$D_id))
  resTab <- create_response(all_elements, n_donors)
  # print("************************ response table was generated ************************")
  dir.create(paste0(path_out, "/"), showWarnings = FALSE, recursive = TRUE)
  
  fwrite(resTab, paste0(path_out, "/", save_name,
                        ".tsv"), 
         sep = "\t", row.names = F) 
  # print("************************file saved************************")
  
}

# select_cohort <- function(path_donorInfo, 
#                           cohort, exclude_lymph_melanoma = TRUE,
#                           exclude_hyper_mutated = TRUE){
#   
#   donorInfo <- fread(path_donorInfo)
#   
#   if (exclude_lymph_melanoma) {
#     exceptions <- c("Skin-Melanoma", "SKCM-US",
#                     "Lymph-NOS", 
#                     "Lymph-CLL", "CLLE-ES",
#                     "Lymph-BNHL", "MALY-DE", "DLBC-US")
#     donorInfo <- donorInfo[-which(donorInfo$cohort1 %in% exceptions),] # 2279 donors
#   }
#   
#   if (exclude_hyper_mutated) {
#     donorInfo <- donorInfo[which(!donorInfo$HyperMut_donor),] 
#   }
#   
#   if (unique(!cohort %in% c('Pan_Cancer', 'Pancan-no-skin-melanoma-lymph'))) {
#     
#     donorInfo <- donorInfo[which(donorInfo$cohort1 %in% c(cohort)),]
#     
#   } 
#   
#   donorInfo <- donorInfo[,c("D_id","freq" )]
#   colnames(donorInfo) <- c('donor_id', 'totalMutNrs')
#   donorInfo
#   
# }


library(readxl)


get_project_code <- function(histology) {
  
  df <- read_xlsx('../extdata/rawInput/PCAWG_supp/pcawg_specimen_histology_August2016_v9.xlsx')
  df <- df[, c('project_code', 'histology_abbreviation')]
  df <- as.data.frame(df[!duplicated(df),])
  
  # Filter the data frame based on the given histology_abbreviation
  result <- df[df$histology_abbreviation == histology, "project_code"]
  
  # Return unique project codes (to handle multiple matches)
  unique(na.omit(result))
}

# ################################## RUN ##########################################
# study = 'simulated'
# exclude_hyper_mutated = FALSE
# cohorts <- c(
#   'Pancan-no-skin-melanoma-lymph', 'Pan_Cancer',"Liver-HCC",
#              "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
#              "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
#              "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC", "Biliary-AdenoCA",
#   "Breast-AdenoCa", "Eso-AdenoCa", "CNS-GBM",
#   "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
#   "Breast-LobularCa", "Thy-AdenoCA", "Myeloid-MPN",
#              "Bone-Leiomyo",   
#              "Lymph-NOS", "CNS-Medullo", "Myeloid-AML", "CNS-Oligo", "Cervix-SCC",
#              "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Epith", "Bone-Osteosarc",
#              "Cervix-AdenoCA", "Breast-DCIS", "Bone-Cart", "Myeloid-MDS" ) 
# 
# 
# cohorts <- c('Pancan-no-skin-melanoma-lymph', 'Pan_Cancer', 'Skin-Melanoma', 'Biliary-AdenoCA', 'CNS-GBM', 'ColoRect-AdenoCA',  'Eso-AdenoCa',
#              'Head-SCC', 'Lung-AdenoCA', 'Lung-SCC',
#                'Stomach-AdenoCA',  'Uterus-AdenoCA',
#              "Lymph-BNHL")
# 
# 
# if (study == 'observed') {
#   path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
#   path_to_gr='../extdata/procInput/mut_PCAWG/all_samples.RData'
#   base_dir <- '../extdata/procInput/BMRs/observed/'
#   
#   
# } else {
#   path_donorInfo <- '../extdata/procInput/iDriverInputs/simulated/Sanger/donorInfo.tsv'
#   path_to_gr='../extdata/procInput/simulated/sanger/all_sanger_gr_removedDups.RData'
#   base_dir <- '../extdata/procInput/BMRs_2024/simulated/Sanger_withHyperMuts/'
#   cohorts = cohorts[which(!cohorts %in% c("Bone-Leiomyo", "Lymph-NOS", "Bone-Cart"))]
# }
# 
# load(path_to_gr)
# 
# 
# path_beds <- c('../../make_features/external/database/bins/proccessed/PCAWG_callable.bed6',
#                '../../make_features/external/database/bins/proccessed/callable_intergenic_intervals_wo_pcawg.bed6')
# 
# 
# save_names = c('test_y', 'train_y')
# 
# for (cohort in cohorts) {
#   
#   print(cohort)
#   path_out <- paste0(base_dir, cohort, '/')
#   exclude_lymph_melanoma = ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
#   
#   
#   if (study == 'simulated') {
#     # "Bone-Leiomyo"     "Lymph-NOS"        "Bone-Cart"
#     if (cohort %in% c('Pancan-no-skin-melanoma-lymph', 'Pan_Cancer')) {
#       project_code = cohort
#     } else {
#       if (cohort  == 'Breast-AdenoCa'){
#         cohort = 'Breast-AdenoCA'
#       } else if (cohort == 'Breast-LobularCa'){
#         cohort = 'Breast-LobularCA'
#       } else if (cohort == "Eso-AdenoCa"){
#         cohort = "Eso-AdenoCA"
#       }
#       project_code = get_project_code(cohort) 
#       
#       
#     }
#     
#   } else {
#     
#     project_code = cohort
#     
#   }
#   print(project_code)
#   donor_info <- select_cohort(path_donorInfo,
#                               cohort, exclude_lymph_melanoma,
#                               exclude_hyper_mutated)
#   
#   if (study == 'observed') {
#     gr_cohort <- gr[which(gr$D_id %in% donor_info$donor_id)]
#   } else if (study == 'simulated') {
#     gr_cohort <- gr[which(gr$D_id %in% donor_info$sampleID)]
#   }
#   
#   
#   for(i in 1:length(path_beds)){
#     save_responseTable(gr_cohort, path_beds[i], path_out, save_names[i])
#   }
#   print("************************files saved************************")
# }
# 
# 
# 
