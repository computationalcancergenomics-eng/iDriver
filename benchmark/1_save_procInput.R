rm(list = ls())

library(data.table)
library(openxlsx)
library(dplyr)

source("benchmark/functions_process_PCAWG_results.R")
source("benchmark/functions_annotate_elements.R")

path_save_procPCAWG_res <- "../extdata/procInput/PCAWG_results/"
path_rawPCAWG_res <- "../extdata/rawInput/PCAWG_12methods_output/p-values/observed/"
path_genomic_intervals <- "../extdata/rawInput/genomic_intervals_lists/"
path_to_PCAWG_supp <- "../extdata/rawInput/PCAWG_nonCoding_supp/41586_2020_1965_MOESM4_ESM.xlsx"
path_proccessedGE <- "../extdata/procInput/"

path_to_cgc_v101 <- '../extdata/rawInput/v101_Census_allSun Jan 26 14_06_42 2025.csv'
path_CGC_new <- "../extdata/rawInput/v101_Census_allSun Jan 26 14_06_42 2025.csv"
path_oncoKB <- "../extdata/rawInput/cancerGeneList_oncoKB_Jan28_2025.tsv"
path_test_y <- "../extdata/procInput/BMRs/SNVs/Pancan-no-skin-melanoma-lymph/SNVs_test_y.tsv"
path_genhancer_PCAWG <- "../extdata/rawInput/map.enhancer.gene.txt.gz"

########### 1) save the processed PCAWG results as a pRes object ###########
save_pRes(path_rawPCAWG_res, path_save_procPCAWG_res)

########### 2) save the PCAWG_IDs annotated by in_CGC, in_PCAWG, type_of_element,.. ###########
comp_ann_PCAWG_ID <- Complete_annotated_PCAWG_ID(path_to_PCAWG_supp, path_test_y,
                                                 path_CGC_new, path_oncoKB, 
                                                 path_genhancer_PCAWG)

write.csv(comp_ann_PCAWG_ID, file = paste0(path_proccessedGE, 
                                                 "ann_PCAWG_ID_complement_2025.csv"),
          row.names = F)

