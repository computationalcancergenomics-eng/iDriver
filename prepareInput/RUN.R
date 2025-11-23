rm(list = ls())

library(data.table)
library(dplyr)
library(rtracklayer)
library(readxl)

source("prepareInput/functional_impact/functions.R")
source("prepareInput/functions.R")


############## save MAF file of ICGC and TCGA ##################
path_raw_tcga <- "../extdata/rawInput/controlled/final_consensus_passonly.snv_mnv_indel.tcga.controlled.maf.gz"
path_raw_icgc <- "../extdata/rawInput/download?fn=%2FPCAWG%2Fconsensus_snv_indel%2Ffinal_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz"
path_save <- "../extdata/procInput/mut_PCAWG/"

create_mut_df(path_raw_tcga, path_raw_icgc, path_save)

############## add CADD score to the complete MAF file ##################
path_complete_MAF <- "../extdata/procInput/mut_PCAWG/mutationFile_all_withCohorts.tsv"

path_SNV_offLine <- "../extdata/procInput/CADD_scores/final/offLine/SNV/"
path_indel_offLine <- "../extdata/procInput/CADD_scores/final/offLine/InDel/"
path_onLineSet <- "../extdata/procInput/CADD_scores/final/online.tsv"


reannotated_MAF <- add_CADD_2MAF(path_SNV_offLine, path_indel_offLine,
                                 path_onLineSet, path_complete_MAF)

###### just include the white list samples (just retain included donors) ###########

path_pcawgSupp <- "../extdata/rawInput/PCAWG_supp/pcawg_donor_clinical_August2016_v9.xlsx"
reannotated_MAF <- retain_included_Donors(reannotated_MAF, path_pcawgSupp) # dim 47223205       15
reannotated_MAF <- reannotated_MAF[!duplicated(reannotated_MAF),] # 46894619

####################### rank CADD scores #########################
CADDannotated_maf <- add_ranked_scores(reannotated_MAF)

####################### Ceate GRanges from MAF ###################
gr <- MutFile2GRanges(CADDannotated_maf)

gr <- filter_samples(path_pcawgSupp, gr, exclude_lymph_melanoma = FALSE,
                     exclude_hyper_mutated = FALSE)
save(gr, file = "../extdata/procInput/mut_PCAWG/withCADD/all_samples.RData")

