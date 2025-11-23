rm(list = ls())

library(stringr)
library(RColorBrewer)
library(rtracklayer)
library(gridExtra)
library(grid)

source('tidy/functions_iDriver.R')
source('benchmark/plot_functions.R')

################################################################################
path_hits <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/novelHits.tsv'
save_dir <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/lolipops/'

path_gr <- '../extdata/procInput/mut_PCAWG/all_samples.RData'
path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
path_to_GEs = '../extdata/procInput/callable_testGE.RData' 

load(path_gr)
load(path_to_GEs)



novelHits <- fread(path_hits)
novelHits <- as.data.frame(novelHits[which(!novelHits$cohort %in% c('Pan_Cancer', 'Lymph-BNHL', "Pancan-no-skin-melanoma-lymph"))])

exclude_hyper_mutated <- T

for (cohort in unique(novelHits$cohort)) {
  
  exclude_lymph_melanoma <- ifelse(cohort ==  "Pancan-no-skin-melanoma-lymph", T, F)
  
  donorInfo <- select_cohort(path_donorInfo, cohort, exclude_lymph_melanoma,
                            exclude_hyper_mutated)
  
  print(cohort)
  gr_cohort = gr[which(gr$D_id %in% donorInfo$donor_id)]
  cohort_hits <- novelHits[which(novelHits$cohort == cohort),]
  mutation <- generate_mut_df(gr = gr_cohort, testGEs = testGE, sigHits = cohort_hits$PCAWG_ID)
  
  save_lolipops(mutation, testGE, path_out = save_dir, save_name = cohort)
}


save_multiple_lolipops(novelHits, testGE, save_dir)
