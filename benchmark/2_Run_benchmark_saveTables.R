rm(list = ls())
library(pROC)
library(PRROC)
library(data.table)
library(dplyr)
source('benchmark/functions_benchmark.R')


based_on <- 'in_CGC_new'
path_ann_PCAWG_ID <- ('../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv')
path_save <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/'
# path_save <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_withHyperMuts/'


sig_definition_method <- 'fdr'
sig_definition_threshold <- .1
baseMethod <- 'iDriver'
elements <- c('NC', 'CDS')

cohorts <- c(
  "Pan_Cancer", "Pancan-no-skin-melanoma-lymph" ,
  "Biliary-AdenoCA",  "Bladder-TCC", "Bone-Leiomyo",  "Bone-Osteosarc",
  "Breast-AdenoCa",   "Cervix-SCC",  "CNS-GBM",   "CNS-Medullo",
    "CNS-PiloAstro",    "ColoRect-AdenoCA",
  "Eso-AdenoCa",  "Head-SCC",   "Kidney-RCC",
  "Liver-HCC", "Lung-AdenoCA", "Lung-SCC", "Lymph-BNHL",
  "Lymph-CLL" ,  "Myeloid-MPN", "Ovary-AdenoCA",
  "Panc-AdenoCA", "Panc-Endocrine", "Prost-AdenoCA", "Skin-Melanoma",
  "Stomach-AdenoCA" , "Thy-AdenoCA"  ,    "Uterus-AdenoCA"
) # "CNS-Oligo", "Kidney-ChRCC"

# cohorts=c("Pan_Cancer", "Pancan-no-skin-melanoma-lymph" , "ColoRect-AdenoCA","CNS-GBM", 
#           "Head-SCC", "Skin-Melanoma", "Lung-AdenoCA", "Lung-SCC",
#           "Stomach-AdenoCA", "Uterus-AdenoCA" ,  "Biliary-AdenoCA",
#           "Eso-AdenoCa", "Lymph-BNHL"
#           
# )


# Original paths with placeholder for cohort
PATHs_newResults_allCohorts <- c('../extdata/output_release2.0/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/observed_{cohort}.tsv')
# PATHs_newResults_allCohorts <- c('../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_{cohort}.tsv')

all_paths <- lapply(cohorts, function(cohort) {
  gsub("\\{cohort\\}", cohort, PATHs_newResults_allCohorts)
})

names(all_paths) <- cohorts

newRESULTS <- c('iDriver')

for (tissue in cohorts) {
  print(tissue)
  if (tissue == 'Pan_Cancer') {
    cohort = 'Pancan'
  } else {
    cohort = tissue
  }
  print(cohort)
  path_pRes <- paste0('../extdata/procInput/PCAWG_results/', cohort, '.RData')
  
  PATHs_newResults <- all_paths[[tissue]]
  
  for (element in elements) {
    save_Measures(path_save, newRESULTS, PATHs_newResults, path_pRes,
                  sig_definition_method, base_method = baseMethod,
                  sig_definition_threshold, path_ann_PCAWG_ID, tissue,
                  based_on, element,
                  selected_methods = NA,
                  compareRank_new = list(T, 'iDriver'))
  }
  
}

