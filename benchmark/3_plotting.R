rm(list = ls())

library(data.table)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)  
library(dplyr)
source('benchmark/functions_benchmark.R')
source('benchmark/plot_functions.R')

################# save methods color for the legend ####################
path_save <- "../extdata/output_release2.0/Benchmark/method_legend2.png"
save_legend_methodColors(path_save)

############################ Fig3a_c: AUPR grouped barplot  ################
nonHyperCancers <- c(
  "Biliary-AdenoCA",  "Bladder-TCC",# "Bone-Leiomyo",
  "Bone-Osteosarc", "Pancan-no-skin-melanoma-lymph",
  "Breast-AdenoCa",   "CNS-GBM",   "CNS-Medullo",
  #            "CNS-Oligo",   "Cervix-SCC", 
  "CNS-PiloAstro",   "ColoRect-AdenoCA",
  "Eso-AdenoCa",  "Head-SCC" ,  "Kidney-RCC",   "Pan_Cancer",
  "Liver-HCC", "Lung-AdenoCA", "Lung-SCC", "Lymph-BNHL",
  "Lymph-CLL",  "Myeloid-MPN", "Ovary-AdenoCA",
  "Panc-AdenoCA",     "Panc-Endocrine", "Prost-AdenoCA",
  "Skin-Melanoma",
  "Stomach-AdenoCA",  "Thy-AdenoCA", "Uterus-AdenoCA"
) #  for "Cervix-SCC" there is just one method to compare with 
hyperCancers <- c()

path_save_tables_rmHypers <- "../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/tables/{cancer}_table_GoldStd_basedon_in_CGC_new_{elem}_fdr.csv"

path_save_tables_withHypers <- NA

save_name <- ''
path_save <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/plots/'
path_save <- '../extdata/output_release2.0/Figures/MainFigures/'


for (measure in c( 'AUPR')) { # 'AUC', 'F1'
  
  for (element_type in c('CDS', 'NC')) {
    
    save_grouped_barPlotAllCohorts(save_name, path_save, hyperCancers, nonHyperCancers,
                                   path_save_tables_withHypers,
                                   path_save_tables_rmHypers, measure, k = 5,
                                   element_type)
    
  }
}


############################ Fig3b_d: AUPR line plot  ################
PATH_newRes <- '../extdata/output_release2.0/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/observed_Pancan-no-skin-melanoma-lymph.tsv'
path_save_plot <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/plots/'

for (element in c('CDS', 'NC')) {
  save_AUPR_singleCohort(n = 600, PATH_newRes, 'iDriver',
                         element, path_save_plot, based_on = 'in_CGC_new',
                         selected_methods_plots = NA,
                         'Pancan-no-skin-melanoma-lymph', 
                         base_method = 'iDriver',
                         path_ann_PCAWG_ID = '../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv',
                         selected_methods_tables = NA
  )
}

############################ Fig3e_f: top-k methods heatmap ################
rm(list = ls())
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)

source('benchmark/plot_functions.R')
source('benchmark/functions_benchmark.R')

path_save_tables_rmHypers <- "../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/tables/{cancer}_table_GoldStd_basedon_in_CGC_new_{elem}_fdr.csv"
path_saveFig <- '../extdata/output_release2.0/Figures/'

save_heatmap_kTop(path_save_tables_rmHypers, path_saveFig)

############################ SuppFig (cancer type line plots ... rmHyper) ############################
rm(list = ls())

source('benchmark/plot_functions.R')
source('benchmark/functions_benchmark.R')


PATHs_newResults_allCohorts <- '../extdata/output_release2.0/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/observed_{cohort}.tsv'
# path_save_plot <- '../extdata/output_release2.0/Figures/linePlots_allcohorts_iDriverBest/'


nonHyperCancerss <- c(
  "Bladder-TCC", "Bone-Osteosarc", "Breast-AdenoCa", #"Bone-Leiomyo"
  "CNS-Medullo", "CNS-PiloAstro", "Kidney-RCC", "Liver-HCC",
  "Lymph-CLL", "Myeloid-MPN", "Ovary-AdenoCA" ,
  "Panc-AdenoCA", "Panc-Endocrine", "Prost-AdenoCA", "Thy-AdenoCA" ) #

save_AUPR_allCohorts(n = 150, PATHs_newResults_allCohorts, newRESULTS = 'iDriver',
                     path_save_plot =  '../extdata/output_release2.0/Figures/Supplementary/linePlots_allcohorts_noHyperCohorts/', 
                     based_on = 'in_CGC_new',
                     selected_methods_plots = NA,
                     cohorts = nonHyperCancerss, 
                     base_method = 'iDriver',
                     path_ann_PCAWG_ID = '../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv',
                     selected_methods_tables = NA)

hyper_rm_Cancers <- c( "ColoRect-AdenoCA","CNS-GBM", 
                  "Head-SCC", "Skin-Melanoma", "Lung-AdenoCA", "Lung-SCC",
                  "Stomach-AdenoCA", "Uterus-AdenoCA" ,  "Biliary-AdenoCA",
                  "Eso-AdenoCa", "Lymph-BNHL")
save_AUPR_allCohorts(n = 150, PATHs_newResults_allCohorts, newRESULTS = 'iDriver',
                     path_save_plot =  '../extdata/output_release2.0/Figures/Supplementary/linePlots_allcohorts_rmHyperCohorts/', 
                     based_on = 'in_CGC_new',
                     selected_methods_plots = NA,
                     cohorts = hyper_rm_Cancers, 
                     base_method = 'iDriver',
                     path_ann_PCAWG_ID = '../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv',
                     selected_methods_tables = NA)

PATHs_newResults_allCohorts <- '../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_{cohort}.tsv'
hyper_Cancers <- c( "ColoRect-AdenoCA","CNS-GBM", 
                       "Head-SCC", "Skin-Melanoma", "Lung-AdenoCA", "Lung-SCC",
                       "Stomach-AdenoCA", "Uterus-AdenoCA" ,  "Biliary-AdenoCA",
                       "Eso-AdenoCa", "Lymph-BNHL")
save_AUPR_allCohorts(n = 150, PATHs_newResults_allCohorts, newRESULTS = 'iDriver',
                     path_save_plot =  '../extdata/output_release2.0/Figures/Supplementary/linePlots_allcohorts_HyperCohorts/', 
                     based_on = 'in_CGC_new',
                     selected_methods_plots = NA,
                     cohorts = hyper_Cancers, 
                     base_method = 'iDriver',
                     path_ann_PCAWG_ID = '../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv',
                     selected_methods_tables = NA)



cohorts <- c(#"Bone-Leiomyo", "CNS-PiloAstro",  "Myeloid-MPN",
  "Biliary-AdenoCA",  "Bladder-TCC", 
  "Bone-Osteosarc", 
  "Breast-AdenoCa",   "CNS-GBM",   "CNS-Medullo",
  # "CNS-Oligo",   "Cervix-SCC", 
  "ColoRect-AdenoCA",
  "Eso-AdenoCa",  "Head-SCC" ,  "Kidney-RCC", 
  "Liver-HCC", "Lung-AdenoCA", "Lung-SCC", "Lymph-BNHL",
  "Lymph-CLL",  
  "Ovary-AdenoCA",
  "Panc-AdenoCA",     "Panc-Endocrine",   "Prost-AdenoCA",
  "Skin-Melanoma",
  "Stomach-AdenoCA",  "Thy-AdenoCA", "Uterus-AdenoCA"
) 

############################ SuppFig (Pancancers line plots ... withHyper) ############################
rm(list = ls())
source('benchmark/plot_functions.R')
source('benchmark/functions_benchmark.R')

PATH_newResults <- c('../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pancan-no-skin-melanoma-lymph.tsv',
                     '../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv' )
path_save_plot <- '../extdata/output_release2.0/Figures/withHyper_plots/'

tissues <- c('Pancan-no-skin-melanoma-lymph', 'Pan_Cancer')
for (element in c('CDS', 'NC')) {
  for (i in 1:length(tissues)) {
    PATH_newResult = PATH_newResults[i]
    tissue = tissues[i]
    save_AUPR_singleCohort(n = 600, PATH_newResult, newRESULT = 'iDriver',
                           element, path_save_plot, based_on = 'in_CGC_new',
                           selected_methods_plots = NA,
                           tissue, 
                           base_method = 'iDriver',
                           path_ann_PCAWG_ID = '../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv',
                           selected_methods_tables = NA
    )
  }
}
############################ SuppFig (withHyper) ############################
nonHyperCancers <- c("Pan_Cancer", "Pancan-no-skin-melanoma-lymph" , "ColoRect-AdenoCA","CNS-GBM", 
                     "Head-SCC", "Skin-Melanoma", "Lung-AdenoCA", "Lung-SCC",
                     "Stomach-AdenoCA", "Uterus-AdenoCA" ,  "Biliary-AdenoCA",
                     "Eso-AdenoCa", "Lymph-BNHL")
hyperCancers <- c()
path_save_tables_rmHypers <-"../extdata/output_release2.0/Benchmark/pooled_catQuartiles_withHyperMuts/tables/{cancer}_table_GoldStd_basedon_in_CGC_new_{elem}_fdr.csv"
path_save_tables_withHypers <- NA
save_name <- ''
path_save <- '../extdata/output_release2.0/Figures/Supplementary/withHyper/'

for (measure in c( 'AUPR')) { # 'AUC', 'F1'
  
  for (element_type in c('CDS', 'NC')) {
    save_kTop_methods(save_name, path_save, hyperCancers, nonHyperCancers,
                      path_save_tables_withHypers,
                      path_save_tables_rmHypers, measure, k = 2, element_type)
    
    save_grouped_barPlotAllCohorts(save_name, path_save, hyperCancers, nonHyperCancers,
                                   path_save_tables_withHypers,
                                   path_save_tables_rmHypers, measure, k = 5,
                                   element_type)
    
  }
}

path_saveFig <- '../extdata/output_release2.0/Figures/Supplementary/withHyper/'

save_heatmap_kTop(path_save_tables_rmHypers, path_saveFig) # Should change some lines in the prepare_data_kTop() from functions_benchmark.R

################# supp Fig OnlyHyper #####################
rm(list = ls())
source('OnlyHyperMuts/compare_onlyHyper_subsamples.R')
nTop <- 50
PATHs <- c('../extdata/output_release2.0/observed/onlyHyperMuts/update_eMET/observed_withHyperUpdate_cancerSpRho.tsv',
           '../extdata/output_release2.0/subsamples69/subsamples_subsample_1.tsv',
           '../extdata/output_release2.0/subsamples69/subsamples_subsample_2.tsv',
           '../extdata/output_release2.0/subsamples69/subsamples_subsample_3.tsv',
           '../extdata/output_release2.0/subsamples69/subsamples_subsample_4.tsv',
           '../extdata/output_release2.0/subsamples69/subsamples_subsample_5.tsv')

newRES <- c('Hypermutated ', "Subsample 1", "Subsample 2", "Subsample 3", "Subsample 4", "Subsample 5")

for(elem in c('CDS', 'non_coding')){
  save_onlyHyperPlo(elem, nTop, PATHs, newRES, 
                    path_save = '../extdata/output_release2.0/Figures/Supplementary/')
  
}

############################ SuppFig (scatter perElem cadd scores) ############################
rm(list = ls())
library(ggplot2)
library(data.table)
library(dplyr)
library(scales)
library(ggrepel)
source('benchmark/plot_functions.R')

path_res <- '../extdata/output_release2.0/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/observed_Pancan-no-skin-melanoma-lymph.tsv'
path_sigHits <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/allSigHits.tsv'
path_save <- '../extdata/output_release2.0/Figures/MainFigures/'

save_scores_scatterPlots(path_res, path_sigHits, path_save, 'CDS',save_name = 'Fig1')
save_scores_scatterPlots(path_res, path_sigHits, path_save, 'NC',save_name = 'Fig1')

############################ SuppFig (nMuts per cohort) ############################
rm(list = ls())

source('benchmark/plot_functions.R')
library(dplyr)
library(ggplot2)
library(data.table)

path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
path_save <- '../extdata/output_release2.0/Figures/MainFigures/'


nMuts_perCohort_stackedPlot(path_donorInfo, path_save, save_name = 'Fig1a')
############################ Fig (Top 65 mutated genes accross samples with scores (stack)) ############################
rm(list = ls())

library(dplyr)
library(data.table)
library(ggplot2)
library(ggtext)  # for element_markdown()

source('benchmark/plot_functions.R')
source('benchmark/functions_annotate_elements.R')
path_perDonorOnsStats <- '../../../../Projects/bahari_work/pancan_star_perDonorObsStats_MutatedDonorElem.tsv'
path_testY <- '../extdata/procInput/BMRs_2024/observed/Pancan-no-skin-melanoma-lymph/test_y.tsv'
path_save <- '../extdata/output_release2.0/Figures/MainFigures/'

# plot_sumScores_stacked(path_perDonorOnsStats, path_testY, path_save, 'SigHits')
# plot_sumScores_stacked(path_perDonorOnsStats, path_testY, path_save, 'nonDrivers')
plot_sumScores_stacked(path_perDonorOnsStats, path_testY, path_save)
plot_sumScores_stacked_v2(path_perDonorOnsStats, path_testY, path_save)

################### nFPs heatmap ######################
rm(list = ls())
library(circlize)
library(ggtext)
library(ggplot2)
library(data.table)

source('benchmark/plot_functions.R')

path_all_methods_nFPs <- '../extdata/output_release2.0/simulated/all_methods_nFPs/nFPs_vsPCAWG_dpIdrBase_codingNCoding.tsv'
path_save <- '../extdata/output_release2.0/Figures/MainFigures/'
save_nFPs_heatmap(path_all_methods_nFPs, elem = 'Coding',path_save)
save_nFPs_heatmap(path_all_methods_nFPs, elem = 'Non-coding',path_save)
save_nFPs_heatmap(path_all_methods_nFPs = '../extdata/output_release2.0/simulated/all_methods_nFPs/nFPs_vsPCAWG.tsv', elem = '',path_save)

################### stacked Bar Ablation ######################
rm(list = ls())
library(data.table)
library(ggplot2)

source('benchmark/plot_functions.R')

path_tables <- '../extdata/output_release2.0/Benchmark/iDriver_variations/finalTabs/'
path_save <- '../extdata/output_release2.0/Figures/MainFigures/'

save_ablation_stackedPlot(path_tables, path_save)



################### upset Ablation ######################
rm(list = ls())

library(data.table)
library(UpSetR)
library(dplyr)

source('benchmark/plot_functions.R')

path_files <- c('../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/allSigHits.tsv',
                '../extdata/output_release2.0/Benchmark/populationLevel_catQuartiles_rmHyperMuts/SigHits/allSigHits_popLvl.tsv',
                '../extdata/output_release2.0/Benchmark/only_rec_SigHits_allSigHits.tsv')
path_save <- '../extdata/output_release2.0/Figures/MainFigures/'


save_upsetAblations(path_files, path_save, save_name = '', Just_cds = T)
save_upsetAblations(path_files, path_save, save_name = 'allElems', Just_cds = F)

################### compare FI only methods (CADD) ################### 
rm(list = ls())

library(data.table)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)  
source('benchmark/functions_benchmark.R')
source('benchmark/plot_functions.R')

PATH_newResults <- c('../extdata/output_release1.0/observed/burdenFree_rmHypermutated_Rho0/observed_observed_Pancan-no-skin-melanoma-lymph_betaRho0.tsv',
                     # '../extdata/output_release2.0/observed/Burden_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pancan-no-skin-melanoma-lymph.tsv',
                     '../extdata/output_release1.0_tmp_testing/observed/compCADD/output.tsv',
                     '../extdata/output_release1.0_tmp_testing/observed/compCADD/pcawg_oncodriveFML_pancannoLnoM.tsv')

newRES <- c('iDriver', 'compositeDriverCADD', "oncodriveFML_cadd")
tissue <- 'Pancan-no-skin-melanoma-lymph'
path_save_plot <- '../extdata/output_release1.0/Benchmark/iDriver_onlyFI_plots/'#'../extdata/output_release2.0/Benchmark/iDriver_onlyFI_plots/'


for (element in c('CDS', 'NC')) {
  
  PATH_newRes = PATH_newResults
  save_AUPR_singleCohort_woPCAWGmeyhods(n = 600, PATH_newRes, newRES, 
                                        element, path_save_plot, based_on = 'in_CGC_new',
                                        tissue, base_method = 'iDriver',
                                        path_ann_PCAWG_ID = '../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv',
                                        include_legend = T
  )
}
#################### iDriver variations (rmHyper--AUPR line plots) #################### 
rm(list = ls())

library(data.table)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)  
source('benchmark/functions_benchmark.R')
source('benchmark/plot_functions.R')

PATH_newResults <- c('../extdata/output_release2.0/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/observed_Pancan-no-skin-melanoma-lymph.tsv',
                     '../extdata/output_release2.0/observed/Burden_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pancan-no-skin-melanoma-lymph.tsv',
                     '../extdata/output_release2.0/observed/population_level/pooled_betaRhoPancan_lengthCat_quartilesElemType/observed_Pancan-no-skin-melanoma-lymph.tsv',
                     '../extdata/output_release2.0/observed/Rho0/pooled_betaRho0_lengthCat_quartilesElemType/observed_Pancan-no-skin-melanoma-lymph.tsv',
                     '../extdata/output_release2.0/observed/FI_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pancan-no-skin-melanoma-lymph.tsv')

newRES <- c('iDriver', "iDriver (FI-Only)", "iDriver (Population-level)",
            "iDriver (Uncorrelated)", "iDriver (Recurrence-Only)")

path_save_plot <- '../extdata/output_release2.0/Benchmark/iDriver_variations/LinePlots/'


for (element in c('CDS', 'NC')) {
  
  PATH_newRes = PATH_newResults
  save_AUPR_singleCohort_woPCAWGmeyhods(n = 600, PATH_newRes, newRES, 
                                        element, path_save_plot, based_on = 'in_CGC_new',
                                        'Pancan-no-skin-melanoma-lymph', base_method = 'iDriver',
                                        path_ann_PCAWG_ID = '../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv',
                                        include_legend = T
  )
}

###################
library(data.table)
library(dplyr)
source('benchmark/functions_benchmark.R')
source('benchmark/plot_functions.R')

cancer_types <- c(#"Bone-Leiomyo", "CNS-PiloAstro",  "Myeloid-MPN",
  "Biliary-AdenoCA",  "Bladder-TCC", 
  "Bone-Osteosarc", 
  "Breast-AdenoCa",   "CNS-GBM",   "CNS-Medullo",
  # "CNS-Oligo",   "Cervix-SCC", 
  "ColoRect-AdenoCA",
  "Eso-AdenoCa",  "Head-SCC" ,  "Kidney-RCC", 
  "Liver-HCC", "Lung-AdenoCA", "Lung-SCC", "Lymph-BNHL",
  "Lymph-CLL",  
  "Ovary-AdenoCA",
  "Panc-AdenoCA",     "Panc-Endocrine",   "Prost-AdenoCA",
  "Skin-Melanoma",
  "Stomach-AdenoCA",  "Thy-AdenoCA", "Uterus-AdenoCA"
) 

path_save_tables_rmHypers <- "../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/tables/{cancer}_table_GoldStd_basedon_in_CGC_new_{elem}_fdr.csv"



get_topCohorts <- function(path_save_tables_rmHypers, cancer_types){
  methods <- extract_topK_methods(cancer_types, measure = 'AUPR', k = 1,
                                  path_save_tables_rmHypers, element_type = 'CDS')
  methods <- methods[methods$method == 'iDriver',]
  
  methods_nc <- extract_topK_methods(cancer_types, measure = 'AUPR', k = 1,
                                     path_save_tables_rmHypers, element_type = 'NC')
  methods_nc <- methods_nc[methods_nc$method == 'iDriver',]
  both_top <- intersect(methods$cancer_type, methods_nc$cancer_type)
  nc_top <- methods_nc$cancer_type[!methods_nc$cancer_type %in% both_top]
  
  list('nc_top' = nc_top,
       'both_top' = both_top)
}

get_topCohorts(path_save_tables_rmHypers, cancer_types)