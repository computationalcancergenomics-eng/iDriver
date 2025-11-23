rm(list = ls())
library(data.table)
library(ggplot2)
library(data.table)

source('tidy/functions_iDriver.R')
source('benchmark/functions_expression.R')
source('benchmark/plot_functions.R')

################################## Inputs ########################################
path_sample_sheet <- '../extdata/rawInput/pcawg_sample_sheet.tsv'
path_IDs_maf <- '../extdata/sampleID_donorID.tsv'
path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
path_CN <- '../extdata/rawInput/PCAWG_data/gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt.gz'
path_FPKM_UQ <- '../extdata/rawInput/PCAWG_data/expression/tophat_star_fpkm_uq.v2_aliquot_gl_ncg.tsv.gz'
path_mutData <- '../extdata/procInput/iDriverInputs/mutData.tsv'



path_sigHits <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/novelHits_expressionInput.tsv'
path_save_table <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/expressionTables/'
output_pathPlot = '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/expressionPlots/'

################################# save sigHits table #################################################

save_expressionTable_sigHits(path_sigHits, path_mutData, path_sample_sheet, path_IDs_maf,
                                         path_donorInfo, path_FPKM_UQ, path_CN)

############################################### plotting #############################################

path_sigHits <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/expressionTables/drivers_expression.tsv'

save_all_expression_boxplots(path_sigHits, path_mutData,
                             path_sample_sheet, path_IDs_maf,
                             path_donorInfo, path_CN, path_FPKM_UQ, 
                             output_pathPlot, save_name = "expression_boxplots")


