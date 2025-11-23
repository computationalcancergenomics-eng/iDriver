#################### non hypers ########################
rm(list = ls())

library(data.table)
library(dplyr)

source('benchmark/functions_benchmark.R')

sig_threshold = .1

cohorts <- c("Biliary-AdenoCA",  "Bladder-TCC", "Bone-Leiomyo",
             "Bone-Osteosarc", "Pancan-no-skin-melanoma-lymph",
             "Breast-AdenoCa",   "Cervix-SCC",  "CNS-GBM",   "CNS-Medullo",
              "CNS-Oligo",  
             "CNS-PiloAstro",   "ColoRect-AdenoCA",
             "Eso-AdenoCa",  "Head-SCC" ,  "Kidney-RCC",   "Pan_Cancer",
             "Liver-HCC", "Lung-AdenoCA", "Lung-SCC", "Lymph-BNHL",
             "Lymph-CLL",  "Myeloid-MPN", "Ovary-AdenoCA",
             "Panc-AdenoCA",     "Panc-Endocrine",   "Prost-AdenoCA",
             "Skin-Melanoma",
             "Stomach-AdenoCA",  "Thy-AdenoCA", "Uterus-AdenoCA", "Kidney-ChRCC") #


# cohorts <-c(
#   "Pan_Cancer", "Pancan-no-skin-melanoma-lymph" , "Cervix-SCC",
#   "Biliary-AdenoCA",  "Bladder-TCC", "Bone-Leiomyo",  "Bone-Osteosarc",
#   "Breast-AdenoCa",   "CNS-GBM",   "CNS-Medullo",
#   "CNS-PiloAstro",    "ColoRect-AdenoCA",
#   "Eso-AdenoCa",  "Head-SCC",   "Kidney-RCC",
#   "Liver-HCC", "Lung-AdenoCA", "Lung-SCC", "Lymph-BNHL",
#   "Lymph-CLL" ,  "Myeloid-MPN", "Ovary-AdenoCA",
#   "Panc-AdenoCA", "Panc-Endocrine", "Prost-AdenoCA", "Skin-Melanoma",
#   "Stomach-AdenoCA" , "Thy-AdenoCA"  ,    "Uterus-AdenoCA"
# )

hits <- c()
for (cohort in cohorts) {
  cols <- c("PCAWG_ID",   "sum_scores", "pElement",   "nPred" ,     "mu0"   ,     "sd0"   ,     "sum_nMut" ,
            "obs_stat" ,  "alpha" ,     "beta" ,      "rho" ,
            "n4var"  , "length_category",
            "pvals" ,"qvals" , "cohort")
  # path_file <- gsub('\\{cohort\\}', cohort, '../extdata/output_release2.0/observed/FI_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_{cohort}.tsv')
  
  # path_file <- gsub('\\{cohort\\}', cohort, '../extdata/output_release1.0/observed/pooled/betaRho/observed_betaRho_{cohort}.tsv')
  path_file <- gsub('\\{cohort\\}', cohort, '../extdata/output_release2.0/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/observed_{cohort}.tsv')
  hit_cohort <- get_sigHits_cohort(path_file, sig_threshold)
  hit_cohort$cohort <- cohort
  hit_cohort <- data.frame(hit_cohort)
  hit_cohort <- hit_cohort[,which(colnames(hit_cohort) %in% cols)]
  
  hits <- rbind(hits, hit_cohort)
}

ann <- fread('../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv')
hits <- left_join(hits, ann, by = c('PCAWG_ID' = 'PCAWG_IDs'))


dir.create('../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/', recursive = T, showWarnings = F)
fwrite(hits, file = '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/allSigHits.tsv', sep = '\t')

sum(hits$cohort == 'Lymph-BNHL')
sum(hits$cohort == 'Pan_Cancer')

sum(hits$in_CGC_new | hits$in_pcawg | hits$in_oncoKB | hits$in_CGC_literature)

newHits <- hits[(!(hits$in_CGC_new | hits$in_pcawg | hits$in_oncoKB | hits$in_CGC_literature)),]
fwrite(newHits, '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/novelHits.tsv', sep = '\t')
newHits1 <- newHits[which(!newHits$cohort %in% c("Pan_Cancer", "Pancan-no-skin-melanoma-lymph", "Lymph-BNHL")),]

CDS_new <- newHits1[newHits1$type_of_element == 'gc19_pc.cds',]
newHits1[newHits1$type_of_element == 'gc19_pc.5utr',]
newHits1[newHits1$type_of_element == 'gc19_pc.promCore',]
newHits1[newHits1$type_of_element == 'gc19_pc.3utr',]
nc <- newHits1[newHits1$type_of_element != 'gc19_pc.cds',]

newHits2 <- newHits[which(newHits$cohort %in% c("Pancan-no-skin-melanoma-lymph")),]
table(newHits2$type_of_element)
cds <- newHits2[newHits2$type_of_element == 'gc19_pc.cds',]


hitspop <- fread('../extdata/output_release2.0/Benchmark/populationLevel_catQuartiles_rmHyperMuts/SigHits/allSigHits_popLvl.tsv')
hitspop <- hitspop[,c('PCAWG_ID', 'cohort')]
hitspop$pop_levelHit <- TRUE
hits <- left_join(hits, hitspop, by = c('PCAWG_ID', 'cohort'))
hits$pop_levelHit <- ifelse(is.na(hits$pop_levelHit) , F, T)

fwrite(hits, file = '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/allSigHits_popLvlAdded.tsv', sep = '\t')

hits_rec <- fread('../extdata/output_release2.0/Benchmark/only_rec_SigHits_allSigHits.tsv')
hits_rec <- hits_rec[,c('PCAWG_ID', 'cohort')]
hits_rec$recOnly <- TRUE
hits <- fread('../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/allSigHits_popLvlAdded.tsv')
hits <- left_join(hits, hits_rec, by = c('PCAWG_ID', 'cohort'))
hits$recOnly <- ifelse(is.na(hits$recOnly) , F, T)

##################### with hyperMuts #######################
rm(list = ls())

library(data.table)
library(dplyr)

source('benchmark/functions_benchmark.R')

sig_threshold = .1

# List all files in the directory
files <- list.files("../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/",
                    pattern = "^observed_.*\\.tsv$", full.names = FALSE)

# Extract cohort names
cohorts <- sub("^observed_(.*)\\.tsv$", "\\1", files)



hits <- c()
for (cohort in cohorts) {
  cols <- c("PCAWG_ID",   "sum_scores", "pElement",   "nPred" ,     "mu0"   ,     "sd0"   ,     "sum_nMut" ,
            "obs_stat" ,  "alpha" ,     "beta" ,      "rho" ,
            "n4var"  ,
            "pvals" ,"qvals" , "cohort")
  
  path_file <- gsub('\\{cohort\\}', cohort, '../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_{cohort}.tsv')
  hit_cohort <- get_sigHits_cohort(path_file, sig_threshold)
  hit_cohort$cohort <- cohort
  hit_cohort <- data.frame(hit_cohort)
  hit_cohort <- hit_cohort[,which(colnames(hit_cohort) %in% cols)]
  
  hits <- rbind(hits, hit_cohort)
}

ann <- fread('../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv')
hits <- left_join(hits, ann, by = c('PCAWG_ID' = 'PCAWG_IDs'))


dir.create('../extdata/output_release2.0/Benchmark/pooled_catQuartiles_withHyperMuts/SigHits/', recursive = T, showWarnings = F)
fwrite(hits, file = '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_withHyperMuts/SigHits/allSigHits_wHyper.tsv', sep = '\t')
###################### Only hypers ######################
rm(list = ls())

library(data.table)
library(dplyr)

source('benchmark/functions_benchmark.R')

path_file <- '../extdata/output_release2.0/observed/onlyHyperMuts/update_eMET/observed_withHyperUpdate_cancerSpRho.tsv'
hits <- get_sigHits_cohort(path_file, sig_threshold = .1)
dir.create('../extdata/output_release2.0/Benchmark/onlyHypers/SigHits/', recursive = T, showWarnings = F)
fwrite(hits, file = '../extdata/output_release2.0/Benchmark/onlyHypers/SigHits/allSigHits_Hyper.tsv', sep = '\t')

############################################
# rm(list = ls())
# library(ggplot2)
# library(dplyr)
# 
# 
# hits <- fread('../extdata/output_release1.0/Benchmark/pooled_betaRho_rmHyperMuts/SigHits/allSigHits_betaRho.tsv')
# hits <- hits[which(hits$qvals <= .1),]
# 
# hits_splitted <- split_cod_nonCod(hits)
# 
# 
# df <- hits_splitted$CDS
# 
# 
# 
# 
# # Convert p-values to -log10(p-value)
# df$log_pval <- -log10(df$pvals)
# 
# # Generate box plot
# ggplot(df, aes(x = cohort, y = log_pval)) +
#   geom_boxplot(outlier.shape = NA) +  # Box plot without showing outlier points
#   geom_jitter(aes(color = as.factor(in_CGC_new), label = geneSymbol), width = 0.2, size = 2) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +  # Red for CGC_new genes
#   theme_bw() +
#   labs(x = "Cancer Type", y = "-log10(p-value)", color = "In CGC_new") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_text(aes(label = geneSymbol), vjust = -0.5, size = 3, check_overlap = TRUE)  # Show gene symbols
# 
# 
# ########################################################33
# library(ggplot2)
# library(ggrepel)
# library(ggbreak)  # For breaking the y-axis
# 
# # Assuming df is your dataset
# df$log_pval <- -log10(df$pvals)
# 
# # df <- df[which(!df$cohort %in% c("Pancan-no-skin-melanoma-lymph", "Skin-Melanoma" , "Lymph-BNHL",  "Eso-AdenoCa" ))]
# 
# # Create the plot
# # ggplot(df, aes(x = cohort, y = log_pval)) +
# #   geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
# #   geom_jitter(aes(color = as.factor(in_CGC_new)), width = 0.2, size = 2) +
# #   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +  # Red for CGC_new genes
# #   
# #   # Add gene labels using ggrepel to avoid overlap
# #   geom_text_repel(aes(label = geneSymbol, color = as.factor(in_CGC_new)), 
# #                   size = 3, max.overlaps = Inf, box.padding = 0.4) +
# #   
# #   # # Break the y-axis for better visibility
# #   # scale_y_break(c(30, 180)) +  # Breaking the y-axis between 30 and 180
# #   # scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +  # Avoid cutting off labels
# #   
# #   # Theme and labels
# #   theme_bw() +
# #   labs(x = "Cancer Type", y = "-log10(p-value)", color = "In CGC_new") +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #######################################################
# 
# df <- hits_splitted$non_coding
# 
# 
# library(ggplot2)
# library(ggrepel)
# library(ggbreak)
# 
# # Assuming df contains: cohort, element_type, pvals, gene_symbol, in_CGC_new
# df$log_pval <- -log10(df$pvals)
# 
# ggplot(df, aes(x = cohort, y = log_pval)) +
#   geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
#   geom_jitter(aes(color = as.factor(in_CGC_new)), width = 0.2, size = 2) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +  # Red for CGC_new genes
#   
#   # Add gene labels using ggrepel to avoid overlap
#   geom_text_repel(aes(label = geneSymbol, color = as.factor(in_CGC_new)), 
#                   size = 3, max.overlaps = Inf, box.padding = 0.4) +
#   
#   # # Break the y-axis for better visibility
#   # scale_y_break(c(30, 180)) +  # Adjust the break point as needed
#   # scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +  # Avoid cutting off labels
#   
#   # Multi-panel plot for each element type
#   facet_wrap(~type_of_element, scales = "free_y") +  # Separate panels by element_type
#   
#   # Theme and labels
#   theme_bw() +
#   labs(x = "Cancer Type", y = "-log10(p-value)", color = "In CGC_new") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ################################################3
# 
# 
# library(ggplot2)
# library(ggrepel)
# library(dplyr)
# 
# # Assuming your dataframe is df and it has the necessary columns
# df$log_pval <- -log10(df$pvals)
# 
# # Create a custom color column based on the conditions
# df$color <- with(df, ifelse(in_CGC_new == TRUE, "red", 
#                             ifelse(in_pcawg == TRUE, "blue", 
#                                    ifelse(in_oncoKB == TRUE, "green", "black"))))
# 
# # Create the plot
# ggplot(df, aes(x = cohort, y = log_pval)) +
#   geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
#   geom_jitter(aes(color = color), width = 0.2, size = 2) +  # Color points based on new column
#   
#   # Add gene labels using ggrepel to avoid overlap
#   geom_text_repel(aes(label = geneSymbol, color = color), 
#                   size = 3, max.overlaps = Inf, box.padding = 0.4) +
#   
#   # Multi-panel plot for each element type
#   facet_wrap(~type_of_element, scales = "free_y") +  # Separate panels by element_type
#   
#   # Theme and labels
#   theme_bw() +
#   labs(x = "Cancer Type", y = "-log10(p-value)", color = "Gene Category") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   
#   # Color scale for points and labels
#   scale_color_identity()
