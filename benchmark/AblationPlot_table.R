rm(list = ls())
library(pROC)
library(PRROC)
library(data.table)
library(dplyr)
source('benchmark/functions_benchmark.R')

based_on <- 'in_CGC_new'
path_ann_PCAWG_ID <- ('../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv')

path_save <- '../extdata/output_release2.0/Benchmark/iDriver_variations/'
sig_definition_method <- 'fdr'
sig_definition_threshold <- .1
baseMethod <- 'iDriver'
elements <- c('NC', 'CDS')

cohorts <- c( 
  "Biliary-AdenoCA",  "Bladder-TCC",  "Bone-Osteosarc",
  "Breast-AdenoCa",  "Cervix-SCC" ,  "CNS-GBM",  "CNS-Medullo", #"CNS-Oligo" ,
  "CNS-PiloAstro" , "ColoRect-AdenoCA" ,
  "Eso-AdenoCa",  "Head-SCC",
  "Kidney-RCC", "Kidney-ChRCC",  "Liver-HCC",  "Lung-AdenoCA",
  "Lung-SCC",  "Lymph-BNHL",  "Lymph-CLL","Myeloid-MPN",
  "Ovary-AdenoCA",
  "Pan_Cancer" ,
  "Panc-AdenoCA",  "Panc-Endocrine",
  "Prost-AdenoCA", "Stomach-AdenoCA",
  "Skin-Melanoma","Thy-AdenoCA",  "Uterus-AdenoCA",
  "Pancan-no-skin-melanoma-lymph" 
)

# ---- paths to new results (same for all cohorts except cohort name) ----
get_paths <- function(cohort) {
  c(
    paste0('../extdata/output_release2.0/observed/FI_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_', cohort, '.tsv'),
    paste0('../extdata/output_release2.0/observed/Burden_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_', cohort, '.tsv'),
    paste0('../extdata/output_release2.0/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/observed_', cohort, '.tsv'),
    paste0('../extdata/output_release2.0/observed/Rho0/pooled_betaRho0_lengthCat_quartilesElemType/observed_', cohort, '.tsv'),
    paste0('../extdata/output_release2.0/observed/population_level/pooled_betaRhoPancan_lengthCat_quartilesElemType/observed_', cohort, '.tsv')
  )
}

newRESULTS <- c('Only recurrence', 'Only FI', 'iDriver', 'Uncorrelated scores', 'Population-level')

# ---- loop over cohorts and elements ----
for (cohort in cohorts) {
  message("Processing cohort: ", cohort)
  if (cohort == 'Pan_Cancer') {
    tissue = 'Pancan'
  } else {
    tissue = cohort
  }
  PATHs_newResults <- get_paths(cohort)
  # path_pRes <- paste0('../extdata/procInput/PCAWG_results/', tissue, '.RData')
  
  for (element in elements) {
    save_Measures2(
      path_save, newRESULTS, PATHs_newResults,
      sig_definition_method, base_method = baseMethod,
      sig_definition_threshold, path_ann_PCAWG_ID, cohort,
      based_on, element,
      selected_methods = NA,
      compareRank_new = list(TRUE, 'iDriver')
    )
  }
}

#########################################################################################
rm(list = ls())
elements <- c('NC', 'CDS')

cohorts <- c( 
  "Biliary-AdenoCA",  "Bladder-TCC",  "Bone-Osteosarc",
  "Breast-AdenoCa",  "Cervix-SCC" ,  "CNS-GBM",  "CNS-Medullo", #"CNS-Oligo" ,
  "CNS-PiloAstro" , "ColoRect-AdenoCA" ,
  "Eso-AdenoCa",  "Head-SCC",
  "Kidney-RCC", "Kidney-ChRCC",  "Liver-HCC",  "Lung-AdenoCA",
  "Lung-SCC",  "Lymph-BNHL",  "Lymph-CLL","Myeloid-MPN",
  "Ovary-AdenoCA",
  "Pan_Cancer" ,
  "Panc-AdenoCA",  "Panc-Endocrine",
  "Prost-AdenoCA", "Stomach-AdenoCA",
  "Skin-Melanoma","Thy-AdenoCA",  "Uterus-AdenoCA",
  "Pancan-no-skin-melanoma-lymph" 
)


path_nFPs <- list('../extdata/output_release2.0/simulated/allCohorts_nFPs/BurdenFree_lengthCat_quartilesElemType_perCohort.tsv',
                  '../extdata/output_release2.0/simulated/allCohorts_nFPs/FIfree_lengthCat_quartilesElemType_perCohort.tsv',
                  '../extdata/output_release2.0/simulated/allCohorts_nFPs/PopulationLevel_lengthCat_quartilesElemType_perCohort.tsv',
                  '../extdata/output_release2.0/simulated/allCohorts_nFPs/Rho0_lengthCat_quartilesElemType_perCohort.tsv',
                  '../extdata/output_release2.0/simulated/allCohorts_nFPs/rmHyperMuts_betaRho_lengthCat_quartilesElemType_perCohort.tsv')

iDriver_variation <- c("Only FI", "Only recurrence", "Population-level", "Uncorrelated scores", "iDriver")

nFPs <- lapply(path_nFPs, fread)
names(nFPs) <- iDriver_variation


nFPs <- lapply(names(nFPs), function(nm) {
  df <- nFPs[[nm]]
  df[, method := nm]   # add method column
  return(df)
})

path_save <- '../extdata/output_release2.0/Benchmark/iDriver_variations/finalTabs/'
dir.create(path = path_save, recursive = T, showWarnings = F)
nFPs <- data.frame(rbindlist(nFPs))
path_table_raw <- '../extdata/output_release2.0/Benchmark/iDriver_variations/tables/{cohort}_table_GoldStd_basedon_in_CGC_new_{elementTypes}_fdr.csv'

for (cohort in cohorts) {
  print(cohort)
  cancer_nFPs <- nFPs[which(nFPs$cohort == cohort),]
  
  
  path_table_cancer <- gsub("\\{cohort\\}", cohort, path_table_raw)
  
  for (elementType in elements) {
    print(elementType)
    path_table <- gsub("\\{elementTypes\\}", elementType, path_table_cancer)
    df <- fread(path_table)
    df_elem <- df[which(df$method %in% c("Uncorrelated scores", "iDriver", "Population-level",
                                         "Only FI", "Only recurrence")),]
    
    
    cancer_nFPs_elem <- cancer_nFPs[, c('method', paste0('nFP_', elementType))]
    colnames(cancer_nFPs_elem) <- c('method', 'nFPs')
    df_elem <- left_join(df_elem, cancer_nFPs_elem, by = 'method')
    df_elem <- df_elem[, c("method", "nTPs", "nHits" , "nHits_row",  "nHits_desired", "PRECISIONs",
                           "Recalls", "F1",  "AUC",  "AUPR",   "n_elemnts", "nFPs")]
    
    colnames(df_elem) <- c("method", "nTPs", "nHits" , "nHits_row",  "nHits_desired", "Precision",
                           "Recall", "F1",  "AUC",  "AUPR",   "n_elemnts", "nFPs")
    
    fwrite(df_elem, paste0(path_save,
           cohort, '_', elementType, '.tsv'), sep = '\t')
  }
}


###################################################################################
rm(list = ls())
library(ggplot2)
library(reshape2)

cohorts <- c( 
  "Pan_Cancer" , "Pancan-no-skin-melanoma-lymph" ,
  "Eso-AdenoCa", "Liver-HCC",  "Skin-Melanoma",
  "Panc-AdenoCA"
)

elemType <- 'CDS'

performance_df <- paste0('../extdata/output_release2.0/Benchmark/iDriver_variations/finalTabs/{cohort}_', elemType,'.tsv')

dfs <- c()
for (cohort in cohorts) {
  df <-  fread(gsub("\\{cohort\\}", cohort, performance_df))
  df$cohort <- cohort
  dfs <- rbind(dfs, df)
}


df_long <- melt(dfs,
                id.vars = c("cohort", "method"),
                measure.vars = c("Precision", "Recall", "F1", "AUC", "AUPR"),
                variable.name = "Metric",
                value.name = "Value")


ggplot(df_long, aes(x = method, y = Value, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~cohort, strip.position = "bottom", ncol = 3) +
  labs(x = "Method", y = "Value", title = "Performance of Methods across Cancer Types") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.placement = "outside")


# Make sure method order is preserved (if needed)
df_long$method <- factor(df_long$method, 
                         levels = c("Uncorrelated scores", "iDriver", "Population-level", 
                                    "Only recurrence", "Only FI"))

# Line plot
ggplot(df_long, aes(x = method, y = Value, group = Metric, color = Metric)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~cohort, strip.position = "bottom", ncol = 3) +
  labs(x = "Method", y = "Performance Value", 
       title = "Performance of Methods across Cancer Types") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        strip.background = element_rect(color = "black", fill = "grey90", size = 0.8),
        
        strip.placement = "outside")




# Plot grouped bar plot
ggplot(df_long, aes(x = Metric, y = Value, fill = method)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~cohort, strip.position = "bottom", ncol = 3) +
  labs(x = "Metric", y = "Value") +
       # , title = "Performance Metrics of Methods across Cancer Types") +
  theme_bw(base_size = 14) +
    theme_classic() + 
  coord_flip() +
    
    
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        # legend.position = "right",
        # legend.title = element_blank(), 
        axis.line = element_line(colour = "black"),
        # strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        strip.background = element_rect(color = "black", fill = "grey90", size = 0.8),
        strip.placement = "outside")

################################## final ################################## 
rm(list = ls())
library(data.table)



prepare_data_perElemType1 <- function(path_tables, elem_type){
  files <- list.files(path_tables, full.names = T)
  
  cds_files <- files[grepl('_CDS', files)]
  NC_files <- files[grepl('_NC', files)]
  
  
  if (elem_type == 'CDS') {
    cancer_types <- cds_files[
      !grepl("Pan_Cancer|Pancan-no-skin-melanoma-lymph", cds_files)
      ]
    
    PanCancer <- cds_files[
      grepl("Pancan-no-skin-melanoma-lymph", cds_files)
      ]
  } else if (elem_type == 'NC'){
    cancer_types <- NC_files[
      !grepl("Pan_Cancer|Pancan-no-skin-melanoma-lymph", NC_files)
      ]
    
    PanCancer <-  NC_files[
      grepl("Pancan-no-skin-melanoma-lymph", NC_files)
      ]
  }
  
  
  PanCancer <- fread(PanCancer)
  
  x <- lapply(cancer_types, fread)
  df <- do.call(rbind, x)
  
  setDT(df)
  
  df_sums <- df[, .(
    nTPs = sum(nTPs, na.rm = TRUE),
    nHits = sum(nHits, na.rm = TRUE)
  ), by = method]
  
  list('Pancancer' = PanCancer, 'Cancer-specific' = df_sums)
  
}

prepare_data_perElemType2 <- function(df_sums, elem_type, cancer_category){
  df_sums <- df_sums[df_sums$method %in% c("iDriver", "Population-level", "Only recurrence", "Only FI"),]
  df_sums$notReported <- df_sums$nHits - df_sums$nTPs
  df_elem <- df_sums[,c("method", "nTPs", "notReported")]
  df_elem$elem_type <- elem_type
  df_elem$cancer_category <- cancer_category
  df_elem
  
}

path_tables <- '../extdata/output_release2.0/Benchmark/iDriver_variations/finalTabs/'


df <- c()
for (elem_type in c('CDS', 'NC')) {
  cancerList_df <- prepare_data_perElemType1(path_tables, elem_type)
  df_elem <- rbind(prepare_data_perElemType2(cancerList_df[['Cancer-specific']], elem_type, cancer_category = 'Cancer-specific'),
                  prepare_data_perElemType2(cancerList_df[['Pancancer']], elem_type, cancer_category = 'Pancancer'))
  df <- rbind(df, df_elem)
}

df$cancer_category <- ifelse(df$cancer_category == 'Pancancer', 'Pancancer*', df$cancer_category)
library(ggplot2)
library(data.table)

setDT(df)

# reshape into long format
df_long <- melt(
  df,
  id.vars = c("method", "elem_type", "cancer_category"),
  measure.vars = c("nTPs", "notReported"),
  variable.name = "type",
  value.name = "count"
)

# ensure factor ordering
df_long[, method := factor(method, levels = c("iDriver", "Population-level", "Only recurrence", "Only FI"))]
df_long[, type := factor(type, levels = c("notReported", "nTPs"))]  # nTPs first → bottom

# ggplot(df_long, aes(x = method, y = count, fill = type)) +
#   geom_bar(stat = "identity", position = "stack") +
#   facet_grid(cancer_category , scales = "free_x", space = "free_x") +
#   labs(x = "Method", y = "Count", fill = "") +
#   theme_bw(base_size = 14) +
#   theme(panel.grid = element_blank(), 
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.background = element_rect(fill = "grey90"),
#     panel.spacing = unit(1, "lines")
#   ) +
#   scale_fill_manual(
#     values = c("nTPs" = "#d7301f", "notReported" = 'grey'),
#     breaks = c("nTPs", "notReported")  # order in legend
#   )


library(ggplot2)

ggplot(df_long[df_long$elem_type == "CDS" & df_long$count > 0, ],
       aes(x = method, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4) +
  facet_grid(~cancer_category ,
             scales = "free_x", space = "free_x") +
  labs(x = "Method", y = "Count", fill = "") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    panel.spacing = unit(1, "lines")
  ) +
  scale_fill_manual(
    values = c("nTPs" ="#d7301f", #,  "#c55466"
               "notReported" = "grey"),
    breaks = c("nTPs", "notReported")
  )

ggplot(df_long[df_long$count > 0, ], 
       aes(x = method, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # <- add border
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4) +
  facet_grid(~cancer_category ,
             scales = "free_x", space = "free_x") +
  labs(x = "Method", y = "Count", fill = "") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    panel.spacing = unit(1, "lines")
  ) +
  scale_fill_manual(
    values = c("nTPs" = "#d7301f", "notReported" = "grey"),
    breaks = c("nTPs", "notReported")
  )

######################### sigHits (iDriver vs Only rec) ##############################
id <- fread('../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/allSigHits.tsv')
id <- id[!id$cohort %in% c('Pan_Cancer', 'Lymph-BNHL',  'Pancan-no-skin-melanoma-lymph'),]
rec <- fread('../extdata/output_release2.0/Benchmark/cohorts_onlyRec.tsv')
just_in_id <- id[which(!id$PCAWG_ID %in% rec$PCAWG_ID),]

hits <- just_in_id
newHits <- hits[(!(hits$in_CGC_new | hits$in_pcawg | hits$in_oncoKB )),]
drivers <- hits[((hits$in_CGC_new | hits$in_pcawg | hits$in_oncoKB )),]
notable <- drivers[(drivers$in_CGC_new  | drivers$in_oncoKB | (drivers$in_CGC_new  & drivers$in_oncoKB)),]
