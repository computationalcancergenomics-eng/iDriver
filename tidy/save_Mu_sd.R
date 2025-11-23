rm(list = ls())
library(data.table)
library(dplyr)
library(nlme)

source('tidy/functions_iDriver.R')

#################################### element-type Mu_sd ############################################
retriev_element_specific_scoreParams <- function(df, excludeDrivers){
  elements <- c("gc19_pc.cds", "enhancers", "gc19_pc.3utr",
                "gc19_pc.5utr", "gc19_pc.promCore", "gc19_pc.ss",
                "lncrna.ncrna", "lncrna.promCore")
  
  if (excludeDrivers) {
    dr <- fread('../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv')
    
    nonDrivers <- unique(dr[!(dr$in_CGC | dr$in_CGC_literature | dr$in_CGC_new | dr$in_oncoKB | dr$in_pcawg), ]$PCAWG_IDs)
    df <- data.frame(df)
    df <- df[which(df$PCAWG_ID %in% nonDrivers),]
  }
  

  params <- c()
  for (i in 1:length(elements)) {

    elem <- df[which(as.character(df$GenomicElement) == elements[i]),]
    scores <- elem$score
    elem_scores <- as.numeric(unlist(scores))
    Mu_0 <- mean(elem_scores, na.rm=T)
    sd_0 <- sd(elem_scores, na.rm=T)

    param <- c(Mu_0, sd_0)
    params <- rbind(params, param)
  }

  params <- data.frame(params)
  colnames(params) <- c('Mu', 'sd')
  params$GenomicElement <- elements
  rownames(params) <- c()

  params
}

# retriev_element_specific_scoreParams <- function(df){ 
#   elements <- c("gc19_pc.cds", "enhancers", "gc19_pc.3utr",
#                 "gc19_pc.5utr", "gc19_pc.promCore",
#                 "gc19_pc.ss", "lncrna.ncrna", "lncrna.promCore"
#                 )
#   GenomicElements <- c()
#   params <- c()
#   for (i in 1:length(elements)) {
#     print(elements[i])
#     elem <- df[which(as.character(df$GenomicElement) == elements[i]),]
#     
#     
#     elem$score2 <- as.numeric(elem$score)
#     elem <- elem[which(!is.na(elem$score2)),]
#     elem <- as.data.table(elem)
#     
#     # Filter PCAWG_IDs that appear at least twice
#     elem_filtered <- elem[, if (.N >= 5) .SD, by = PCAWG_ID]
#     
#     # # Check variance of score2 by PCAWG_ID
#     # var_by_group <- elem_filtered[, .(var_score = var(score2, na.rm = TRUE)), by = PCAWG_ID]
#     # idx = which(var_by_group$var_score != 0 )
#     # elem_filtered = elem_filtered[idx,]
#     
#     if (nrow(elem_filtered) != 0) {
#       
#       model <- lme(score2 ~ 1, random = ~1 | PCAWG_ID, data = elem_filtered)
#       
#       mu <- fixed.effects(model)
#       tau2 <- as.numeric(VarCorr(model)[1,1])  # might return both random and residual
#       sigma2 <- model$sigma^2
#       
#       param <- c(mu, tau2, sigma2)
#       params <- rbind(params, param)
#       GenomicElements = c(GenomicElements, elements[i])
#     }
#     
#   }
#   
#   params <- data.frame(params)
#   colnames(params) <- c('mu', 'tau2', 'sigma2')
#   params$GenomicElement <- GenomicElements
#   rownames(params) <- c()
#   
#   params  
# }


calc_elemType_observedScore_params <- function(study, path_observed_scores,
                                               cohort, score_method,
                                               exclude_hyper_mutated = FALSE){
  
  if(study == "simulated") {
    path_MutData <- "../extdata/procInput/iDriverInputs/simulated/Sanger/mutData.tsv"
    path_donorInfo <- '../extdata/procInput/iDriverInputs/simulated/Sanger/donorInfo.tsv'
    save_dir = '../extdata/procInput/Mu_sd/elemType/simulated/'
    excludeDrivers = F
  } else if(study == "observed") {
    path_MutData <- "../extdata/procInput/iDriverInputs/mutData.tsv"
    path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
    save_dir = '../extdata/procInput/Mu_sd/elemType/observed/'
    excludeDrivers = T
  } else if(study == "subsamples") {
    path_MutData <- "../extdata/procInput/iDriverInputs/mutData.tsv"
    path_donorInfo <- '../extdata/procInput/iDriverInputs/subsamples_dInfo.tsv'
    save_dir = '../extdata/procInput/Mu_sd/elemType/subsamples69/'
    excludeDrivers = T
  }
  
  exclude_lymph_melanoma <- ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
  
  donorInfo <- select_cohort(path_donorInfo,
                             cohort,
                             exclude_lymph_melanoma,
                             exclude_hyper_mutated)
  
  mutData = create_per_cohort_MutData(donorInfo, path_MutData, path_observed_scores)
  
  elemTypeParams <- retriev_element_specific_scoreParams(mutData, excludeDrivers)
  
  x <- left_join(mutData, elemTypeParams, by = 'GenomicElement')
  x <- x[,c('PCAWG_ID', 'Mu', 'sd')]
  x= x[!duplicated(x),]
  
  dir.create(save_dir, showWarnings = F, recursive = T)
  path_save <- paste0(save_dir, cohort, '_', score_method,'_elemTypeParams.tsv')
  
  fwrite(x, file = path_save, sep = '\t')
}


cohorts = c( 
  'Pancan-no-skin-melanoma-lymph',
  'Pan_Cancer', "Liver-HCC",
  "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
  "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
  "Panc-Endocrine", "Head-SCC",
  "Breast-AdenoCa", "Biliary-AdenoCA", "Eso-AdenoCa", "CNS-GBM",
  "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
  "Breast-LobularCa",
  "Thy-AdenoCA", "Myeloid-MPN", "Bone-Leiomyo",    
  "Lymph-NOS", "CNS-Medullo", "Myeloid-AML", "CNS-Oligo", "Cervix-SCC",
  "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Epith", "Bone-Osteosarc",
  "Cervix-AdenoCA", "Breast-DCIS", "Bone-Cart",
  "Stomach-AdenoCA", "Skin-Melanoma" #error convergence
  # , "Myeloid-MDS" all elems just have 1 muts
  )

path_observed_scores = '../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv'
score_method <- 'CADD6'
# path_observed_scores = '../extdata/procInput/iDriverInputs/cadd_scores.tsv'


for (cohort in cohorts) {
  print(paste0('Calculating element-type Params bsed on ', score_method,'scores for ', cohort))

  calc_elemType_observedScore_params(study = 'simulated',
                                     path_observed_scores,
                                     cohort, score_method)

  print('=====================================================')
}

############################################# element-type Mu_sd (SNV-nonSNV) ############################################
calc_elemType_observedScore_params_SNVnonSNV <- function(study, var_type, path_observed_scores,
                                               cohort, score_method,
                                               exclude_hyper_mutated = FALSE){
  
  # if(study == "simulated") {
  #   path_MutData <- "../extdata/procInput/iDriverInputs/simulated/Sanger/mutData.tsv"
  #   path_donorInfo <- '../extdata/procInput/iDriverInputs/simulated/Sanger/donorInfo.tsv'
  #   save_dir = '../extdata/procInput/Mu_sd/elemType/simulated/'
  #   excludeDrivers = F
  # } else {
  if(study == "observed"){
    path_MutData <- paste0("../extdata/procInput/iDriverInputs/splittedSNV_nonSNVs/observed/mutData/", var_type, "s_mutData.tsv")
    path_donorInfo <- paste0("../extdata/procInput/iDriverInputs/splittedSNV_nonSNVs/observed/donorInfo/obs_", var_type, "s_donorInfo.tsv")
    save_dir = paste0('../extdata/procInput/iDriverInputs/splittedSNV_nonSNVs/observed/Mu_sd/', score_method, '/', var_type, '/element_type/')
    excludeDrivers = T
  }
  
  exclude_lymph_melanoma <- ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
  
  donorInfo <- select_cohort(path_donorInfo,
                             cohort,
                             exclude_lymph_melanoma,
                             exclude_hyper_mutated)
  
  mutData = as.data.frame(create_per_cohort_MutData(donorInfo, path_MutData, path_observed_scores))
  
  elemTypeParams <- retriev_element_specific_scoreParams(mutData, excludeDrivers)
  
  x <- left_join(mutData, elemTypeParams, by = 'GenomicElement')
  x <- x[,c('PCAWG_ID', 'Mu', 'sd')]
  x= x[!duplicated(x),]
  
  dir.create(save_dir, showWarnings = F, recursive = T)
  path_save <- paste0(save_dir, cohort, '.tsv')
  
  fwrite(x, file = path_save, sep = '\t')
}


path_observed_scores = c('../extdata/procInput/iDriverInputs/cadd_scores.tsv')

score_method <- c('CADD6')
for(i in 1:length(score_method)){
  for (var_type in c('SNV', 'nonSNV')) {
    for (cohort in cohorts) {
      
      print(paste0('Calculating element-type Mu-Sd bsed on ', score_method[i],' scores for ', cohort))
      
      calc_elemType_observedScore_params_SNVnonSNV(study = 'observed', var_type,
                                         path_observed_scores[i],
                                         cohort, score_method[i], exclude_hyper_mutated = T)
      
      print('=====================================================')
    }
  }
}



# ########################################### element Mu-sd (random scores) ##############################################
# rm(list = ls())
# 
# library(data.table)
# library(seqminer)
# library(doParallel)
# library(foreach)
# 
# generate_random_numbers <- function(start, end, count) {
#   return(sample(start:end, count, replace=TRUE))
# }
# 
# 
# retrieve_perElem_random_mutScores <- function(elemID, path_all_cadd_scores, 
#                                               pcawg_bed6, nMut = 1000){
#   
#   elem <- pcawg_bed6[pcawg_bed6$V4 == elemID, ]
#   
#   
#   df <- data.frame(
#     'chr' = elem$V1,
#     'start' = elem$V2,
#     'end' = elem$V3)
#   
#   
#   colnames(df) <- c('chr', 'start', 'end')
#   
#   # Generate random numbers for each range in the dataframe
#   all_random_numbers <- c()
#   for (i in 1:nrow(df)) {
#     start <- df$start[i]
#     end <- df$end[i]
#     if (start == end) {
#       random_numbers = start
#     } else {
#       if (nMut < nrow(df)) {
#         count = 1
#       } else {
#         count <- nMut %/% nrow(df)  # Ensure total count is distributed among ranges
#       }
#       
#       random_numbers <- generate_random_numbers(start, end, count)
#     }
#     
#     
#     all_random_numbers <- c(all_random_numbers, random_numbers)
#   }
#   
#   pos <- all_random_numbers[!is.na(all_random_numbers)]
#   
#   
#   vcf_data <- data.frame('chr' = rep(unique(df$chr), length(pos)),
#                          'pos' = pos)
#   
#   result <- data.frame(chr = vcf_data$chr, pos = vcf_data$pos)
#   
#   result$tbix <- paste0(gsub("chr", "", result$chr), ":", 
#                         result$pos,  
#                         "-", 
#                         result$pos,
#                         sep = "")
#   result <- result[!grepl("NA", result$tbix),]
#   result2 = tabix.read.table(path_all_cadd_scores, result$tbix )
#   # elem_scores <- result2[sample(1:nrow(result2), nMut, replace=FALSE),]
#   elem_scores = result2
#   scores <- paste0(elem_scores$RawScore, collapse = ",")
#   
#   df <- data.frame('PCAWG_ID' = elemID,
#                    'nMut' = nrow(result),
#                    'Mu' = mean(elem_scores$RawScore),
#                    'sd' = sd(elem_scores$RawScore),
#                    'scores' = scores)
#   
#   df
# }
# 
# 
# path_all_cadd_scores <- "../extdata/rawInput/cadd/whole_genome_SNVs.tsv.gz"
# set.seed(123)
# nMutation = 50
# 
# pcawg_bed6 <- fread('../../make_features/external/database/bins/proccessed/PCAWG_callable.bed6')
# pcawg_bed6$V2 <- pcawg_bed6$V2 + 1
# 
# path_save <- '../extdata/procInput/CADD_scores/possible_muts_perElem/'
# dir.create(path_save, showWarnings = F, recursive = T)
# elemIDs <- unique(pcawg_bed6$V4)
# 
# # Define the chunk size
# chunk_size <- 5000
# 
# # Get the total number of chunks
# num_chunks <- ceiling(length(elemIDs) / chunk_size)
# 
# # Loop over each chunk
# for (chunk_index in 1:num_chunks) {
#   # Determine the start and end indices for the current chunk
#   start_index <- (chunk_index - 1) * chunk_size + 1
#   end_index <- min(chunk_index * chunk_size, length(elemIDs))
#   
#   # Extract the current chunk of elemIDs
#   current_elemIDs <- elemIDs[start_index:end_index]
#   
#   # Define the filename for the current chunk
#   filename <- paste0(path_save, "chunk_", chunk_index, "_results.csv")
#   
#   if (!file.exists(filename)) {
#     # Define the foreach loop to run in parallel for the current chunk
#     registerDoParallel(cores = 30) 
#     results <- foreach(i = 1:length(current_elemIDs), .combine = rbind) %dopar% {
#       elemID <- current_elemIDs[i]
#       
#       s <- retrieve_perElem_random_mutScores(elemID, path_all_cadd_scores, pcawg_bed6, nMut = nMutation)
#       if ((nrow(s) != 1 | ncol(s) != 5)) {
#         print(dim(s))
#         print(elemID)
#       }
#       if (i %% 100 == 0) {
#         print(paste("Processed", i, "elements"))
#       }
#       s
#     }
#     
#     # Filter out NULL results (for empty elements)
#     chunk_results <- results[!is.null(results), ]
#     
#     # Clean up parallel processing
#     stopImplicitCluster()
#     
#     # Save the results for the current chunk to a file
#     write.csv(chunk_results, file = filename, row.names = FALSE)
#   } else {
#     # Print a message if the file already exists
#     cat("Chunk", chunk_index, "already exists - Skipping processing.\n")
#   }
#   
# }
# 
# 
# length(elemIDs)
# all_scores <- do.call(rbind, lapply( list.files(path_save, full.names = T), fread))
# length(unique(all_scores$PCAWG_ID))
# fwrite(all_scores[, 1:4], file = '../extdata/procInput/CADD_scores/randomScores50base.tsv', sep = '\t')
# 
# ########################################### element context Mu-sd #########################
# rm(list = ls())
# elem_type = 'gc19_pc.promCore'
# 
# load(paste0('../extdata/procInput/CADD_context/all_192Cntxt/', elem_type, '_192CADDcntxt.RData'))
# 
# unified_CADD_contxt <- CADD_context
# rm(list = c('CADD_context'))
# 
# # Load necessary packages
# library(dplyr)
# 
# # Define the computeMuSD function
# computeMuSD <- function(context) {
#   f = context$freq / sum(context$freq)
#   mu = sum(f * context$mu_raw)
#   moment2 = context$mu_raw^2 + context$var_raw
#   indexes = which(!is.na(context$var_raw))
#   f2 = context$freq[indexes] / sum(context$freq[indexes])
#   moment2 = moment2[indexes]
#   sd = sqrt(sum(f2 * moment2) - mu^2)
#   list(mu = mu, sd = sd)
# }
# 
# # Initialize an empty list to store the results
# 
# 
# # Assuming unified_CADD_contxt is the combined dataframe from all CADD_contxt files
# # Get unique PCAWG_IDs
# unique_PCAWG_IDs <- unique(unified_CADD_contxt$PCAWG_ID)
# print(length(unique_PCAWG_IDs))
# 
# t0 <- Sys.time()
# results <- list()
# # Iterate over each PCAWG_ID and perform the calculations
# for (i in 1: length(unique_PCAWG_IDs)) {
#   if (i %% 200 == 0) {
#     print(paste0(i, ' out of ', length(unique_PCAWG_IDs)))
#   }
#   elem_name= unique_PCAWG_IDs[i]
#   possibleCntxt <- unified_CADD_contxt[which(unified_CADD_contxt$PCAWG_ID == elem_name),]
#   
#   # Calculate mean, variance, and frequency for each triCntxt in possibleCntxt
#   stats <- possibleCntxt %>%
#     group_by(triCntxt) %>%
#     summarise(
#       mu_raw = mean(RawScore, na.rm = TRUE),
#       var_raw = var(RawScore, na.rm = TRUE),
#       freq = n()
#     )
#   
#   # Compute Mu and SD using the computeMuSD function
#   stat_elem <- computeMuSD(stats)
#   
#   # Store the results in the list
#   results[[elem_name]] <- data.frame(
#     PCAWG_ID = elem_name,
#     Mu = stat_elem$mu,
#     sd = stat_elem$sd
#   )
# }
# 
# 
# print(Sys.time() - t0)
# 
# # Combine all results into a single dataframe
# final_result <- do.call(rbind, results)
# 
# row.names(final_result) <- c()
# 
# # Save the result to a new RData file
# dir.create('../extdata/procInput/Mu_sd/tmp_perElem/', recursive = T, showWarnings = F)
# fwrite(final_result, file = paste0('../extdata/procInput/Mu_sd/tmp_perElem/', elem_type, 'mu_sd.tsv'), sep = '\t')
# 
# # load('../extdata/procInput/CADD_context/all_192Cntxt/gc19_pc.promCore_192CADDcntxt.RData')
# # 
# # unified_CADD_contxt <- CADD_context
# 
# 
# # possibleCntxt <- CADD_context[which(CADD_context$PCAWG_ID == elem_name),]
# # 
# # # Calculate mean and standard deviation for each triCntx in possibleCntxt
# # stats <- possibleCntxt %>%
# #   group_by(triCntxt) %>%
# #   summarise(
# #     mu_raw = mean(RawScore, na.rm = TRUE),
# #     var_raw = var(RawScore, na.rm = TRUE),
# #     freq = n()
# #   )
# # 
# # computeMuSD<-function(context) {
# #   f = context$freq/ sum(context$freq)
# #   mu = sum(f * context$mu_raw)
# #   moment2 = context$mu_raw^2 + context$var_raw
# #   indexes = which(!is.na(context$var_raw))
# #   f2 = context$freq[indexes]/ sum(context$freq[indexes])
# #   moment2 = moment2[indexes]
# #   sd = sqrt(sum(f2*moment2) - mu^2)
# #   list(mu = mu, sd=sd)
# # }
# # 
# # 
# # stat_elem <- computeMuSD(stats)
# # 
# # 
# # data.frame(cbind('PCAWG_ID'= elem_name,
# #            'Mu' = stat_elem$mu,
# #            'sd' = stat_elem$sd))
# 
# ######################## element-type context Mu_sd ##################
# elems <- c( "gc19_pc.cds", "gc19_pc.promCore", "enhancers", "gc19_pc.3utr", "gc19_pc.5utr",
#             #"lncrna.ncrna" , "lncrna.promCore",
#             "gc19_pc.ss" )
# 
# # Initialize an empty list to store the loaded data
# data_list <- list()
# 
# # Load each RData file and store the CADD_contxt object in the list
# for (elem_type in elems) {
#   file = paste0('../extdata/procInput/CADD_context/all_192Cntxt/', elem_type,'_192CADDcntxt.RData')
#   load(file)
#   
#   stats <- CADD_context %>%
#     group_by(triCntxt) %>%
#     summarise(
#       mu_raw = mean(RawScore, na.rm = TRUE),
#       var_raw = var(RawScore, na.rm = TRUE),
#       freq = n()
#     )
#   
#   # Compute Mu and SD using the computeMuSD function
#   stat_elemType <- computeMuSD(stats)
#   
#   
#   data_list[[elem_type]] <- data.frame(
#     GenomicElement = elem_type,
#     Mu = stat_elemType$mu,
#     sd = stat_elemType$sd
#   )
# }
# 
# # Combine all data frames in the list into a single data frame
# unified_CADD_contxt <- do.call(rbind, data_list)
# 
