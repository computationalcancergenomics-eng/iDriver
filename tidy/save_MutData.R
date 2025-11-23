rm(list = ls())
map_muts <- function(gr, path_to_GEs, path_to_callable, element,callable_GE_train, 
                     count_blocks=TRUE, test=NULL, train=NULL){
  
  if (element == "test") {
    
    load(path_to_GEs)
    
    
    testGE_gr <- unlist(testGE)
    
    g <- names(testGE_gr)
    s <- strsplit(g, "[::]")
    GenomicElement <- unlist(lapply(s, function(x){x[1]}))
    
    lenElement <- width(testGE_gr)
    
    mcols(testGE_gr) <- DataFrame(mcols(testGE_gr), GenomicElement, lenElement)
    
    callable_GE <- relist(testGE_gr, testGE)
    lenElement_new <- sum(width(callable_GE))
    
  } else if (element == "train") {
    
    #callable_GE <- make_callable_trainGE(path_to_GEs, path_to_callable)
    callable_GE <- callable_GE_train
    lenElement_new <- sum(width(callable_GE))
  }
  
  
  
  ov <- findOverlaps(gr, callable_GE, ignore.strand = TRUE)
  
  # if we run this function for each element, length(unique(queryHits(ov))) is nMut and
  # length(unique(subjectHits(ov))) is number of mutated blocks. >>> ov <- findOverlaps(gr, testGE$`gc19_pc.cds::gencode::ARID1A::ENSG00000117713.13`, ignore.strand = TRUE)
  
  idx_gr <- queryHits(ov)
  idx_callable_GE <- subjectHits(ov)
  
  gr_mappedGE <- gr[idx_gr]
  
  GE_Mutated <- callable_GE[idx_callable_GE]  
  
  if (element == "test") {
    mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE), 
                                    GenomicElement=unique(mcols(GE_Mutated, level = "within")[,"GenomicElement"]), 
                                    elemenntLength = sum(mcols(GE_Mutated, level = "within")[,"lenElement"]),
                                    name = names(GE_Mutated))
  } else if (element == "train") {
    mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE), 
                                    name = names(GE_Mutated))
  }
  
  idx_completely_black <- which(lenElement_new == 0)
  
  
  
  if(length(idx_completely_black) != 0){
    lenElement_new <- lenElement_new[-idx_completely_black]
    callable_GE <- callable_GE[-idx_completely_black]
  }
  
  if (count_blocks) {
    
    GE_blocks <- unlist(callable_GE)
    
    ov_block <- findOverlaps(gr, GE_blocks, ignore.strand = TRUE)
    idx_gr_block <- queryHits(ov_block)
    idx_GE_block <- subjectHits(ov_block)
    
    gr_blocks_mappedGE <- gr[idx_gr_block]
    GE_blocks_mutated <- GE_blocks[idx_GE_block]
    
    mcols(gr_blocks_mappedGE) <- DataFrame(mcols(gr_blocks_mappedGE), 
                                           binID=names(GE_blocks_mutated),
                                           block_ID=paste0(names(GE_blocks_mutated), "**",1:length(gr_blocks_mappedGE)),
                                           block_length=width(GE_blocks_mutated))
    
    list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE, 
         lenAllTests = lenElement_new, gr_blocks = gr_blocks_mappedGE)
  } else {
    list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE, 
         lenAllTests = lenElement_new)
  }
  
}

get_ref_context <- function(genome, chr, position, context_range = 10){
  getSeq(genome, chr, start = position - context_range, end = position + context_range)
}

replace_second_char <- function(original_string, new_char, var_type) {
  if (!var_type %in% c('SNP', 'SNV')) {
    return(original_string)  # Return original if string is indel, mnv, ...
  }
  first_part <- substr(original_string, 1, 1)
  last_part <- substr(original_string, 3, nchar(original_string))
  new_string <- paste0(first_part, new_char, last_part)
  return(new_string)
}


generate_cohort_DonorID <- function(gr){
  
  info <- data.frame(table(gr$cohort, gr$D_id))
  colnames(info) <- c("cohort", "D_id", "freq")
  info <- info[info$freq !=0,]
  
  # Grouping by 'D_id' and aggregating the cohorts
  grouped_df <- info %>%
    group_by(D_id) %>%
    summarise(
      cohorts = paste(unique(cohort), collapse = ", "),
      freq = sum(freq)
    ) %>%
    ungroup()
  
  # Splitting the cohorts into two columns
  grouped_df <- grouped_df %>%
    separate(cohorts, into = c("cohort1", "cohort2"), sep = ", ", fill = "right", extra = "merge")
  
  grouped_df <- grouped_df[, c("D_id", "cohort1", "cohort2")]
  grouped_df
}



generate_Donor_samplesInfo <- function(path_sample_sheet, path_mapDonInfo,
                                       path_mapSampInfo){
  
  sampleInfo <- fread(path_mapSampInfo)
  pcawg_supp <- as.data.frame(fread(path_sample_sheet))
  pcawg_supp <- pcawg_supp[,c("icgc_donor_id", "aliquot_id", 
                              "dcc_specimen_type", "donor_wgs_exclusion_white_gray")]
  pcawg_supp <- pcawg_supp[which(pcawg_supp$aliquot_id %in% sampleInfo$D_id),]
  pcawg_supp <- pcawg_supp[!(duplicated(pcawg_supp)),]
  
  donorInfo2 <- left_join(sampleInfo, pcawg_supp, by = c('D_id' = 'aliquot_id'))
  
  donorInfo_pcawgObs <- fread(path_mapDonInfo)
  
  donorInfo3 <- left_join(donorInfo2, donorInfo_pcawgObs, by = c('icgc_donor_id' = 'D_id'))
  
  donorInfo3 <- donorInfo3[which(donorInfo3$donor_wgs_exclusion_white_gray == "Whitelist"),]
  unique_donorInfo3 <- donorInfo3[!duplicated(donorInfo3$icgc_donor_id),]
  
  unique_donorInfo3 <- unique_donorInfo3[, c("D_id",  "cohort1.x", 
                                             "icgc_donor_id", 
                                             "cohort1.y", "cohort2.y")]
  colnames(unique_donorInfo3) <- c("sampleID",  "project_code", 
                                   "D_id", 
                                   "cohort1", "cohort2")
  unique_donorInfo3
}


prepareMutDataSanger <- function(path_gr, path_to_GEs, path_to_callable, 
                                 path_Donor_samplesInfo){
  
  load(path_gr)
  info <- fread(path_Donor_samplesInfo) 
  
  gr <- gr[which(gr$D_id %in% info$sampleID)]
  
  gr$var_ID <- paste0('M_', 1:length(gr))
  
  all_elements <- map_muts(gr, path_to_GEs, path_to_callable, element = 'test',
                           callable_GE_train = '', count_blocks = F)
  
  obs_Muts <- all_elements$MAF_GE_mapped
  
  genome <- BSgenome.Hsapiens.UCSC.hg19
  
  if (is.null(obs_Muts$Ref)) {
    
    obs_Muts$Ref <- as.character(get_ref_context(genome, as.character(seqnames(obs_Muts)), start(obs_Muts), context_range = 1))
    # names(obs_Muts$Ref) <- c()
  } else {
    obs_Muts$Ref <- toupper(substr(obs_Muts$Ref_cntxt, 10, 12))
  }
  
  obs_Muts$Alt <- mapply(replace_second_char, obs_Muts$Ref, obs_Muts$alt_al, obs_Muts$var_type)
  
  obs_Muts$Ref <- ifelse(!obs_Muts$var_type %in% c('SNP', 'SNV'), NA, obs_Muts$Ref)
  obs_Muts$Alt <- ifelse(!obs_Muts$var_type %in% c('SNP', 'SNV'), NA, obs_Muts$Alt)
  
  obs_Muts$triCntx <- paste0(obs_Muts$Ref, ' > ', obs_Muts$Alt)
  
  df <- as.data.frame(DataFrame(obs_Muts))
  
  
  mutData <- df[,c('D_id' , 'name', 'Ref','Alt', 'triCntx', 'var_type', 
                   'GenomicElement', 'X.start', "var_ID", "elemenntLength" )]
  
  mutData2 <- left_join(mutData, info, by = c('D_id' = 'sampleID'))
  
  colnames(mutData2) <- c("sampleID",  "PCAWG_ID", 'Ref','Alt', 'triCntx', 
                          'var_type', 'GenomicElement', 'position',
                          "var_ID", "elemenntLength", "project_code", "donor_id",
                          "cohort1", "cohort2")
  
  mutData2 <- mutData2[,c("sampleID",  "PCAWG_ID", 'Ref','Alt', 'triCntx', 
                          'var_type', 'GenomicElement', 'position',
                          "var_ID", "elemenntLength")]
  
  donor_totMutTab <- data.frame(table(gr$D_id))
  colnames(donor_totMutTab) <- c("sampleID", "donor_totMut")
  
  info2 <- data.frame(table(gr$cohort, gr$D_id))
  info2 <- info2[info2$Freq !=0,]
  colnames(info2) <- c("project_code", "sampleID", "freq")
  
  info <- left_join(info2, info, by = c("project_code", "sampleID"))
  
  info$HyperMut_donor <- (info$freq/3000 >30)
  
  cadd_scores <- df[,c("var_ID", "CADD_RawScore")]
  colnames(cadd_scores) <- c('var_ID', 'score')
  
  list('mutData' = mutData2,
       'donorInfo' = info,
       'CADD_scores' = cadd_scores)
}



################################## Run #######################################


library(data.table)
library(rtracklayer)
library(readxl)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(tidyr)

# # 1-save the information of cohorts and D_ids based on observed PCAWG data
# # and sampleIDs and project_codes based on simulated data
# 
# load('../extdata/procInput/mut_PCAWG/all_samples.RData')
# info <- generate_cohort_DonorID(gr)
# fwrite(info, '../extdata/procInput/map_cohort_donID.tsv', sep = '\t')
# 
# load('../extdata/procInput/simulated/sanger/all_sanger_gr_removedDups.RData')
# info2 <- generate_cohort_DonorID(gr)
# fwrite(info2, '../extdata/procInput/map_ProjCode_sampID.tsv', sep = '\t')
# 
# 
# path_mapDonInfo <- '../extdata/procInput/map_cohort_donID.tsv'
# path_mapSampInfo <- '../extdata/procInput/map_ProjCode_sampID.tsv'
# path_sample_sheet <- '../extdata/rawInput/pcawg_sample_sheet.tsv'
# 
# Donor_samplesInfo <- generate_Donor_samplesInfo(path_sample_sheet, path_mapDonInfo,
#                                        path_mapSampInfo)
# fwrite(Donor_samplesInfo, '../extdata/procInput/DonorSampleInfo.tsv', sep = '\t')

######## save sanger sim input data #############
path_to_GEs = "../extdata/procInput/callable_testGE.RData"
path_to_callable = "../extdata/procInput/callable_regenerated.bed"
path_gr_sanger <- '../extdata/procInput/simulated/sanger/all_sanger_gr_removedDups.RData'
path_Donor_samplesInfo <- '../extdata/procInput/DonorSampleInfo.tsv'

mut_list <- prepareMutDataSanger(path_gr_sanger, path_to_GEs, path_to_callable,
                           path_Donor_samplesInfo)

dir.create('../extdata/procInput/iDriverInputs/simulated/Sanger/', recursive = T, showWarnings = F)
fwrite(mut_list$mutData, file = '../extdata/procInput/iDriverInputs/simulated/Sanger/mutData.tsv', sep = '\t')
fwrite(mut_list$donorInfo, file = '../extdata/procInput/iDriverInputs/simulated/Sanger/donorInfo.tsv', sep = '\t')
fwrite(mut_list$CADD_scores, file = '../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv', sep = '\t')
