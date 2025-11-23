
create_mut_df <- function(path_raw_tcga, path_raw_icgc, path_save){
  
  fname <- paste0(path_save, "mutationFile_all_rmDups.tsv")
  
  if(sum(file.exists(fname)) != length(fname)){
    
    maf_tcga <- fread(path_raw_tcga)    #dim >>> 29504368       43
    maf_icgc <- fread(path_raw_icgc)    #dim >>>  23159591       43
    
    maf <- rbind(maf_tcga, maf_icgc)  # dim >>>  52663959       43
    
    df=maf[,c( "Chromosome",  "Start_position", "End_position",
               "Reference_Allele", "Tumor_Seq_Allele2",
               "Donor_ID", "ref_context", "Variant_Type",
               "Project_Code", "Variant_Classification")]   # >>> dim: 52663959  9 without Variant_Classification
    
    fwrite(df, file = fname, sep = "\t")
  } else {
    print("mutation file exists")
  }
}

add_CADD_2MAF <- function(path_SNV_offLine, path_indel_offLine, 
                            path_onLineSet, path_complete_MAF){
  
  all_scores <- aggregate_Scores(path_SNV_offLine, path_indel_offLine, path_onLineSet)
  all_scores <- all_scores[,-c("ID")]
  all_scores <- all_scores[!duplicated(all_scores),]
  
  all_MAF <- fread(path_complete_MAF)
  
  all_MAF <- define_cadd_vcf_cols(all_MAF)
  MAF_annotated=left_join(all_MAF, all_scores, by = c("Chrom", "Pos", "Ref", "Alt"))
  MAF_annotated
}

retain_included_Donors <- function(reannotated_MAF, path_pcawgSupp){
  
  info <- read_xlsx(path_pcawgSupp)
  included <- as.data.frame(info[which(info$donor_wgs_included_excluded == "Included"),
                                 "icgc_donor_id"], row.names = c())
  reannotated_MAF <- reannotated_MAF[which(reannotated_MAF$Donor_ID %in% included[,1]),]
  reannotated_MAF
}


add_ranked_scores <- function(df){
  
  df <- df[order(as.numeric(df$RawScore)),]
  
  Rank_portion <- (1:nrow(df))/nrow(df)
  df$Rank_portion <- ifelse(is.na(df$RawScore), NA, Rank_portion)
  df
  
}

################## create GRanges ###############
########## MAF_Obj2GRanges function
# Arguments: @CADDannotated_maf
MutFile2GRanges <- function(CADDannotated_maf){
  
  gr <- GenomicRanges::GRanges(seqnames = paste0("chr", CADDannotated_maf$Chrom), 
                               #strand = CADDannotated_maf$Strand,
                               IRanges::IRanges(start = CADDannotated_maf$Start_position,
                                                end = CADDannotated_maf$End_position),
                               ref_al = CADDannotated_maf$Reference_Allele,
                               alt_al = CADDannotated_maf$Tumor_Seq_Allele2,
                               D_id = CADDannotated_maf$Donor_ID,
                               cohort = CADDannotated_maf$Project_Code, 
                               Ref_cntxt = CADDannotated_maf$ref_context,
                               var_type = CADDannotated_maf$Variant_Type, 
                               variant_class = CADDannotated_maf$Variant_Classification,
                               CADD_RawScore = CADDannotated_maf$RawScore,
                               CADD_PHRED = CADDannotated_maf$PHRED,
                               CADD_ranked = CADDannotated_maf$Rank_portion)
  
  strand(gr) <- Rle(rep("+", nrow(CADDannotated_maf)))
  
  gr
}


filter_samples <- function(path_pcawgSupp, gr,
                           exclude_lymph_melanoma = TRUE,
                           exclude_hyper_mutated = TRUE){
  
  info <- read_xlsx(path_pcawgSupp)
  
  included <- as.data.frame(info[which(info$donor_wgs_included_excluded == "Included"),
                                 "icgc_donor_id"], row.names = c())
  
  gr <- gr[which(gr$D_id %in% included[,1])]
  
  if (exclude_lymph_melanoma) {
    exceptions <- c("Skin-Melanoma", "Lymph-NOS", "Lymph-CLL", "Lymph-BNHL")
    gr <- gr[-which(gr$cohort %in% exceptions)]
  }
  
  if (exclude_hyper_mutated) {
    df <- data.frame(table(gr$cohort, gr$D_id))
    colnames(df) <- c("cohort", "D_id", "freq")
    df <- df[df$freq !=0,]
    
    nonHyperMut_donors <- unique(df[-which(df$freq/3000 >30), "D_id"])
    
    gr <- gr[which(gr$D_id %in% nonHyperMut_donors)] 
  }
  
  gr
}