###### TO DO !!!!!!!!!!!! ##################
# write functions such that if a variant is occurred in the PAR of the Y chromosome convert it to its X coordinate

create_mut_df_caddScores <- function(input_maf_path, path_save){
  
  var_types <- c("SNVs", "MNVs", "InDels")
  
  fnames <- paste0(path_save, "mutationFile_", var_types, ".tsv")
  
  if(sum(file.exists(fnames)) != length(fnames)){
    
    
    maf=fread(input_maf_path)
    
    df=maf[,c( "Chromosome",  "Start_position", "End_position",
               "Reference_Allele", "Tumor_Seq_Allele2",
               "Donor_ID", "ref_context", "Variant_Type")]
    
    SNVs <- df[which(df$Variant_Type == "SNP"),] # dim >>> 48272127        8
    InDels <- df[which((df$Variant_Type == "DEL") | 
                         (df$Variant_Type == "INS")),] # dim  3970387       8
    
    MNVs <- df[which(!((df$Variant_Type == "DEL") | 
                         (df$Variant_Type == "INS") | 
                         (df$Variant_Type == "SNP"))),] # dim  421445      8
    
    #extract unique
    SNVs <- SNVs[!duplicated(SNVs),] # dim >>>>>>  47978999        8
    InDels <- InDels[!duplicated(InDels),] # dim >>>>>>>  3931623       8
    
    for (i in 1:length(var_types)) {
      fwrite(mget(var_types[i])[[1]], fnames[i], sep = "\t")
    }
  } else {
    print("mutation files exist")
  }
}


create_pcawg_vcf <- function(path_mut_pcawg, varType){
  
  pcawg_vcf <- fread(path_mut_pcawg)
  pcawg_vcf$ID <- paste0("ID_", 1:nrow(pcawg_vcf))
  
  if(varType == "SNV" | varType ==  "MNV"){
    
    Chrom <- pcawg_vcf$Chromosome
    Pos <- pcawg_vcf$Start_position
    ID <- pcawg_vcf$ID
    Ref <- pcawg_vcf$Reference_Allele
    Alt <- pcawg_vcf$Tumor_Seq_Allele2
    pcawg_vcf <- data.frame(cbind(Chrom, Pos, ID, Ref, Alt))
    
  } else if(varType == "InDel"){
    
    # del <- pcawg_vcf[which((pcawg_vcf$Variant_Type == "DEL") | (pcawg_vcf$Variant_Type == "INS") |
    #                          (pcawg_vcf$Variant_Type == "InDel") ),]
    # Chrom <- del$Chromosome
    # Pos <- (del$Start_position) -1
    # ID <- del$ID
    # Ref <- toupper(gsub("-", "", paste0(substr(del$ref_context, 10,10), del$Reference_Allele)))
    # Alt <- toupper(gsub("-", "", paste0(substr(del$ref_context, 10,10), del$Tumor_Seq_Allele2)))
    # 
    # pcawg_vcf <- data.frame(cbind(Chrom, Pos, ID, Ref, Alt))
    
    
    del <- pcawg_vcf[which(pcawg_vcf$Variant_Type == "DEL"),]
    
    gr_del <- GenomicRanges::GRanges(seqnames = paste0("chr", del$Chrom),
                                     IRanges::IRanges(start = del$Start_position -1,
                                                      width = nchar(del$Reference_Allele) +1),
                                     ref_al = del$Reference_Allele,
                                     alt_al = del$Tumor_Seq_Allele2,
                                     id = del$ID)
    DELs <- BSgenomeViews(Hsapiens, gr_del)
    
    del_chr <- as.character(seqnames(DELs))
    del_pos <- start(DELs)
    del_ref <- as.character(DNAStringSet(DELs))
    del_alt <- substr(del_ref, 1, 1) 
    del_id <- mcols(DELs)$id
    
    DEL_vcf <- data.frame(cbind(del_chr, del_pos, del_id, del_ref, del_alt))
    colnames(DEL_vcf) <- c("Chrom", "Pos", "ID", "Ref", "Alt")
    
    ins <- pcawg_vcf[which(pcawg_vcf$Variant_Type == "INS"),]
    
    gr_ins <- GenomicRanges::GRanges(seqnames = paste0("chr", ins$Chrom),
                                     IRanges::IRanges(start = ins$Start_position -1,
                                                      width = nchar(ins$Reference_Allele)),
                                     ref_al = ins$Reference_Allele,
                                     alt_al = ins$Tumor_Seq_Allele2,
                                     id = ins$ID)
    
    INSs <- BSgenomeViews(Hsapiens, gr_ins)
    
    ins_chr <- as.character(seqnames(INSs))
    ins_pos <- start(INSs)
    ins_ref <- as.character(DNAStringSet(INSs))
    ins_alt <- paste0(ins_ref, mcols(INSs)$alt_al) 
    ins_id <- mcols(INSs)$id
    
    INS_vcf <- data.frame(cbind(ins_chr, ins_pos, ins_id, ins_ref, ins_alt))
    colnames(INS_vcf) <- c("Chrom", "Pos", "ID", "Ref", "Alt")
    
    
    pcawg_vcf <- rbind(INS_vcf, DEL_vcf)
    
  } 
  
  pcawg_vcf$Pos <- as.integer(pcawg_vcf$Pos)
  pcawg_vcf
  
}


retrieve_mut_scores <- function(path_mut_pcawg, varType,
                                path_mut_cadd, path_save_Scores,
                                n_cores){
  
  dir.create(path_save_Scores, showWarnings = F, recursive = T)
  mut_pcawg <- create_pcawg_vcf(path_mut_pcawg, varType)
  
  mut_pcawg$tbix <- paste0(gsub("chr", "", mut_pcawg$Chrom), ":", 
                           mut_pcawg$Pos,  
                                     "-", 
                           mut_pcawg$Pos,
                                     sep = "")
  t <- proc.time()
  
  #var_pos_unique <- unique(mut_pcawg$tbix)
  
  n <- 1:nrow(mut_pcawg)
  chunks <- split(n, ceiling(seq_along(n)/100000))
  nRow_last_chunck <- nrow(mut_pcawg) %% 100000
  
  registerDoParallel(n_cores)
  
  tmp_dir <- paste0("../extdata/tmp/", varType, "/")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  chunk_scores <- foreach(i = 1:length(chunks), .combine = 'rbind') %dopar% {
    if(i %% 10 == 0){
      print(paste("chunk", i, "out of", length(chunks)))
    }
    chunkIndex = chunks[[i]]
    
    checking_RData <- paste0(tmp_dir, "chunck_",i, ".RData")
    
    if (!file.exists(checking_RData)) {
      
      save(i, file = checking_RData)
      output_file <- paste0(path_save_Scores, "chunck_",i,
                                  "_", varType, "scores.tsv")
      
      if (!file.exists(output_file)) {
        
        chunk = mut_pcawg[chunkIndex, ]
        chunk$Chrom = as.character(chunk$Chrom)
        chunk$Pos = as.integer(chunk$Pos)
        
        mut_scores <- tabix.read.table(path_mut_cadd, chunk[, "tbix"])
        mut_scores <- mut_scores[!(duplicated(mut_scores)),]
        
        mut_scores$Chrom = as.character(mut_scores$Chrom)
        mut_scores$Pos = as.integer(mut_scores$Pos)
        
        x <- left_join(mut_scores, chunk, by =c("Chrom", "Pos", "Ref",  "Alt"))
        pcawg_caddScores <- x[!(is.na(x$ID)),]
        pcawg_caddScores=pcawg_caddScores[,c(1:7)]
        
        # if (sum(dim(pcawg_caddScores) == c(100000, 7)) != 2 |
        #     sum(dim(pcawg_caddScores) == c(nRow_last_chunck, 7)) != 2) {
        #   
        #   stop(paste0("----- dimension error for the ", i, "th chunk!!! -----"))
        # }
        # 
        # if (sum(is.na(pcawg_caddScores)) != 0) {
        #   stop(paste0("----- NA's introduced for the ", i, "th chunk!!! -----"))
        # }
        
        fwrite(pcawg_caddScores, file = output_file, 
               sep = "\t")
        
        print(proc.time() - t) 
        print(sprintf("--------chunk %s finished--------", i))
      }
    }
  }
  
  print("********************************************************")
  print("Total time:")
  print(proc.time() - t) 
}


get_final_offLine_scores <- function(path_offLine_scores, pcawg_vcf){
  
  fnames <- list.files(path_offLine_scores, full.names = T)
  offLine_scores <- lapply(fnames, fread)
  offLine_scores <- do.call(rbind, offLine_scores)
  offLine_scores <- offLine_scores[!duplicated(offLine_scores),]
  
  x <- left_join(offLine_scores, pcawg_vcf, by =c("Chrom", "Pos", "Ref",  "Alt"))
  pcawg_caddScores <- x[which(x$ID.x == x$ID.y),]
  pcawg_caddScores
  
}

save_onLine_set <- function(path_offLine_scores, path_mnv_pcawg,
                            path_indel_pcawg, path_save_onLine_Set){
  
  indel_vcf <- create_pcawg_vcf(path_indel_pcawg, "InDel")
  
  pcawg_caddScores <- get_final_offLine_scores(path_offLine_scores, indel_vcf)
  
  
  for_onLine <- indel_vcf[which(!(indel_vcf$ID %in% pcawg_caddScores$ID.x)),]
  MNVs <- create_pcawg_vcf(path_mnv_pcawg, "MNV")
  
  MNVs$ID <- paste0(MNVs$ID, "*")
  onLine_set <- rbind(for_onLine, MNVs)
  
  fwrite(onLine_set, file = paste0(path_save_onLine_Set, "CADD_onLineSet.tsv"),
         sep = "\t")
  
}


aggregate_Scores <- function(path_SNV_offLine,
                             path_indel_offLine, path_onLineSet_final){
  
  SNVs <- lapply(list.files(path_SNV_offLine, full.names = T), fread)
  SNVs <- do.call(rbind, SNVs)
  SNVs$ID <- paste0(SNVs$ID, "_snv")
  
  if(path_indel_offLine != ""){
    
    InDels <- lapply(list.files(path_indel_offLine, full.names = T), fread)
    InDels <- do.call(rbind, InDels)
    
  } else InDels <- data.frame(cbind("Chrom"=c(), "Pos" = c(),
                                    "Ref"=c(), "Alt"=c(),
                                    "RawScore"=c(), "PHRED"=c(), 
                                    "ID" = c()))
  
  
  onLineSet <- fread(path_onLineSet_final)
  all_cadd_scores <- rbind(SNVs, InDels, onLineSet)
  
  all_cadd_scores
}

add_scores_2MAF <- function(path_SNV_offLine, path_indel_offLine, 
                            path_onLineSet, input_maf_path){
  
  all_scores <- aggregate_Scores(path_SNV_offLine, path_indel_offLine, path_onLineSet)
  all_scores <- all_scores[,-c("ID")]
  all_scores <- all_scores[!duplicated(all_scores),]
  
  all_MAF <- fread(input_maf_path)
  
  all_MAF <- define_cadd_vcf_cols(all_MAF)
  
  MAF_annotated=left_join(all_MAF, all_scores, by = c("Chrom", "Pos", "Ref", "Alt"))
  MAF_annotated
  
}


define_cadd_vcf_cols <- function(maf){
  colnames(maf) <- c("Chrom", "Start_position", "End_position",
                     "Reference_Allele", "Tumor_Seq_Allele2", "Donor_ID",
                     "ref_context", "Variant_Type", "Project_Code",
                     "Variant_Classification")
  
  maf$Pos <- ifelse(maf$Variant_Type == "INS" | maf$Variant_Type == "DEL",
                    maf$Start_position -1, maf$Start_position)
  
  maf$Ref <- ifelse(maf$Variant_Type == "INS" | maf$Variant_Type == "DEL",
                    
                    toupper(gsub("-", "", paste0(substr(maf$ref_context, 10,10), 
                                                 maf$Reference_Allele))), 
                    maf$Reference_Allele)
  
  
  maf$Alt <- ifelse(maf$Variant_Type == "INS" | maf$Variant_Type == "DEL",
                    
                    toupper(gsub("-", "", paste0(substr(maf$ref_context, 10,10), 
                                                 maf$Tumor_Seq_Allele2))), 
                    maf$Tumor_Seq_Allele2)
  
  maf
}

add_varType_sangerSim <- function(Sanger_sim){
  
  mnvs <- Sanger_sim[which(((nchar(Sanger_sim$ref)) != 1) &
                             ((nchar(Sanger_sim$mut)) != 1)),]
  
  mnvs$Variant_Type <- rep("MNV", nrow(mnvs))
  
  
  snvs <- Sanger_sim[which(((nchar(Sanger_sim$ref)) == 1) & ((nchar(Sanger_sim$mut)) == 1) &
                             !(Sanger_sim$ref == "-" | Sanger_sim$mut == "-")),]
  
  snvs$Variant_Type <- rep("SNV", nrow(snvs))
  
  Insertions <- Sanger_sim[which(Sanger_sim$ref == "-"),]
  Deletions <- Sanger_sim[which(Sanger_sim$mut == "-"),]
  
  indels <- rbind(Insertions, Deletions)
  indels$Variant_Type <- c(rep("INS", nrow(Insertions)), rep("DEL", nrow(Deletions)))
  
  SangerSim_varType <- rbind(mnvs, snvs, indels)
  
  return(SangerSim_varType)
}

define_cadd_vcf_cols2 <- function(Sanger_sim){
  
  maf <- add_varType_sangerSim(Sanger_sim)
  
  # colnames(maf) <- c("Donor_ID", "Chrom", "Pos",
  #                    "Ref", "Alt",
  #                    "Project_Code", "Variant_Type")
  # 
  
  snv_mnv <- maf[which(maf$Variant_Type == "SNV" | maf$Variant_Type == "MNV"),]
  snv_mnv <- data.frame(cbind(snv_mnv, "Chrom"= snv_mnv$chr,
                   "Pos" = snv_mnv$pos,
                   "Ref" = snv_mnv$ref,
                   "Alt" = snv_mnv$mut,
                   "Donor_ID" = snv_mnv$sampleID,
                   "Project_Code" = snv_mnv$Project_Code,
                   "Variant_Type" = snv_mnv$Variant_Type)) 

  del <- maf[which(maf$Variant_Type == "DEL"),]
    
    gr_del <- GenomicRanges::GRanges(seqnames = paste0("chr", del$chr),
                                     IRanges::IRanges(start = del$pos -1,
                                                      width = nchar(del$ref) +1),
                                     ref_al = del$ref,
                                     alt_al = del$mut)
    DELs <- BSgenomeViews(Hsapiens, gr_del)
    
    DEL_vcf <- data.frame(cbind("Chrom"= gsub("chr", "", as.character(seqnames(DELs))),
                                "Pos" = start(DELs),
                                "Ref" = as.character(DNAStringSet(DELs)),
                                "Alt" = substr(as.character(DNAStringSet(DELs)), 1, 1),
                                "Donor_ID" = del$sampleID,
                                "Project_Code" = del$Project_Code,
                                "Variant_Type" = del$Variant_Type))
    DEL_vcf <- cbind(del, DEL_vcf)
    
    ins <- maf[which(maf$Variant_Type == "INS"),]
    
    gr_ins <- GenomicRanges::GRanges(seqnames = paste0("chr", ins$chr),
                                     IRanges::IRanges(start = ins$pos -1,
                                                      width = nchar(ins$ref)),
                                     ref_al = ins$ref,
                                     alt_al = ins$mut)
    
    INSs <- BSgenomeViews(Hsapiens, gr_ins)
    
    INS_vcf <- data.frame(cbind("Chrom"= gsub("chr", "", as.character(seqnames(INSs))),
                                "Pos" = start(INSs), 
                                "Ref" = as.character(DNAStringSet(INSs)),
                                "Alt" = paste0(as.character(DNAStringSet(INSs)), mcols(INSs)$alt_al), 
                                "Donor_ID" = ins$sampleID,
                                "Project_Code" = ins$Project_Code,
                                "Variant_Type" = ins$Variant_Type))

    INS_vcf <- cbind(ins, INS_vcf)
    sim_vcf <- rbind(snv_mnv, INS_vcf, DEL_vcf)

    return(sim_vcf)
}