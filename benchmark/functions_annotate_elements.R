concate_bed_files <- function(path_genomic_intervals_list){
  GEs_bedFiles <- list.files(path_genomic_intervals_list, full.names = T)
  GEs <- lapply(GEs_bedFiles, fread)
  x <- do.call(rbind, GEs)
  x
}

define_element_type <- function(binID_vector){
  
  s <- strsplit(binID_vector, "[::]")
  GenomicElement <- unlist(lapply(s, function(x){x[1]}))
  GenomicElement
  
}



extract_cancerGenes <- function(path_to_PCAWG_supp){
  
  cancer_genes <- read.xlsx(path_to_PCAWG_supp, sheet = "Table 7 List of cancer genes")
  
  CGC_v80 <- na.omit(cancer_genes$`603.genes.from.a.curated.version.of.the.Cancer.Gene.Census.(v80)`)
  
  CGC_and_literature <- cancer_genes$`757.genes:.union.of.603.CGC.genes.and.369.genes.identified.in.cancer.exome.studies.(Lawrence.et.al.2014,.Martincorena.et.al.2017).`
  
  list("CGC_v80" = CGC_v80, "CGC_and_literature" = CGC_and_literature)
  
}


get_newCGCs <- function(path_to_PCAWG_supp, path_to_cgc_v101){
  cancer_genes <- read.xlsx(path_to_PCAWG_supp, sheet = "Table 7 List of cancer genes")
  
  CGC_v80 <- na.omit(cancer_genes$`603.genes.from.a.curated.version.of.the.Cancer.Gene.Census.(v80)`)
  
  
  CGC_v101 <- fread(path_to_cgc_v101)
  CGC_v101 <- CGC_v101$`Gene Symbol`
  
  
  newly_added_genes <- CGC_v101[which(!CGC_v101 %in% CGC_v80)]
  newly_added_genes
}

extract_HugoSymbol <- function(binID_vector){
  
  s <- strsplit(binID_vector, "[::]")
  GenomicElement <- unlist(lapply(s, function(x){x[5]}))
  GenomicElement
  
}

extract_PCAWG_drivers <- function(path_to_PCAWG_supp){
  
  
  PCAWG_drivers_coding <- read.xlsx(path_to_PCAWG_supp, sheet = "Table 4 Protein-coding point mu")
  
  PCAWG_drivers_nonCoding <- read.xlsx(path_to_PCAWG_supp, sheet = "Table 5 Non-coding point mut dr")
  
  PCAWG_drivers <- rbind(PCAWG_drivers_coding, PCAWG_drivers_nonCoding)
  PCAWG_drivers <- PCAWG_drivers[which(PCAWG_drivers$`Pre-filter.q-value` <.1),]
  
  PCAWG_drivers 
}

PCAWG_ID_annotations <- function(PCAWG_IDs, path_to_PCAWG_supp, path_CGC_new, 
                                 path_oncoKB){
  x=data.frame(cbind("PCAWG_IDs" = PCAWG_IDs,
                     "type_of_element" = define_element_type(PCAWG_IDs)))
  
  CGC_new <- fread(path_CGC_new)$'Gene Symbol'
  cancer_genes <- extract_cancerGenes(path_to_PCAWG_supp)
  CGC_goldstandard <- cancer_genes[["CGC_v80"]]
  hugoSymbol <- extract_HugoSymbol(PCAWG_IDs)
  x$in_CGC <- hugoSymbol %in% CGC_goldstandard
  x$in_CGC_literature <- hugoSymbol %in% cancer_genes[["CGC_and_literature"]]
  x$in_CGC_new <- hugoSymbol %in% CGC_new
  
  oncoKB <- fread(path_oncoKB)
  oncoKB <- oncoKB$`Hugo Symbol`
  x$in_oncoKB <- hugoSymbol %in% oncoKB
  
  PCAWG_drivers <- extract_PCAWG_drivers(path_to_PCAWG_supp)
  PCAWG_drivers <- PCAWG_drivers[,c("ID", "tissue")]
  names(PCAWG_drivers) <- c("ID", "tissue")
  in_pcawg <- PCAWG_IDs %in% PCAWG_drivers$ID
  x$in_pcawg <- in_pcawg
  y = PCAWG_drivers %>% group_by(ID)
  ID_groups <- group_split(y)
  IDs <- unlist(lapply(ID_groups, function(s) {s$ID[1]}))
  tissues <- unlist(lapply(ID_groups, function(s) {paste0(s$tissue, collapse = ",")}))
  df <- data.frame(cbind("PCAWG_IDs" = IDs, "tissues" = tissues))
  x_inpcawg <- x[in_pcawg,]
  x_inpcawg <- left_join(x_inpcawg, df, by = "PCAWG_IDs")
  x_not_inpcawg <- x[!in_pcawg,]
  x_not_inpcawg$tissues <- rep(NA, nrow(x_not_inpcawg))
  x <- rbind(x_inpcawg, x_not_inpcawg)
  x
}


save_annotated_PCAWG_IDs <- function(path_out, all_beds, path_to_PCAWG_supp, 
                                     path_test_y, path_CGC_new, 
                                     path_oncoKB){
  test_y <- fread(path_test_y)
  # PCAWG_IDs <- unique(c(all_beds$V4,test_y$binID)) 
  PCAWG_IDs <- test_y$binID
  ann_PCAWG_ID <- PCAWG_ID_annotations(PCAWG_IDs, path_to_PCAWG_supp, path_CGC_new, 
                                       path_oncoKB)
  
  write.csv(ann_PCAWG_ID, 
            file = paste0(path_out, "annotated_PCAWG_IDs.csv"),
            row.names = FALSE)
}

# tmpComplete_annotated_PCAWG_ID <- function(ann_PCAWG_ID, path_test_y){
#   
#   test_y <- fread(path_test_y)
#   colnames(test_y) <- c("PCAWG_IDs", "length", "nMut", "nSample", "N")
#   
#   comp_ann_PCAWG_ID <- full_join(test_y, ann_PCAWG_ID)
#   comp_ann_PCAWG_ID
# }

extract_ensembleIDs <- function(PCAWG_IDs){
  
  x <- unlist(lapply(strsplit(PCAWG_IDs, "::"), function(s){
    unlist(lapply(s[4], function(z){
      unlist(lapply(strsplit(z, "\\."), function(b){
        b[1]
      }))
    }))
  }))
  
  return(x)
}

Complete_annotated_PCAWG_ID <- function(path_to_PCAWG_supp, path_test_y,
                                        path_CGC_new, path_oncoKB, path_genhancer_PCAWG){
  test_y <- fread(path_test_y)
  colnames(test_y) <- ifelse(colnames(test_y) == 'binID', "PCAWG_IDs", colnames(test_y)) 
  test_y <- test_y[,c("PCAWG_IDs", "length")]
  
  ann_PCAWG_ID <- PCAWG_ID_annotations(test_y$PCAWG_IDs, path_to_PCAWG_supp, 
                                       path_CGC_new, path_oncoKB)
  enhancers <- ann_PCAWG_ID[which(ann_PCAWG_ID$type_of_element == "enhancers"),]
  enhancers$enhancer_pos <- unlist(lapply(strsplit(enhancers$PCAWG_IDs, "::"), 
                                          function(s){
                                            s[2]
                                          }))
  CDSs <- ann_PCAWG_ID[which(ann_PCAWG_ID$type_of_element == "gc19_pc.cds"),]
  CDSs$ENSG_ids <- extract_ensembleIDs(CDSs$PCAWG_IDs)
  gene_enhancer <- fread(path_genhancer_PCAWG, 
                         header = F, col.names = c("enhancer_pos", "ENSG_id"))
  ensemblIDs_enhancers <- strsplit(gene_enhancer$ENSG_id, ";")
  gene_enhancer_long <- data.frame(enhancer_pos = rep(gene_enhancer$enhancer_pos, 
                                                      unlist(lapply(ensemblIDs_enhancers, 
                                                                    length))),
                                   ENSG_ids = unlist(ensemblIDs_enhancers))
  CDS_enhancer <- left_join(CDSs, gene_enhancer_long, by = "ENSG_ids")
  enhancers2 <- left_join(enhancers, CDS_enhancer, by = c("enhancer_pos"))
  enhancers2$geneSymbol <- extract_HugoSymbol(enhancers2$PCAWG_IDs.y)
  
  enhancers2 <- enhancers2[,c("PCAWG_IDs.x", "type_of_element.x",
                              "in_CGC.y", "in_CGC_literature.y", 
                              "in_CGC_new.y", "in_oncoKB.y", "in_pcawg.x", 
                              "tissues.x", "geneSymbol")]
  colnames(enhancers2) <- c(colnames(ann_PCAWG_ID), "geneSymbol")
  
  
  df <- enhancers2 %>% group_by(PCAWG_IDs) %>% 
    summarise(sum_inCGC = sum((in_CGC)),
              sum_inCGC_litr = sum(in_CGC_literature),
              sum_inPCAWG = sum(in_pcawg),
              sum_inCGC_new = sum(in_CGC_new),
              sum_inoncoKB = sum(in_oncoKB),
              type_of_element = type_of_element,
              tissues = unique(tissues),
              geneSymbol = paste(unique(geneSymbol), collapse = ','),
              .groups = "keep")
  df <- df[!duplicated(df),]
  df$in_CGC <- ifelse(df$sum_inCGC == 0 | is.na(df$sum_inCGC), FALSE, TRUE )
  df$in_CGC_literature <- ifelse(df$sum_inCGC_litr == 0 | is.na(df$sum_inCGC_litr),
                                 FALSE, TRUE )
  df$in_pcawg <- ifelse(df$sum_inPCAWG == 0 | is.na(df$sum_inPCAWG),
                        FALSE, TRUE )
  df$in_CGC_new <- ifelse(df$sum_inCGC_new == 0 | is.na(df$sum_inCGC_new),
                          FALSE, TRUE )
  df$in_oncoKB <- ifelse(df$sum_inoncoKB == 0 | is.na(df$sum_inoncoKB),
                         FALSE, TRUE )
  df <- df[,c(colnames(ann_PCAWG_ID), "geneSymbol")]
  df_nonEnhancer <- ann_PCAWG_ID[which(ann_PCAWG_ID$type_of_element != "enhancers"),]
  df_nonEnhancer$geneSymbol <- extract_HugoSymbol(df_nonEnhancer$PCAWG_IDs)
  ann_PCAWG_ID <- rbind(df_nonEnhancer, df)
  comp_ann_PCAWG_ID <- full_join(test_y, ann_PCAWG_ID)
  comp_ann_PCAWG_ID$in_CGC_afterPCAWG <- comp_ann_PCAWG_ID$in_CGC_new & !comp_ann_PCAWG_ID$in_CGC
  comp_ann_PCAWG_ID
}

