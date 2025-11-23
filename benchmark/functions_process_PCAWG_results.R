

readPvalues <- function(fname) {
  pval_col_names = tolower(c("p-value", "pval", "p", "raw-p", "raw_p","p_value", "pmax", "pCVmax", "combined.p"))
  qval_col_names =  tolower(c("q-value", "q", "qval"))
  id_col_names = tolower(c("id", "name", "element_ID", "gene", "PCAWG_id", "id.element"))
  other_col_names = tolower(c("Chrom", "Start", "End", "chr", "shortname", "pcvmid", "pcvmin", "qcv",
                              "pcv", "pcl", "pfn",
                              "mis_mle" , "non_mle" , "mis_low", "non_low", "mis_high", "non_high",
                              "wmis", "wtrunc", "npat", "npat_exp", "freq", "freq_excess",
                              "obsexp_subs",   "obsexp_indels"))
  all_col_names = c(pval_col_names, qval_col_names, id_col_names, other_col_names)
  
  
  dat = as.data.frame(fread(fname))
  if(nrow(dat) == 0) {
    return(NULL)
  }


  if(sum(tolower(colnames(dat)) %in% pval_col_names) != 1 )  {
    print(head(dat))
    stop(paste("[pvalue] error in reading ", fname, "---", paste(colnames(dat), collapse = ' ')))
  }
  
      
  if(sum(tolower(colnames(dat)) %in% tolower(id_col_names)) != 1) {
    if( tolower("PCAWG_id") %in% tolower(colnames(dat))) {
      id_col_names = tolower("PCAWG_id")
    } else {
      print(head(dat))
      stop(paste("[id] error in reading ", fname, "---", paste(colnames(dat), collapse = ' ')))
    }
  }
  
  pval_index = which(tolower(colnames(dat)) %in% pval_col_names)
  qval_index = which(tolower(colnames(dat)) %in% qval_col_names)
  id_col_index = which(tolower(colnames(dat)) %in% id_col_names)
  
  if(length(qval_index) == 0) {
    q_values = rep(NA, nrow(dat)) 
    if(ncol(dat) > 2) {
      tmp = setdiff(tolower(colnames(dat)), all_col_names)
      if(length(tmp) > 0 ) {
        print(head(dat))
        print(tmp)
        print(fname)
        print("----------------")
        stop()
        
      }
    }
  } else {
    q_values = dat[, qval_index]
  }
  
  dat = dat[, c(id_col_index, pval_index)]
  dat$q_value = q_values
  colnames(dat) = c("PCAWG_IDs",  "p_value", "q_value")
  indexes = which(!is.na(dat$p_value))
  if(length(indexes) > 0) {
    dat = dat[indexes, ]
  }
  
  dat$fdr = p.adjust(dat$p_value, "fdr")
  dat
}

createSinglePRes <- function(path) {
  fnames = list.files(path, full.names = T)
  zero_fnames = fnames[which(file.size(fnames) == 0)]
  if(length(zero_fnames) > 0 ) {
    warning(paste("Zero file size:", paste(zero_fnames, collapse='\n')))  
  }
  
  fnames = fnames[which(file.size(fnames)> 0)]
  
  tmpList <-  lapply(fnames, readPvalues)
  do.call(rbind, tmpList)
}

createPRes <- function(tissue_path) {
  methodNames = list.files(tissue_path, full.names = F)
  paths = paste0(tissue_path, "/", methodNames)
  tmpList <-  lapply(paths, createSinglePRes)
  class(tmpList) = "pRes"
  names(tmpList) = methodNames
  
  tmpList
}

save_pRes <- function(tissue_PCAWG_res, path_procPCAWG_res){
  
  tissue_names <- list.files(tissue_PCAWG_res)
  tissue_paths = list.files(tissue_PCAWG_res, full.names = T)
  dir.create(path_procPCAWG_res, showWarnings = F, recursive = T) 
  
  for( k in 1:length(tissue_paths)) {
    tissue_path = tissue_paths[k]
    
    pRes = createPRes(tissue_path)
    print(tissue_path)
    print(object.size(pRes))
    
   
   save(pRes, file = paste0(path_procPCAWG_res, tissue_names[k], ".RData"))
   
  }
  
}
