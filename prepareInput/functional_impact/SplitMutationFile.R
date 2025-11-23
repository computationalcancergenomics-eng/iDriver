##cat mutationFile_all.tsv | awk '{print $1,$2,".",$4,$5}' | sed -n '1,20 p' | gzip > mutationFile_all.1_20.vcf.gz
out_path <- "../extdata/procInput/simulated/sanger/onlineSet_splitted"

dir.create(out_path, recursive=TRUE, showWarnings = FALSE)
reverse = FALSE

fname = "splitMuts.sh"

sink(fname)

# wiritn the rhead
cat("#!/usr/bin/env bash\n \n")


split_2Range <- function(n, thr){ # n: total number of rows     & thr: threshold for which we want to split the data frame(in this case up to 100K variant can be submitted to CADD website)
  
  num <- ceiling(n/thr)   # by splitting the original file we obtain num files.
  
  y <- ceiling(n/num)
  x <- ceiling(n/num)
  
  vec_ends <- c()
  for(i in 1:num){ 
    vec_ends <- c(vec_ends, x)
    x <- y+x
  }
  
  vec_starts <- c(1, vec_ends[1:length(vec_ends)-1]+1)
  vec_ends <- c(vec_ends[1:length(vec_ends)-1] , n)
  
  df <- cbind(vec_starts, vec_ends)
  df <- data.frame(df)
  df
}

x <- split_2Range(3959761, 100000) # nrow(df_onlinePCAWG) >>> 2763622  #nrow(indel_mnv_simulated) >>> 2634414  #nrow(sanger_tmp_onLineSet_cadd_annotated) >>> 3959761

path_save="../extdata/procInput/simulated/sanger/onlineSet_splitted/"
dir.create(path_save, showWarnings = F, recursive = T)

for(i in 1:nrow(x)) {
  
  st = x$vec_starts[i]
  end = x$vec_ends[i]
  
  cat(paste0("cat ../extdata/procInput/simulated/sanger/tmp_onLineSet_cadd_annotated.tsv | sed -n '", st,",", end,
             " p' | gzip > ../extdata/procInput/simulated/sanger/onlineSet_splitted/MNV_InDel_onLineSet.", 
             st, "_", end, ".vcf.gz \n"))
  
  cat("echo ----------------------- \n \n")
}
