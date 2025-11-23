library(mvtnorm)
library(doParallel)
library(Rcpp)
library(dplyr)
library(tidyr)
library(VGAM)
library(PoissonBinomial)
library(data.table)
library(openxlsx)
library(metap)
library(nlme)

select_cohort <- function(path_donorInfo, 
                          cohort, exclude_lymph_melanoma = TRUE,
                          exclude_hyper_mutated = TRUE,
                          only_hyperMut_pancancer = FALSE){
  
  donorInfo <- fread(path_donorInfo)
  
  if (only_hyperMut_pancancer) {
    donorInfo <- donorInfo[which(donorInfo$HyperMut_donor),]
  } else {
    if (exclude_lymph_melanoma) {
      exceptions <- c("Skin-Melanoma", "SKCM-US",
                      "Lymph-NOS", 
                      "Lymph-CLL", "CLLE-ES",
                      "Lymph-BNHL", "MALY-DE", "DLBC-US")
      donorInfo <- donorInfo[-which(donorInfo$cohort1 %in% exceptions),] # 2279 donors
    }
    
    if (exclude_hyper_mutated) {
      donorInfo <- donorInfo[which(!donorInfo$HyperMut_donor),] 
    }
    
    if (!cohort %in% c('Pan_Cancer', 'Pancan-no-skin-melanoma-lymph')) {
      if (grepl('subsample', cohort)) {
        donorInfo <- donorInfo[which(donorInfo$subsample_id == cohort),]
      } else {
        donorInfo <- donorInfo[which(donorInfo$cohort1 == cohort),]
      }
      
    } 
  }
  
  
  
  if ('sampleID' %in% colnames(donorInfo)) {
    donorInfo <- donorInfo[,c('sampleID', "D_id","freq" )]
    colnames(donorInfo) <- c('sampleID', 'donor_id', 'totalMutNrs')
  } else {
    donorInfo <- donorInfo[,c("D_id","freq" )]
    colnames(donorInfo) <- c('donor_id', 'totalMutNrs')
  }
  
  donorInfo
  
  
  
}


selectDonorInfo <- function(donorInfo, cohort, hyper_filter, driver_filter) {
  # Check if the cohort is valid
  if (cohort != "Pan-Cancer" && !(cohort %in% donorInfo$cohort1)) {
    stop("Cohort not found in donorInfo$cohort1")
  }
  
  # Filter by cohort if not Pan-Cancer
  if (cohort != "Pan-Cancer") {
    donorInfo <- donorInfo[donorInfo$cohort1 == cohort, ]
  }
  
  # Filter by hypermutation status
  if (hyper_filter == "only_hyper") {
    donorInfo <- donorInfo[donorInfo$HyperMut_donor == TRUE, ]
  } else if (hyper_filter == "only_none_hyper") {
    donorInfo <- donorInfo[donorInfo$HyperMut_donor == FALSE, ]
  } # else no filter
  
  # Filter by driver status
  cohort_cgc <- paste0(cohort, "_CGC")
  cohort_pcawg <- paste0(cohort, "_PCAWG")
  
  # Some cohorts may not have a corresponding CGC or PCAWG column
  has_cgc <- cohort_cgc %in% colnames(donorInfo)
  has_pcawg <- cohort_pcawg %in% colnames(donorInfo)
  
  if (driver_filter == "only_driver") {
    if (has_cgc && has_pcawg) {
      donorInfo <- donorInfo[donorInfo[[cohort_cgc]] == 1 | donorInfo[[cohort_pcawg]] == 1, ]
    } else if (has_cgc) {
      donorInfo <- donorInfo[donorInfo[[cohort_cgc]] == 1, ]
    } else if (has_pcawg) {
      donorInfo <- donorInfo[donorInfo[[cohort_pcawg]] == 1, ]
    } else {
      warning("No CGC or PCAWG columns for this cohort — skipping driver filtering")
    }
  } else if (driver_filter == "only_none_driver") {
    if (has_cgc && has_pcawg) {
      donorInfo <- donorInfo[donorInfo[[cohort_cgc]] == 0 & donorInfo[[cohort_pcawg]] == 0, ]
    } else if (has_cgc) {
      donorInfo <- donorInfo[donorInfo[[cohort_cgc]] == 0, ]
    } else if (has_pcawg) {
      donorInfo <- donorInfo[donorInfo[[cohort_pcawg]] == 0, ]
    } else {
      warning("No CGC or PCAWG columns for this cohort — skipping driver filtering")
    }
  } # else no filter
  
  cols = c("D_id", "cohort1", "cohort2", "freq", "HyperMut_donor", "in_mut_data")
  
  if(has_pcawg) {
    donorInfo$in_pcawg = donorInfo[, cohort_pcawg]
    cols = c(cols, "in_pcawg")
  }
  
  if(has_cgc) {
    donorInfo$in_cgc = donorInfo[, cohort_cgc]
    cols = c(cols, "in_cgc")
  }
  
  
  donorInfo = donorInfo[, cols]
  donorInfo
}



add_methodSpecific_params <- function(mutDat_long, path_Mu_sd, path_pElem){
  
  # Add pElement to mutDat_long
  predRates <- fread(path_pElem)
  mutDat_long <- data.frame(mutDat_long)
  
  predRates <- predRates[, c('PCAWG_ID', 'pElement', 'nPred')]
  mutDat_long <- left_join(mutDat_long, predRates, by = 'PCAWG_ID')
  # mutDat_long$pElement <- (mutDat_long$predRate * mutDat_long$elemenntLength * N)/all_cohort_mutsNr
  
  # Add Mu_sd to mutDat_long
  mu_sd <- fread(path_Mu_sd)
  mutDat_long <- left_join(mutDat_long, mu_sd, by = 'PCAWG_ID')
  
  mutDat_long
}


create_per_cohort_MutData <- function(donorInfo, path_MutData, path_observed_scores){
  
  # create the mutData long format for the cohort of interest and scores of interest
  mutData <- fread(path_MutData)
  
  if ('donor_id' %in% colnames(mutData)) {
    mutData <- left_join(donorInfo, mutData, by = 'donor_id')
  } else {
    mutData <- left_join(donorInfo, mutData, by = "sampleID")
  }
  
  Scores <- fread(path_observed_scores)
  Scores <- Scores[!duplicated(Scores),]
  mutData <- left_join(mutData, Scores, by = 'var_ID')
  mutData
}

create_IndvLvl_geneInfo <- function(donorInfo, path_MutData, path_observed_scores, 
                                    path_Mu_sd, path_pElem){
  
  mutData = create_per_cohort_MutData(donorInfo, path_MutData, path_observed_scores)
  
  # add Method specific parameters to mutData
  completeMutDat <- add_methodSpecific_params(mutData, path_Mu_sd, path_pElem)
  completeMutDat <- completeMutDat[!grepl('lncrna.ncrna', completeMutDat$PCAWG_ID),]
  completeMutDat <- completeMutDat[!grepl('lncrna.promCore', completeMutDat$PCAWG_ID),]
  completeMutDat <- completeMutDat[!grepl('gc19_pc.ss', completeMutDat$PCAWG_ID),]
  
  print(paste0('# NAs for pElem that were excluded from analysis: ', sum(is.na(completeMutDat$pElement))))
  completeMutDat <- completeMutDat[which(!is.na(completeMutDat$pElement)),]
  
  completeMutDat <- completeMutDat[!duplicated(completeMutDat),]
  
  gene_groups <- completeMutDat %>% group_by(PCAWG_ID)
  gene_groups <- group_split(gene_groups)
  
  IDs <- unlist(lapply(gene_groups, function(s){
    unique(s$PCAWG_ID)
  }))
  
  names(gene_groups) <- IDs
  gene_groups
}


Compute_perDonor_obsStat <- function(n, df, mu, sd, pElement){
  df$pElement <- as.numeric(pElement)
  df$mu0 <- as.numeric(mu)
  df$sd0 <- as.numeric(sd)
  df$obs_stat <- apply(df, 1, applyCalcObservedStats, n= n)
  df
}

applyCalcObservedStats <- function(n, row) {
  nMut <- as.numeric(row["nMut"])
  sumScore <- as.numeric(row["sumScore"])
  pElement <- as.numeric(row["pElement"])
  donor_totMut <- as.numeric(row["totalMutNrs"])
  mu0 <- as.numeric(row["mu0"])
  sd0 <- as.numeric(row["sd0"])
  return (calcObservedStats(n, nMut, sumScore, donor_totMut,
                            pElement, mu0, sd0))
}

calcObservedStats <- function(n, n_obs, score_obs, totalMutNr, pElement, mu, sd) {
  if(n_obs == 0) {
    return(0)
  }
  p0_prec <- return_p0_prec(n,totalMutNr, pElement)
  p0 = p0_prec$p0
  p_rec = p0_prec$p_rec
  Fs = pnorm(score_obs, (1:n)*mu, sqrt(1:n)*sd)
  sum(p_rec*Fs)/(1-p0)
}

return_p0_prec <- function(n, totalMutNr, pElement){ 
  
  p0 = dbinom(0, totalMutNr, pElement)
  tmp = dbinom(1:(n-1), totalMutNr, pElement)
  p_n = pbinom(n, totalMutNr, pElement, lower.tail = FALSE)
  list(p0=p0, p_rec=c(tmp, p_n))
}


IH_CDF <- function(x, n, rho) {
  sigma = sqrt( (1/12) * (n + n * (n - 1) * rho))
  return(pnorm(x, mean = n/2, sd = sigma, lower.tail = TRUE))
}

compute_pval_exact <- function(t, pElement, totMutDon, nPred, rho = 0){
  success_probs = 1-dbinom(0, totMutDon, pElement)
  
  prob = 0 
  if(floor(t/2) > length(totMutDon)) {
    return (0)
  }
  
  for(k in max(0, floor(t-1)): max(10, min(floor(6*max(t, nPred)), length(success_probs)))) {
    prob = prob + exp(dpbinom(k, success_probs, log=T) + log(1-IH_CDF(t, k, rho)))
  }
  
  prob
}





# compute_beta_tail_prob <- function(t, alpha, beta, K, log = FALSE) {
#   # Compute mean and variance for the sum of K Beta distributions
#   mu <- K * alpha / (alpha + beta)
#   sigma_sq <- K * (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
#   sigma <- sqrt(sigma_sq)
#   
#   # Approximate the tail probability using the normal distribution
#   p_approx <- 1 - pnorm(t, mean = mu, sd = sigma)
#   
#   # Return probability in log scale if requested
#   if (log) {
#     return(log(p_approx))
#   } else {
#     return(p_approx)
#   }
# }



compute_beta_tail_prob <- function(t, alpha, beta, K, rho = 0, log = FALSE) {
  # Mean of one Beta(alpha, beta)
  mu_single <- alpha / (alpha + beta)
  
  # Variance of one Beta(alpha, beta)
  var_single <- (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
  
  # Total variance of sum with correlation rho among K variables
  # Var(Sum) = K * var + K*(K - 1)*rho*var
  sigma_sq <- var_single * (K + K * (K - 1) * rho)
  sigma <- sqrt(sigma_sq)
  
  # Mean of the sum
  mu <- K * mu_single
  
  # Normal approximation to tail probability
  p_approx <- 1 - pnorm(t, mean = mu, sd = sigma)
  
  # Return result
  if (log) {
    return(log(p_approx))
  } else {
    return(p_approx)
  }
}



compute_pval_exact_betaTail <- function(t, pElement, nPred, totMutDon, alpha, beta, rho = 0){
  success_probs = 1-dbinom(0, totMutDon, pElement)

  prob = 0
  if(floor(t/2) > length(totMutDon)) {
    return (0)
  }
  
  for(k in max(0, floor(t-1)): max(10, min(floor(6*max(t, nPred)), length(success_probs)))) {
  # for(k in max(0, floor(t-1)): max(10, length(success_probs))) {
    prob = prob + exp(dpbinom(k, success_probs, log=T) +
                        log(compute_beta_tail_prob(t, alpha, beta, k, rho, log = FALSE)))
  }

  prob
}

create_complete_donors_gene_df <- function(donorInfo, s){
  
  if(sum(is.na(s$score)) != 0){
    warning(paste0("At least one of the variants doesn't have score for ", unique(s$PCAWG_ID)))
  }
  
  genePatient <- s %>% group_by(donor_id, PCAWG_ID) %>% 
    summarise(sumScores = sum(score, na.rm = T),
              nMut=n(), .groups = "keep") %>% 
    arrange(desc(nMut))
  
  df <- left_join(donorInfo, genePatient, by = c("donor_id"))
  df$sumScore <- ifelse(is.na(df$sumScore), 0, df$sumScore)
  df$nMut <- ifelse(is.na(df$nMut), 0, df$nMut)
  df$PCAWG_ID <- rep(unique(s$PCAWG_ID), nrow(df))
  df
}






compute_rho <- function(cohort_obsStats){
  
  x <- as.data.table(cohort_obsStats %>%
                       group_by(PCAWG_ID) %>%
                       filter(n() >= 2) %>%
                       ungroup())
  
  # Check variance of obsStat by PCAWG_ID
  var_by_group <- x[, .(var_score = var(obs_stat, na.rm = TRUE)), by = PCAWG_ID]
  
  idx = which(var_by_group$var_score != 0 )
  x = x[idx,]
  n4var = nrow(x)
  
  if (n4var != 0) {
    # Fit random intercept model
    model <- lme(obs_stat ~ 1, random = ~1 | PCAWG_ID, data = x)
    
    # Extract variance components
    var_between <- as.numeric(VarCorr(model)[1, "Variance"])  # between-element variance
    var_within  <- as.numeric(VarCorr(model)[2, "Variance"])  # within-element residual
    
    # Compute ICC (intra-class correlation)
    rho <- var_between / (var_between + var_within)
    
    list('rho' = rho, 'n4var' = n4var)
    
  } else {
    NA
  }
  
  
  
}



##################### FI-free #########################################
Compute_perDonor_obsStat_FIfree <- function(n, df, pElement){
  df$pElement <- as.numeric(pElement)
  
  df$obs_stat <- apply(df, 1, applyCalcObservedStats_FIfree, n= n)
  df
}

applyCalcObservedStats_FIfree <- function(n, row) {
  nMut <- as.numeric(row["nMut"])
  sumScore <- as.numeric(row["sumScore"])
  pElement <- as.numeric(row["pElement"])
  donor_totMut <- as.numeric(row["totalMutNrs"])
  
  return (calcObservedStats_FIfree(n, nMut, sumScore, donor_totMut,
                            pElement))
}

calcObservedStats_FIfree <- function(n, n_obs, score_obs, totalMutNr, pElement) {
  if(n_obs == 0) {
    return(0)
  }
  p0_prec <- return_p0_prec(n,totalMutNr, pElement)
  p0 = p0_prec$p0
  p_rec = p0_prec$p_rec
  # Fs = pnorm(score_obs, (1:n)*mu, sqrt(1:n)*sd)
  # sum(p_rec*Fs)/(1-p0)
  sum(p_rec)/(1-p0)
}

get_cohortElemType_obsStats <- function(n, elemTypeIndvLvl, donorInfo, n_cores, iDriver_variation){
  
  dr <- fread('../extdata/procInput/ann_PCAWG_ID_complement_2025_cancerTypes.csv')
  
  nonDrivers <- dr[!(dr$in_CGC | dr$in_CGC_literature | dr$in_CGC_new | dr$in_oncoKB | dr$in_pcawg), ]
  
  registerDoParallel(n_cores)
  
  elemTypeIndvLvl <- elemTypeIndvLvl[which(names(elemTypeIndvLvl) %in% nonDrivers$PCAWG_IDs)]
  t0 <- Sys.time()
  
  cohort_obsStats <- foreach(i = 1:length(elemTypeIndvLvl), .combine = rbind) %dopar% {
    s <- elemTypeIndvLvl[[i]]
    PCAWG_ID <- unique(s$PCAWG_ID)
    elemenntLength <- unique(s$elemenntLength)
    
    df <- create_complete_donors_gene_df(donorInfo, s)
    pElement = as.numeric(unique(s$pElement))
    nPred = as.numeric(unique(s$nPred))
    
    if (iDriver_variation == 'FI-free') {
      df <- Compute_perDonor_obsStat_FIfree(n, df, pElement)
      
    } else if (iDriver_variation %in% c('Original', 'Burden-free')) {
      
      mu = unique(s$Mu)
      sd = unique(s$sd)
      df <- Compute_perDonor_obsStat(n, df, mu, sd, pElement)
      
    } else {
      cat('"iDriver_variation" argument should be one of the following: "Original", "Burden-free", "Burden-free (Rho = 0)", or "FI-free".')
    }
    
    
    nonzero_donors <- df[df$nMut > 0, ]
    elem_obsStats <- nonzero_donors[, c("PCAWG_ID", "obs_stat")]
    elem_obsStats$elemenntLength <- rep(elemenntLength, nrow(elem_obsStats))
    
    return(elem_obsStats)
  }
  
  # print(Sys.time() - t0)
  
  
  cohort_obsStats
  
}



##################### burden-free #########################################
compute_pval_exact_burdenFree <- function(t, pElement, sum_nMut, rho = 0){
  # success_probs = 1-dbinom(0, totMutDon, pElement)
  # 
  # prob = 0 
  # if(floor(t/2) > length(totMutDon)) {
  #   return (0)
  # }
  # 
  # for(k in max(0, floor(t-1)): max(10, min(floor(6*t), length(success_probs)))) {
  #   prob = prob + exp(dpbinom(k, success_probs, log=T) + log(1-IH_CDF(t, k, rho)))
  # }
  
  prob = 1-IH_CDF(t, sum_nMut, rho)
  prob
}


compute_pval_exact_betaTail_burdenFree <- function(t, pElement, nPred, sum_nMut, alpha, beta, rho = 0){
  # success_probs = 1-dbinom(0, totMutDon, pElement)
  # 
  # prob = 0
  # if(floor(t/2) > length(totMutDon)) {
  #   return (0)
  # }
  # 
  # for(k in max(0, floor(t-1)): max(10, min(floor(6*max(t, nPred)), length(success_probs)))) {
  #   # for(k in max(0, floor(t-1)): max(10, length(success_probs))) {
  #   prob = prob + exp(dpbinom(k, success_probs, log=T) +
  #                       log(compute_beta_tail_prob(t, alpha, beta, k, rho, log = FALSE)))
  # }
  
  prob = compute_beta_tail_prob(t, alpha, beta, sum_nMut, rho, log = FALSE)
  prob
}


compute_pval_helper <- function(rho_list, cohort_obsStats, elemTypeIndvLvl, n, donorInfo,
                                n_cores, iDriver_variation = 'Original'){
  
  print(paste0('Started calculating p-values with iDriver ', iDriver_variation, ' settings'))
  
  if (length(rho_list) == 2) {
    
    rho = rho_list$rho
    n4var = rho_list$n4var
    
    if (sum(c('alpha', 'beta') %in% colnames(cohort_obsStats)) == 2) {
      
      alpha <- unique(cohort_obsStats$alpha)
      beta <- unique(cohort_obsStats$beta)
      print(alpha)
      print(beta)
      
    } else {
      
      mean_x <- mean(cohort_obsStats$obs_stat)
      var_x  <- var(cohort_obsStats$obs_stat)
      
      alpha <- mean_x * ((mean_x * (1 - mean_x) / var_x) - 1)
      beta <- (1 - mean_x) * ((mean_x * (1 - mean_x) / var_x) - 1)
      
    }
    
    if (n4var != 0 & rho != 0) {
      print(n4var)
      print(rho)
    }
    
    
    N <- 1:length(elemTypeIndvLvl)
    chunks <- split(N, ceiling(seq_along(N)/5))
    
    registerDoParallel(n_cores)
    
    
    t1 = proc.time()
    pvals_needed <- foreach(i = 1:length(chunks),
                            .combine = 'rbind') %dopar% {
                              if(i %% 100 == 0){
                                print(paste("chunk", i, "out of", length(chunks)))
                              }
                              x <- lapply(elemTypeIndvLvl[chunks[[i]]], 
                                          function(s, n, alpha, beta, rho, n4var){
                                            # function(s, n, rho){
                                            elem <- unique(s$PCAWG_ID)
                                            df <- create_complete_donors_gene_df(donorInfo, s)
                                            pElement = as.numeric(unique(s$pElement))
                                            nPred = as.numeric(unique(s$nPred))
                                            mu = unique(s$Mu)
                                            sd = unique(s$sd)
                                            
                                            if (iDriver_variation == 'FI-free') {
                                              #### for FIfree there isnt any mu and sd ####
                                              df <- Compute_perDonor_obsStat_FIfree(n, df,
                                                                                    pElement)
                                            } else if (iDriver_variation %in% c('Original', 'Burden-free', "Burden-free (Rho = 0)")) {
                                              df <- Compute_perDonor_obsStat(n, df,
                                                                             mu, sd, 
                                                                             pElement)
                                            } else {
                                              cat('"iDriver_variation" argument should be one of the following: "Original", "Burden-free", "Burden-free (Rho = 0)", or "FI-free".')
                                            }
                                            
                                            
                                            
                                            
                                            t <- sum(df$obs_stat)
                                            
                                            totMutDon <- df$totalMutNrs
                                            sum_nMut = sum(df$nMut)
                                            sum_nSample = length(unique(s$donor_id))  # equals to nrow(df[which(df$nMut != 0),])
                                            
                                            if (iDriver_variation %in%  c('Burden-free', "Burden-free (Rho = 0)")){
                                              ##### for burden-free #####
                                              pvals <- compute_pval_exact_betaTail_burdenFree(t, pElement, nPred, sum_nSample, alpha, beta, rho)
                                              
                                            } else if (iDriver_variation %in% c('Original', 'FI-free')) {
                                              ##### for orig and FI-free #####
                                              pvals <- compute_pval_exact_betaTail(t, pElement, nPred, totMutDon, alpha, beta, rho)
                                              
                                            }
                                            
                                            
                                            data.frame(cbind("PCAWG_ID" = elem, 
                                                             "sum_scores" = ifelse(is.na(sum(df$sumScore)), NA, sum(df$sumScore)),
                                                             "pElement" = pElement,
                                                             "nPred" =nPred,
                                                             "mu0" = mu,
                                                             "sd0" = sd,
                                                             "sum_nMut" = sum(df$nMut),
                                                             "sum_nSample" = sum_nSample,
                                                             "obs_stat" = sum(df$obs_stat),
                                                             'alpha' = alpha,
                                                             'beta' = beta,
                                                             "rho" = rho,
                                                             "n4var" = n4var,
                                                             "pvals" = pvals))
                                          }, n=n, alpha = alpha, beta = beta, rho = rho, n4var)
                              # }, n=n, rho = rho)
                              
                              proc.time()-t1
                              
                              
                              p_vals <- do.call(rbind, x)
                              p_vals
                            }
    
    
    
  } else {
    pvals_needed <- c()
  }
  
  pvals_needed
}



# compute_iDriver_pVals <- function(n, donorInfo, indvLvl_geneInfo,
#                                   n_cores){
#   all_pvals_sh <- c()
#   all_pvals_ln <- c()
#   
#   elem_types <- c("gc19_pc.cds", "enhancers",
#                   # "lncrna.ncrna","lncrna.promCore",  "gc19_pc.ss"
#                   "gc19_pc.3utr", "gc19_pc.promCore", "gc19_pc.5utr")
#   
#   for (elem_type in elem_types) {
#     print(elem_type)
#     
#     if (sum(grepl(elem_type, names(indvLvl_geneInfo))) != 0) {
#       
#       elemTypeIndvLvl <- indvLvl_geneInfo[grepl(elem_type, names(indvLvl_geneInfo))]
#       cohort_obsStats <- get_cohortElemType_obsStats(n, elemTypeIndvLvl, donorInfo, n_cores)
#       
#       length_distr <- get_length_distr(elemTypeIndvLvl)
#       
#       cohort_obsStats_shorts <- cohort_obsStats[which(cohort_obsStats$elemenntLength < length_distr),]
#       cohort_obsStats_longs <- cohort_obsStats[which(cohort_obsStats$elemenntLength >= length_distr),]
# 
#       rho_list_shorts <- compute_rho(cohort_obsStats_shorts)
#       rho_list_longs <- compute_rho(cohort_obsStats_longs)
#       
#       # rho_list <- compute_rho(cohort_obsStats)
#       
#       # pvals_needed <- compute_pval_helper(rho_list, cohort_obsStats, elemTypeIndvLvl, n, donorInfo)
#       
#       idx_short <- unlist(lapply(elemTypeIndvLvl, function(s){
#         unique(s$elemenntLength) < length_distr
#       }))
#       
#       idx_long <- unlist(lapply(elemTypeIndvLvl, function(s){
#         unique(s$elemenntLength) >= length_distr
#       }))
#       
#       short_elems <- elemTypeIndvLvl[idx_short]
#       long_elems <- elemTypeIndvLvl[idx_long]
#       
#       pvals_needed_sh <- compute_pval_helper(rho_list_shorts, cohort_obsStats_shorts, 
#                                              short_elems, n, donorInfo)
#       
#       
#       pvals_needed_ln <- compute_pval_helper(rho_list_longs, cohort_obsStats_longs, 
#                                              long_elems, n, donorInfo)
#       
#       
#     }
#     all_pvals_sh <- rbind(all_pvals_sh, pvals_needed_sh)
#     all_pvals_ln <- rbind(all_pvals_ln, pvals_needed_ln)
#   }
#   all_pvals <- rbind(all_pvals_ln, all_pvals_sh)
#   all_pvals
# }
# 
# 
# get_length_distr <- function(elemTypeIndvLvl){
#   elem_lengths <- unlist(lapply(elemTypeIndvLvl, function(s){
#     unique(s$elemenntLength)
#   }))
#   
#   
#   qs <- summary(elem_lengths)
#   
#   med_lengths <- qs[3]
#   
#   med_lengths
#   
# }


# define_breaks_param_estimate <- function(all_elems_theSame = TRUE){
#   
#   if (all_elems_theSame) {
#     breaks_NC <- c(0, 100, 500, 1000, 1500, 2000, Inf)
#     labels_NC <- c("0-100 bp", "100-500 bp", "500-1000 bp", "1000-1500 bp",
#                 "1500-2000 bp", ">2000 bp")
#     
#     
#     
#     list('CDS'= list('breaks' = breaks,
#                      'labels' = labels),
#          'NC' = list('breaks' = breaks,
#                              'labels' = labels))
#   } else {
#     breaks_NC <- c(0, 100, 500, 1000, 1500, 2000, Inf)
#     labels_NC <- c("0-100 bp", "100-500 bp", "500-1000 bp", "1000-1500 bp",
#                    "1500-2000 bp", ">2000 bp")
#     
#     breaks_CDS <- c(0, 100, 500, 1000, 1500, 2000, Inf)
#     labels_CDS <- c("0-100 bp", "100-500 bp", "500-1000 bp", "1000-1500 bp",
#                    "1500-2000 bp", ">2000 bp")
#     
#     
#     list('CDS'= list('breaks' = breaks,
#                      'labels' = labels),
#          'NC' = list('breaks' = breaks,
#                      'labels' = labels))
#   }
#   
#   
# }

define_breaks_param_estimate <- function(elem_type) {
  
  elemtype_breaks <- list(
   
    'gc19_pc.cds' = list(
      breaks = c(0, 781, 1294, 2094, Inf),
      labels = c("≤781", "781–1294", "1294–2094", ">2094")
    ),
    'enhancers' = list(
      breaks = c(0, 199, 313, 490, Inf),
      labels = c("≤199", "199–313", "313–490", ">490")
    ),
    'gc19_pc.3utr' = list(
      breaks = c(0, 397, 1044, 2266, Inf),
      labels = c("≤397", "397–1044", "1044–2266", ">2266")
    ),
    'gc19_pc.5utr' = list(
      breaks = c(0, 129, 302, 595, Inf),
      labels = c("≤129", "129–302", "302–595", ">595")
    ),
    'gc19_pc.promCore' = list(
      breaks = c(0, 378, 572,930, Inf),
      labels = c("≤378", "378–572", "572–930", ">930")
    ),
    'gc19_pc.ss' = list(
      breaks = c(0, 104, 208, 364, Inf),
      labels = c("≤104", "104–208", "208–364", ">364")
    )
  )
  
  return(elemtype_breaks[[elem_type]])
}


compute_iDriver_pVals <- function(n, donorInfo, indvLvl_geneInfo, n_cores, iDriver_variation = 'Original') {
  all_pvals <- list()

  elem_types <- c("gc19_pc.cds", "enhancers",
                  # "lncrna.ncrna","lncrna.promCore",  "gc19_pc.ss"
                  "gc19_pc.3utr", "gc19_pc.promCore", "gc19_pc.5utr")



  for (elem_type in elem_types) {
    print(elem_type)

    # Define the breaks and labels for element length categories
    breaks_params_elemType <- define_breaks_param_estimate(elem_type)

    breaks <- breaks_params_elemType$breaks
    labels <- breaks_params_elemType$labels

    if (sum(grepl(elem_type, names(indvLvl_geneInfo))) != 0) {

      elemTypeIndvLvl <- indvLvl_geneInfo[grepl(elem_type, names(indvLvl_geneInfo))]


      cohort_obsStats <- get_cohortElemType_obsStats(n, elemTypeIndvLvl, donorInfo, n_cores, iDriver_variation)

      # Bin cohort_obsStats by element length
      cohort_obsStats$length_category <- cut(cohort_obsStats$elemenntLength,
                                             breaks = breaks, labels = labels, right = FALSE)

      for (cat_label in labels) {
        cat_obsStats <- cohort_obsStats[cohort_obsStats$length_category == cat_label, ]
        if (nrow(cat_obsStats) == 0) next

        if (iDriver_variation == 'Burden-free (Rho = 0)') {
          rho_list <- list('rho' = 0, 'n4var' = NA)
        } else  if (iDriver_variation %in% c("Original", "FI-free", "Burden-free", 'Burden-free (Rho = 0)'))  {
          if (var(cat_obsStats$obs_stat) < .000001 ) {

            if (var(cohort_obsStats$obs_stat) < .000001) {
              get_n4var <- function(cohort_obsStats){
                x <- as.data.table(cohort_obsStats %>%
                                     group_by(PCAWG_ID) %>%
                                     filter(n() >= 2) %>%
                                     ungroup())

                # Check variance of obsStat by PCAWG_ID
                var_by_group <- x[, .(var_score = var(obs_stat, na.rm = TRUE)), by = PCAWG_ID]

                idx = which(var_by_group$var_score != 0 )
                x = x[idx,]
                n4var = nrow(x)
                n4var
              }


              rho_list <- list('rho' = 0, 'n4var' = paste0('Failed to calculate rho, n4var = ', get_n4var(cat_obsStats)))

            } else {
              rho_list <- compute_rho(cohort_obsStats)
              rho_list <- list('rho' = rho_list$rho, 'n4var' = paste0('Failed to calculate rho, n4var = ', rho_list$n4var))
            }
          } else {
            rho_list <- compute_rho(cat_obsStats)
          }

        } else {
          cat('"iDriver_variation" argument should be one of the following: "Original", "Burden-free", "Burden-free (Rho = 0)", or "FI-free".')
        }

        idx <- unlist(lapply(elemTypeIndvLvl, function(s){
          unique_cat <- cut(unique(s$elemenntLength), breaks = breaks, labels = labels, right = FALSE)
          as.character(unique_cat) == cat_label
        }))

        cat_elems <- elemTypeIndvLvl[idx]

        pvals_cat <- compute_pval_helper(rho_list, cat_obsStats,
                                         cat_elems, n, donorInfo, n_cores, iDriver_variation)

        pvals_cat$length_category <- cat_label
        all_pvals[[paste(elem_type, cat_label, sep = "_")]] <- pvals_cat
      }
    }
  }

  # Combine all p-values into one data frame
  final_pvals <- do.call(rbind, all_pvals)
  return(final_pvals)
}

# #computing iDriver p-values without categorizing elements by length for burden-free or any desired pattern
# compute_iDriver_pVals <- function(n, donorInfo, indvLvl_geneInfo, n_cores, iDriver_variation = 'Original') {
#   all_pvals <- list()
#   
#   elem_types <- c("gc19_pc.cds", "enhancers",
#                   # "lncrna.ncrna","lncrna.promCore",  "gc19_pc.ss"
#                   "gc19_pc.3utr", "gc19_pc.promCore", "gc19_pc.5utr")
#   
#   
#   per_length_category <- ifelse(iDriver_variation == 'Burden-free', F, T)
#   
#   for (elem_type in elem_types) {
#     print(elem_type)
#     
#     # Define the breaks and labels for element length categories
#     breaks_params_elemType <- define_breaks_param_estimate(elem_type)
#     
#     breaks <- breaks_params_elemType$breaks
#     labels <- breaks_params_elemType$labels
#     
#     if (sum(grepl(elem_type, names(indvLvl_geneInfo))) != 0) {
#       
#       elemTypeIndvLvl <- indvLvl_geneInfo[grepl(elem_type, names(indvLvl_geneInfo))]
#       
#       
#       cohort_obsStats <- get_cohortElemType_obsStats(n, elemTypeIndvLvl, donorInfo, n_cores, iDriver_variation)
#       
#       # Bin cohort_obsStats by element length
#       cohort_obsStats$length_category <- cut(cohort_obsStats$elemenntLength, 
#                                              breaks = breaks, labels = labels, right = FALSE)
#       
#       if (per_length_category){
#         for (cat_label in labels) {
#           cat_obsStats <- cohort_obsStats[cohort_obsStats$length_category == cat_label, ]
#           if (nrow(cat_obsStats) == 0) next
#           
#           if (iDriver_variation == 'Burden-free') {
#             rho_list <- list('rho' = 0, 'n4var' = NA)
#           } else  if (iDriver_variation %in% c("Original", "FI-free"))  {
#             if (var(cat_obsStats$obs_stat) < .000001 ) {
#               
#               if (var(cohort_obsStats$obs_stat) < .000001) {
#                 get_n4var <- function(cohort_obsStats){
#                   x <- as.data.table(cohort_obsStats %>%
#                                        group_by(PCAWG_ID) %>%
#                                        filter(n() >= 2) %>%
#                                        ungroup())
#                   
#                   # Check variance of obsStat by PCAWG_ID
#                   var_by_group <- x[, .(var_score = var(obs_stat, na.rm = TRUE)), by = PCAWG_ID]
#                   
#                   idx = which(var_by_group$var_score != 0 )
#                   x = x[idx,]
#                   n4var = nrow(x)
#                   n4var
#                 }
#                 
#                 
#                 rho_list <- list('rho' = 0, 'n4var' = paste0('Failed to calculate rho, n4var = ', get_n4var(cat_obsStats)))
#                 
#               } else {
#                 rho_list <- compute_rho(cohort_obsStats)
#                 rho_list <- list('rho' = rho_list$rho, 'n4var' = paste0('Failed to calculate rho per length category (used elemtype rho), n4var = ', rho_list$n4var))
#               } 
#             } else {
#               rho_list <- compute_rho(cat_obsStats)
#             }
#             
#           } else {
#             cat('"iDriver_variation" argument should be one of the following: "Original", "Burden-free", or "FI-free".')
#           }
#           
#           idx <- unlist(lapply(elemTypeIndvLvl, function(s){
#             unique_cat <- cut(unique(s$elemenntLength), breaks = breaks, labels = labels, right = FALSE)
#             as.character(unique_cat) == cat_label
#           }))
#           
#           cat_elems <- elemTypeIndvLvl[idx]
#           
#           pvals_cat <- compute_pval_helper(rho_list, cat_obsStats, 
#                                            cat_elems, n, donorInfo, n_cores, iDriver_variation)
#           
#           pvals_cat$length_category <- cat_label
#           all_pvals[[paste(elem_type, cat_label, sep = "_")]] <- pvals_cat
#         }
#         
#         
#         
#       } else {
#         print('parameters are estimated per element-type without categorizing based on length')
#         if (iDriver_variation == 'Burden-free') {
#           rho_list <- list('rho' = 0, 'n4var' = NA)
#           
#         }
#         
#         pvals_elemType <- compute_pval_helper(rho_list, cohort_obsStats, 
#                                            elemTypeIndvLvl, n, donorInfo, n_cores, iDriver_variation)
#         all_pvals[[elem_type]] <- pvals_elemType
#       }
#     }
#   }
#   # Combine all p-values into one data frame
#   final_pvals <- do.call(rbind, all_pvals)
#   
#   return(final_pvals)
# }

        #                              filter(n() >= 2) %>%
        #                              ungroup())
        #         
        #         # Check variance of obsStat by PCAWG_ID
        #         var_by_group <- x[, .(var_score = var(obs_stat, na.rm = TRUE)), by = PCAWG_ID]
        #         
        #         idx = which(var_by_group$var_score != 0 )
        #         x = x[idx,]
        #         n4var = nrow(x)
        #         n4var
        #       }
        #       
        #       
        #       rho_list <- list('rho' = 0, 'n4var' = paste0('Failed to calculate rho, n4var = ', get_n4var(cat_obsStats)))
        #       
        #     } else {
        #       rho_list <- compute_rho(cohort_obsStats)
        #       rho_list <- list('rho' = rho_list$rho, 'n4var' = paste0('Failed to calculate rho, n4var = ', rho_list$n4var))
        #     }
        #   } else {
        #     rho_list <- compute_rho(cat_obsStats)
        #   }
        # rho_list <- compute_rho(cat_obsStats)
        
        
        # } else {
        #   cat('"iDriver_variation" argument should be one of the following: "Original", "Burden-free", or "FI-free".')
        # }
#         rho_list <- list('rho' = unique(cat_obsStats$rho),
#                          'n4var' = unique(cat_obsStats$n4var))
#         
#         
#         
#         idx <- unlist(lapply(elemTypeIndvLvl, function(s){
#           unique_cat <- cut(unique(s$elemenntLength), breaks = breaks, labels = labels, right = FALSE)
#           as.character(unique_cat) == cat_label
#         }))
#         
#         cat_elems <- elemTypeIndvLvl[idx]
#         
#         pvals_cat <- compute_pval_helper(rho_list, cat_obsStats,
#                                          cat_elems, n, donorInfo, n_cores, iDriver_variation)
#         
#         pvals_cat$length_category <- cat_label
#         all_pvals[[paste(elem_type, cat_label, sep = "_")]] <- pvals_cat
#       }
#     }
#   }
#   
#   # Combine all p-values into one data frame
#   final_pvals <- do.call(rbind, all_pvals)
#   return(final_pvals)
# }


compute_iDriver_pVals_preCalacParams <- function(n, donorInfo, indvLvl_geneInfo, n_cores,
                                                 precalcParams, iDriver_variation = 'Original') {
  all_pvals <- list()
  
  elem_types <- c("gc19_pc.cds", "enhancers",
                  # "lncrna.ncrna","lncrna.promCore",  "gc19_pc.ss"
                  "gc19_pc.3utr", "gc19_pc.promCore", "gc19_pc.5utr")
  
  
  
  for (elem_type in elem_types) {
    print(elem_type)
    
    # Define the breaks and labels for element length categories
    breaks_params_elemType <- define_breaks_param_estimate(elem_type)
    
    breaks <- breaks_params_elemType$breaks
    labels <- breaks_params_elemType$labels
    
    if (sum(grepl(elem_type, names(indvLvl_geneInfo))) != 0) {
      
      elemTypeIndvLvl <- indvLvl_geneInfo[grepl(elem_type, names(indvLvl_geneInfo))]
      cohort_obsStats <- precalcParams[which(precalcParams$PCAWG_ID %in% names(elemTypeIndvLvl)),]
      
      
      
      for (cat_label in labels) {
        cat_obsStats <- cohort_obsStats[cohort_obsStats$length_category == cat_label, ]
        if (nrow(cat_obsStats) == 0) next
        
        rho_list <- list('rho' = unique(cat_obsStats$rho),
                         'n4var' = unique(cat_obsStats$n4var))
        
        
        
        idx <- unlist(lapply(elemTypeIndvLvl, function(s){
          unique_cat <- cut(unique(s$elemenntLength), breaks = breaks, labels = labels, right = FALSE)
          as.character(unique_cat) == cat_label
        }))
        
        cat_elems <- elemTypeIndvLvl[idx]
        
        pvals_cat <- compute_pval_helper(rho_list, cat_obsStats,
                                         cat_elems, n, donorInfo, n_cores, iDriver_variation)
        
        pvals_cat$length_category <- cat_label
        all_pvals[[paste(elem_type, cat_label, sep = "_")]] <- pvals_cat
      }
    }
  }
  
  # Combine all p-values into one data frame
  final_pvals <- do.call(rbind, all_pvals)
  return(final_pvals)
}




iDriver <- function(study, path_observed_scores, path_pElem, 
                    path_Mu_sd, save_name, n_cores = 40, cohort= 'Pan_Cancer',
                    save_dir, exclude_lymph_melanoma = TRUE,
                    exclude_hyper_mutated = TRUE,
                    n = 4, iDriver_variation= 'Original',
                    cancer_specific_params = NULL, only_hyperMut_pancancer = F ){
  
  if(study == "simulated") {
    path_MutData <- "../extdata/procInput/iDriverInputs/simulated/Sanger/mutData.tsv"
    path_donorInfo <- '../extdata/procInput/iDriverInputs/simulated/Sanger/donorInfo.tsv'
    # save_dir = '../extdata/output/individual_level/efficientCodes/typeIerror/'
  } else if (study == 'observed') {
    path_MutData <- "../extdata/procInput/iDriverInputs/mutData.tsv"
    path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
    
  }else if (study == 'subsamples'){
    path_donorInfo <- '../extdata/procInput/iDriverInputs/subsamples_dInfo.tsv'
    path_MutData <- "../extdata/procInput/iDriverInputs/mutData.tsv"
    
  }
  
  if (cohort == 'onlyHyperMuts') {
    donorInfo <- select_cohort(path_donorInfo = path_donorInfo,
                               only_hyperMut_pancancer = TRUE)
  } else {
    donorInfo <- select_cohort(path_donorInfo,
                               cohort,
                               exclude_lymph_melanoma,
                               exclude_hyper_mutated)
  }
  
  if (iDriver_variation == 'population_level') {
    donorInfo$totalMutNrs <- as.integer(sum(donorInfo$totalMutNrs)/nrow(donorInfo))
    iDriver_variation = 'Original'
    print(paste0('Running population-level iDriver with ', iDriver_variation, ' settings...'))
  }
  
  indvLvl_geneInfo <- create_IndvLvl_geneInfo(donorInfo, path_MutData,
                                              path_observed_scores,
                                              path_Mu_sd, path_pElem)
  print(length(indvLvl_geneInfo))
  
  included_IDs <- fread(path_Mu_sd)
  included_IDs <- included_IDs[!grepl('lncrna.ncrna', included_IDs$PCAWG_ID),]
  included_IDs <- included_IDs[!grepl('lncrna.promCore', included_IDs$PCAWG_ID),]
  
  indvLvl_geneInfo <- indvLvl_geneInfo[which(names(indvLvl_geneInfo) %in% included_IDs$PCAWG_ID)]
  print(length(indvLvl_geneInfo))
  if (is.null(cancer_specific_params)) {
    
    print('Calculating p-values using cancer-specific parameters...')
    
    p_vals <- compute_iDriver_pVals(n, donorInfo, indvLvl_geneInfo, n_cores, iDriver_variation)
  } else {
    
    print('We are using precalculated params for this cohort:')
    precalcParams <- fread(cancer_specific_params)
    precalcParams <- precalcParams[,c('PCAWG_ID', 'alpha', 'beta', 'rho', "length_category", "n4var")]
    
    p_vals <- compute_iDriver_pVals_preCalacParams(n, donorInfo, indvLvl_geneInfo, n_cores, precalcParams, iDriver_variation)
    
  }
  
  
  if (study == "simulated") {
    p_vals$fdr <- p.adjust(p_vals$pvals, method = 'fdr')
    
    print(sum(p_vals$fdr <.1))
    print(sum(p_vals$fdr <.05))
  }
  
  p_vals$pvals <- as.numeric(p_vals$pvals)
  p_vals <- p_vals[order(p_vals$pvals),]
  
  path_save = paste0(save_dir, study, "_", save_name, ".tsv")
  dir.create(save_dir, showWarnings = F, recursive = T)
  fwrite(p_vals, file = path_save, sep = '\t')
  
}



iDriver_vartype <- function(study, path_observed_scores, 
                            save_name, is_SNV, n_cores,
                            cohort,
                            exclude_lymph_melanoma = TRUE,
                            exclude_hyper_mutated = TRUE,
                            n = 4, cancer_specific_params = NULL){
  
  varType <- ifelse(is_SNV, 'SNV', 'nonSNV')
  save_name <- paste0(save_name, '_', varType)
  
  if(study == "simulated") {
    path_MutData <- "../extdata/procInput/iDriverInputs/simulated/Sanger/mutData.tsv"
    path_donorInfo <- '../extdata/procInput/iDriverInputs/simulated/Sanger/donorInfo.tsv'
    save_dir = paste0('../extdata/output_release2.0/simulated/SNV_nonSNV/', cohort, '/tmp/')
    
  } else if(study == "observed") {
    path_MutData <- paste0('../extdata/procInput/SNV_nonSNV/mutData/', varType, '_mutData.tsv')
    path_donorInfo <- paste0('../extdata/procInput/SNV_nonSNV/donorInfo/obs_', varType, '_donorInfo.tsv')
    # save_dir = paste0('../extdata/output_release2.0/', study, '/SNV_nonSNV/gbm_pElems/', cohort, '/tmp/')
    save_dir = paste0('../extdata/output_release2.0/', study, '/SNV_nonSNV/eMET_pElems/', cohort, '/tmp/')
    path_Mu_sd <- paste0('../extdata/procInput/SNV_nonSNV/Mu_sd/', varType, '/element_type/', cohort, '.tsv')
    path_pElem <- paste0('../extdata/procInput/SNV_nonSNV/BMRs/observed/', cohort, '/', varType, '/pElems/eMET_orig_pElmens.tsv') #Original pElems
    # path_pElem <- paste0('../extdata/procInput/SNV_nonSNV/BMRs/observed/', cohort, '/', varType, '/pElems/GBM_orig_pElmens.tsv') #Original pElems
  }
  
  donorInfo <- select_cohort(path_donorInfo,
                             cohort,
                             exclude_lymph_melanoma,
                             exclude_hyper_mutated)
  
  indvLvl_geneInfo <- create_IndvLvl_geneInfo(donorInfo, path_MutData,
                                              path_observed_scores,
                                              path_Mu_sd, path_pElem)
  print(length(indvLvl_geneInfo))
  
  included_IDs <- fread(path_Mu_sd)
  included_IDs <- included_IDs[!grepl('lncrna.ncrna', included_IDs$PCAWG_ID),]
  included_IDs <- included_IDs[!grepl('lncrna.promCore', included_IDs$PCAWG_ID),]
  
  indvLvl_geneInfo <- indvLvl_geneInfo[which(names(indvLvl_geneInfo) %in% included_IDs$PCAWG_ID)]
  print(length(indvLvl_geneInfo))
  
  if (is.null(cancer_specific_params)) {
    
    print('Calculating p-values using cancer-specific parameters...')
    
    p_vals <- compute_iDriver_pVals(n, donorInfo, indvLvl_geneInfo, 
                                    n_cores, iDriver_variation = 'Original')
  } else {
    
    print('We are using precalculated params for this cohort:')
    precalcParams <- fread(cancer_specific_params)
    precalcParams <- precalcParams[,c('PCAWG_ID', 'alpha', 'beta', 'rho', "length_category", "n4var")]
    
    p_vals <- compute_iDriver_pVals_preCalacParams(n, donorInfo, indvLvl_geneInfo,
                                                   n_cores, precalcParams, 
                                                   iDriver_variation = 'Original')
    
  }
  
  if (study == "simulated") {
    p_vals$fdr <- p.adjust(p_vals$pvals, method = 'fdr')
    
    print(sum(p_vals$fdr <.1))
    print(sum(p_vals$fdr <.05))
  }
  p_vals <- p_vals[order(p_vals$pvals),]
  path_save = paste0(save_dir, study, "_", save_name, ".tsv")
  dir.create(save_dir, showWarnings = F, recursive = T)
  fwrite(p_vals, file = path_save, sep = '\t')
  
  
}
