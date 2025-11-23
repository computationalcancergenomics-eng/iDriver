rm(list=ls())

source("tidy/functions_iDriver.R")


# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
study <- args[1]
path_observed_scores <- args[2]
path_pElem <- args[3]
path_Mu_sd <- args[4]
save_name <- args[5]
n_cores <- as.numeric(args[6])
cohort <- args[7]
save_dir <- args[8]
exclude_lymph_melanoma <- as.logical(args[9])
exclude_hyper_mutated <- as.logical(args[10])
n <- as.numeric(args[11])
iDriver_variation=args[12]
if (length(args) >= 13 && args[13] != "NULL") {
  params <- args[13]
} else {
  params <- NULL
}

# Run the iDriver function with the specified parameters
iDriver(study, path_observed_scores, path_pElem, path_Mu_sd, save_name,
        n_cores, cohort, save_dir, exclude_lymph_melanoma,
        exclude_hyper_mutated, n, iDriver_variation, params)


print(' ************** Job Done **************')


# Run the R script with the arguments
# nohup Rscript tidy/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/iDriverInputs/pElems/Pancan-no-skin-melanoma-lymph/sangerSim/originalPred.tsv" "../extdata/procInput/Mu_sd/element_type_sangerParams_alpha1.5.tsv" "PR2_1.5alpha"  40 "Pan_Cancer" TRUE TRUE 4 > ../nohup_iDriver &