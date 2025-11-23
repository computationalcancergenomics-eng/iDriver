#!/bin/bash
 
# List of cohorts

cohorts=( #"Pan_Cancer" 
   # "Pancan-no-skin-melanoma-lymph"
  "Biliary-AdenoCA"  "Bladder-TCC" "Bone-Leiomyo"  "Bone-Osteosarc"  "Lymph-BNHL"
  "Breast-AdenoCa"   "Cervix-SCC"  "CNS-GBM"
  "CNS-Medullo" "CNS-Oligo"  "CNS-PiloAstro"    "ColoRect-AdenoCA"
  "Eso-AdenoCa"  "Head-SCC"   "Kidney-RCC" "Liver-HCC"
   "Lung-AdenoCA" "Lung-SCC"
  "Lymph-CLL"        "Myeloid-MPN"      "Ovary-AdenoCA"
  "Panc-AdenoCA"     "Panc-Endocrine"   "Prost-AdenoCA"    "Skin-Melanoma"
  "Stomach-AdenoCA"  "Thy-AdenoCA"      "Uterus-AdenoCA" "Kidney-ChRCC"

  )

# cohorts=( "Lung-AdenoCA"
#   "Pancan-no-skin-melanoma-lymph" "ColoRect-AdenoCA" "CNS-GBM"
#   "Head-SCC"         "Skin-Melanoma"         "Lung-SCC"
#   "Stomach-AdenoCA"  "Uterus-AdenoCA"   "Biliary-AdenoCA"  "Eso-AdenoCa" "Lymph-BNHL"
# )

###### Rho = 0 ###### 
# template_command='Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/observed/Rho0/pooled_betaRho0_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/procInput/iDriverInputs/Pan_CancerRho0FILE.tsv"'
# template_command='Rscript tidy/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/simulated/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 75 "{cohort}" "../extdata/output_release2.0/simulated/Rho0/betaRho0_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/procInput/iDriverInputs/Pan_CancerRho0FILE.tsv"'

###### population-level ###### 
# template_command='Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/observed/population_level/pooled_betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "population_level" "../extdata/output_release2.0/observed/population_level/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'
# template_command='Rscript tidy/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/simulated/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/simulated/population_level_rmHyperMuts/betaRhoPancan_ElemType/" {flag} TRUE 4 "population_level" "../extdata/output_release2.0/observed/population_level/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'

###### Burden-free ###### 
# template_command='Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed_withHyperMutated/Pan_Cancer/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output_release2.0/observed/Burden_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "Burden-free"'
# template_command='Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/observed/Burden_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Burden-free" "../extdata/output_release2.0/observed/Burden_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'
# template_command='Rscript tidy/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/simulated/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/simulated/Burden_free/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Burden-free" "../extdata/output_release2.0/observed/Burden_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'

###### FI-free (version 1) ######
# template_command='Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed_withHyperMutated/Pan_Cancer/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output_release2.0/observed/FI_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "FI-free"'
# template_command='Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/observed/FI_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "FI-free" "../extdata/output_release2.0/observed/FI_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'
# template_command='Rscript tidy/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/simulated/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 75 "{cohort}" "../extdata/output_release2.0/simulated/FI_free/rmHyperMuts_betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "FI-free" "../extdata/output_release2.0/observed/FI_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'



###### Original (withHyperMuts) ######
# Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed_withHyperMutated/Pan_Cancer/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "Original"
# template_command='Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed_withHyperMutated/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" {flag} FALSE 4 "Original" "../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'



###### Original (rmHyperMuts) ######
template_command='Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'
# template_command='Rscript tidy/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/simulated/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 75 "{cohort}" "../extdata/output_release2.0/simulated/rmHyperMuts_betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'

###### Original (rmHyperMuts_GBM pElems) ######
# Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed_withHyperMutated/Pan_Cancer/pElems/GBM_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output_release2.0/observed/pooled_withHyperMuts/gbm_pElems/" FALSE FALSE 4 "Original"
# template_command='Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/GBM_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 75 "{cohort}" "../extdata/output_release2.0/observed/pooled/gbm_pElems/" {flag} TRUE 4 "Original" "../extdata/output_release2.0/observed/pooled_withHyperMuts/gbm_pElems/observed_Pan_Cancer.tsv"'


# template_command='Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/observed/pooled/forAverageTis/" {flag} TRUE 4 "Original" "../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'



# Loop over cohorts
for cohort in "${cohorts[@]}"
do
  echo "Running analysis for cohort: $cohort"
 
  # Determine logical flag for excluding Lymphoma and Melanoma
  if [[ "$cohort" == "Pancan-no-skin-melanoma-lymph" ]]; then
    second_last_flag="TRUE"
  else
    second_last_flag="FALSE"
  fi
 
  # Replace placeholders in the template
  modified_command="${template_command//\{cohort\}/$cohort}"
  modified_command="${modified_command//\{flag\}/$second_last_flag}"
 
  # Execute the command
  eval $modified_command
done