<p align="left">
  <img src="logo.png" alt="Project Logo" width="200">
</p>



## Welcome to iDriver

iDriver is a probabilistic graphical model that identifies positively selected coding and non-coding cancer driver elements by jointly integrating mutation recurrence and functional impact at the individual-patient level.


### 1. Installation

```bash
git clone https://github.com/computationalcancergenomics-eng/iDriver.git 
```

### 2. Required data files

iDriver requires five data files as input:


- Somatic mutations (SNVs & indels) from a cancer cohort
- Donor information table
- The mutation probability of each element for a cancer cohort
- Functional impact (FI) scores for all variants in a cancer cohort
- The mean and standard deviation of FI scores for each element type in a cancer cohort

#### a. Somatic mutations


#### b.  Donor information table


#### c.  Mutations probabilities


#### d. FI scores


#### e. Mean-SD table


### 3. Run iDriver

```bash
Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"
```


