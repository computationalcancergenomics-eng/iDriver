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

  First of all, please provide unique mutation identifiers to each mutation.

#### a. Somatic mutations
This is a table in which each row represents a single somatic variant. This file provides the essential genomic and contextual information used for driver discovery, and should have the following columns:

- donor_id — Unique patient identifier, showing the identifier of the patients in which the mutation is occoured
- PCAWG_ID — Full genomic element identifier
- Ref — Reference allele
- Alt — Alternative allele
- triCntx — Trinucleotide context of the mutation (e.g., CAG > CGG)
- var_type — Variant type (e.g., SNP, INDEL)
- GenomicElement — Genomic element category (e.g., gc19_pc.cds, gc19_pc.promCore, enhancers)
- position — Genomic coordinate (1-based)
- var_ID — User-assigned unique mutation identifier
- elemenntLength — Length (bp) of the element containing the mutation

Put this mutation data file in the following destination:
"../extdata/procInput/iDriverInputs/mutData.tsv"

#### b.  Donor information table
This file provides the essential information for each patient.

Required Columns:
- D_id — Unique patient identifier, showing the identifier of the patients (e.g., DO9226)
- cohort1 — Cancer project identifier (e.g., ColoRect-AdenoCA) 
- freq — Total number of mutations for the patient
- HyperMut_donor — Is the patient a hypermutated patient or not (TRUE/FALSE)

Put this donor information file in the following destination:
"../extdata/procInput/iDriverInputs/donorInfo.tsv"

#### c.  Mutations probabilities
This file provides the predicted number of mutations for each element and the corresponding probability of observing mutation in that element in one patient for a cancer cohort.


Required Columns:
- PCAWG_ID — Full genomic element identifier
- pElement — The probability of observing mutation in each element based on total number of mutations observed in a cancer cohort
- nPred — Estimated number of mutations for the element in a cancer cohort of N patients ( For estimating user can provide the estimated number of mutations from tools that can compute the Background mutation rate (BMR), e.g. DriverPower, eMET, ... . We used eMET for this section.

... to be completed with cmd scripts of save_pElems.R ...

#### d. FI scores


#### e. Mean-SD table


### 3. Run iDriver

```bash
Rscript tidy/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"
```


