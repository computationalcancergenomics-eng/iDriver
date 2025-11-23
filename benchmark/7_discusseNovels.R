rm(list = ls())

library(data.table)
source('tidy/functions_iDriver.R')


path_sigHits <- '../extdata/output_release2.0/Benchmark/pooled_catQuartiles_rmHyperMuts/SigHits/allSigHits.tsv'
path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
path_MutData <- '../extdata/procInput/iDriverInputs/mutData.tsv'
bed <- fread('../extdata/rawInput/PCAWG_test_genomic_elements.bed12.gz')
bed <- bed[, c('V1', 'V4')]
colnames(bed) <- c('chr', 'PCAWG_ID')

hit_id <- "gc19_pc.cds::gencode::ZIC1::ENSG00000152977.5" 
cohort <- "CNS-Medullo"
exclude_hyper_mutated <- T
path_observed_scores <- '../extdata/procInput/iDriverInputs/cadd_scores.tsv'
path_Mu_sd <- paste0("../extdata/procInput/Mu_sd/elemType/observed/", cohort, "_CADD6_elemTypeParams.tsv")
path_pElem <- paste0('../extdata/procInput/BMRs_2024/observed/', cohort,'/pElems/eMET_orig_pElmens.tsv')

hits <- fread(path_sigHits)
hits <- hits[which(!hits$cohort %in% c("Pan_Cancer", "Pancan-no-skin-melanoma-lymph", "Lymph-BNHL")),]

newHits <- hits[(!(hits$in_CGC_new | hits$in_pcawg | hits$in_oncoKB)),]


newHits[newHits$PCAWG_ID == hit_id,]


exclude_lymph_melanoma <- ifelse(cohort == "Pancan-no-skin-melanoma-lymph", T, F)
donorInfo <- select_cohort(path_donorInfo,
                           cohort,
                           exclude_lymph_melanoma,
                           exclude_hyper_mutated = T)

indvLvl_geneInfo <- create_IndvLvl_geneInfo(donorInfo, path_MutData,
                                            path_observed_scores,
                                            path_Mu_sd, path_pElem)



geneInfo <- indvLvl_geneInfo[[hit_id]]
geneInfo <- left_join(geneInfo, bed)
paste0(geneInfo$chr, ':', geneInfo$position)

summary(donorInfo$totalMutNrs)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
87     780    1258    1587    2001    5956

donor_id totalMutNrs
<char>       <int>
  1:  DO35755         325
2:  DO48893        3976
3:  DO48969         644

geneInfo$var_ID
[1] "M_46755482" "M_46867373" "M_46845339"
funseq <- fread('../extdata/procInput/iDriverInputs/Scores/funseq_obsSNV_scores.tsv')
x <- funseq[which(funseq$var_ID %in% geneInfo$var_ID),]
3.48;Yes;VA=1:ZIC1:ENSG00000152977.5:+:nonsynonymous:1/1:ZIC1-001:ENST00000282928.4:1344_1036_346_S->R;REG;.;.;.;.;.;.;.;ZIC1;.;3;.;.;.


3.89;Yes;VA=1:ZIC1:ENSG00000152977.5:+:nonsynonymous:1/1:ZIC1-001:ENST00000282928.4:1344_853_285_A->T;REG;.;.;.;.;.;.;.;ZIC1;.;3;.;.;.


alpha-missense scores
M_46755482 >>> 0.9971
M_46845339 >>> 0.9244

iDriver
0.006260713 0.09923231


compositeDriver
gc19_pc.cds::gencode::ZIC1::ENSG00000152977.5 1.67686
sum(compositeDriver$`p-value` > 1) # 33
##########################################################3
zic1 <- mc3[which(mc3$Chromosome == '3'),]
zic1 <- zic1[which(zic1$Start_position %in% c(147130358, 147131233, 147128752)),]
path_funseqScores <- '../../../../Projects/bahari_work/funseq2/hg19_score.funseq216.bed.bgz'

library(seqminer)
positions <- c(147130358, 147131233, 147128752)
tbix_data <- data.frame(cbind('IDs' = c('DO35755', 'DO48893', 'DO48969'), 
                              'tbix' = paste0('3', ":", positions, "-", positions)))
scores <- tabix.read.table(path_funseqScores, tbix_data$tbix)
