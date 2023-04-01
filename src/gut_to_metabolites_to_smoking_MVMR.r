
##ieu ID
### Age Of Smoking Initiation           ieu-b-24
### cigarettes per day                  ieu-b-25
### lifetime smoking
### smoking initiation                  ieu-b-4877
### Smoking cessation                   bbj-a-82

###Tryptophan色氨酸                     met-a-304
###tyrosine酪氨酸                       met-a-325  (met-c-938)2016  (met-d-Tyr)
###phenylalanine苯丙氨酸                met-a-308  (met-c-919)2016  (met-d-Phe)
###glutamate谷氨酸盐                    met-a-466
###glycine甘氨酸                        met-a-468
###SCFA-valerate丙酸                    met-a-575

## phylum.Actinobacteria                ebi-a-GCST90017110
## class.Actinobacteria                 ebi-a-GCST90016908  
## order.Bifidobacteriales              ebi-a-GCST90017093
## family.Bifidobacteriaceae            ebi-a-GCST90016929
## genus.Bifidobacterium                ebi-a-GCST90016970
## genus.Peptococcus                    ebi-a-GCST90017042

#########MVMR adjusted for six neurotransmitter-associated or bacterial metabolites simultaneously######
library(data.table)
library(TwoSampleMR)
library(openxlsx)
library(MendelianRandomization)

lifetime <- fread('E:/science/paper/gut_and_smoke/data/2019.10.02 Lifetime Smoking GWAS Data Sheet 1.txt')
head(lifetime)
lifetime$phenotype <- 'lifetime'

gut <- c('ebi-a-GCST90017110', 'ebi-a-GCST90016908', 'ebi-a-GCST90017093', 'ebi-a-GCST90016929', 'ebi-a-GCST90016970', 'ebi-a-GCST90017042')

table <- data.frame(gut=gut, M1=c(rep('met-a-304', time=length(gut))), 
                             M2=c(rep('met-a-325', time=length(gut))),
                             M3=c(rep('met-a-308', time=length(gut))),
                             M4=c(rep('met-a-466', time=length(gut))),
                             M5=c(rep('met-a-468', time=length(gut))),
                             M6=c(rep('met-a-575', time=length(gut))),
                             smoking=c(rep('ieu-b-24',time=length(gut)),
                                       rep('ieu-b-25',time=length(gut))))


for (i in 1:length(table[,1])){

 id_exposure <- c(table[i,1], table[i,2], table[i,3], table[i,4], table[i,5], table[i,6], table[i,7])
 id_outcome <- table[i,8]

 exposure_mvdat <- mv_extract_exposures(id_exposure, 
                                       clump_r2 = 0.01,
                                       clump_kb = 250,
                                       pval_threshold = 1e-06)
 outcome_mvdat <- extract_outcome_data(exposure_mvdat$SNP, id_outcome)
 #outcome_mvdat <- format_data(dat=lifetime,
                             #type="outcome",
                             #snps=exposure_mvdat$SNP,
                             #header=TRUE,
                             #phenotype_col='phenotype',
                             #snp_col='SNP',
                             #beta_col='BETA',
                             #se_col='SE',
                             #effect_allele_col='EFFECT_ALLELE',
                             #other_allele_col='OTHER_ALLELE',
                             #eaf_col='EAF',
                             #pval_col='P')
 mvdat <- mv_harmonise_data(exposure_mvdat, outcome_mvdat)
 exposure <- as.data.frame(mvdat[["exposure_beta"]])

 MRMVInputObject <- mr_mvinput(bx=cbind (exposure[,1],
                                        exposure[,2],
                                        exposure[,3],
                                        exposure[,4],
                                        exposure[,5],
                                        exposure[,6],
                                        exposure[,7]),
                              bxse=cbind (exposure[,1],
                                          exposure[,2],
                                          exposure[,3],
                                          exposure[,4],
                                          exposure[,5],
                                          exposure[,6],
                                          exposure[,7]),
                              by=mvdat[["outcome_beta"]],
                              byse=mvdat[["outcome_se"]])

 ##run MVMR-IVW
 IVW <- mr_mvivw(MRMVInputObject)
 IVW
 IVW@Heter.Stat

 ##run MVMR-Egger
 Egger <- mr_mvegger(MRMVInputObject)
 Egger
 Egger@Pvalue.Int

 #summarize results
 ans = list()
  
   ans$exposure = table[i,1]
   ans$outcome = table[i,8]
   #ans$outcome = Lifetime

   ans$MVMR_IVW_beta = IVW@Estimate[1]
   ans$MVMR_IVW_se = IVW@StdError[1]
   ans$MVMR_IVW_p = IVW@Pvalue[1]
   ans$MVMR_IVW_heterogeneity =  IVW@Heter.Stat[2]

   ans$MVMR_Egger_beta = Egger@Estimate[1]
   ans$MVMR_Egger_se = Egger@StdError.Est[1]
   ans$MVMR_Egger_p = Egger@Pvalue.Est[1]
   ans$egger_p_intercept = Egger@Pvalue.Int
   ans$Q_pval = Egger@Heter.Stat[2]

 #output results
  if(i==1){
   results <- as.data.frame(ans)
 }else{
   results <- rbind(results, data.frame(ans))
 }

}
write.xlsx(results, 'E:/science/paper/gut_and_smoke/data/results_MVMR_1.xlsx')

#########MVMR adjusted for six neurotransmitter-associated or bacterial metabolites separately######
library(data.table)
library(TwoSampleMR)
library(openxlsx)
library(MendelianRandomization)

lifetime <- fread('E:/science/paper/gut_and_smoke/data/2019.10.02 Lifetime Smoking GWAS Data Sheet 1.txt')
head(lifetime)
lifetime$phenotype <- 'lifetime'

gut <- c('ebi-a-GCST90016908','ebi-a-GCST90017042')
table <- data.frame(gut=gut, M=c(rep('met-a-304', time=length(gut)), 
                                 rep('met-a-325', time=length(gut)),
                                 rep('met-a-308', time=length(gut)),
                                 rep('met-a-466', time=length(gut)),
                                 rep('met-a-468', time=length(gut)),
                                 rep('met-a-575', time=length(gut))))

table <- table[order(table$gut),]

for (i in 1:length(table[,1])){

 id_exposure <- c(table[i,1], table[i,2])

 exposure_mvdat <- mv_extract_exposures(id_exposure, 
                                       clump_r2 = 0.01,
                                       clump_kb = 250,
                                       pval_threshold = 1e-06)
outcome_mvdat <- format_data(dat=lifetime,
                             type="outcome",
                             snps=exposure_mvdat$SNP,
                             header=TRUE,
                             phenotype_col='phenotype',
                             snp_col='SNP',
                             beta_col='BETA',
                             se_col='SE',
                             effect_allele_col='EFFECT_ALLELE',
                             other_allele_col='OTHER_ALLELE',
                             eaf_col='EAF',
                             pval_col='P')
 mvdat <- mv_harmonise_data(exposure_mvdat, outcome_mvdat)
 exposure <- as.data.frame(mvdat[["exposure_beta"]])

 MRMVInputObject <- mr_mvinput(bx=cbind (exposure[,1],
                                        exposure[,2]),
                              bxse=cbind (exposure[,1],
                                          exposure[,2]),
                              by=mvdat[["outcome_beta"]],
                              byse=mvdat[["outcome_se"]])

 ##run MVMR-IVW
 IVW <- mr_mvivw(MRMVInputObject)
 IVW
 IVW@Heter.Stat

 ##run MVMR-Egger
 Egger <- mr_mvegger(MRMVInputObject)
 Egger
 Egger@Pvalue.Int

 #summarize results
 ans = list()

   ans$exposure = table[i,1]
   ans$outcome = Lifetime

   ans$MVMR_IVW_beta = IVW@Estimate[1]
   ans$MVMR_IVW_se = IVW@StdError[1]
   ans$MVMR_IVW_p = IVW@Pvalue[1]
   ans$MVMR_IVW_heterogeneity =  IVW@Heter.Stat[2]

   ans$MVMR_Egger_beta = Egger@Estimate[1]
   ans$MVMR_Egger_se = Egger@StdError.Est[1]
   ans$MVMR_Egger_p = Egger@Pvalue.Est[1]
   ans$egger_p_intercept = Egger@Pvalue.Int
   ans$Q_pval = Egger@Heter.Stat[2]

 #output results
 if(i==1){
   results <- as.data.frame(ans)
 }else{
   results <- rbind(results, data.frame(ans)) 
 }

}
write.xlsx(results, 'E:/science/paper/gut_and_smoke/data/results_MVMR_2.xlsx')