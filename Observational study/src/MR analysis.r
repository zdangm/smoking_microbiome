
library(dplyr)
library(tidyr)
library(optparse)

option_list = list(
  make_option("--taxa_id", action="store", default=NA, type='integer')
)
opt = parse_args(OptionParser(option_list=option_list))

library(TwoSampleMR)
library(readr)
library(readxl)
library(data.table)
library(stringr)
library(ieugwasr)

exposure_trait <- c('g__Actinomyces','g__Atopobium','g__Haemophilus','g__Lachnospira','g__Pseudopropionibacterium','g__Turicibacter')[opt$taxa_id]
BBJ_trait <- read_excel('~/BBJ/BBJ_trait.xlsx') %>% as.data.frame()
all_trait <- BBJ_trait$Trait
trait_conti <- BBJ_trait$Trait[which(BBJ_trait$Status=='conti')]
outcome_trait <- BBJ_trait$Trait
#outcome_trait <- c('A02B','AIH','AP','CG','Cho','ChP','Cir','CL','Co','CP','GBP','GERD','GP','GU','Ile','InH','PeD','UC')
#outcome_trait <- c('A02B','Cho','COPD','CRP','GU','Ile','LuC','MI','N02A','NEU','PAD','PrC','SAP','TBil','UAP','WBC')

exposure_source <- paste0('~/Smoking_microbe/',exposure_trait,'_MAF_0.05.txt')
result_path <- "~/projects/smoking_microbiome_BBJ"

#select exposure SNP:beta(OR):se:effect_allele:other_allele:(or eaf):pval colunm index from exposure raw data
exposure_select<-c(10,6,7,3,4,5,8)
exposure_colnames <- c("exp_SNP","exp_beta","exp_se","exp_effect_allele","exp_other_allele","exp_eaf","exp_pval")
outcome_colnames <- c("out_SNP","out_beta","out_se","out_effect_allele","out_other_allele","out_eaf","out_pval")

######## set p value to filter snps and other parameters
p_select <- 5e-6
clump_kb <- 1000
clump_r2 <- 0.01

################### exposure raw data input
exposure_raw_data<-as.data.frame(fread(exposure_source))
print("the number of snps in exposure_raw_data:")
nrow(exposure_raw_data)

########## select exposure SNP,beta,se,effect_allele,other_allele,eaf,pval coloum
exposure_raw_data<-exposure_raw_data[,exposure_select]
colnames(exposure_raw_data)<- exposure_colnames
exposure_raw_data$exp_effect_allele <- str_to_upper(as.character(exposure_raw_data$exp_effect_allele))
exposure_raw_data$exp_other_allele <- str_to_upper(as.character(exposure_raw_data$exp_other_allele))

######### filter snp without RSID, beta, se value
exposure_raw_data <- exposure_raw_data[complete.cases(exposure_raw_data[,c(1:3)]),]
nrow(exposure_raw_data)

############# select variant with pvalue < p_select
exposure_filter_pvalue <- exposure_raw_data[exposure_raw_data$exp_pval<p_select,]

library(foreach)  
library(doParallel) 

id = 1

mr_func = function(id){

  if (outcome_trait[id] %in% trait_conti){
  outcome_select<-c(2,12,13,6,7,8,15)
  }else{
  outcome_select<-c(4,11,12,6,5,8,14)
  }
  outcome_source <- paste0('~/BBJ/rawdata/hum0197.v3.BBJ.',outcome_trait[id],'.v1/GWASsummary_',outcome_trait[id],'_Japanese_SakaueKanai2020.auto.txt.gz')
  
    ans_NA = list()

    ans_NA$exposure = exposure_trait[id]
    ans_NA$outcome = outcome_trait

    ans_NA$nsnp = NA

    ans_NA$ivw_beta = NA
    ans_NA$ivw_se = NA
    ans_NA$ivw_p = NA

    ans_NA$egger_beta = NA
    ans_NA$egger_se = NA
    ans_NA$egger_p = NA
    ans_NA$egger_p_intercept = NA
    ans_NA$Q_pval = NA

    ans_NA$wme_beta = NA
    ans_NA$wme_se = NA
    ans_NA$wme_p = NA

    ans_NA$presso_p_global = NA
    ans_NA$presso_beta_raw = NA
    ans_NA$presso_p_raw = NA
    ans_NA$presso_beta_corrected = NA
    ans_NA$presso_p_corrected = NA
  
  ################### outcome raw data input
  outcome_raw_data<-as.data.frame(fread(outcome_source))
  print("the number of snps in outcome_raw_data:")
  nrow(outcome_raw_data)

  ########## select outcome SNP,beta,se,effect_allele,other_allele,eaf,pval coloum
  outcome_raw_data<-outcome_raw_data[,outcome_select]
  colnames(outcome_raw_data)<- outcome_colnames
  outcome_raw_data$out_effect_allele <- str_to_upper(as.character(outcome_raw_data$out_effect_allele))
  outcome_raw_data$out_other_allele <- str_to_upper(as.character(outcome_raw_data$out_other_allele))

  ######### filter snp without RSID, beta, se value
  outcome_raw_data <- outcome_raw_data[complete.cases(outcome_raw_data[,c(1:3)]),]
  nrow(outcome_raw_data)

  ############# select exposure variants with SNP code which also were contained in outcome
  merge_exposure_and_outcome <- merge(exposure_filter_pvalue,outcome_raw_data,by.x="exp_SNP",by.y="out_SNP")
  df <- merge_exposure_and_outcome

  if (nrow(df)>=1){

     ########## build exposure data frame
        exposure_df <- data.frame(
          SNP = df$exp_SNP,
          beta = df$exp_beta,
          se = df$exp_se,
          effect_allele = df$exp_effect_allele,
          other_allele=df$exp_other_allele,
          eaf=df$exp_eaf,
          pval=df$exp_pval,
          Phenotype=exposure_trait
        )

     ########## build outcome data frame
        outcome_df <- data.frame(
          SNP = df$exp_SNP,
          beta = df$out_beta,
          se = df$out_se,
          effect_allele = df$out_effect_allele,
          other_allele=df$out_other_allele,
          eaf=df$out_eaf,
          pval=df$out_pval,
          Phenotype=outcome_trait[id]
        )

      exposure_dat <- format_data(exposure_df, type="exposure")
      outcome_dat <- format_data(outcome_df, type="outcome")
      print("the number of snps after mergeing exposue_df and outcome_df:")
      nrow(exposure_dat)
      nrow(outcome_dat)
      exposure_dat <- exposure_dat[complete.cases(exposure_dat[,c('effect_allele.exposure','other_allele.exposure')]),]
      outcome_dat <- outcome_dat[complete.cases(outcome_dat[,c('effect_allele.outcome','other_allele.outcome')]),]

      ############ harmonise
        dat <- harmonise_data(
          exposure_dat = exposure_dat, 
          outcome_dat = outcome_dat,
          action=2
        )
      print("the number of snps after harmonise:")
      length(which(dat$mr_keep=="TRUE"))
      dat <-dat[dat$mr_keep=="TRUE",]

     if (nrow(dat)>=1){

        #Clumping
        #dat <- clump_data(dat, clump_kb=clump_kb, clump_r2=clump_r2)
        #print("the number of snps after clumping:")
        #nrow(dat)

        dat_clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$exposure),clump_kb=clump_kb, clump_r2=clump_r2,
                              pop='EAS',
                              plink_bin = genetics.binaRies::get_plink_binary(),
                              bfile = "~/projects/addi_scan/other/reference/EAS/EAS")
        dat <- dat[which(dat$SNP %in% dat_clump$rsid),]
        print("the number of snps after clumping:")
        nrow(dat)
  
        if (nrow(dat)>=1){
 
          if (nrow(dat)<=2){

          #run MR
          ans_regular = mr(dat)

          #summarize results
          ans = list()

          ans$exposure = as.character(exposure_trait)
          ans$outcome = as.character(outcome_trait[id])
          ans$nsnp = ans_regular[1,'nsnp']

          ans$ivw_beta = ans_regular$b
          ans$ivw_se = ans_regular$se
          ans$ivw_p = ans_regular$pval

          ans$egger_beta = NA
          ans$egger_se = NA
          ans$egger_p = NA
          ans$egger_p_intercept = NA
          ans$Q_pval = NA

          ans$wme_beta = NA
          ans$wme_se = NA
          ans$wme_p = NA

          ans$presso_p_global = NA
          ans$presso_beta_raw = NA
          ans$presso_p_raw = NA
          ans$presso_beta_corrected = NA
          ans$presso_p_corrected = NA

          }else if (nrow(dat)==3){
         
          #run MR
          ans_regular = mr(dat)
  
          ans_egger = mr_egger_regression(b_exp = dat$beta.exposure,
                                          se_exp = dat$se.exposure,
                                          b_out = dat$beta.outcome,
                                          se_out = dat$se.outcome)
        
          #summarize results
          ans = list() 

          ans$exposure = as.character(exposure_trait)
          ans$outcome = as.character(outcome_trait[id])
          ans$nsnp = ans_regular[1,'nsnp']

          row_index_ivw = which(ans_regular[,'method'] == 'Inverse variance weighted')
          ans$ivw_beta = ans_regular[row_index_ivw,'b']
          ans$ivw_se = ans_regular[row_index_ivw,'se']
          ans$ivw_p = ans_regular[row_index_ivw,'pval']

          row_index_egger = which(ans_regular[,'method'] == 'MR Egger')
          ans$egger_beta = ans_regular[row_index_egger,'b']
          ans$egger_se = ans_regular[row_index_egger,'se']
          ans$egger_p = ans_regular[row_index_egger,'pval']
          ans$egger_p_intercept = ans_egger$pval_i
          ans$Q_pval = ans_egger$Q_pval

          row_index_wme = which(ans_regular[,'method'] == 'Weighted median')
          ans$wme_beta = ans_regular[row_index_wme,'b']
          ans$wme_se = ans_regular[row_index_wme,'se']
          ans$wme_p = ans_regular[row_index_wme,'pval']
    
          ans$presso_p_global = NA
          ans$presso_beta_raw = NA
          ans$presso_p_raw = NA
          ans$presso_beta_corrected = NA
          ans$presso_p_corrected = NA
            
          }else if (nrow(dat)>3){

          #run MR
          ans_regular = mr(dat)
  
          ans_egger = mr_egger_regression(b_exp = dat$beta.exposure,
                                          se_exp = dat$se.exposure,
                                          b_out = dat$beta.outcome,
                                          se_out = dat$se.outcome)
    
          ans_presso<-run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)

          #summarize results
          ans = list()

          ans$exposure = as.character(exposure_trait)
          ans$outcome = as.character(outcome_trait[id])
          ans$nsnp = ans_regular[1,'nsnp']

          row_index_ivw = which(ans_regular[,'method'] == 'Inverse variance weighted')
          ans$ivw_beta = ans_regular[row_index_ivw,'b']
          ans$ivw_se = ans_regular[row_index_ivw,'se']
          ans$ivw_p = ans_regular[row_index_ivw,'pval']

          row_index_egger = which(ans_regular[,'method'] == 'MR Egger')
          ans$egger_beta = ans_regular[row_index_egger,'b']
          ans$egger_se = ans_regular[row_index_egger,'se']
          ans$egger_p = ans_regular[row_index_egger,'pval']
          ans$egger_p_intercept = ans_egger$pval_i
          ans$Q_pval = ans_egger$Q_pval

          row_index_wme = which(ans_regular[,'method'] == 'Weighted median')
          ans$wme_beta = ans_regular[row_index_wme,'b']
          ans$wme_se = ans_regular[row_index_wme,'se']
          ans$wme_p = ans_regular[row_index_wme,'pval']
    
          ans$presso_p_global = ans_presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
          ans$presso_beta_raw = ans_presso[[1]]$`Main MR results`[1,'Causal Estimate']
          ans$presso_p_raw = ans_presso[[1]]$`Main MR results`[1,'P-value']
          ans$presso_beta_corrected = ans_presso[[1]]$`Main MR results`[2,'Causal Estimate']
          ans$presso_p_corrected = ans_presso[[1]]$`Main MR results`[2,'P-value']

          }
        }else{
        ans = ans_NA
        }
      }else{
      ans = ans_NA
      }
    }else{
    ans = ans_NA
    }    

  return(ans)

}


registerDoParallel(cores = detectCores()/10)

trials <- length(outcome_trait)
result <- foreach(i = 1:trials, .combine=rbind) %dopar% {
  
  mr_ans = mr_func(i)
  cat(paste0('INFO ',i,' completed \n'))
  unlist(mr_ans)
}

stopImplicitCluster()

result <- as.data.frame(result)
result <- merge(BBJ_trait, result, by.x='Trait', by.y='outcome', all.x=F,all.y=T)
result$BH <- p.adjust(result$ivw_p, method='BH',n=length(result$ivw_p))
result$bonferroni <- p.adjust(result$ivw_p, method='bonferroni',n=length(result$ivw_p))
write.table(result,paste0(result_path,"/",exposure_trait,"_to_BBJ_",p_select,".txt"), quote =F,sep = '\t',row.names = F)

