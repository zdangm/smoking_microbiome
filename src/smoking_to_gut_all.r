
args<-as.numeric(commandArgs(TRUE))

library(TwoSampleMR)
library(readr)
library(data.table)
library(stringr)
library(MRcML)

bac_id_files<-list.files("/data/coxvgi/zhoud2/projects/gut_addiction/gwas/")
bac_id<-bac_id_files[startsWith(bac_id_files,"genus")]
bac_id<-str_remove(bac_id,".summary.txt.gz")
bac_id_files<-list.files("/data/coxvgi/zhoud2/projects/gut_addiction/gwas/")
bac_id_genus<-bac_id_files[startsWith(bac_id_files,"genus")]
bac_id_phylum<-bac_id_files[startsWith(bac_id_files,"phylum")]
bac_id_class<-bac_id_files[startsWith(bac_id_files,"class")]
bac_id_order<-bac_id_files[startsWith(bac_id_files,"order")]
bac_id_family<-bac_id_files[startsWith(bac_id_files,"family")]
bac_id<-c(bac_id_genus,bac_id_phylum,bac_id_class,bac_id_order,bac_id_family)
bac_id<-str_remove(bac_id,".summary.txt.gz")


exposure_trait<-c("CigarettesPerDay","AgeOfInitiation","SmokingInitiation","SmokingCessation","Lifetime")[args]
outcome_trait<-"gut"
result_path<-"/data/coxvgi/zhoud2/projects/gut_addiction/mr/results/" 
exposure_source<-paste0("/data/coxvgi/zhoud2/data/pheno/smoking/2019/",exposure_trait,".txt.gz")


#select exposure SNP:beta(OR):se:effect_allele:other_allele:(or eaf):pval colunm index from exposure raw data
exp_OR_input<-"no"

if(exposure_trait == 'Lifetime'){
  exp_col_select<-c(1,8,9,4,5,10)
}else{
  exp_col_select<-c(3,9:10,5,4,8)
}

if(exposure_trait == 'Lifetime'){
  sample_size = 462690
}else if(exposure_trait == 'CigarettesPerDay'){
  sample_size = 337334
}else if(exposure_trait == 'AgeOfInitiation'){
  sample_size = 341427
}else if(exposure_trait == 'SmokingInitiation'){
  sample_size = 1232091
}else if(exposure_trait == 'SmokingCessation'){
  sample_size = 547219
}


#select outcome SNP:beta(OR):se:effect_allele:other_allele:(or eaf):pval colunm index from outcome raw data
out_OR_input<-"no"
out_col_select<-c(4,7:8,6,5,10) 
######## set p value to filter snps
p_select<- 1e-6 #5e-8
######## p value in exposure data need to adjust or not, with "yes" or "no" two options
p_adjust<-"no"
clump_kb<-250
clump_r2<-0.01
get_pos_col_exposure_sep<-":"


################### exposure raw data input
exposure_raw_data<-fread(exposure_source)
exposure_raw_data<-as.data.frame(exposure_raw_data)
#data.frame(exposure_col=colnames(exposure_raw_data),index=1:ncol(exposure_raw_data))
print("the number of snps in exposure_raw_data:")
nrow(exposure_raw_data)

########## select exposure SNP,beta,se,effect_allele,other_allele,eaf,pval coloum
exposure_raw_data<-exposure_raw_data[,exp_col_select]
print("select exposure SNP:beta:se:effect_allele:other_allele:(or eaf):pval colunm from exposure raw data")
colnames(exposure_raw_data)
exposure_colnames_without_eaf<-c("exp_SNP","exp_beta","exp_se","exp_effect_allele","exp_other_allele","exp_pval")
exposure_colnames_with_eaf<-c("exp_SNP","exp_beta","exp_se","exp_effect_allele","exp_other_allele","exp_eaf","exp_pval")
if(length(exp_col_select)==6){
  colnames(exposure_raw_data)<-exposure_colnames_without_eaf} else{
    colnames(exposure_raw_data)<-exposure_colnames_with_eaf
  }
exposure_raw_data$exp_effect_allele<-str_to_upper(as.character(exposure_raw_data$exp_effect_allele))
exposure_raw_data$exp_other_allele<-str_to_upper(as.character(exposure_raw_data$exp_other_allele))
if(exp_OR_input=="yes")
{
  exposure_raw_data$exp_beta<-log(exposure_raw_data$exp_beta, 2.718)
}

######### filter snp without beta, se value
exposure_raw_data<-exposure_raw_data[complete.cases(exposure_raw_data[,c(2:3)]),]
############# select variant with pvalue <p_select
if(p_adjust=="yes")
{
  exposure_filter_pvalue = exposure_raw_data[p.adjust(exposure_raw_data$exp_pval,method = 'BH')<p_select,]
}else
{
  exposure_filter_pvalue = exposure_raw_data[exposure_raw_data$exp_pval<p_select,]
}
if(length(which(exposure_filter_pvalue$exp_SNP =="."))>0)
{
  exposure_filter_pvalue<-exposure_filter_pvalue[-1*which(exposure_filter_pvalue$exp_SNP =="."),]
}
print("the number of snps with p_value< p select:")
nrow(exposure_filter_pvalue)



library(foreach)  
library(doParallel)   


bac_i = 1

mr_func = function(bac_i){
  
  
  outcome_source<-paste0("/data/coxvgi/zhoud2/projects/gut_addiction/gwas/",bac_id[bac_i],".summary.txt.gz")
  outcome_fread_sep<-ifelse(!str_detect(outcome_source,".csv"),"\t",",")
  
  outcome_trait_bac<-paste0("gut_addiction","_",bac_id[bac_i])
  
  #outcome raw data input
  outcome_raw_data<-fread(outcome_source,sep =outcome_fread_sep,header = T,stringsAsFactors=F)
  outcome_raw_data<-as.data.frame(outcome_raw_data)
  data.frame(outcome_col=colnames(outcome_raw_data),index=1:ncol(outcome_raw_data))
  print("the number of snps in outcome_raw_data:")
  nrow(outcome_raw_data)
  
  ########## select outcome SNP,beta,se,effect_allele,other_allele,eaf,pval coloum
  outcome_raw_data<-outcome_raw_data[,out_col_select]
  print("select outcome SNP:beta:se:effect_allele:other_allele:(or eaf):pval colunm from outcome raw data")
  colnames(outcome_raw_data)
  outcome_colnames_without_eaf<-c("out_SNP","out_beta","out_se","out_effect_allele","out_other_allele","out_pval")
  colnames(outcome_raw_data)<-outcome_colnames_without_eaf
  outcome_raw_data$out_effect_allele<-str_to_upper(as.character(outcome_raw_data$out_effect_allele))
  outcome_raw_data$out_other_allele<-str_to_upper(as.character(outcome_raw_data$out_other_allele))
  if(out_OR_input=="yes")
  {
    outcome_raw_data$out_beta<-log(outcome_raw_data$out_beta, 2.718)
  }
  
  ### select exposure variants with SNP code which also were contained in outcome
  merge_exposure_and_outcome<-merge(exposure_filter_pvalue,outcome_raw_data,by.x="exp_SNP",by.y="out_SNP")
  
  #df<-merge_exposure_and_outcome[,c(2:4,11:13,16:21)]
  df<-merge_exposure_and_outcome
  #df$Beta = ifelse(df$ALT_biovu == df$EA, df$Beta, df$Beta*-1) #flip allele

  ########## build exposure data frame
    exposure_df <- data.frame(
      SNP = df$exp_SNP,
      beta = df$exp_beta,
      se = df$exp_se,
      effect_allele = df$exp_effect_allele,
      other_allele=df$exp_other_allele,
      #eaf=df$exp_eaf,
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
      #eaf=df$out_eaf,
      pval=df$out_pval,
      Phenotype=outcome_trait_bac
    )
    
  exposure_dat <- format_data(exposure_df, type="exposure")
  print("the number of snps after mergeing exposue_df and outcome_df:")
  nrow(exposure_dat)
  
  outcome_dat <- format_data(outcome_df, type="outcome")
  
  ############ harmonise
  dat <- harmonise_data(
    exposure_dat = exposure_dat, 
    outcome_dat = outcome_dat,
    action=2
  )
  print("the number of snps after harmonise:")
  length(which(dat$mr_keep=="TRUE"))
  dat <-dat[dat$mr_keep=="TRUE",]
  
  #Clumping
  exposure_df <- data.frame(
    SNP = dat$SNP,
    beta = dat$beta.exposure,
    se = dat$se.exposure,
    effect_allele = dat$effect_allele.exposure,
    other_allele=dat$other_allele.exposure,
    pval=dat$pval.exposure,
    eaf=dat$eaf.exposure,
    Phenotype=exposure_trait
  )
  exposure_dat <- format_data(exposure_df, type="exposure")
  exposure_dat <- clump_data(exposure_dat,clump_kb=clump_kb,clump_r2 = clump_r2)
  print("the number of snps after clumping:")
  nrow(exposure_dat)
  
  ###### extract outcome data
  dat <- dat[dat$SNP %in% exposure_dat$SNP,]
  print("the number of mr_input snps:")
  length(which(dat$mr_keep=="TRUE"))
  
  
  #run MR
  ans_regular = mr(dat)
  
  ans_egger = mr_egger_regression(b_exp = dat$beta.exposure,
                                  se_exp = dat$se.exposure,
                                  b_out = dat$beta.outcome,
                                  se_out = dat$se.outcome)

  #run MRcML
  ans_cML = mr_cML(dat$beta.exposure,
                   dat$beta.outcome,
                   dat$se.exposure,
                   dat$se.outcome,
                   n=sample_size,
                   random_start=100,
                   random_seed=2023)
  
  ans_presso<-run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
  
  #summarize results
  
  ans = list()
  
  ans$exposure = as.character(ans_regular[1,'exposure'])
  ans$outcome = as.character(ans_regular[1,'outcome'])
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

  ans$cML_beta = ans_cML[["MA_BIC_theta"]]
  ans$cML_se = ans_cML[["MA_BIC_se"]]
  ans$cML_p = ans_cML[["MA_BIC_p"]]

  ans$presso_p_global = ans_presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
  ans$presso_beta_raw = ans_presso[[1]]$`Main MR results`[1,'Causal Estimate']
  ans$presso_p_raw = ans_presso[[1]]$`Main MR results`[1,'P-value']
  ans$presso_beta_corrected = ans_presso[[1]]$`Main MR results`[2,'Causal Estimate']
  ans$presso_p_corrected = ans_presso[[1]]$`Main MR results`[2,'P-value']

  return(ans)
  
}


registerDoParallel(cores = detectCores())

trials <- length(bac_id)

result <- foreach(i = 1:trials, .combine=rbind) %dopar% {
  
  mr_ans = mr_func(i)
  cat(paste0('INFO ',i,' completed \n'))
  unlist(mr_ans)
}

stopImplicitCluster()


write.table(result,paste0(result_path,"/",exposure_trait,"_to_",outcome_trait,".txt"), quote =F,sep = '\t',row.names = F)





