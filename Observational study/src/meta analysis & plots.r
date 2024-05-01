
############## Collation of MaAsLin results ############## (take ZMSC 16S for example)
library(data.table)
library(readxl)
library(ggvenn)
library(ggplot2)
library(openxlsx)

## mian results
select <- c('','_adjdrink')[1]

#Current smoker
ZMSC_raw_current <- fread(paste0('~/ZMSC/results/MaAsLin2_16S_smoke_male',select,'/all_results.tsv'))
ZMSC_raw_current$feature <- gsub("\\.", "_", ZMSC_raw_current$feature)
ZMSC_tidy_current <- ZMSC_raw_current[which(ZMSC_raw_current$metadata=='smoking_status'),][,c(1,3,4,5,8:9)]
colnames(ZMSC_tidy_current) <- c('Taxa','Smoking_status','Estimate','SE','P','Padj_BH')

#past smoker
ZMSC_raw_past <- fread(paste0('~/ZMSC/results/MaAsLin2_16S_smoke_male_pastsmoker',select,'/all_results.tsv'))
ZMSC_raw_past$feature <- gsub("\\.", "_", ZMSC_raw_past$feature)
ZMSC_tidy_past <- ZMSC_raw_past[which(ZMSC_raw_past$metadata=='smoking_status'),][,c(1,3,4,5,8:9)]
colnames(ZMSC_tidy_past) <- c('Taxa','Smoking_status','Estimate','SE','P','Padj_BH')

#Current+past
smoking <- c('Current smoker v.s. Non-smoker', 'Former smoker v.s. Non-smoker', 'Former smoker v.s. Current smoker')

ZMSC_tidy <- rbind(ZMSC_tidy_current, ZMSC_tidy_past)
ZMSC_tidy <- ZMSC_tidy[order(ZMSC_tidy$P),]
for (i in 1:length(smoking)){
  ZMSC_tidy$Padj_BH[which(ZMSC_tidy$Smoking_status==smoking[i])] <- p.adjust(ZMSC_tidy$P[which(ZMSC_tidy$Smoking_status==smoking[i])], 
                                                                    method='BH',
                                                                    n=length(ZMSC_tidy$P[which(ZMSC_tidy$Smoking_status==smoking[i])]))
  ZMSC_tidy$Padj_Bonferroni[which(ZMSC_tidy$Smoking_status==smoking[i])] <- p.adjust(ZMSC_tidy$P[which(ZMSC_tidy$Smoking_status==smoking[i])], 
                                                                            method='bonferroni',
                                                                            n=length(ZMSC_tidy$P[which(ZMSC_tidy$Smoking_status==smoking[i])]))
}

write.xlsx(ZMSC_tidy, file=paste0('~/Results/V1/MaAsLin_16S_ZMSC',select,'.xlsx'))

############## Collation of alpha/beta results ############## (take 16S for example)
  #### data input ####
microbiota <- c('16S','metagenome','pathway')[1]

load(paste0('~/ZMSC/results/alpha_Richness_',microbiota,'.Rdata')) #(take ZMSC for example)
ZMSC_Richness <- Richness                                          #(the same for Shannon, Simpson, Pielou, PCoA)

  #### plot patch ####
library(ggplot2)
library(patchwork)

pdf(paste0('~/Results/V1/alpha_beta_',microbiota,'_results.pdf'), width = 16, height = 10)
ZMSC_Pielou + GNHS_Pielou + ZMSC_PCoA + GNHS_PCoA +
GGMP_Pielou + SRRSHS_Pielou + GGMP_PCoA + SRRSHS_PCoA +
  plot_layout(ncol=4, widths =c(1,1,1.5,1.5), guides = 'collect') +
  plot_annotation() &
  theme(legend.position='bottom',
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
dev.off()


############## 16S Meta analysis ############## (take 16S Current-smokers v.s. Non-smokers for example)
  ##### Venn plot ####
library(data.table)
library(readxl)
library(ggvenn)
library(ggplot2)
library(openxlsx)

ZMSC_raw <- fread('~/ZMSC/results/MaAsLin2_16S_smoke_male/all_results.tsv')  #(the same for GNHS, GGMP, SRRSHS)
ZMSC_raw$feature <- gsub("\\.", "_", ZMSC_raw$feature)
ZMSC_raw <- ZMSC_raw[which(ZMSC_raw$value=='Current smoker'),][,c(1,3,4,5,8)]

venn_data <- list(GGMP = GGMP_raw$feature[which(GGMP_raw$pval<0.05)],
                  ZMSC = ZMSC_raw$feature[which(ZMSC_raw$pval<0.05)],
                  GNHS = GNHS_raw$feature[which(GNHS_raw$pval<0.05)], 
                  SRRSHS=SRRSHS_raw$feature[which(SRRSHS_raw$pval<0.05)])
pdf(paste0("~/Meta/Venn_16S_sign.pdf"), width = 8, height = 6)
ggvenn(venn_data, fill_color = c('#247BA0', '#DB7093', '#4CAF50','#F2B134'),
       fill_alpha = 0.45, text_size=6.5, set_name_size=6.5) #8-10
dev.off()

  ##### meta analysis ####
##data processing
ZMSC_raw$cohort <- 'ZMSC'
GNHS_raw$cohort <- 'GNHS'
GGMP_raw$cohort <- 'GGMP'
SRRSHS_raw$cohort <- 'SirRR'

cohort_all <- merge(ZMSC_raw, GNHS_raw, by='feature', all=T)
cohort_all <- merge(cohort_all, GGMP_raw, by='feature', all=T)
cohort_all <- merge(cohort_all, SRRSHS_raw, by='feature', all=T)

ZMSC_meta <- cohort_all[,c(1:6)]
colnames(ZMSC_meta) <- c('Taxa','value','coef','stderr','pval','cohort')  #(the same for GNHS, GGMP, SRRSHS)
ZMSC_meta$cohort <- 'ZMSC'
ZMSC_meta$No.sample <- 383
ZMSC_meta$No.smoker <- 161

cohort_meta <- rbind(ZMSC_meta, GNHS_meta, GGMP_meta, SRRSHS_meta)
cohort_meta <- cohort_meta[order(cohort_meta$Taxa),]

##analysis
library(meta)
library(ggplot2)
library(forestplot)
library(readxl)
library(haven)

meta <- metagen(data=cohort_meta, 
                studlab=paste0(cohort_meta$cohort,'  ', cohort_meta$No.smoker,'/', cohort_meta$No.sample,'  current_smokers'),
                #cluster=Taxa,
                coef, seTE=stderr,
                sm='β',
                subgroup = Taxa)

meta_common <- cbind(meta[["TE.common.w"]],meta[["seTE.common.w"]], meta[["pval.common.w"]], meta[["I2.w"]], meta[["tau2.w"]], meta[["pval.Q.w"]]) %>% as.data.frame()
colnames(meta_common) <- c('Estimate.fix', 'SE.fix', 'Pval', 'Heterogeneity.I2','Heterogeneity.tau2','Heterogeneity.Pval')
meta_common$Padj_BH <- p.adjust(meta_common$Pval, method='BH', n=length(meta_common$Pval))
meta_common$Padj_BF <- p.adjust(meta_common$Pval, method='bonferroni', n=length(meta_common$Pval))
meta_common[which(meta_common$Padj_BF<0.05),]

write.xlsx(meta_common, '~/Meta/Meta_16S_results.xlsx')

##forest plot
cohort_sign <- cohort_meta[which(cohort_meta$Taxa %in% meta_common$Taxa[which(meta_common$Padj_BH<0.05)][c(1,2,5,8,6)]),] 
meta_sign <- metagen(data=cohort_sign, 
                     studlab=paste0(cohort_sign$cohort,'  ', cohort_sign$No.smoker,'/', cohort_sign$No.sample,'  current_smokers'),
                     #cluster=Taxa,
                     coef, seTE=stderr,
                     sm='β',
                     subgroup = Taxa)
pdf(paste0("~/Meta/Forest_16S_smoker_BF.pdf"), width = 10, height = 10)
forest(meta_sign, overall = F, random = F, overall.hetstat = NULL, col.diamond.common='#770000', col.subgroup='#770000', print.stat = T) 
dev.off()


############## ZMSC regression analysis for health-related outcomes ##############
  #### data input ####
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library(openxlsx)

setwd("~/ZMSC/raw")
load('Demo_smoke_male_16S_tran.Rdata')
Demo_tran <- Demo
rm(Demo)

## bioche data
bioche <- read_excel('2020_all.xlsx') %>% as.data.frame()

bioche_name <- c('Systolic blood pressure(SBP)','Diastolic blood pressure(DBP)',
                 'Carotid artery thickness(left)', 'Carotid plaque count(left)','Carotid artery thickness(right)','Carotid plaque count(right)',
                 'Alanine transaminase(ALT)','Aspertate aminotransferase(AST)','AST/ALT',
                 'Total cholesterol','High density lipoprotein cholesterol(HDL)','Low density lipoprotein cholesterol(LDL)','Triglyceride',
                 'Fasting blood glucose','Glycated haemoglobin',
                 'Erythrocyte count','Hematokrit(HCT)','Haemoglobin','Mean corpuscular hemoglobin(MCH)','Mean corpuscular hemoglobin concentration(MCHC)','Mean corpuscular volume(MCV)',
                 'White blood cell count')

colnames(bioche) <- c('ID','Sample_ID','Name', c(bioche_name))

## dise data
load('Dise.RData')
Demo_tran$`Type 2 diabetes` <- ifelse(Demo_tran$身份证号码 %in% 糖尿病_case$AAE135[糖尿病_case$OPDATE<'2020-06-30 06:00:00'],
                                 'yes','no') %>% as.factor() # (take T2D for example)

dise_name <- c("Stroke","Coronary heart disease",
               "Peptic ulcer","Chronic gastritis","Chronic enteritis","Cholecystitis","Reflux disease",
               "Gastroenteritis","Gastrointestinal polyps","Constipation","Gastrointestinal hemorrhage",
               "Chronic lung disease","Chronic kidney disease","Type 2 diabetes")

Demo <- merge(Demo_tran, bioche[,c('Sample_ID', bioche_name)], by='Sample_ID', all.x=T, all.y=F)

  #### taxa > bioche regression ####
taxa <- c('g__Actinomyces','g__Haemophilus','g__Turicibacter','g__Lachnospira')

tmp_j <- list()
for (j in 1:length(taxa)){
  tmp_i <- list()
  for (i in 1:length(bioche_name)){
      fit <- glm(Demo[,bioche_name[i]]~log(Demo[,taxa[j]])+age+BMI, data=Demo)
      results <- as.data.frame(summary(fit)$coefficients)
      results$Beta <- results$Estimate
      results$Lower <- results$Estimate - 1.96*results$`Std. Error`      
      results$Upper <- results$Estimate + 1.96*results$`Std. Error`
      #results$OR <- exp(results$Estimate)                                     
      #results$lower <- exp(results$Estimate - 1.96*results$`Std. Error`)      
      #results$upper <- exp(results$Estimate + 1.96*results$`Std. Error`)
      tmp_i[[i]] <- results[2,] 
    }
  results_glm_taxa_bioche <- do.call(rbind, tmp_i)
  results_glm_taxa_bioche$Taxa <- taxa[j]
  results_glm_taxa_bioche$Phenotyps <- bioche_name
  results_glm_taxa_bioche$`Median±SD` <- paste0(sprintf("%0.2f",apply(Demo[,bioche_name],2,function(x){median(x,na.rm=T)})),
                                         "±", 
                                         sprintf("%0.2f",apply(Demo[,bioche_name],2,function(x){sd(x,na.rm=T)})))
  tmp_j[[j]] <- results_glm_taxa_bioche
}
results_glm_taxa_bioche <- do.call(rbind, tmp_j)
results_glm_taxa_bioche <- results_glm_taxa_bioche[,c(8:10,1,2,4)]
colnames(results_glm_taxa_bioche) <- c('Taxa','Phenotypes','Info.phenotypes','Estimate','SE','P')
results_glm_taxa_bioche$status <- 'continuous'

  #### taxa > dise regression ####
tmp_j <- list()
for (j in 1:length(taxa)){
  tmp_i <- list()
  for (i in 1:length(dise_name)){
    fit <- glm(Demo[,dise_name[i]]~log(Demo[,taxa[j]])+age+BMI, data=Demo, family=binomial)
    results <- as.data.frame(summary(fit)$coefficients)
    #results$beta <- results$Estimate
    #results$lower <- results$Estimate - 1.96*results$`Std. Error`      
    #results$upper <- results$Estimate + 1.96*results$`Std. Error`
    results$OR <- exp(results$Estimate)                                     
    results$Lower <- exp(results$Estimate - 1.96*results$`Std. Error`)      
    results$Upper <- exp(results$Estimate + 1.96*results$`Std. Error`)
    tmp_i[[i]] <- results[2,]
  }
  results_glm_taxa_dise <- do.call(rbind, tmp_i)
  results_glm_taxa_dise$Taxa <- taxa[j]
  results_glm_taxa_dise$Phenotyps <- dise_name
  results_glm_taxa_dise$case <- apply(Demo[,dise_name],2,function(x){sum(x=='yes')})
  tmp_j[[j]] <- results_glm_taxa_dise
}
results_glm_taxa_dise <- do.call(rbind, tmp_j)
results_glm_taxa_dise <- results_glm_taxa_dise[,c(8:10,1:2,4)]
colnames(results_glm_taxa_dise) <- c('Taxa','Phenotypes','Info.phenotypes','Estimate','SE','P')
results_glm_taxa_dise$status <- 'binary'

  #### bioche + dise  results ####
results_glm_taxa_trait <- rbind(results_glm_taxa_bioche, results_glm_taxa_dise)

  #### bioche + dise plot ####
## observation 
library(readxl)
library(ggplot2)

df <- results_glm_taxa_trait    #read_excel('~/ZMSC/results/glm_taxa_trait.xlsx')
df$sig <- ifelse(df$Padj_BH<0.05, '#',
                 ifelse(df$Padj_BH>=0.05 & df$P<0.01, '**',
                        ifelse(df$P>0.01 & df$P<0.05, '*', '')))

df$Phenotypes <- factor(df$Phenotypes, levels = c(dise_name, bioche_name))

pdf(paste0("~/Meta/taxa_phenotypes_zhoushan.pdf"), width = 6.6, height = 11)
ggplot(df,aes(x=Taxa, y=Phenotypes)) +
  geom_tile(color="#C2C2C2",fill="white",size=0.3)+
  geom_point(aes(fill=Estimate), pch=22, color="white",size=7)+
  scale_fill_gradient2(name = 'effect',
                       limit = c(min(df$Estimate),max(df$Estimate)),
                       low='#66C7B4',mid='#F2F2F2',high='#F4511E', 
                       midpoint = 0)+
  geom_text(aes(label=df$sig), color="black", size=5) +
  labs(x=NULL,y=NULL,title=NULL)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=11,color="black"),
        axis.text.x = element_text(angle=45,hjust=1),
        axis.text.y= element_text(color=c(rep('#C6403D', length(dise_name)),rep('#246BAE', length(bioche_name)))))+
  guides(fill=guide_colorbar(barheight = 30))
dev.off()

## Mendelian randomization
library(ggplot2)
library(patchwork)
library(readxl)

Actinomyces <- read_excel('~/BBJ/results/Taxa_BBJ.xlsx', sheet = 1) # (The same for Atopobium, Haemophilus, Lachnospira, Turicibacter)

results_MR_taxa_trait <- rbind(Actinomyces,Atopobium,Haemophilus,Lachnospira,Turicibacter)
results_MR_taxa_trait$Status <- ifelse(results_MR_taxa_trait$Status=='bin', 'binary', 'continuous')
results_MR_taxa_trait$sam <- ifelse(results_MR_taxa_trait$No.cases!='-', results_MR_taxa_trait$No.cases, results_MR_taxa_trait$No.samples) %>% as.numeric()

results_MR_taxa_trait <- results_MR_taxa_trait[which(results_MR_taxa_trait$sam >1800),]
results_MR_taxa_trait <- results_MR_taxa_trait[which(results_MR_taxa_trait$Phenotypes!='Breast cancer'),]
results_MR_taxa_trait <- results_MR_taxa_trait[which(results_MR_taxa_trait$Phenotypes!='Ovarian cyst'),]
results_MR_taxa_trait <- results_MR_taxa_trait[which(results_MR_taxa_trait$Phenotypes!='Uterine fibroid'),]
taxa <- unique(results_MR_taxa_trait$Taxa)

df <- results_MR_taxa_trait     #read_excel('~/ZMSC/results/MR_taxa_trait.xlsx')
df$sig <- ifelse(df$IVW_padj_BH<0.05, '#',
                 ifelse(df$IVW_padj_BH>=0.05 & df$IVW_p<0.01, '**',
                        ifelse(df$IVW_p>0.01 & df$IVW_p<0.05, '*', '')))
bioche_name <- df$Phenotypes[which(df$Status=='continuous')] %>% unique()
dise_name <- df$Phenotypes[which(df$Status=='binary')] %>% unique()

df$Phenotypes <- factor(df$Phenotypes, levels = c(dise_name, bioche_name))

pdf(paste0("~/Meta/taxa_phenotypes_MR.pdf"), width =15, height = 15.5)
p1 <- ggplot(df[which(df$Status=='continuous'),],aes(x=Taxa, y=Phenotypes)) +
  geom_tile(color="#C2C2C2",fill="white",size=0.3)+
  geom_point(aes(fill=IVW_beta), pch=22, color="white",size=7.5)+
  scale_fill_gradient2(name = 'effect',
                       limit = c(min(df$IVW_beta),max(df$IVW_beta)),
                       low='#66C7B4',mid='#F2F2F2',high='#F4511E', 
                       midpoint = 0)+
  geom_text(aes(label=sig), color="black", size=7.5) +
  labs(x=NULL,y=NULL,title=NULL)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=13,color="black"),
        axis.text.x = element_text(angle=45,hjust=1),
        axis.text.y= element_text(color=c('#246BAE')))+
  guides(fill=guide_colorbar(barheight = 30))
p2 <- ggplot(df[which(df$Status=='binary'),],aes(x=Taxa, y=Phenotypes)) +
  geom_tile(color="#C2C2C2",fill="white",size=0.3)+
  geom_point(aes(fill=IVW_beta), pch=22, color="white",size=6.5)+
  scale_fill_gradient2(name = 'effect',
                       limit = c(min(df$IVW_beta),max(df$IVW_beta)),
                       low='#66C7B4',mid='#F2F2F2',high='#F4511E', 
                       midpoint = 0)+
  geom_text(aes(label=sig), color="black", size=7.5) +
  labs(x=NULL,y=NULL,title=NULL)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=13,color="black"),
        axis.text.x = element_text(angle=45,hjust=1),
        axis.text.y= element_text(color=c('#C6403D')))+
  guides(fill=guide_colorbar(barheight = 30))

p1 + p2 + plot_layout(ncol=2, widths=c(1,1.3), guides = 'collect') +
  #plot_annotation(tag_levels = 'a') &
  theme(legend.position='right')
dev.off()

############## ZMSC regression analysis for metabolites（2015）##############
  #### data input ####
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library(openxlsx)

## metabolites
setwd("~/ZMSC/raw")
load('Demo_smoke_male_16S_tran.Rdata')
Demo_tran <- Demo
rm(Demo)

load('~/ZMSC/raw/smoke_male_2015.Rdata')
Demo <- merge(Demo_tran, Demo_2015, by='身份证号码', all.x=T, all.y=F)

BloodID <- read_excel('E:/science/paper/Metabolites/raw/ALL_ID_LH.xlsx') %>% as.data.frame()
metab_info <- read_excel('E:/science/paper/Metabolites/raw/LH_phenoID_list.xlsx') %>% as.data.frame()
metabolite <- read_excel('E:/science/paper/Metabolites/raw/LH_metabdata_BLOODID.xlsx') %>% as.data.frame()
str(metabolite)
rownames(metabolite) <- metabolite$Index
metabolite <- metabolite[,-1]
metabolite <- t(metabolite) %>% as.data.frame()
metab_list <- colnames(metabolite)

##log transform $ z transform
metabolite_log <- apply(metabolite, 2, function(x) log(x)) %>% as.data.frame()
metabolite_log <- apply(metabolite_log, 2, function(x) ifelse(abs((x - mean(x)) / sd(x)) > 3, NA, x))  %>% as.data.frame()
metabolite_standardized <- apply(metabolite_log, 2, function(x) (x - mean(x,na.rm = TRUE)) / sd(x,na.rm = TRUE))  %>% as.data.frame()

Demo <- merge(Demo, BloodID[,c(1,2)], by.x='编号',by.y='QID', all.x=T, all.y=F)
Demo <- merge(Demo,metabolite_standardized,by='BLOODID', all.x=T,all.y=F)

taxa <- c('g__Actinomyces','g__Haemophilus','g__Turicibacter','g__Lachnospira')

  #### smoke_2015 > metab_2015 regression ####
tmp <- list()
for (i in 1:length(metab_list)){
  fit <- glm(Demo[,metab_list[i]]~smoking_status_2015+age_2015+BMI_2015, data=Demo)
  results <- as.data.frame(summary(fit)$coefficients)
  results$beta <- results$Estimate
  results$Lower <- results$Estimate - 1.96*results$`Std. Error`      
  results$Upper <- results$Estimate + 1.96*results$`Std. Error`
  #results$OR <- exp(results$Estimate)                                     
  #results$lower <- exp(results$Estimate - 1.96*results$`Std. Error`)      
  #results$upper <- exp(results$Estimate + 1.96*results$`Std. Error`)
  tmp[[i]] <- results[2:3,]
}
results_glm_smoke_metab <- do.call(rbind, tmp)
results_glm_smoke_metab$smoking_status <- rep(c('Current smoker', 'Former smoker'), length=2)
results_glm_smoke_metab$metabolite[which(results_glm_smoke_metab$smoking_status=='Former smoker')] <- metab_list
results_glm_smoke_metab$metabolite[which(results_glm_smoke_metab$smoking_status=='Current smoker')] <- metab_list

results_glm_smoke_metab <- merge(results_glm_smoke_metab, metab_info, by.x='metabolite',by.y='Index',all.x=T,all.y=F)
results_glm_smoke_metab <- results_glm_smoke_metab[,c(9,14,10,2:3,5,11,16,18,20:31)]
colnames(results_glm_smoke_metab)[1:6] <- c('Smoking','Metabolite','Info.metablite','Estimate','SE','P')

write.xlsx(results_glm_smoke_metab, file='~/ZMSC/results/glm_smoke_metab.xlsx')

  #### smoke_2015 > metab_2015 plot ####
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(RColorBrewer)
library(ggrepel)

df <- read_xlsx('~/ZMSC/results/glm_smoke_metab.xlsx')
df <- df[df$Smoking=='Current smoker',]

class <- unique(df$`Class I`)
color <- c(brewer.pal(12, "Set3")[c(1,3:12)], brewer.pal(12, "Paired")[1:12])
df$pos <- NA
df$color <- NA
for (i in 1:length(class)){
  if (i==1){
    df$pos[which(df$`Class I`==class[i])] <- 1:length(df$`Class I`[which(df$`Class I`==class[i])]) 
    df$color[which(df$`Class I`==class[i])] <- color[i]
  }else if (i>1){
    df$pos[which(df$`Class I`==class[i])] <- (max(df$pos, na.rm=T)+1):(max(df$pos, na.rm=T)+length(df$`Class I`[which(df$`Class I`==class[i])]))
    df$color[which(df$`Class I`==class[i])] <- color[i]
  }
}
df$label <- ifelse(df$Padj_BH<0.05, df$Metabolite, NA)
df$shape <- ifelse(df$Estimate>0, '+', '-')
df$size <-ifelse(is.na(df$label)==TRUE,'-','+')
pos_break <- df %>% group_by(`Class I`) %>%
       summarise(pos_break = (max(pos)+min(pos))/2)

pdf(paste0("~/Results/V1/smoke_metab_zhoushan.pdf"), width =10.5, height = 5)
ggplot(df, aes(x=pos, y=-log(P), label = label)) +
  scale_x_continuous(breaks=pos_break$pos_break, labels = pos_break$`Class I`) +  #x axis
  #scale_y_continuous(trans='log10') +  #y axis log
  #scale_y_break(breaks=c(40,120), space=0.1, scales=0.5) +
  geom_point(data = df, aes(x = pos, y = -log(P),
                            color = color,
                            shape = shape,
                            size = size,
                            fill= color)) + #points
  scale_color_manual(breaks=unique(df$color),
                     values=color) +
  scale_shape_manual(breaks=c('+','-'),
                     values=c(24,25))+ #24,25
  scale_fill_manual(breaks=unique(df$color),
                    values=color)+ #24,25
  scale_size_manual(breaks=c('+','-'),
                    values=c(2,1))+
  labs(x = "Metabolites", y = "-log(P)") +
  geom_hline(yintercept=-log(1.594524e-04),linetype=2, size=0.50, col = "red") +
  theme_bw() +  #rm background
  theme(panel.grid =element_blank(), #rm grids
        legend.position= "none", #set legend
        axis.text.x = element_text(angle=35,hjust=1),
        axis.text = element_text(size=6,color="black"))+ 
  geom_label_repel( #Non overlapped labels
    size=3,
    #color = df$color,
    #nudge_x=5e4, #shift to the right
    segment.alpha = 0.2,  #transparent of segment
    min.segment.length = 1,
    segment.size = 0.2,
    fill=rgb(255, 255, 255, 210, maxColorValue=255),
    fontface='italic'
  )
dev.off()

  #### metab_2015 > taxa_2020 regression ####
tmp_j <- list()
for (j in 1:length(taxa)){
  tmp_i <- list()
  for (i in 1:length(metab_list)){
    fit <- glm(log(Demo[,taxa[j]])~Demo[,metab_list[i]]+age_2015+BMI_2015, data=Demo)
    results <- as.data.frame(summary(fit)$coefficients)
    results$beta <- results$Estimate
    results$lower <- results$Estimate - 1.96*results$`Std. Error`      
    results$upper <- results$Estimate + 1.96*results$`Std. Error`
    #results$OR <- exp(results$Estimate)                                     
    #results$lower <- exp(results$Estimate - 1.96*results$`Std. Error`)      
    #results$upper <- exp(results$Estimate + 1.96*results$`Std. Error`)
    tmp_i[[i]] <- results[2,]
  }
  results_glm_metab_taxa <- do.call(rbind, tmp_i)
  results_glm_metab_taxa$taxa <- taxa[j]
  results_glm_metab_taxa$Index <- metab_list
  results_glm_metab_taxa$Info.metabolite <- paste0(sprintf("%0.2e",apply(Demo[,c(metab_list)],2,function(x){median(x,na.rm=T)})),
                                         "±",
                                         sprintf("%0.2e",apply(Demo[,c(metab_list)],2,function(x){sd(x,na.rm=T)})))
  results_glm_metab_taxa$Padj_BH <- p.adjust(results_glm_metab_taxa$`Pr(>|t|)`,method='BH',n=length(results_glm_metab_taxa$`Pr(>|t|)`))
  results_glm_metab_taxa <- merge(results_glm_metab_taxa, metab_info, by='Index', all.x=T, all.y=F)
  tmp_j[[j]] <- results_glm_metab_taxa
}
results_glm_metab_taxa <- do.call(rbind, tmp_j)
colnames(results_glm_metab_taxa)[1:6] <- c('Metabolite','Info.metablite','Taxa','Estimate','SE','P')
write.xlsx(results_glm_metab_taxa, file='~/ZMSC/results/glm_metab_taxa.xlsx')

metab_sign <- c('MEDN0074','MEDN1970','MEDP0364','MEDP1243','MW0008522')
tmp_j <- list()
for (j in 1:length(taxa)){
  tmp_i <- list()
  for (i in 1:length(metab_sign)){
    fit <- glm(log(Demo[,taxa[j]])~Demo[,metab_sign[i]]+age_2015+BMI_2015, data=Demo)
    results <- as.data.frame(summary(fit)$coefficients)
    results$beta <- results$Estimate
    results$lower <- results$Estimate - 1.96*results$`Std. Error`      
    results$upper <- results$Estimate + 1.96*results$`Std. Error`
    tmp_i[[i]] <- results[2,]
  }
  results_glm_metab_taxa <- do.call(rbind, tmp_i)
  results_glm_metab_taxa$taxa <- taxa[j]
  results_glm_metab_taxa$Index <- metab_sign
  results_glm_metab_taxa$Info.metabolite <- paste0(sprintf("%0.2e",apply(Demo[,c(metab_sign)],2,function(x){median(x,na.rm=T)})),
                                                   "±",
                                                   sprintf("%0.2e",apply(Demo[,c(metab_sign)],2,function(x){sd(x,na.rm=T)})))
  results_glm_metab_taxa$Padj_BH <- p.adjust(results_glm_metab_taxa$`Pr(>|t|)`,method='BH',n=length(results_glm_metab_taxa$`Pr(>|t|)`))
  results_glm_metab_taxa <- merge(results_glm_metab_taxa, metab_info, by='Index', all.x=T, all.y=F)
  tmp_j[[j]] <- results_glm_metab_taxa
}
results_glm_metab_taxa_sign <- do.call(rbind, tmp_j)
colnames(results_glm_metab_taxa_sign)[1:6] <- c('Metabolite','Info.metablite','Taxa','Estimate','SE','P')
write.xlsx(results_glm_metab_taxa_sign, file='~/ZMSC/results/glm_metab_taxa_sign.xlsx')

  #### metab_2015 > taxa_2020 plot ####
library(ggplot2)

df <- read_xlsx('~/ZMSC/results/glm_metab_taxa_sign.xlsx')
df <- df[,1:7]

df$sig <- ifelse(df$Padj_BH<0.05, '#',
                 ifelse(df$Padj_BH>=0.05 & df$P<0.01, '**',
                        ifelse(df$P>0.01 & df$P<0.05, '*', '')))
df$Metabolite[which(df$Metabolite=='γ-Glu-Cys')] <- '	gamma-Glu-Cys'

pdf(paste0("~/Meta/metab_taxa_zhoushan.pdf"), width =4.5, height = 4.5)
ggplot(df,aes(x=Taxa, y=Metabolite)) +
  geom_tile(color="#C2C2C2",fill="white",size=0.3)+
  geom_point(aes(fill=Estimate), pch=22, color="white",size=7.5)+
  scale_fill_gradient2(name = 'effect',
                       limit = c(min(df$Estimate),max(df$Estimate)),
                       low='#66C7B4',mid='#F2F2F2',high='#F4511E', 
                       midpoint = 0)+
  geom_text(aes(label=sig), color="black", size=5) +
  labs(x=NULL,y=NULL,title=NULL)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10,color="black"),
        axis.text.x = element_text(angle=45,hjust=1))+
  guides(fill=guide_colorbar(barheight = 10))
dev.off()


############## five genus for three group in each study ##############
  #### input data ####
library(data.table)
library(readxl)
library(ggvenn)
library(ggplot2)
library(openxlsx)

load('~/ZMSC/raw/Demo_smoke_male_16S_tran.Rdata')  # (the same for GNHS,GGMP)
ZMSC <- Demo
colnames(ZMSC)[24:156] <- sub("^(g__)(.*)$", "\\2", colnames(ZMSC)[24:156])
colnames(ZMSC)[24:156] <- paste0('g_', colnames(ZMSC)[24:156])
rm(Demo)

  #### box plot ####
library(ggplot2)
library(ggsignif)
library(customLayout)

taxa <- c('g_Actinomyces','g_Turicibacter','g_Lachnospira','g_Haemophilus','g_Atopobium')
dataset <- c('ZMSC','GNHS','GGMP')
tag=data.frame(j1=letters[1:5], j2=letters[6:10], j3=letters[11:15]) %>% t()

j=1
i=2

#raw
for (j in 1:length(dataset)){
  for (i in 1:length(taxa)){
    df <- eval(parse(text=dataset[j]))
    if (is.numeric(df[[taxa[i]]])){
      plot <- ggplot(data=df, aes_string(x='smoking_status',y=df[[taxa[i]]], fill='smoking_status')) +
        geom_violin(width=0.55) +
        geom_boxplot(width=0.4, color="#3f3f3f", alpha=0.3) +
        scale_y_continuous(trans='log10') + 
        scale_fill_manual(breaks=c('Non-smoker', 'Current smoker', 'Former smoker'),
                          values=c('#66C7B4','#FFC107','#7C9FB0'))+
        labs(x = "", y = paste0(taxa[i])) +
        geom_signif(comparisons = list(c('Non-smoker','Current smoker'),
                                       c('Non-smoker','Former smoker'),
                                       c('Current smoker','Former smoker')),
                    map_signif_level=F, textsize=4, na.rm=T,
                    test = "wilcox.test",
                    y_position = c(log(max(df[[taxa[i]]])*7),
                                   log(max(df[[taxa[i]]])*16),
                                   log(max(df[[taxa[i]]])*10.5)), extend_line=0) +
        theme_bw()
      eval(parse(text = paste0('p',tag[j,i], ' = plot')))
    }else{
      plot <- ggplot() + theme_bw() + labs(x = "", y = paste0(taxa[i],' (residual ajdusted for age and BMI)'))
      eval(parse(text = paste0('p',tag[j,i], ' = plot')))
    }
  }
}

  #### patch plot ####
library(patchwork)

pdf(paste0('~/Results/V1/taxa_raw_among_three_groups.pdf'), width = 18, height = 13)
pa + pb + pc + pd + pe +
pf + pg + ph + pi + pj +
pk + pl + pm + pn + po + 
  plot_layout(ncol=5, widths =c(1,1,1,1,1), guides = 'collect') +
  plot_annotation() & #tag_levels = 'a'
  theme(legend.position='bottom',
        legend.key.size = unit(1.2, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
dev.off()

############## sankey diagram for Actinomyces ##############
  #### plot ####
library(networkD3)
library(tidyr)

load(paste0("~/ZMSC/raw/physeq_smoke_metagenome_without0.01.Rdata"))
otu_table <- as.data.frame(physeq@otu_table)
tax_table <- as.data.frame(physeq@tax_table)
Actinomyces <- otu_table[grep('Actinomyces', rownames(otu_table)),]
Actinomyces <- as.data.frame(rowSums(Actinomyces)) 
df <- data.frame(Genus=rep('g__Actinomyces ↑',length(Actinomyces$`rowSums(Actinomyces)`)),
                 Species=rownames(Actinomyces),
                 Value=Actinomyces$`rowSums(Actinomyces)`)
df$Species <- c("s__Actinomyces_cardiffensis","s__Actinomyces_georgiae",
                "**s__Actinomyces_graevenitzii ↑",
                "s__Actinomyces_hongkongensis",    
                "**s__Actinomyces_johnsonii ↓","**s__Actinomyces_massiliensis ↓","**s__Actinomyces_naeslundii ↓",
                "s__Actinomyces_odontolyticus",    
                "s__Actinomyces_oris",
                "**s__Actinomyces_sp_HMSC035G02 ↑","**s__Actinomyces_sp_HPA0247 ↑",
                "s__Actinomyces_sp_ICM47",         
                "**s__Actinomyces_sp_S6_Spd3 ↑",
                "s__Actinomyces_sp_oral_taxon_170","s__Actinomyces_sp_oral_taxon_180",
                "**s__Actinomyces_sp_oral_taxon_181 ↑",
                "s__Actinomyces_sp_oral_taxon_414","s__Actinomyces_sp_oral_taxon_448","s__Actinomyces_sp_oral_taxon_897","s__Actinomyces_turicensis")

                 
node <- data.frame(name=c(as.character(df$Genus), as.character(df$Species)) %>% unique())
df$IDgenus=match(df$Genus, node$name)-1
df$IDspecies=match(df$Species, node$name)-1

sankeyNetwork(Links = df, Nodes = node,
              Source = "IDgenus", Target = "IDspecies",
              Value = "Value", NodeID = "name",
              sinksRight=FALSE,
              nodeWidth=45, fontSize=15, nodePadding=12)
