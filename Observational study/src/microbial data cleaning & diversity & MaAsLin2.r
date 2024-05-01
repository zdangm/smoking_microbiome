
library(readxl)
library(openxlsx)
library(data.table)
library(phyloseq)
library(qiime2R)
library(dplyr)
library(tidyr)
library(GUniFrac)
library(vegan)
library(ggplot2)
library(UpSetR)
library(MicrobiotaProcess)

######## phyloseq data for ZMSC/GNHS/GGMP/SRRSHS######## (take ZMSC 16S for example)
Cohort <- 'ZMSC'
  ####input data####
ASV_all <- fread('~/ZMSC/16S/taxonomy_g_level.tsv') %>% as.data.frame()
rownames(ASV_all) <- ASV_all$`#OTU ID`
ASV_all <- ASV_all[,-1]
ASV_all = apply(ASV_all, 2, function(x){return(x/sum(x))}) %>% as.data.frame()

load('~/ZMSC/raw/smoke_male_2020.Rdata')
group_all <- smoke_male
rm(smoke_male)
rownames(group_all) <- group_all$Sample_ID

taxnomy_all <- rownames(ASV_all) %>% as.data.frame()
colnames(taxnomy_all) <- 'Taxon'
rownames(taxnomy_all) <- taxnomy_all$Taxon
taxnomy_all <- separate(taxnomy_all, Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")
taxnomy_all <- subset(taxnomy_all,!is.na(Genus)) #remove those microbes that are not annotated to genus level
taxnomy_all <- subset(taxnomy_all,Genus!='__') 
taxnomy_all <- subset(taxnomy_all, taxnomy_all$Family!='f__Mitochondria') #remove Mitochondria
taxnomy_all <- subset(taxnomy_all, taxnomy_all$Family!='f__Chloroplast') #remove Chloroplast
taxnomy_all <- subset(taxnomy_all, taxnomy_all$Genus!='g__uncultured') #remove genera with uncultured annotations

  ####data cleaning####
ASV_input <- ASV_all[rowMeans(ASV_all)>=0.0001,]  #remove ASVs with average relative abundance <= 0.01%
ASV_input <- ASV_input[rowMeans(ASV_input>0) >= 0.1,] #remove ASVs with prevalence <= 10%
ASV_input <- ASV_input[,which(colnames(ASV_input) %in% group_all$Sample_ID)]
ASV_input <- ASV_input[which(rownames(ASV_input) %in% rownames(taxnomy_all)),]

taxnomy_input <- taxnomy_all[which(rownames(taxnomy_all) %in% rownames(ASV_input)),]
identical(rownames(ASV_input), rownames(taxnomy_input)) #Check for alignment
rownames(ASV_input) <- taxnomy_input$Genus
rownames(taxnomy_input) <- taxnomy_input$Genus

group_input <- group_all[which(rownames(group_all) %in% colnames(ASV_input)),]
group_input <- group_input[order(group_input$Sample_ID,decreasing = F),]
identical(colnames(ASV_input), rownames(group_input)) #Check for alignment

  ####merge phyloseq####
taxnomy_input <- as.matrix(taxnomy_input)
ASV = otu_table(ASV_input, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(taxnomy_input)
Group = sample_data(group_input[order(group_input$Sample_ID),])

identical(rownames(Group), colnames(ASV_input)) #Check for alignment
physeq <- phyloseq(ASV, TAX, Group)

save(physeq, file=paste0("~/",Cohort,"/raw/physeq_smoke_16S.Rdata"))


######## Analysis for ZMSC/GNHS/GGMP/SRRSHS########
Cohort <- c('ZMSC','GNHS','GGMP','SRRSHS')[1]
microbiota <- c('16S','metagenome_without0.01','metagenome','pathway')[1]

load(paste0("~/",Cohort,"/raw/physeq_smoke_",microbiota,".Rdata"))

  #### Demo info ####
library(tableone)
tmp <- cbind(physeq@sam_data$Sample_ID, physeq@sam_data[,c('age','BMI','drinking_status','smoking_status')]) #,'cigarettes_per_day'
allVars <- c('age','BMI','drinking_status','cigarettes_per_day')
fvars <- c('drinking_status')
Svytab <- CreateTableOne(vars=allVars, factorVars=fvars, strata="smoking_status", data=tmp, addOverall=T)
Svytab

  #### Taxa info ####
taxa_info <- data.frame(Taxa = rownames(physeq@otu_table),
                        Prevalence = rowMeans(physeq@otu_table>0))
tmp <- summary(t(physeq@otu_table)) %>% as.data.frame()
tmp <- separate(tmp, Freq, into = c("Q", "Abundance"), sep = ":")
tmp$Abundance <- as.numeric(tmp$Abundance)
taxa_info$Abundance.25th <- sprintf("%0.2e",tmp$Abundance[which(tmp$Q=='1st Qu.')])
taxa_info$Abundance.50th <- sprintf("%0.2e",tmp$Abundance[which(tmp$Q=='Median ')])
taxa_info$Abundance.75th <- sprintf("%0.2e",tmp$Abundance[which(tmp$Q=='3rd Qu.')])
taxa_info$Abundance.Mean <- sprintf("%0.2e",tmp$Abundance[which(tmp$Q=='Mean   ')])
taxa_info$Abundance.SD <- sprintf("%0.2e",apply(t(physeq@otu_table),2,sd))

if (microbiota=='16S'){
  taxa_info <- merge(taxa_info, physeq@tax_table, by.x='Taxa', by.y='Genus', all.x=T, all.y=F)
  taxa_info$Taxa <- gsub("\\.", "_", taxa_info$Taxa)
}else if (microbiota=='metagenome' | microbiota=='metagenome_without0.01'){
  taxa_info <- merge(taxa_info, physeq@tax_table, by.x='Taxa', by.y='Species', all.x=T, all.y=F)
  taxa_info$Taxa <- gsub("\\.", "_", taxa_info$Taxa)
  taxa_info$Genus <- gsub("\\.", "_", taxa_info$Genus)
}

if (Cohort %in% c('ZMSC','GNHS','GGMP') & microbiota=='16S'){
  taxa_info <- taxa_info[,c(9:12,1:7)] 
  colnames(taxa_info)[colnames(taxa_info) == "Taxa"] <- "Genus"
}else if (Cohort=='SRRSHS' & microbiota=='16S'){
  taxa_info <- taxa_info[,c(10:13,1:7)]
  colnames(taxa_info)[colnames(taxa_info) == "Taxa"] <- "Genus"
}else if (Cohort %in% c('ZMSC','GNHS') & microbiota=='metagenome'){
  taxa_info <- taxa_info[,c(9:13,1:7)] 
  colnames(taxa_info)[colnames(taxa_info) == "Taxa"] <- "Species"
}else if (Cohort %in% c('ZMSC','GNHS') & microbiota=='metagenome_without0.01'){
  taxa_info <- taxa_info[,c(9:13,1:7)] 
  colnames(taxa_info)[colnames(taxa_info) == "Taxa"] <- "Species"
}

write.xlsx(taxa_info, paste0("~/",Cohort,"/results/physeq_smoke_",microbiota,"_taxa_info.xlsx"))

  #### Venn plot ####
library(UpSetR)
upsetda <- get_upset(physeq, factorNames="smoking_status")
upsetda

if (microbiota=='16S'){
  label='No.genus'
}else if (microbiota=='metagenome'){
  label='No.species'
}else if (microbiota=='pathway'){
  label='No.pathways'
}

pdf(paste0("~/",Cohort,"/results/Genus_",microbiota,"_smoke_male_group.pdf"), width = 6, height = 5) #pathway 8 5 #16S 6 5
upset(upsetda, 
      matrix.color = '#777777',main.bar.color = '#777777',
      sets.x.label = label,
      point.size = 3, line.size = 2,
      text.scale = 2.2, set_size.angles = 0,
      sets=c("Non-smoker","Current smoker", "Former smoker"), 
      sets.bar.color = c('#66C7B4','#FFC107','#7C9FB0'),
      order.by = "freq", empty.intersections = "on")
dev.off()

  #### α diversity ####
library(vegan)
library(ggplot2)
library(ggsignif)
library(customLayout)

asv_data <- t(physeq@otu_table) %>% as.data.frame()
identical(physeq@sam_data$Sample_ID, rownames(asv_data))
Demo <- physeq@sam_data
Demo$richness <- specnumber(asv_data)
Demo$shannon <- diversity(asv_data, index = "shannon")
Demo$simpson <- diversity(asv_data, index = "simpson")
Demo$Pielou <- Demo$shannon/log(Demo$richness,exp(1))
Demo <- cbind(Demo, asv_data)

Pielou <- ggplot(data=Demo, aes(x=smoking_status,y=Pielou, fill=smoking_status)) +  #(take Pielou for example)
  geom_violin(width=0.5) +
  geom_boxplot(width=0.35, color="#3f3f3f", alpha=0.3) +
  geom_jitter(color="grey", size=0.6, alpha=0.9) +
  scale_fill_manual(breaks=c('Non-smoker', 'Current smoker', 'Former smoker'),
                    values=c('#66C7B4','#FFC107','#7C9FB0'))+
  labs(x = "", y = "Pielou") +
  geom_signif(comparisons = list(c('Non-smoker','Current smoker'),c('Non-smoker','Former smoker'),c('Current smoker','Former smoker')),
              map_signif_level=F, textsize=4, na.rm=T,
              test = "wilcox.test",
              y_position = c(max(Demo$Pielou)*1.01, max(Demo$Pielou)*1.15, max(Demo$Pielou)*1.06), extend_line=0) +
  theme_bw()

  #### β diversity ####
library(GUniFrac)
library(vegan)
set.seed(2023)
adonis_result <- adonis2(Demo[,c(rownames(physeq@otu_table))]~Demo$smoking_status+ 
                           Demo$age+Demo$BMI, 
                         data=Demo, distance = 'bray', permutations = 999, na.action=na.omit)
adonis_result
 #PCoA作图
bray_dis <- vegdist(Demo[,c(rownames(physeq@otu_table))], method = 'bray')    
pcoa <- cmdscale(bray_dis, k = 2, eig = TRUE)

pcoa_exp <- pcoa$eig/sum(pcoa$eig)
site <- scores(pcoa)

site <- data.frame(pcoa$point)[1:2]
site$Sample_ID <- rownames(site)
identical(site$Sample_ID, Demo$Sample_ID)
site$group <- Demo$smoking_status

#pdf(paste0("~/",Cohort,"/results/beta_",microbiota,"_smoke_male.pdf"), width = 6, height = 5)
PCoA <- ggplot(data = site, aes(X1, X2)) +
  geom_point(aes(color = group)) +
  stat_ellipse(aes(fill = group, color = group), linetype = 2, size= 1.3, geom = 'polygon', 
               level = 0.99, alpha = 0.02, show.legend = FALSE) +  
  scale_color_manual(breaks=c('Non-smoker', 'Current smoker', 'Former smoker'),
                     values=c('#66C7B4','#FFC107','#7C9FB0')) +
  scale_fill_manual(breaks=c('Non-smoker', 'Current smoker', 'Former smoker'),
                    values=c('#66C7B4','#FFC107','#7C9FB0')) +
  theme_bw() +
  #theme(legend.position= "none") +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.1) +
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.1) +
  labs(x = paste('PC1(', round(100*pcoa_exp[1], 2), '%)'), 
       y = paste('PC2(', round(100*pcoa_exp[2], 2), '%)'),
       title = paste0('R squared: ',round(adonis_result$R2[1],4),'; P-value: ',round(adonis_result$`Pr(>F)`[1],3),' (PERMANOVA)'))
#dev.off()

save(Pielou, file=paste0('~/',Cohort,'/results/alpha_Pielou_',microbiota,'.Rdata'))
save(PCoA, file=paste0('~/',Cohort,'/results/beta_PCoA_',microbiota,'.Rdata'))

  #### Rdata save ####
save(Demo, file=paste0('~/',Cohort,'/raw/Demo_smoke_male_',microbiota,'_raw.Rdata'))
for (i in which(match(colnames(Demo), rownames(physeq@otu_table))==1) : ncol(Demo)){
  Demo[[i]] <- ifelse(Demo[[i]] == 0, min(Demo[[i]][which(Demo[[i]]!= 0)])/10, Demo[[i]])
}
save(Demo, file=paste0('~/',Cohort,'/raw/Demo_smoke_male_',microbiota,'_tran.Rdata'))

  #### MaAsLin2 ####
library(Maaslin2)
MaAsLin = Maaslin2( # (take mian analysis for example)
  input_data = Demo[,c('Sample_ID', rownames(physeq@otu_table))], 
  input_metadata = Demo[,c("Sample_ID",'smoking_status','age', 'BMI')], 
  output = paste0('~/',Cohort,'/results/MaAsLin2_',microbiota,'_smoke_male'), 
  min_abundance = 0.0,
  min_prevalence = 0.1,
  min_variance = 0.0,
  normalization = "NONE", #TSS, CLR, CSS, NONE, TMM
  transform = "LOG",         #LOG, LOGIT, AST, NONE
  analysis_method = "LM",   #LM, CPLM, NEGBIN, ZINB
  max_significance = 0.25,
  fixed_effects = c('smoking_status','age', 'BMI'),
  correction = "BH",
  standardize = TRUE,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 20,
  reference = c('smoking_status,Non-smoker'))

