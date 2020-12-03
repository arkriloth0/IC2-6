#install packages
install.packages("BiocManager")
BiocManager::install(version='devel')
BiocManager::install("depmap")
BiocManager::install("remotes")
BiocManager::install("uclouvain-cbio/depmap")

#load libraries
library(tidyverse)
library(data.table)
library(ExperimentHub)
library(depmap)

#load depmap proteomic data
prot <- depmap_proteomic()

#load sample info
sample_info <- fread("sample_info.csv")

#load list of breast cancer cell lines
BreastCA_lines <- fread("Breast cell lines.csv", header = FALSE, col.names = "DepMap_ID")


#Get copy number of ZNF703
ZNF_CN <- fread("CCLE_gene_cn.csv", select = c("V1","ZNF703 (80139)"))
colnames(ZNF_CN) <- c("DepMap_ID","ZNF703_copy_number")

#write.csv(ZNF_CN,file="ZNF_CN.csv")

#Get absolute copy number of ZNF703
ZNF_absCN <- fread("Copy_Number_(Absolute).csv", select = c("V1","ZNF703"),col.names = c("DepMap_ID","abs_CN"))
colnames(ZNF_absCN) <- c("DepMap_ID","ZNF703_abs_CN")

#write.csv(ZNF_CN,file="ZNF_absCN.csv")

#Get mRNA expression of ZNF703
ZNF_expression <- fread("CCLE_expression.csv", select = c("V1","ZNF703 (80139)"))
colnames(ZNF_expression) <- c("DepMap_ID","ZNF703_mRNA_expression")

#write.csv(ZNF_expression,file="ZNF_expression.csv")

#get protein expression of ZNF703
ZNF_prot <- prot %>% filter(entrez_id == "80139")%>%
  select(depmap_id,protein_expression)
colnames(ZNF_prot)<- c("DepMap_ID","ZNF703_protein_expression")

#write.csv(ZNF_prot,file="ZNF_prot.csv")

#get CRISPR dependency data of ZNF703
ZNF_dep <- fread("Achilles_gene_effect.csv", select = c("DepMap_ID","ZNF703 (80139)"))
colnames(ZNF_dep) <- c("DepMap_ID","ZNF703_CERES")

#write.csv(ZNF_dep,file="ZNF_dep.csv")

#get CRISPR dependency p value of ZNF703
ZNF_pdep <- fread("Achilles_gene_dependency.csv", select = c("DepMap_ID","ZNF703 (80139)"))
colnames(ZNF_pdep) <- c("DepMap_ID","ZNF703_p_dep")

#write.csv(ZNF_dep,file="ZNF_pdep.csv")

#Get expression of ZNF503
ZNF503_expression <- fread("CCLE_expression.csv", select = c("V1","ZNF503 (84858)"))
colnames(ZNF503_expression) <- c("DepMap_ID","ZNF503_expression")

#write.csv(ZNF503_expression,file="ZNF503_expression.csv")

#joining all ZNF703 tables
ZNF_all <- ZNF_CN %>%
  full_join(ZNF_absCN) %>%
  full_join(ZNF_expression) %>%
  full_join(ZNF_prot) %>%
  #full_join(ZNF503_expression) %>%
  full_join(ZNF_dep) %>%
  full_join(ZNF_pdep)

#write.csv(ZNF_all, file="ZNF_all.csv")

#Adding sample info/annotation
#ZNF_all <- full_join(ZNF_all,sample_info) %>%
#  .[, c(1,8:31,2:7)]

#write.csv(ZNF_all, file="ZNF_all.csv")

#selecting only breast cancer cell lines
#ZNF_breast <- semi_join(ZNF_all,BreastCA_lines)
#write.csv(ZNF_breast, file="ZNF_breast.csv")


##PLPP5

#Get copy number of PLPP5
plpp5_CN <- fread("CCLE_gene_cn.csv", select = c("V1","PLPP5 (84513)"))
colnames(plpp5_CN) <- c("DepMap_ID","PLPP5_copy_number")

#Get absolute copy number of PLPP5
plpp5_absCN <- fread("Copy_Number_(Absolute).csv", select = c("V1","PLPP5"),col.names = c("DepMap_ID","abs_CN"))
colnames(plpp5_absCN) <- c("DepMap_ID","PLPP5_abs_CN")

#Get mRNA expression of PLPP5
plpp5_expression <- fread("CCLE_expression.csv", select = c("V1","PLPP5 (84513)"))
colnames(plpp5_expression) <- c("DepMap_ID","PLPP5_mRNA_expression")

#get protein expression of PLPP5
plpp5_prot <- prot %>% filter(entrez_id == "84513")%>%
  select(depmap_id,protein_expression)
colnames(plpp5_prot) <- c("DepMap_ID","PLPP5_protein_expression")

#get CRISPR dependency data of PLPP5
plpp5_dep <- fread("Achilles_gene_effect.csv", select = c("DepMap_ID","PLPP5 (84513)"))
colnames(plpp5_dep) <- c("DepMap_ID","PLPP5_CERES")

#get CRISPR dependency p value of PLPP5
plpp5_pdep <- fread("Achilles_gene_dependency.csv", select = c("DepMap_ID","PLPP5 (84513)"))
colnames(plpp5_pdep) <- c("DepMap_ID","PLPP5_p_dep")


#joining all PLPP5 tables
plpp5_all <- plpp5_CN %>%
  full_join(plpp5_absCN) %>%
  full_join(plpp5_expression) %>%
  full_join(plpp5_prot) %>%
  full_join(plpp5_dep) %>%
  full_join(plpp5_pdep)

##FGFR1

#Get copy number of FGFR1
fgfr1_CN <- fread("CCLE_gene_cn.csv", select = c("V1","FGFR1 (2260)"))
colnames(fgfr1_CN) <- c("DepMap_ID","FGFR1_copy_number")

#Get absolute copy number of FGFR1
fgfr1_absCN <- fread("Copy_Number_(Absolute).csv", select = c("V1","FGFR1"),col.names = c("DepMap_ID","abs_CN"))
colnames(fgfr1_absCN) <- c("DepMap_ID","FGFR1_abs_CN")

#Get mRNA expression of FGFR1
fgfr1_expression <- fread("CCLE_expression.csv", select = c("V1","FGFR1 (2260)"))
colnames(fgfr1_expression) <- c("DepMap_ID","FGFR1_mRNA_expression")

#get protein expression of FGFR1
fgfr1_prot <- prot %>% filter(entrez_id == "2260")%>%
  select(depmap_id,protein_expression)
colnames(fgfr1_prot) <- c("DepMap_ID","FGFR1_protein_expression")

#get CRISPR dependency data of FGFR1
fgfr1_dep <- fread("Achilles_gene_effect.csv", select = c("DepMap_ID","FGFR1 (2260)"))
colnames(fgfr1_dep) <- c("DepMap_ID","FGFR1_CERES")

#get CRISPR dependency p value of FGFR1
fgfr1_pdep <- fread("Achilles_gene_dependency.csv", select = c("DepMap_ID","FGFR1 (2260)"))
colnames(fgfr1_pdep) <- c("DepMap_ID","FGFR1_p_dep")

#joining all FGFR1 tables
fgfr1_all <- fgfr1_CN %>%
  full_join(fgfr1_absCN) %>%
  full_join(fgfr1_expression) %>%
  full_join(fgfr1_prot) %>%
  full_join(fgfr1_dep) %>%
  full_join(fgfr1_pdep)

#joining all 8p12 tables
chr8p12_all <- ZNF_all %>%
  full_join(plpp5_all) %>%
  full_join(fgfr1_all)

#Adding sample info/annotation
chr8p12_all <- full_join(sample_info,chr8p12_all)

write.csv(chr8p12_all, file="chr8p12_all.csv")

#selecting only breast cancer cell lines
chr8p12_breast <- semi_join(chr8p12_all,BreastCA_lines)
write.csv(chr8p12_breast, file="chr8p12_breast.csv")

##other stuff
dep <- fread("Achilles_gene_effect.csv") 
dep_SUM52PE <- filter(dep, dep$DepMap_ID == "ACH-001396") %>%
  pivot_longer(cols = 2:18120, names_to = "gene", values_to = "CERES")
mean(dep_SUM52PE$CERES)

## create ExperimentHub query object and display web service
eh <- ExperimentHub()
query(eh, "depmap")
d <- display(eh)


#           title             
#EH3797 | crispr_20Q3       
#EH3798 | copyNumber_20Q3   
#EH3799 | TPM_20Q3          
#EH3800 | mutationCalls_20Q3
#EH3801 | metadata_20Q3   
