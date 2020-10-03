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

#Get copy number of ZNF703
ZNF_CN <- fread("CCLE_gene_cn.csv", select = c("V1","ZNF703 (80139)"))
colnames(ZNF_CN) <- c("DepMap_ID","copy_number")

write.csv(ZNF_CN,file="ZNF_CN.csv")

#Get absolute copy number of ZNF703
ZNF_absCN <- fread("Copy_Number_(Absolute).csv", select = c("V1","ZNF703"))
colnames(ZNF_absCN) <- c("DepMap_ID","abs_CN")

write.csv(ZNF_CN,file="ZNF_absCN.csv")

#Get mRNA expression of ZNF703
ZNF_expression <- fread("CCLE_expression.csv", select = c("V1","ZNF703 (80139)"))
colnames(ZNF_expression) <- c("DepMap_ID","mRNA_expression")

write.csv(ZNF_expression,file="ZNF_expression.csv")

#get protein expression of ZNF703
prot <- depmap_proteomic()
ZNF_prot <- prot %>% filter(gene_name == "ZNF703")%>%
  select(depmap_id,protein_expression)
colnames(ZNF_prot)[1] <- "DepMap_ID"

write.csv(ZNF_prot,file="ZNF_prot.csv")

#get CRISPR dependency data of ZNF703
ZNF_dep <- fread("Achilles_gene_effect.csv", select = c("DepMap_ID","ZNF703 (80139)"))
colnames(ZNF_dep) <- c("DepMap_ID","CERES")

write.csv(ZNF_dep,file="ZNF_dep.csv")

#get CRISPR dependency p value of ZNF 703
ZNF_pdep <- fread("Achilles_gene_dependency.csv", select = c("DepMap_ID","ZNF703 (80139)"))
colnames(ZNF_pdep) <- c("DepMap_ID","p_dep")

write.csv(ZNF_dep,file="ZNF_pdep.csv")

#Get expression of ZNF503
ZNF503_expression <- fread("CCLE_expression.csv", select = c("V1","ZNF503 (84858)"))
colnames(ZNF503_expression) <- c("DepMap_ID","ZNF503_expression")

write.csv(ZNF503_expression,file="ZNF503_expression.csv")

#joining all tables + annotation
sample_info <- fread("sample_info.csv")
ZNF_all <- ZNF_CN %>%
  full_join(ZNF_absCN) %>%
  full_join(ZNF_expression) %>%
  full_join(ZNF_prot) %>%
  #full_join(ZNF503_expression) %>%
  full_join(ZNF_dep) %>%
  full_join(ZNF_pdep) %>%
  full_join(sample_info)

ZNF_all <- ZNF_all[, c(1,8:31,2:7)]
write.csv(ZNF_all, file="ZNF_all.csv")

BreastCA_lines <- fread("Breast cell lines.csv", header = FALSE, col.names = "DepMap_ID")

ZNF_breast <- semi_join(ZNF_all,BreastCA_lines)
write.csv(ZNF_breast, file="ZNF_breast.csv")

ggplot(ZNF_all, aes(x=abs_CN,y=expression)) + geom_point()


dep <- fread("Achilles_gene_effect.csv") 
dep_SUM52PE <- filter(dep, dep$DepMap_ID == "ACH-001396") %>%
  pivot_longer(cols = 2:18120, names_to = "gene", values_to = "CERES")
mean(dep_SUM52PE$CERES)
                     
## create ExperimentHub query object and display web service
eh <- ExperimentHub()
query(eh, "depmap")
d <- display(eh)

