setwd("~/Documents/CASB185")

########################
####INSTALL PACKAGES####
########################
BiocManager::install("RRHO")
install.packages("dplyr")
BiocManager::install("biomaRt")

########################
#####LOAD LIBRARIES#####
########################
library(stats)
library (readr)
library(dplyr)
library(tidyverse)
library(RRHO)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)

########################
######READ IN DATA######
########################
load(file="processed_data.rdata")
load(file="gene_annotation.rdata")
load(file="LongTerm_pariedVoom_results.rdata")
load(file="fgsea_results.rdata")
load(file="hypermatAll.rdata")
load(file="genedat.rdata")
hcs_expression <- readRDS("/Users/faithyuan/Downloads/hcs_expression_data.rds")
brainspan_expression <- readRDS("/Users/faithyuan/Downloads/brainspan_expression_data.rds")
new_asd_sfari_all <- read.csv("SFARI-Gene_genes_01-11-2022release_03-06-2022export.csv")
mouse_asd_sfari_all <- read.csv("SFARI-Gene_animal-genes_01-11-2022release_03-07-2022export.csv")


View(new_asd_sfari_all)
View(datExpr_reg_batch)
View(hcs_expression)
View(brainspan_expression)

########################
#####HUMAN ASD DATA#####
########################

##Updated dataset with 100+ new genes included
View(new_asd_sfari_all)

##Filter dataset to only include syndromic or if there is a gene-score higher 2 or higher
new_asd_sfari<- filter(new_asd_sfari_all, (new_asd_sfari_all$gene.score <= "2" | new_asd_sfari_all$gene.score == "S") | new_asd_sfari_all$syndromic == "1")
View(new_asd_sfari)

##Find intersection of ASD genes and HCS dataset
new_expressed_sfari <- intersect(new_asd_sfari$ensembl.id,rownames(datExpr_reg_batch))

##Use intersected ensembles list to find gene expressions for specific list
new_sfari_expr <- datExpr_reg_batch[new_expressed_sfari,]

##Find intersection of ASD genes and brain span dataset
new_brain_expressed_sfari <-intersect(new_asd_sfari$ensembl.id,rownames(brainspan_expression))
new_brain_sfari_expr <- brainspan_expression[new_brain_expressed_sfari,]

##RRHO on ASD/HCS gene expressions and ASD/brain span gene expressions
new_rrho_sfari <- RRHO(new_sfari_expr,new_brain_sfari_expr, alternative = "two.sided")
image(new_rrho_sfari$hypermat)

##Heatmaps of both gene expressions

pal <- rev(colorRampPalette(brewer.pal("RdBu", n = 11)[c(1,1,2,3,4,6,8,9,10,11,11)])(100))
pheatmap(new_sfari_expr, color = pal, scale = "row", )
pheatmap(new_brain_sfari_expr, color = pal, scale = "row", )


########################
#####MOUSE ASD DATA#####
########################

##Filter mouse data with the same parameters as before
mouse_asd_sfari<- filter(mouse_asd_sfari_all, (mouse_asd_sfari_all$human.gene.score <= "2" | mouse_asd_sfari_all$human.gene.score == "S") | mouse_asd_sfari_all$human.syndromic == "1")
View(mouse_asd_sfari$gene.symbol)

##Read in converted mouse genes to human ensembles text
mouse_ensemble <- read.table("input.txt")
View(mouse_ensemble)

##Change to character datatype
temp4 <- as.character(mouse_ensemble)
View(temp4)

##Compare the mouse gene names with ASD/HCS gene list
temp5 <- datExpr_reg_batch[temp4,]


##Compare the mouse gene names with ASD/HCS gene list
temp3 <- intersect(mouse_asd_sfari$gene.symbol, new_asd_sfari$gene.symbol)
temp6<- getBM(mouse_asd_sfari$gene.name)
temp7<-useMart(biomart="ENSEMBL_MART_ENSEMBL", data= mouse_asd_sfari$gene.name, )
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), values = mouse_asd_sfari$gene.name, mart=mart)
