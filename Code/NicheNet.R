
library(Seurat)
library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)


ligand_target_matrix <- readRDS('ligand_target_matrix_nsga2r_final.rds')
lr_network <- readRDS('lr_network_human_21122021.rds')
weighted_networks <- readRDS('weighted_networks_nsga2r_final.rds')
weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

bcell <- readRDS('B_obj_annotation.rds')
Idents(bcell) <- "Treatment"
bcell <- subset(bcell,idents="Pre_treatment")

Idents(bcell) <- "B_main.type"
list_expressed_genes_sender = get_expressed_genes(c("AtM"), bcell, pct = 0.10)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


#receiver：CD8
tcell <- readRDS("T_obj_annotation_final_Response.rds")
Idents(tcell) <- "Treatment"
tcell <- subset(tcell,idents="Pre_treatment")

Idents(tcell) <- "T_main.type"
expressed_genes_receiver = get_expressed_genes(c('CD8_Tcm','CD8_Tex_GZMK','CD8_Tem','CD8_Tex_CXCL13','CD8_HSP'), tcell, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## geneset_oi
tcell <- readRDS("T_obj_annotation_final_Response.rds")
Idents(tcell) <- "Treatment"
tcell <- subset(tcell,idents="Pre_treatment")

Idents(tcell) <- "T_main.type"
DE_table_receiver <-  FindMarkers(object = tcell,
                                  ident.1 = "CD8_Tex_CXCL13",ident.2 = c('CD8_Tcm','CD8_Tex_GZMK','CD8_Tem','CD8_HSP'),only.pos=T,
                                  min.pct = 0.25) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]


bcell <- readRDS('B_obj_annotation.rds')
Idents(bcell) <- "Treatment"
bcell <- subset(bcell,idents="Pre_treatment")

Idents(bcell) <- "B_main.type"
obj_use <- subset(bcell,idents="AtM")

source("NiCheNet_nicheFun.R")
res = nicheFun(background_expressed_genes,
               expressed_genes_sender,
               expressed_genes_receiver,
               geneset_oi,
               ligand_target_matrix,
               lr_network, 
               weighted_networks,
               weighted_networks_lr,
               obj_use)

pdf("AtM_CD8_Tex_CXCL13.pdf",width = 13,height = 8)
res$combined
dev.off()




