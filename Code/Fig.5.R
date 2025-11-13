# ------------------------------------------------------------------
# Fig. 5
# ------------------------------------------------------------------

#-----Fig. 5b-----####

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(circlize)

B_obj_integrated <- readRDS("B_obj_integrated.rds")

DefaultAssay(B_obj_integrated) <- "RNA"
scale_data <- GetAssayData(B_obj_integrated, slot = "scale.data")

cell_types <- Idents(B_obj_integrated)

ABC_gene_list_combined <- unique(unlist(ABC_gene_list))
valid_genes <- ABC_gene_list_combined[ABC_gene_list_combined %in% rownames(scale_data)]

avg_expression <- tapply(
  1:ncol(scale_data),
  cell_types,
  function(cells) rowMeans(scale_data[valid_genes, cells, drop = FALSE])
)

avg_expression_matrix <- do.call(cbind, avg_expression)

gene_groups <- unlist(lapply(names(ABC_gene_list), function(group) {
  rep(group, length(ABC_gene_list[[group]]))
}))

names(gene_groups) <- unlist(ABC_gene_list)

valid_gene_groups <- gene_groups[valid_genes]

library(ComplexHeatmap)

row_annotation <- rowAnnotation(
  Group = factor(valid_gene_groups, levels = unique(valid_gene_groups)),
  col = list(Group = c(
    "ABC_marker" = "#98C1B6",
    "Activation" = "#EDA9B8",
    "Chemokine_receptors" = "#CA6651",
    "Inflammatory_cytokine_receptors" = "#EED392",
    "MHCII" = "#98CBCB",
    "Transcription" = "#7DABCA"
  )),
  annotation_legend_param = list(title = "Gene Groups"),
  annotation_width = unit(1, "mm")
)

col_fun <- colorRamp2(
  breaks = c(min(avg_expression_matrix), 0, max(avg_expression_matrix)),
  colors = c("#377594", "white", "#9F1626") # 定义颜色范围
)


pdf("ABC_genelist_expr_heatmap.pdf",width = 5,height = 8)
Heatmap(avg_expression_matrix,
        name = "Expression",
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_split = valid_gene_groups,
        left_annotation = row_annotation,
        heatmap_legend_param = list(title = "Avgrage Scaled Expression"),
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x, y, width = width, height = height, 
                    gp = gpar(col = "black", lwd = 0.5))
        })
dev.off()


ABC_obj <- subset(B_obj_integrated,idents="ABC")
Idents(ABC_obj) <- "Response"

chemokine_receptors <- c("CXCR5")

DefaultAssay(ABC_obj) <- "RNA"
Idents(ABC_obj) <- "Response"

Idents(ABC_obj) <- factor(Idents(ABC_obj),levels = rev(c("Pre_nPCR","Pre_PCR","Post_nPCR","Post_PCR")))

pdf("CXCR5_ABC_expr.pdf",width = 3,height = 2.5)
DotPlot(object = ABC_obj, features = chemokine_receptors) +
  scale_color_gradient2(low = "#377594", mid = "#ffffff", high = "#9F1626")+
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,color = "black"),
        axis.text.y =element_text(color = "black"))+
  xlab('')+ylab('')
dev.off()


#-----Fig. 5c-----####
library(Seurat)
library(sscVis)
library(data.table)
library(grid)
library(cowplot)
library(ggrepel)
library(readr)
library(plyr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
group_res <- readRDS("B_OR.group_OR.rds")

all_OR_res <- group_res
heatdata <- all_OR_res[,c(1,2,5)]
library(reshape2)
heatdata <- as.data.frame(heatdata)

heatdata <- dcast(heatdata, rid ~ cid, value.var = "OR")
rownames(heatdata) <- heatdata[,1]
heatdata <- heatdata[,-1]

library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(0, 1, 5), c("#6DC2E9", "#FFFFFF", "#D64756"))

heatdata <- heatdata[c(4,3,1,2,5),]
heatdata <- heatdata[,c(1,3,2,5,4)]

pdf("OR_heatmap.pdf",width = 5,height = 4)
pheatmap(as.matrix(heatdata),
         color = col_fun,display_numbers = TRUE,
         name='OR',number_color = "black",
         cellwidth = 25,
         cellheight = 25,
         border_color = 'black', 
         # scale = "row", 
         cluster_rows = F,
         cluster_cols = F, 
         legend = TRUE,
         legend_breaks = c(0, 3, 5), 
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 8,
         gaps_col=c(1,1+2,3+2)
         # gaps_row=c(3,3+2)
)
dev.off()


# miloR
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)
library(ggbeeswarm)
library(scales)
library(forcats)
library(data.table)
library(stringr)

sc.epi_milo <- readRDS("sc.epi_milo.rds")

da_results <- annotateNhoods(sc.epi_milo, da_results, coldata_col = "B_main.type")
da_results$B_main.type <- factor(da_results$B_main.type, 
                                 levels = rev(c("Bn","Bmen","ABC","Bgc","Plasma")))

pdf("/plotDAbeeswarm.pdf",width = 3.5,height = 3)
plot_celltype <- plotDAbeeswarm(da_results, group.by = "B_main.type") + 
  scale_color_gradient2(low="#FFB600",                     
                        mid="lightgrey",                
                        high="#77B6BE") +  
  labs(x="", y="Log2 Fold Change") +  
  theme_bw(base_size=10)+
  theme(axis.text.x=element_text(vjust=1,size=10,color = "black"),
        axis.text.y=element_text(vjust=1,size=10,color = "black"))
print(plot_celltype)
dev.off()




#-----Fig. 5i-----####

GO_res <- read.csv("B_pre_nPCR_PCR_GO_res.txt",quote=FALSE,sep="\t",row.names=F)
plot_data <- GO_res[,c(1,3,6,10)]

pos_plot <- plot_data[which(plot_data$Cluster=='Pre_PCR'),]
pos_plot <- pos_plot[which(pos_plot$Description %in% PCR_pathway),]
pos_plot$abs <- -log10(pos_plot$pvalue)
pos_plot <- pos_plot[order(pos_plot$abs,decreasing = F),]
pos_plot$Description <- factor(pos_plot$Description,levels = pos_plot$Description)


neg_plot <- plot_data[which(plot_data$Cluster=='Pre_nPCR'),]
neg_plot <- neg_plot[which(neg_plot$Description %in% nPCR_pathway),]
neg_plot$abs <- log10(neg_plot$pvalue)
neg_plot <- neg_plot[order(neg_plot$abs,decreasing = F),]
neg_plot$Description <- factor(neg_plot$Description,levels = neg_plot$Description)


plot_df <- rbind(neg_plot,pos_plot)

plot_df$Description <- factor(plot_df$Description,levels = plot_df$Description)

mytheme <- theme(
  legend.position = 'none',
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'grey60',size = 1.1),
  axis.text = element_text(size = 12)
)

pos <- plot_df[which(plot_df$abs > 0),]
neg <- plot_df[which(plot_df$abs < 0),]

group_color <- c("Pre_PCR"="#A8541E","Pre_nPCR"="#FC8F62")

pdf('B_pre_nPCR_PCR_GO_res.pdf',width = 6,height = 6)
ggplot(plot_df,
       aes(x =abs, y = Description, fill = Cluster))+
  geom_col()+
  theme_bw()+mytheme+
  geom_text(data = pos,
            aes(x = -0.2, y = Description, label = Description),
            size = 3.5,
            hjust = 1)+
  geom_text(data = neg,
            aes(x = 0.2, y = Description, label = Description),
            size = 3.5,
            hjust = 0)+
  labs(x = 'log(pvalue)', y = ' ', title = 'Enriched GO Biological  Pathway') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14))+
  scale_fill_manual(values = group_color)
dev.off() 


#-----Extended Data Fig. 6a-----####

Bn_genes <- c("TCL1A","IGHD","FCER2","IL4R")
Bmen_genes <- c("CCR7","CD27","AIM2","TNFRSF13B")
ABC_genes <- c("DUSP4","FCRL5","FCRL4","ITGAX","TBX21","CR2")
Bgc_genes <- c("BCL6", "AICDA", "STMN1","MEF2B")
Plasma_genes <- c("MZB1", "PRDM1","IRF4","XBP1","JCHAIN")

DefaultAssay(B_obj_integrated) <- "RNA"
features <- list("Bn"=Bn_genes,
                 "Bmen"=Bmen_genes,
                 "ABC"=ABC_genes,
                 "Bgc"=Bgc_genes,
                 "Plasma"=Plasma_genes)

Idents(B_obj_integrated) <- "B_main.type"
Idents(B_obj_integrated) <- factor(Idents(B_obj_integrated), 
                                   levels = rev(c("Bn","Bmen","ABC","Bgc","Plasma")))

pdf("B_celltype_marker_exp.pdf",width = 7,height = 3)
DotPlot(object = B_obj_integrated, features = features) +
  scale_color_gradient2(low = "#377594", mid = "#ffffff", high = "#9F1626")+
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,color = "black"),
        axis.text.y =element_text(color = "black"),
        legend.position="bottom")+
  xlab('')+ylab('')
dev.off()


#-----Extended Data Fig. 6c-----####
library(scales)
library(ggplot2)
library(ggthemes)

GO_res <- read.csv("ABC_GO_res.txt",stringsAsFactors = F,header = T,sep = "\t")

GO_res$logP <- -log10(GO_res$p.adjust)

library(dplyr)
GO_top20 <- GO_res %>% 
  group_by(Cluster) %>% 
  top_n(n = 20, wt = logP)

plot_data <- as.data.frame(GO_top20)

plot_data <- plot_data %>%
  dplyr::arrange(logP) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))
custom_colors <- colorRampPalette(c("#3daeb7", "#eeeeee", "#db643e"))(100)


pdf("ABC_GO_res.pdf",width = 9,height = 5)
ggplot(plot_data, aes(x = Count, y = Description)) +
  geom_col(aes(fill = -log10(p.adjust)), width = 0.8) +
  scale_fill_gradientn(colors = custom_colors, 
                       values = rescale(seq(0, 6, length.out = 100)), 
                       name = "-Log10(p.adjust)")+
  labs(x = "Number of genes", y = NULL)+
  theme_few()+
  theme(axis.text = element_text(color = "black"),
        axis.text.x =element_text(size=10,color = "black"), 
        axis.text.y=element_text(size=10,color = "black"))+
  ggtitle("ABC")
dev.off()



#-----Extended Data Fig. 6d-----####

ABC_TAAB_score <- read.csv("ABC_TAAB_score.txt",stringsAsFactors = F,header = T,sep = "\t")

pdf("ABC_TAAB_core_TCGA.pdf",width = 5,height = 5)
ggplot(ABC_TAAB_score,
       mapping=aes(x=ABC,y=PMID39047727))+
  geom_point(size=3,color = "#A92226")+
  geom_smooth(method = 'lm',#线性回归
              formula = 'y ~ x',color = "#A92226",fill = "#A92226", 
              alpha = 0.3)+
  labs(x="ABC Signature Score",
       y="TAAB (Zhang et al.) Signature Score")+
  theme_few()+
  theme(axis.text.x =element_text(size=10,color = "black"), 
        axis.text.y=element_text(size=10,color = "black"),
        aspect.ratio=1)+
  stat_cor(method='pearson',label.sep = "\n",size=5)
dev.off()


#-----Extended Data Fig. 7a-----####

all_VJ_FC_res <- read.csv("VJ_FC_res.txt",stringsAsFactors = F,header = T,sep = "\t")
all_VJ_FC_res$New_Tag <- all_VJ_FC_res$Label


all_VJ_FC_res[which(all_VJ_FC_res$Label!="Stable"),]$New_Tag <- paste0(all_VJ_FC_res[which(all_VJ_FC_res$Label!="Stable"),]$Tag,
                                                                       "_",
                                                                       all_VJ_FC_res[which(all_VJ_FC_res$Label!="Stable"),]$Label)

unique_genes <- unique(all_VJ_FC_res$Gene)
set.seed(1) 
random_order <- sample(1:length(unique_genes))
gene_order <- setNames(random_order, unique_genes)

all_VJ_FC_res$GeneOrder <- gene_order[all_VJ_FC_res$Gene]
all_VJ_FC_res$New_logP <- all_VJ_FC_res$logP
all_VJ_FC_res[which(all_VJ_FC_res$LogFC < 0),]$New_logP <- -(all_VJ_FC_res[which(all_VJ_FC_res$LogFC < 0),]$logP)


all_VJ_FC_res$LogFC <- as.numeric(all_VJ_FC_res$LogFC)

all_VJ_FC_res2 <- all_VJ_FC_res
all_VJ_FC_res2$New_logP <- ifelse(all_VJ_FC_res2$New_logP < -10, -7, all_VJ_FC_res2$New_logP)

all_VJ_FC_res2$New_group <- ""
library(dplyr)
update_new_group <- function(df) {
  df <- df %>%
    mutate(New_group = case_when(
      abs(LogFC) > 2 & Pvalue < 0.01 ~ Gene,
      TRUE ~ New_group
    ))
  return(df)
}

all_VJ_FC_res2 <- update_new_group(all_VJ_FC_res2)
head(all_VJ_FC_res2)

group_col <- c("pre_post_Post_treatment_up"="#DB7EB7","pre_post_Pre_treatment_up"="#DB7EB7","T_N_T_up"="#EEA271","T_N_N_up"="#EEA271",
               "T_P_T_up"="#66D8F3","T_P_P_up"="#66D8F3","P_N_P_up"="#E6688E","P_N_N_up"="#E6688E","Stable"="#CCCCCC")

all_VJ_FC_res2$New_Tag <- factor(all_VJ_FC_res2$New_Tag,levels = c("pre_post_Post_treatment_up","pre_post_Pre_treatment_up",
                                                                   "T_N_T_up","T_N_N_up","T_P_T_up","T_P_P_up","P_N_P_up",
                                                                   "P_N_N_up","Stable"))


pdf("all_group_VJ_volcanoplot.pdf",width = 8,height = 6)
ggplot(all_VJ_FC_res2, aes(x = GeneOrder, y = New_logP, size = abs(LogFC), color = New_Tag)) +
  geom_point(alpha = 1) +
  scale_size_continuous(range = c(0.5, 5)) +
  scale_color_manual(values = group_col) +
  geom_text_repel(aes(label = New_group), size = 3, max.overlaps = Inf) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank()
  ) +
  labs(y = "-log10(p-value)", size = "Log2FC", color = "Compare group") 
dev.off()


#-----Extended Data Fig. 7b-----####
library(tidyverse)
library(ggsci)
library(magrittr)
library(ggplot2)
library(ggforce)

sampletype_isotype_per <- read.csv("sampletype_isotype_per.txt",stringsAsFactors = F,header = T,sep = "\t")


# Pre_treatment
pre_data <- sampletype_isotype_per[which(sampletype_isotype_per$Treatment=="Pre_treatment"),]
pre_data <- pre_data[,-1]

df1 <- as.data.frame(table(pre_data$c_gene)) %>% 
  set_colnames( c("c_gene", "freq")) %>% 
  mutate(per = freq/sum(freq)) %>% 
  mutate(label = paste0(round(per,3)*100,'%'))
df1
df1$c_gene <- factor(df1$c_gene,levels = c("IGHM","IGHD","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4"))

pdf("pre_Isotype_per.pdf",width = 5,height = 5)
ggplot(df1, aes(x = 3,y = freq,fill = c_gene)) +
  geom_col(width = 1.5,color = 'white', alpha = 0.8)  + 
  coord_polar(theta = "y")  +
  xlim(c(0.8, 3.8))  +
  scale_fill_manual(values = Isotype_col)+
  theme_void()+
  theme(
    strip.text.x = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  ) + 
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4)+
  ggtitle("Pre_treatment")
dev.off()


#Post_treatment
post_data <- sampletype_isotype_per[which(sampletype_isotype_per$Treatment=="Post_treatment"),]
post_data <- post_data[,-1]

df2 <- as.data.frame(table(post_data$c_gene)) %>% 
  set_colnames( c("c_gene", "freq")) %>% 
  mutate(per = freq/sum(freq)) %>% 
  mutate(label = paste0(round(per,3)*100,'%'))
df2
df2$c_gene <- factor(df2$c_gene,levels = c("IGHM","IGHD","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4"))

pdf("post_Isotype_per.pdf",width = 5,height = 5)
ggplot(df2, aes(x = 3,y = freq,fill = c_gene)) +
  geom_col(width = 1.5,color = 'white', alpha = 0.8)  +
  coord_polar(theta = "y")  +
  xlim(c(0.8, 3.8))  +
  scale_fill_manual(values = Isotype_col)  + 
  theme_void()  +
  theme(
    strip.text.x = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  ) + 
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4)+
  ggtitle("Post_treatment")
dev.off()



#-----Extended Data Fig. 7e-----####
library(Seurat)
library(ggplot2)
library(MASS)

expr <- read.csv("TCGA_processed_expr.txt",stringsAsFactors = F,header = T,sep = "\t")
rownames(expr) <- expr[,1]
expr <- expr[,-1]

ESCC_sample <- read.csv("TCGA_ESCC_sample.txt",stringsAsFactors = F,header = T,sep = "\t")
ESCC_sample$sample <- gsub("-",".",ESCC_sample$sample)
expr <- expr[,which(colnames(expr) %in% ESCC_sample$sample)]

expr_use <- expr[which(rownames(expr) %in% c("CD79A","FCRL4")),]
expr_use <- as.data.frame(t(expr_use))

dens <- kde2d(expr_use$CD79A, expr_use$FCRL4, n = 100, h = c(0.5, 0.5))

data <- expr_use

threshold1 <- 8
threshold2 <- 3
data$group <- ifelse(data$CD79A > threshold1 & data$FCRL4 > threshold2, "High-High",
                     ifelse(data$CD79A > threshold1 & data$FCRL4 <= threshold2, "High-Low",
                            ifelse(data$CD79A <= threshold1 & data$FCRL4 > threshold2, "Low-High", "Low-Low")))

group_counts <- table(data$group)
total_cells <- nrow(data)
group_percentages <- round(group_counts / total_cells * 100, 2)

label_data <- data.frame(
  group = names(group_counts),
  x = c(threshold1 + 1, threshold1 + 1, threshold1 - 1, threshold1 - 1),
  y = c(threshold2 + 3, threshold2 - 1.5, threshold2 + 3, threshold2 - 1.5),
  label = paste0(names(group_counts), "\n", group_percentages, "% (", group_counts, " samples)")
)

myPalette <- colorRampPalette(c("#9ECCD9","#F5E890","#ED3320"))
gradientColors <- myPalette(1000)

pdf("TCGA_CD79A_FCRL4_exp_density.pdf",width = 4.5,height = 4)
ggplot(data, aes(x = CD79A, y = FCRL4)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", bins = 30) +
  scale_fill_gradientn(colors = gradientColors) +
  geom_text(data = label_data, aes(x = x, y = y, label = label), color = "black", size = 3) +
  geom_vline(xintercept = threshold1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = threshold2, linetype = "dashed", color = "red") +
  theme_light() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.text.x =element_text(size=10,color = "black"), 
        axis.text.y=element_text(size=10,color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8)) + 
  coord_fixed() +
  labs(x = paste("CD79A", "Expression"), y = paste("FCRL4", "Expression"))
dev.off()


#-----Extended Data Fig. 7g-----####
cellper_B <- read.csv("cellper_B.txt",stringsAsFactors = F,header = T,sep = "\t")
cellper_T <- read.csv("cellper_T.txt",stringsAsFactors = F,header = T,sep = "\t")

correlation_matrix <- cor(cellper_B, cellper_T, method = "pearson") 

B_obj_integrated <- readRDS("B_obj_annotation.rds")

B_cells <- unique(Idents(B_obj_integrated))
T_cells <- unique(Idents(T_obj_integrated))

cor_matrix <- matrix(NA, nrow = length(B_cells), ncol = length(T_cells))
p_matrix <- matrix(NA, nrow = length(B_cells), ncol = length(T_cells))

for (i in 1:length(B_cells)) {
  for (j in 1:length(T_cells)) {
    cor_test <- cor.test(cellper_B[, i], cellper_T[, j], method = "pearson")
    cor_matrix[i, j] <- cor_test$estimate
    p_matrix[i, j] <- cor_test$p.value
  }
}

cor_matrix_df <- as.data.frame(cor_matrix)
p_matrix_df <- as.data.frame(p_matrix)

rownames(cor_matrix_df) <- B_cells
colnames(cor_matrix_df) <- T_cells
rownames(p_matrix_df) <- B_cells
colnames(p_matrix_df) <- T_cells

cor_df <- cor_matrix_df[,c("CD8_Tex_CXCL13","CD4_Treg_TNFRSF9+")]
cor_df$B_main.type <- rownames(cor_df)
cor_df <- cor_df[,-2]
colnames(cor_df)[1] <- "Cor"

pvalue_df <- p_matrix_df[,c("CD8_Tex_CXCL13","CD4_Treg_TNFRSF9+")]
pvalue_df$B_main.type <- rownames(pvalue_df)
pvalue_df <- pvalue_df[,-2]
colnames(pvalue_df)[1] <- "pvalue"

cor_p_df <- merge(cor_df,pvalue_df,by="B_main.type")

B_color <- c('Bmen'='#86C8EF','Bn'='#C58CBD','Plasma'='#F38185','Bgc'='#F9BA4E','ABC'='#1FAA9F')

cor_p_df <- cor_p_df %>%
  mutate(group = ifelse(Cor > 0, "positive", "negative"))

cor_p_df$B_main.type <- factor(
  cor_p_df$B_main.type, 
  levels = c(
    cor_p_df$B_main.type[cor_p_df$group == "positive"][order(-cor_p_df$Cor[cor_p_df$group == "positive"])],
    cor_p_df$B_main.type[cor_p_df$group == "negative"][order(-cor_p_df$Cor[cor_p_df$group == "negative"])]
  )
)

cor_p_df$significance <- ifelse(cor_p_df$pvalue < 0.05, "*", "")

pdf("B_CD8_CXCL13_cor.pdf",width = 4,height = 3)
ggplot(cor_p_df, aes(x = B_main.type, y = Cor, color = B_main.type)) +
  geom_segment(aes(xend = B_main.type, yend = 0), size = 1, color = "gray") +
  geom_point(size = 5) +
  geom_text(aes(label = significance), vjust = -0.5, size = 6) + 
  scale_color_manual(values = B_color) +
  theme_minimal() +
  labs(x = "", y = "Pearson correlation coefficient") +
  theme(
    axis.text.x = element_text(size = 10, color = "black", angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black", size = 0.5),
    legend.position = "none",
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black")
  )
dev.off()



#-----Extended Data Fig. 7h-----####

B_T_cellper <- read.csv("B_T_cor.txt",stringsAsFactors = F,header = T,sep = "\t")

pdf("ABC_CD8_Tex_CXCL13_cor.pdf",width = 4,height = 4)
ggplot(B_T_cellper, aes(x = CD8_Tex_CXCL13, y = ABC))+
  geom_point(size=3,color = "#74ADD4")+
  geom_smooth(method = 'lm',#线性回归
              formula = 'y ~ x',color = "#C44A53")+
  # labs(x="Intra-tumor heterogeneity"))+
  theme_few()+
  theme(axis.text.x =element_text(size=10,color = "black"), 
        axis.text.y=element_text(size=10,color = "black"),
        aspect.ratio=1)+
  stat_cor(method='pearson',label.sep = "\n",size=5)
dev.off()

