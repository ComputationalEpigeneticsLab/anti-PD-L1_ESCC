# ------------------------------------------------------------------
# Fig. 3
# ------------------------------------------------------------------

#-----Fig. 3b-----####
library(dplyr)
library(Seurat)
library(patchwork)
library(ggrepel)
library(ggthemes)

T_obj_integrated_filted <- readRDS("T_obj_annotation_final.rds")
T_obj_integrated_filted@meta.data$T_main.type <- factor(T_obj_integrated_filted@meta.data$T_main.type,
                                                        levels = c('CD4_Tn','CD4_Tfh','CD4_Th17','CD4_Treg_TNFRSF9+','CD4_Treg_TNFRSF9-',
                                                                   'CD8_Tcm','CD8_Tem','CD8_Tex_GZMK','CD8_Tex_CXCL13','CD8_HSP','MAIT',
                                                                   'NK_CD16','NK_XCL1','T_Pro','T_MALAT1'))

gene_use<-c(
  'CCR7','TCF7','SELL','LEF1',
  'PASK','FAAH2',
  'IL17F','IL17A',
  'FOXP3','IKZF2','TNFRSF9',
  'IL7R','FOS','ANXA1',
  "GZMK","KLRG1","NKG7","CST7","CCL5","GZMA","GZMB",'GNLY','RPF1','FGFBP2','GZMH', 'CXCL13',
  'TIGIT','LAG3','PDCD1','HAVCR2',
  'HSPA1A','HSPA1B','HSPA6','HSPD1',
  'SLC4A10','AQP3','KLRB1','TRAV1-2',
  'FCGR3A','AREG','XCL1',
  'TOP2A','MKI67',
  'MALAT1','NEAT1'
)

gene_cell_exp <- AverageExpression(T_obj_integrated_filted,
                                   features = unique(gene_use),
                                   group.by = 'T_main.type',
                                   slot = 'data',assays = 'RNA') 

gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
colnames(marker_exp)

marker_exp <- marker_exp[,c('CD4_Tn','CD4_Tfh','CD4_Th17','CD4_Treg_TNFRSF9+','CD4_Treg_TNFRSF9-',
                            'CD8_Tcm','CD8_Tem','CD8_Tex_GZMK','CD8_Tex_CXCL13','CD8_HSP','MAIT',
                            'NK_CD16','NK_XCL1','T_Pro','T_MALAT1')]


colnames(marker_exp)
library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(-2, 0, 2), c("#4796B3", "white", "#E53845"))

marker_exp <- marker_exp[,c('CD4_Tn','CD4_Tfh','CD4_Th17','CD4_Treg_TNFRSF9+','CD4_Treg_TNFRSF9-',
                          'CD8_Tcm','CD8_Tem','CD8_Tex_GZMK','CD8_Tex_CXCL13','CD8_HSP','MAIT',
                          'NK_CD16','NK_XCL1','T_Pro','T_MALAT1')]

pdf('marker_expression.pdf',width = 6,height =10)
pheatmap(marker_exp, 
         color = col_fun,
         name='Expression',
         cellwidth = 15,
         cellheight = 10,
         scale = "row", 
         cluster_rows = F, 
         cluster_cols = F, 
         legend = TRUE, 
         legend_breaks = c(-2, 0, 2), 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize = 8,
         gaps_row=c(4,4+2,6+2,8+3,11+3,14+12,26+4,30+4,34+4,38+3,41+2,43+2),
         gaps_col = c(5,5+5,10+1,11+2,13+2)
)
dev.off()


#-----Fig. 3c-----####
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

sample.type_res <- readRDS("T_sample_type_OR.rds")
treat_res <- readRDS("T_treatment_OR.rds")
all_OR_res <- rbind(sample.type_res,treat_res)

heatdata <- all_OR_res[,c(1,2,5)]
library(reshape2)
heatdata <- as.data.frame(heatdata)

heatdata <- dcast(heatdata, rid ~ cid, value.var = "OR")
rownames(heatdata) <- heatdata[,1]
heatdata <- heatdata[,-1]

library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(0, 1, 5), c("#6DC2E9", "#FFFFFF", "#D64756"))
heatdata <- heatdata[c(3,1,2,5,4,7,8,10,9,6,11,12,13,15,14),]


pdf("OR_heatmap.pdf",width = 8,height = 4)
pheatmap(t(heatdata),
         color = col_fun,display_numbers = TRUE,
         name='OR',number_color = "black",
         cellwidth = 25,
         cellheight = 25,
         border_color = 'black', 
         # scale = "row", 
         cluster_rows = F,
         cluster_cols = F, 
         legend = TRUE,
         legend_breaks = c(0, 1.5, 3), 
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 8,
         gaps_col=c(5,5+6,11+2,13+2),
         gaps_row=c(3,3+2)
)
dev.off()



#-----Fig. 3f-----####
library(Seurat)
library(dplyr)
library(ggthemes)
library(tidydr)

tmp <- read.csv("T_obj_tcr.txt",stringsAsFactors = F,header = T,sep = "\t")
tmp <- na.omit(tmp)

tmp_clean <- tmp[,c('T_main.type','clone','UMAP1','UMAP2','cloneSize')]

avg_umap_with_clonesize <- tmp_clean %>%
  group_by(clone, T_main.type) %>%
  summarise(
    avg_UMAP1 = mean(UMAP1, na.rm = TRUE),
    avg_UMAP2 = mean(UMAP2, na.rm = TRUE),
    total_cloneSize = sum(cloneSize, na.rm = TRUE)
  )
avg_umap_with_clonesize <- as.data.frame(avg_umap_with_clonesize)

avg_umap_with_clonesize$cloneSizelog <- log2(avg_umap_with_clonesize$total_cloneSize+1)

avg_umap_with_clonesize <- avg_umap_with_clonesize[order(avg_umap_with_clonesize$cloneSizelog),]

tmp_n2 <- avg_umap_with_clonesize[avg_umap_with_clonesize$cloneSizelog>1,]

T_color <- c("CD8_Tex_CXCL13"="#639EBC","CD8_Tcm"="#D47AAD","CD8_Tem"="#94C54B", "CD8_Tex_GZMK"="#E9866D","MAIT"="#D79ABF",
             "Unknown"="#F1AE7A","CD8_HSP"="#87BEE7",
             
             "CD4_Treg_TNFRSF9+"="#BE589D","CD4_Treg_TNFRSF9-"="#D1A8CC","CD4_Tcm"="#C84853","CD4_Th17"="#EFC384","CD4_Tn"="#5B98C7",
             "CD4_Tfh"="#6F983A",
             "NK_CD16"="#F0BB64","NK_XCL1"="#96B3AE","T_Pro"="#C6A37C","T_MALAT1"="#66C4A9")

library(ggplot2)
library(dplyr)

size_breaks <- c(1, 4, 8, 12, 18)
size_labels <- c("1-4", "5-8", "9-12", "13-18")
size_values <- c(1, 2, 3, 4)

tmp_n2 <- tmp_n2 %>%
  mutate(size_category = cut(cloneSizelog, breaks = size_breaks, labels = size_labels))

pdf("T_cluster_TCR_point.pdf",width = 6,height = 5)
ggplot(tmp_n2, aes(x = avg_UMAP1, y = avg_UMAP2, fill = T_main.type, size = size_category)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_manual(values = T_color) +
  scale_size_manual(values = size_values) +
  theme_dr() +
  theme(panel.grid=element_blank(),aspect.ratio = 1)
dev.off()



#-----Fig. 3g-----####
library(dplyr)
library(Seurat)
library(patchwork)
library(ggrepel)
library(ggthemes)

augur <- readRDS("augur_res.rds")

data <- plot_lollipop(augur)$data
data <- data[order(data$auc, decreasing = TRUE), ]
data$cell_type <- factor(data$cell_type, levels = data$cell_type)

data_order <- data[order(data$auc,decreasing = T),]


pdf("T_main.type_Augur.pdf",width = 6,height = 7)
plot_lollipop(augur) +
  geom_segment(aes(xend = cell_type, yend = 0.5), size = 1) +
  geom_point(size = 6, aes(color = cell_type)) +
  scale_color_manual(values = T_color)+theme(aspect.ratio = 1)+
  theme(axis.text.x =element_text(size=10,color = "black"), 
        axis.text.y=element_text(size=10,color = "black"))
dev.off()



#-----Extended Data Fig. 4d-----####
sample_type_TCR <- read.csv("scatterClonotype_res.rds")

class_col <- c("N unique"="#1AA9BC","N expanded"="#529CA1","T unique"="#EDB7B7","T expanded"="#C3638A","Other shared"="#D8D8D8")

pdf("N_T_shared_TCR.pdf",width = 6,height = 5)
ggplot(sample_type_TCR, aes(x = N.fraction, y = T.fraction, colour = class, size = sum)) +
  geom_point(aes(fill = class), shape = 21, colour = "black") +
  scale_fill_manual(values = class_col) +
  scale_size_continuous(range = c(1, 10), name = "Total n") +
  theme_few()+
  geom_abline(slope = 1, intercept = 0,
              alpha = 0.4, lty = 2) + scale_y_sqrt() + scale_x_sqrt()+
  theme(axis.text.x=element_text(vjust=1,size=10,color = "black"))+
  theme(axis.text.y=element_text(vjust=1,size=10,color = "black"),
        aspect.ratio=1)
dev.off()


class_T_main.type <- read.csv("TCR_res.txt",stringsAsFactors = F,header = T,sep = "\t")

pdf("N_T_shared_TCR_T.maintype_distribution.pdf",width = 4,height = 5)
class_T_main.type %>% 
  ggplot(aes(x = class, fill = T_main.type)) +
  geom_bar(position = position_fill(),width = 0.8,color="black") + 
  scale_fill_manual(values = T_color)+
  theme_classic() +
  labs(y = 'Percent',x="") +
  theme(legend.position="right")+scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x =element_text(size=10,color = "black"), axis.text.y=element_text(size=10,color = "black"))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))
dev.off()



#-----Extended Data Fig. 4f-g-----####
library(ggrepel)

# CD8
CD8_TCR_share <- read.csv("CD8_TCR_share.rds")

pdf('CD8_shared_point.pdf',width = 6,height = 5)
ggplot(CD8_TCR_share, aes(x=logpvalue, y=number,color=pair)) + 
  geom_point(size=5)+
  geom_hline(yintercept= 10,linetype='dashed')+
  geom_vline(xintercept= -log10(0.05),linetype='dashed')+
  theme_few()+
  scale_color_manual(values = c('#FFC6AC','#FF9BC5','#A1CCD1','#A8CD9F','#E78F81',
                                '#cca8e9','#E5E483','#10ddc2','#99ddcc','#D885A3'))+
  geom_text_repel(aes(label=pair),CD8_TCR_share,min.segment.length = 0.1)+
  xlab('-log10 (Pvalue of Correlation)')+
  theme(aspect.ratio = 1,
        axis.text.x =element_text(size=10,color = "black"), 
        axis.text.y=element_text(size=10,color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  ylab('Number of shared clonotypes')+ggtitle('CD8 T cell clusters')+NoLegend()
dev.off()


# CD4
CD4_TCR_share <- read.csv("CD4_TCR_share.rds")

pdf('CD4_shared_point.pdf',width = 6,height = 5)
ggplot(CD4_TCR_share, aes(x=logpvalue, y=number,color=pair)) + 
  geom_point(size=5)+
  geom_hline(yintercept= 10,linetype='dashed')+
  geom_vline(xintercept= -log10(0.05),linetype='dashed')+
  theme_few()+
  scale_color_manual(values = c('#FFC6AC','#FF9BC5','#A1CCD1','#A8CD9F','#E78F81',
                                '#cca8e9','#E5E483','#10ddc2','#99ddcc','#D885A3'))+
  geom_text_repel(aes(label=pair),CD4_TCR_share,min.segment.length = 0.1)+
  xlab('-log10 (Pvalue of Correlation)')+
  theme(aspect.ratio = 1,
        axis.text.x =element_text(size=10,color = "black"), 
        axis.text.y=element_text(size=10,color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  ylab('Number of shared clonotypes')+ggtitle('CD4 T cell clusters')+NoLegend()
dev.off()













