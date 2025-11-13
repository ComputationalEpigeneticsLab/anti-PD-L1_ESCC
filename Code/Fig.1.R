
# ------------------------------------------------------------------
# Fig. 1
# ------------------------------------------------------------------

#-----Fig. 1e-----####
plot_data <- read.csv("Cellratio.txt",stringsAsFactors = F,header = T,sep = "t")

pdf("radarchart.pdf",width = 6,height = 5)
radarchart(plot_data,
           axistype = 1,
           pcol = mycolor_sampletype, #线条颜色
           pfcol = alpha(mycolor_sampletype, 0.1), #填充颜色
           plwd = 2, #线宽
           plty = 1, #线形
           cglcol = 'black', #网格颜色
           # cglty = 1, #网格线型
           axislabcol = 'white', #轴标签颜色
           caxislabels = seq(0,0.5,0.1), #轴标签显示
           cglwd = 0.8, #网格线宽
           vlcex = 1 #分组标签大小
)+theme(aspect.ratio = 1)

legend(x = 1.1, y = 0.5, 
       legend = rownames(plot_data[-c(1,2),]), 
       bty = "n", pch = 20 , 
       col = mycolor_sampletype,
       text.col = "black", 
       cex = 1, pt.cex = 2)
dev.off()


#-----Fig. 1f-----####

library(ggplot2)
library(Seurat)
library(cowplot)
library(Matrix)
library(viridis)

obj_use <- readRDS("sample_list.rds")
cell2_res <- read.csv("st_cell2location_res.csv",row.names = 1)

obj_use <- AddMetaData(obj_use, metadata = cell2_res)

colorscale <- magma(8)

celltype <- c("Bcells","Endothelial","Epithelial","Fibroblasts","Mast","Myeloid","Plasma","T_NK")
colnames(obj_use@meta.data)[5:12] <- c("Bcells","Endothelial","Epithelial","Fibroblasts","Mast","Myeloid","Plasma","T_NK")

for (i in 1:length(celltype)) {
  
  celltype_use <- celltype[i]
  
  pdf(paste0("Visualization/","sample1_",celltype_use,".pdf"))
  p <- SpatialFeaturePlot(obj_use, image.alpha = 0,features=celltype_use, pt.size.factor=1.8,alpha = c(1, 1),stroke=NA)+
    scale_fill_gradientn(colors = colorscale, 
                         limits = c(0,4),
                         oob = scales::squish,
                         na.value = "#FCFDBFFF",
                         name = "fraction",
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black"))+
    theme_void()+ 
    theme(
      axis.ticks=element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(color = "white", fill = NA, size =2),
    )
  print(p)
  dev.off()
  
  cat(i,"\n")
}



#-----Fig. 1g-----####
library(reshape2)

Indir <- "cell_type_fractions_mean/"
file_res <- list.files(Indir,pattern = "csv")

colors = c("#FFF7EC","#EAE9B5","#E3C4A5", "#D3938F", "#C66D7E", "#BB4E70", "#B3476D")
colors1 = colorRampPalette(colors)(50)

factor_res <- read.csv(paste0(Indir,"n_fact14.csv"),stringsAsFactors = F,header = T,sep = ",",row.names = 1)
colnames(factor_res) <- sub("mean_cell_type_factorsfact_", "factor_", colnames(factor_res))

factor_res <- as.data.frame(t(factor_res))
factor_res <- factor_res[c(2,3,4,10,13,14),]
factor_res <- as.matrix(factor_res)
factor_res_long <- melt(factor_res)
colnames(factor_res_long) <- c("Factor","Celltype","value")
factor_res_long$Factor <- factor(factor_res_long$Factor, levels = rev(unique(factor_res_long$Factor)))

pdf(paste0(Indir,"Factor_rev_filnally.pdf"),width = 10,height = 3)
ggplot(factor_res_long, aes(x = Celltype,y = Factor)) +## global aes
  geom_point(aes(fill = value,size =value),shape=21,color = "black")  +
  scale_fill_gradientn(colors = colors1) +
  scale_color_gradientn(colors = colors1) + 
  labs(x='',y='',fill='Fraction',color='Fraction',size='Fraction')+
  scale_size(range = c(0,5), limits = c(0, 1))+
  theme_minimal()+  
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2),        
        panel.grid.minor = element_blank(),        
        panel.background = element_blank(),         
        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
        axis.text.x=element_text(angle=90,hjust=1,color = 'black'), 
        axis.text.y=element_text(hjust=0, color = 'black',size=12))
dev.off()



#-----Extended Data Fig. 1d-----####
library(tidyverse)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggsci)
library(tidyr)
library(dplyr)
library(cowplot)

wide_df <- read.csv("count_df.txt",stringsAsFactors = F,header = T,sep = "\t")

tree <- hclust(vegan::vegdist(t(wide_df), method = 'bray'),
              method = 'average') %>%
  as.phylo()

tree <- groupClade(tree, .node=c(31))

p1 <- ggtree(tree,branch.length="none") + 
  geom_tiplab(color="black",size=2) + xlim(NA, 70)
p1

wide_df$mean = apply(wide_df, 1, mean)
wide_df <- wide_df[order(wide_df$mean),]
wide_df$celltype.main = factor(rownames(wide_df), levels = rownames(wide_df))
wide_df <- wide_df[, -c(ncol(wide_df)-1)]

p1 <- wide_df %>%  
  reshape2::melt(id.vars = 'celltype.main') %>%  
  ggplot(aes(variable, value, fill = celltype.main))+  
  geom_bar(stat = 'identity', position = 'fill')+  
  scale_x_discrete(limits = samples)+  
  scale_fill_manual(values = mycolor)+  
  scale_y_continuous(labels = scales::percent)+  
  theme_classic()+  
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black",angle = 90,hjust = 1),
    axis.line = element_blank())+
  labs(y = 'Percentage')

pdf("sample_celltype_percent_cluster.pdf",width = 7,height = 4)
print(p1)
dev.off()


#-----Extended Data Fig. 1f-----####
library(ggplot2)
library(vegan)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(factoextra)


polt_data <- read.csv("pca_result.txt",stringsAsFactors = F,header = T,sep = "t")
polt_data$Treatment <- factor(polt_data$Treatment,levels = c("Pre_treatment","Post_treatment"))


PCA <- ggplot(polt_data, aes(x = PC1, y = PC2, color = sample.type,shape=Treatment)) +
  geom_point(size = 4) +
  scale_color_manual(values = mycolor_sampletype)+
  
  labs(x = "PC_1 (24.6%)", y = "PC_2 (13.8%)") +
  theme_bw()+
  theme(axis.text.x =element_text(size=10,color = "black"), 
        axis.text.y=element_text(size=10,color = "black"),
        # axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))+
  stat_ellipse(aes(group = Treatment),lwd=1,level = 0.8)

PC1_density <- ggplot(polt_data) +
  geom_density(aes(x=PC1, group=Treatment,fill=Treatment),
               color="black", alpha=0.6,position = 'identity',
               show.legend = F) +
  scale_fill_manual(values=mycolor_treat) +
  scale_linetype_manual(values = c("solid"))+
  scale_y_discrete(expand = c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())

PC2_density <- ggplot(polt_data) +
  geom_density(aes(x=PC2, group=Treatment,fill=Treatment),
               color="black", alpha=0.6,position='identity',
               show.legend = F) +
  scale_fill_manual(values=mycolor_treat) +
  scale_linetype_manual(values = c("solid"))+
  scale_y_discrete(expand = c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  coord_flip()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

library(ggplotify)
library(aplot)
pdf("treatment_sample.type_PCA.pdf",width = 6.7,height = 5,onefile=F)
p2 <- PCA %>% 
  insert_top(PC1_density,height = 0.3) %>% 
  insert_right(PC2_density,width=0.3) %>% 
  as.ggplot()
print(p2)
dev.off()


















