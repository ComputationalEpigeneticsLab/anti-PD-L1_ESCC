# ------------------------------------------------------------------
# Fig. 6
# ------------------------------------------------------------------

#-----Extended Data Fig. 8b-----####

library(tidyr)
library(ggtern)

Mye_obj_integrated <- readRDS("Mye_obj_annotation_final.rds")
metadata <- Mye_obj_integrated@meta.data

celltype_sampletype <- metadata[,c("sample.type","Mye_main.type")]
celltype_sampletype$barcodes <- rownames(celltype_sampletype)
rownames(celltype_sampletype) <- NULL

cell_counts <- celltype_sampletype %>%
  group_by(sample.type, Mye_main.type) %>%
  summarise(count = n()) %>%
  ungroup()
group_totals <- celltype_sampletype %>%
  group_by(sample.type) %>%
  summarise(total = n()) %>%
  ungroup()

percentages <- cell_counts %>%
  left_join(group_totals, by = "sample.type") %>%
  mutate(percentage = count / total * 100)

percentages_per <- as.data.frame(percentages)
percentages_per <- percentages_per[,c(1,2,5)]

wide_df <- pivot_wider(percentages_per, names_from = sample.type, values_from = percentage)
wide_df <- as.data.frame(wide_df)
 
Mye_color <- c('Mon_VCAN'='#D6CFE8','Mast_CTSG'='#F7D49A','Mast_CPA3'='#EA84A3',
               'Macro_APOE'='#2E9AD5','Macro_HLA-DR+'='#EBBABD','Macro_FOLR2'='#BBADFD','Macro_CSF1R'='#77C6CD',
               'Mon_CCL3'='#EE7A75','cDC2'='#EED5BA','Neutrophils'='#84D1AA','Mon_FCGR3A'='#B5B811',
               'pDC'='#D4DB47','tDC'='#D4B16C','cDC1'='#D3D396')

wide_df$Mye_main.type <- factor(wide_df$Mye_main.type, 
                                levels = c('Macro_APOE','Macro_HLA-DR+','Macro_FOLR2','Macro_CSF1R',
                                           'Mon_VCAN','Mon_CCL3','Mon_FCGR3A',
                                           'cDC1','cDC2','tDC','pDC',
                                           'Neutrophils',
                                           'Mast_CTSG','Mast_CPA3'))

pdf("sampletype_celltype_ggtern.pdf",width = 7,height = 6)
ggtern(data = wide_df, aes(x = N, y = P, z = T)) +
  geom_mask() +
  geom_point(aes(color = Mye_main.type), alpha = 1, size = 9) +
  scale_colour_manual(values = Mye_color) +
  geom_text(aes(label = Mye_main.type, color = Mye_main.type), 
            check_overlap = TRUE, vjust = 1, size = 3) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
dev.off()



#-----Extended Data Fig. 9c-----####
library(Seurat)
library(dplyr)

Mye_obj_integrated <- readRDS("Mye_obj_annotation_final.rds")

Idents(Mye_obj_integrated) <- "Mye_main.type"
Mast_obj <- subset(Mye_obj_integrated,idents=c("Mast_CPA3","Mast_CTSG"))
Mast_markers <- FindAllMarkers(Mast_obj,only.pos = T,logfc.threshold = 0.25,min.pct = 0.25,assay = 'RNA',slot = 'data')
Mast_markers.sig <- Mast_markers[which(Mast_markers$p_val_adj<0.05),]

top10_Mast<-Mast_markers.sig %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)

pdf("Mast_DEGs.pdf",width = 11,height = 2)
DotPlot(object = Mast_obj, features = unique(top10_Mast$gene),
) +
  scale_color_gradient2(low = "#377594", mid = "#ffffff", high = "#B55B77")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x =element_text(size=10,color = "black"), 
        axis.text.y=element_text(size=10,color = "black"))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  xlab('')+ylab('')
dev.off()


#-----Extended Data Fig. 9d-----####
library(ggplot2)
library(dplyr)

GO_res <- read.csv('GO_res_mast.txt',row.names = F,col.names = T,sep = '\t',quote = F)

GO_top5 <- GO_res %>% 
  group_by(Cluster) %>% 
  top_n(n = 10, wt = Count)

plot_data <- as.data.frame(GO_top5)
plot_data <- plot_data[,c(1,3,7,10)]

plot_data$xlab <- paste0(plot_data$Cluster,"_",plot_data$Description)

plot_data <- plot_data %>%
  arrange(desc(Count))

Mast_col <- c('Mast_CTSG'='#F7D49A','Mast_CPA3'='#EA84A3')

pdf("Mast_sub_GOBP.pdf",width = 10,height = 5.5)
ggplot(plot_data, aes(x = Count, y = reorder(xlab, -Count), fill = Cluster)) +
  geom_bar(stat = "identity",color = "black") +
  theme_classic()+
  theme(axis.text.x = element_text(size=10,color = "black",angle = 90, hjust = 1),
        axis.text.y=element_text(size=10,color = "black"))+
  scale_fill_manual(values = Mast_col)
dev.off()


