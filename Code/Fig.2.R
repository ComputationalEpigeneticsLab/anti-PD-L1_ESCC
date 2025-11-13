# ------------------------------------------------------------------
# Fig. 2
# ------------------------------------------------------------------

#-----Fig. 2a-----####
library(tidydr)
library(ggplot2)
library(ggrastr)
library(SeuratWrappers)
library(Seurat)

Mal_obj_integrated <- readRDS("Mal_obj_annotation.rds")

pdf("Mal_annotation.pdf",width = 5,height = 5)
p1 <- DimPlot(Mal_obj_integrated, group.by = "Mal_cluster",reduction = "umap",
              label = F,raster=FALSE,pt.size = 0.00001,cols = Mal_color)+
  NoLegend()+theme_dr()+theme(panel.grid=element_blank())+theme(aspect.ratio = 1)
rasterise(p1,dpi = 300)
dev.off()

pdf("Post_PCR.pdf",width = 5,height = 5)
p2 <- DimPlot(subset(Mal_obj_integrated,Response=="PCR"), group.by = "Mal_cluster",reduction = "umap",
              label = F,raster=FALSE,pt.size = 0.00001,cols = Mal_color)+
  NoLegend()+theme_dr()+theme(panel.grid=element_blank())+theme(aspect.ratio = 1)
rasterise(p2,dpi = 300)
dev.off()

pdf("Post_nPCR.pdf",width = 5,height = 5)
p3 <- DimPlot(subset(Mal_obj_integrated,Response=="nPCR"), group.by = "Mal_cluster",reduction = "umap",
              label = F,raster=FALSE,pt.size = 0.00001,cols = Mal_color)+
  NoLegend()+theme_dr()+theme(panel.grid=element_blank())+theme(aspect.ratio = 1)
rasterise(p3,dpi = 300)
dev.off()

pdf("Pre_Treatment.pdf",width = 5,height = 5)
p4 <- DimPlot(subset(Mal_obj_integrated,Response=="Pre_treatment"), group.by = "Mal_cluster",reduction = "umap",
              label = F,raster=FALSE,pt.size = 0.00001,cols = Mal_color)+
  NoLegend()+theme_dr()+theme(panel.grid=element_blank())+theme(aspect.ratio = 1)
rasterise(p4,dpi = 300)
dev.off()



#-----Fig. 2k-----####

# the enrichment result 
library(Seurat)
library(circlize)
library(ComplexHeatmap)

enrichment_res <- read.csv("/Module_enrichement_res_use.txt",stringsAsFactors = F,header = T,sep = "\t")
enrichment_res$logp <- -log10(enrichment_res$p.adjust)

logp <- as.data.frame(enrichment_res$logp)
rownames(logp) <- enrichment_res$Module
colnames(logp) <- "-log10(p.adjust)"

col_fun = colorRamp2(c(1,25, 50), c("#FFFFFF", "#BABABA", "#000000"))

rownames(logp) <- enrichment_res$Cluster
pdf("Module_enrich_num_heatmap.pdf",width = 1.8,height = 6)
pheatmap(as.matrix(logp),
         color = col_fun,
         display_numbers = TRUE,
         name='Gene number',
         number_color = "black",
         number_format = "%.1f",
         # cellwidth = 18,
         # cellheight = 15,
         border_color = 'black', 
         # scale = "row", 
         cluster_rows = F,
         cluster_cols = F, 
         legend = TRUE,
         legend_breaks = c(1,25, 50), 
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 10
)
dev.off()


# the number of module genes
module_gene <- read.csv("module_df_ESCC_example_leiden_base.tsv",stringsAsFactors = F,header = T,sep = "\t")
module_gene <- module_gene[,c(2,3)]
colnames(module_gene) <- c("gene","module")
module_gene$module <- paste0("M",module_gene$module)

module_use <- c('M0','M1','M3','M4','M5','M6','M7','M8','M9','M10','M12','M14','M15','M16','M17','M18','M19','M20')

module_gene <- module_gene[which(module_gene$module %in% module_use),]

module_gene_num <- as.data.frame(table(module_gene$module))
colnames(module_gene_num)[1] <- "module"
rownames(module_gene_num) <- module_gene_num$module

gene_num <- as.data.frame(module_gene_num$Freq)
rownames(gene_num) <- module_gene_num$module
colnames(gene_num) <- "gene_num"

indices <- match(module_use, rownames(gene_num))
gene_num_ordered <- gene_num[indices, , drop = FALSE]


col_fun = colorRamp2(c(20,90, 200), c("#FFFFFF", "#BABABA", "#000000"))

pdf("Module_gene_num_heatmap.pdf",width = 1.8,height = 6)
pheatmap(as.matrix(gene_num_ordered),
         color = col_fun,
         display_numbers = TRUE,
         name='Gene number',
         number_color = "black",
         number_format = "%d",
         # cellwidth = 18,
         # cellheight = 15,
         border_color = 'black', 
         # scale = "row", 
         cluster_rows = F,
         cluster_cols = F, 
         legend = TRUE,
         legend_breaks = c(20,90, 200), 
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 10
)
dev.off()

# comparison of modules in different group
pheatmap(all_group_ave,
         color = col_fun,
         name='Average Expression',number_color = "black",
         cluster_rows = F,
         cluster_cols = F, 
         legend = TRUE,
         legend_breaks = c(-0.4,0, 0.4), 
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 10,
         gaps_col=c(2,2+2,4+4,8+4)
)


#-----Fig. 2l-----####
Mal_obj <- readRDS("Mal_obj_annotation.rds")

module_gene <- read.csv("module_df_ESCC_example_leiden_base.tsv",stringsAsFactors = F,header = T,sep = "\t")
module_gene <- module_gene[,c(2,3)]
colnames(module_gene) <- c("gene","module")
module_gene$module <- paste0("M",module_gene$module)

module_use <- c('M0','M1','M3','M4','M5','M6','M7','M8','M9','M10','M12','M14','M15','M16','M17','M18','M19','M20')

module_gene <- module_gene[which(module_gene$module %in% module_use),]

gcSample <- split(module_gene$gene, module_gene$module)
gcSample <- gcSample[module_use]

gcSample

Mal_obj <- AddModuleScore(object =Mal_obj, features = gcSample,name = names(gcSample),assay = 'RNA')

score_df <- Mal_obj@meta.data[,c(33:50)]
colnames(score_df) <- module_use

correlation_matrix <- cor(score_df, score_df, method = "pearson")

n <- ncol(score_df)
cor_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(colnames(score_df), colnames(score_df)))
p_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(colnames(score_df), colnames(score_df)))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      test_result <- cor.test(score_df[[i]], score_df[[j]], method = "pearson")
      cor_matrix[i, j] <- test_result$estimate
      p_matrix[i, j] <- test_result$p.value
    } else {
      cor_matrix[i, j] <- 1 
      p_matrix[i, j] <- NA
    }
  }
}

cor_matrix[p_matrix > 0.05] <- 0

library(DescTools)
library(ggrastr)

pdf("Module_cor.pdf",width = 10,height = 10)
PlotWeb(m=cor_matrix,col=c("#4682B4","#CD5555"), 
        cex.lab=0.8,las=1,
        args.legend=list(x=-6,y=-4.5,ncol=1,
                         cex=0.6,box.col="grey80"),
        main="")
dev.off()






#-----Extended Data Fig. 2b-----####

PT_obj
add_metadata <- score[-which(score$group=="Normal Epithelial"),] #第一种方法中包含了正常样本细胞的得分
PT_obj <- AddMetaData(PT_obj,metadata = add_metadata)
PT_obj@meta.data <- PT_obj@meta.data[,-27]

pdf("/data/lisi/study/ESCC_New/4-Epithelial/cnv_score_UMAP_method1.pdf",width = 5,height = 5)
FeaturePlot(PT_obj, features = "cnv_score",pt.size = 0.00001,max.cutoff = 'q99') + 
  scale_colour_gradientn(colours = 
                           (rev(c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                                  "#66C2A5", "#5E4FA2"))), 
                         na.value = "transparent", 
                         name = "CNV score",
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  theme_few()+ggtitle("CNV score") + 
  theme_dr()+theme(panel.grid=element_blank())+
  theme(aspect.ratio = 1)
dev.off()



#-----Extended Data Fig. 2c-----####

#----恶性基因的表达####
Mal_color <- c('Normal sample Epithelial'='#5897A5','Malignant Epithelial'='#EA758A','Normal Epithelial'='#67ADDF')

pdf("/data/lisi/study/ESCC_New/4-Epithelial/Mal_annotation_umap.pdf",width = 6,height = 5)
DimPlot(Epi_obj, group.by = "Mal_group",reduction = "umap",label = F,raster=FALSE,pt.size = 0.00001,cols = Mal_color)+
  NoLegend()+theme_dr()+theme(panel.grid=element_blank(),aspect.ratio = 1)+ggtitle('Epithelial (n=56,073 cells)')
dev.off()


library(ggthemes)
library(ggpubr)
library(ggplot2)
Epi_obj@meta.data$Mal_group <- factor(Epi_obj@meta.data$Mal_group,levels = c("Normal sample Epithelial",
                                                                             "Normal Epithelial","Malignant Epithelial"))
gene_use <- c('EPCAM','KRT16','KRT17','KRT18','KRT19','KRT6A','KRT6B')

pdf("/data/lisi/study/ESCC_New/4-Epithelial/Malignant_gene_expr.pdf",width = 12,height = 12)
comparations <- list(c("Malignant Epithelial","Normal sample Epithelial"),c("Malignant Epithelial","Normal Epithelial"))
VlnPlot(Epi_obj,features = gene_use,group.by = 'Mal_group',cols = Mal_color,pt.size = 0)+
  stat_compare_means(comparisons = comparations,method = 'wilcox.test')+
  theme(legend.position = "none")+
  xlab("")+
  theme(
    axis.text.x =element_text(size=10,color = "black"),
    axis.text.y =element_text(size=10,color = "black"))
dev.off()


#-----Extended Data Fig. 3d-----####

library(ggthemes)
library(tidydr)

my_cds <- readRDS("monocle2.rds")

pdf("Mal_Cluster_Pseudotime.pdf",width = 6,height = 5)
plot_cell_trajectory(my_cds, color_by = "Pseudotime",cell_size = 0.8,show_branch_points=F)+ theme_dr()+
  theme(panel.grid=element_blank(),legend.position = "right",aspect.ratio = 1,
  ) + 
  scale_colour_gradient(low = "#FAF5A9", high = "#C52F2D")
dev.off()

Mal_color <- c("Mal_C0"="#F8C87E","Mal_C1"="#9AD0F0","Mal_C2"="#8FC7C6","Mal_C3"="#F196BA")

pdf("Mal_Cluster.pdf",width = 6,height = 5)
plot_cell_trajectory(my_cds, color_by = "Mal_cluster",cell_size = 0.8,show_branch_points=F,backbone_color = "white")+ theme_dr()+
  theme(panel.grid=element_blank(),legend.position = "right",aspect.ratio = 1,
  ) + 
  scale_color_manual(values=Mal_color)
dev.off()

# Proportion of malignant cell subsets
data_df <- t(reducedDimS(my_cds)) %>%
  as.data.frame() %>%
  select_('Component 1' = 1, 'Component 2' = 2) %>%
  rownames_to_column("Cells") %>%
  mutate(pData(my_cds)$State,
         pData(my_cds)$Pseudotime,
         pData(my_cds)$sample,
         pData(my_cds)$Mal_cluster,
         pData(my_cds)$Treatment,
         pData(my_cds)$Response)
colnames(data_df) <- c("cells","Component_1","Component_2","State",
                       "Pseudotime","sample","Mal_cluster","Treatment","Response")

reduced_dim_coords <- reducedDimK(my_cds)
ica_space_df <- Matrix::t(reduced_dim_coords) %>%
  as.data.frame() %>%
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  mutate(sample_name = rownames(.), sample_state = rownames(.))

dp_mst <- minSpanningTree(my_cds)

edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  select_(source = "from", target = "to") %>%
  left_join(ica_space_df %>% 
              select_(source="sample_name", 
                      source_prin_graph_dim_1="prin_graph_dim_1", 
                      source_prin_graph_dim_2="prin_graph_dim_2"), 
            by = "source") %>%
  left_join(ica_space_df %>% 
              select_(target="sample_name", 
                      target_prin_graph_dim_1="prin_graph_dim_1", 
                      target_prin_graph_dim_2="prin_graph_dim_2"), 
            by = "target")


pie_data <- data_df[,c(1,2,3,7)]
x_breaks <- seq(-14, 8, length.out = 11)
y_breaks <- seq(-7, 5, length.out = 11)

pie_data$grid_x <- cut(pie_data$Component_1, breaks = x_breaks, include.lowest = TRUE)
pie_data$grid_y <- cut(pie_data$Component_2, breaks = y_breaks, include.lowest = TRUE)

grid_counts <- table(pie_data$grid_x, pie_data$grid_y, pie_data$Mal_cluster)
grid_counts <- as.data.frame(grid_counts)
colnames(grid_counts) <- c("grid_x", "grid_y", "Mal_cluster", "Count")

grid_counts <- grid_counts %>%
  group_by(grid_x, grid_y) %>%
  mutate(Total = sum(Count), Proportion = Count / Total)

grid_counts$grid_y <- factor(grid_counts$grid_y, levels = rev(levels(grid_counts$grid_y)))

grid_counts <- grid_counts %>%
  mutate(
    x_min = as.numeric(str_extract(grid_x, "-?\\d+\\.?\\d*")),
    x_max = as.numeric(str_extract(grid_x, "(?<=,)-?\\d+\\.?\\d*")),
    y_min = as.numeric(str_extract(grid_y, "-?\\d+\\.?\\d*")),
    y_max = as.numeric(str_extract(grid_y, "(?<=,)-?\\d+\\.?\\d*")),
    x_center = (x_min + x_max) / 2,
    y_center = (y_min + y_max) / 2
  )
picture_data <- as.data.frame(grid_counts)[,c(3,4,5,6,11,12)]

max_total <- max(grid_counts$Total)

pdf("Mal_cluster_pie.pdf",width = 6,height = 5)
ggplot(grid_counts, aes(x = "", y = Proportion, fill = Mal_cluster)) +
  geom_bar(stat = "identity",color = "white", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(grid_y ~ grid_x) +
  scale_fill_manual(values = Mal_color)+
  labs(x = "Grid X",
       y = "Proportion") +
  theme_void() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 1)
dev.off()


#-----Extended Data Fig. 3e-----####
library(ggrepel)
library(ggthemes)

Mal_C0_mark_mat <- read.csv("Mal_C0_pyscenic_res.txt",stringsAsFactors = ,header = T,sep = "\t")

pdf('Mal_C0_tf_rss_V2.pdf',width = 4,height = 5)
ggplot(data = Mal_C0_mark_mat) + 

  geom_point(mapping = aes(x = Order, 
                           y = Mal_C0,
                           color = ifelse(is.na(label_mark), '#807F7F', '#C3423F'),
                           size = ifelse(is.na(label_mark), 1, 1),
                           alpha = ifelse(is.na(label_mark), 1, 1)),
             show.legend = FALSE) + 
  scale_size_identity() +
  scale_alpha_identity() +
  scale_color_identity() +  
  theme_classic() +

  geom_text_repel(data = Mal_C0_mark_mat,
                  aes(x = Order, y = Mal_C0, label = label_mark),
                  size = 4, 
                  color = "#C3423F",
                  max.overlaps = 30,
                  segment.size = 0.4, 
                  min.segment.length = 0.5, 
                  box.padding = 0.2,
                  arrow = arrow(length = unit(0, "npc")), 
                  show.legend = FALSE) +
  labs(x = "Rank", y = "Regulon specificity score") +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"))+
  ggtitle("Mal_C0")
dev.off()


#-----Extended Data Fig. 3g-----####

library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(progeny)
library(Seurat)
library(decoupleR)

Mal_obj_integrated <- readRDS("Mal_obj_annotation.rds")

Idents(Mal_obj_integrated) <- "Mal_cluster"
CellsClusters <- data.frame(Cell = names(Idents(Mal_obj_integrated)),
                            CellType = as.character(Idents(Mal_obj_integrated)),
                            stringsAsFactors = FALSE) 

Mal_obj_integrated <- progeny(Mal_obj_integrated, scale=FALSE, 
                              organism="Human", top=100, perm=1, 
                              return_assay = TRUE)  #"Human" Mouse

Mal_obj_integrated@assays$progeny
Mal_obj_integrated@assays$progeny %>%dim()
Mal_obj_integrated@assays$progeny@data[,1:19]


Mal_obj_integrated <- Seurat::ScaleData(Mal_obj_integrated, assay = "progeny") 

progeny_scores_df <- 
  as.data.frame(t(GetAssayData(Mal_obj_integrated, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 50
myColor = colorRampPalette(c("#63B8CE", "white","#D9533D"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))

pdf("progeny.pdf",width = 5,height = 2)
progeny_hmap = pheatmap(summarized_progeny_scores_df,
                        fontsize=10,fontsize_row = 10, 
                        treeheight_row = 0,treeheight_col = 0,
                        color=myColor, breaks = progenyBreaks, 
                        # main = "PROGENy (100)",
                        angle_col = 90,
                        # treeheight_col = 3,  
                        border_color = NA)
dev.off()





