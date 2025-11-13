# ------------------------------------------------------------------
# Fig. 4
# ------------------------------------------------------------------

#-----Fig. 4a-----####

dm <-readRDS("Tem_Tex_destiny_res.rds")

Tex_col <- c("CD8_Tem"="#94C54B","CD8_Tex_CXCL13"="#639EBC","CD8_Tex_GZMK"="#E9866D")

pdf("destiny_Tem_Tex.pdf",width = 6,height = 5)
ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = cell_type2),size=0.3)+
  theme_dr()+
  scale_color_manual(values = Tex_col)+
  theme(panel.grid=element_blank(),aspect.ratio = 1)
dev.off()



#-----Fig. 4b-----####
library(Seurat)
library(monocle)
library(splines)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(ggrastr)


my_cds <- readRDS("/data/lisi/study/ESCC_New/5-T_NK/CD8/Tem_Tex_monocle2_cds.rds")

pseudotime <- my_cds@phenoData@data$Pseudotime

gene_expression <- exprs(my_cds)
genes <- unique(c('PDCD1', 'CTLA4', 'CXCL13','LAG3'))
genes2 <- intersect(genes, rownames(gene_expression))


p_list <- list()

for(i in genes2) {
  gene_expression_i <- gene_expression[i, ]
  
  fit_data <- data.frame(pseudotime = pseudotime, gene_expression = gene_expression_i)
  
  fit <- lm(gene_expression ~ ns(pseudotime, df = 3), data = fit_data)
  
  new_pseudotime <- seq(min(pseudotime), max(pseudotime), length.out = 100)
  fit_values <- predict(fit, newdata = data.frame(pseudotime = new_pseudotime), interval = "confidence")
  
  fit_curve <- data.frame(
    pseudotime = new_pseudotime,
    fit = fit_values[, "fit"],
    lower = fit_values[, "lwr"],
    upper = fit_values[, "upr"]
  )
  
  p <- ggplot() +
    geom_ribbon(data = fit_curve, aes(x = pseudotime, ymin = lower, ymax = upper), fill = "#9ecae1", alpha = 0.4) +
    geom_line(data = fit_curve, aes(x = pseudotime, y = fit), color = "#223D6C", size = 1) +
    labs(title = i, x = "Pseudotime", y = "Expression") +
    theme_bw() +
    theme(
      aspect.ratio = 1/1,
      axis.title.x = element_text(size = 10),  
      axis.title.y = element_text(size = 10),  
      axis.text.x = element_text(size = 10),  
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 10)
    )
  
  p_list[[i]] <- p
}

pdf("exhauted_expr.pdf", width=8, height=4)
plot_grid(plotlist = p_list, ncol = 2)
dev.off()


plot_df <- read.csv("pheno.data_score.txt",stringsAsFactors = F,header = T,sep = "\t")
pdf("CD8_Tem_Tex_monocle2_exhuastion.pdf",width = 6,height = 5)
p <- ggplot(plot_df, aes(x=Pseudotime, y=exhuastion2,color=T_main.type)) + 
  geom_point(shape=20,size=3,stroke =0,alpha=1)+
  theme_few()+
  # theme(panel.grid=element_blank())+
  theme(aspect.ratio = 1,
        axis.text.x =element_text(size=10,color = "black"),
        axis.text.y =element_text(size=10,color = "black"))+
  scale_color_manual(values = Tex_col)+
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "black")+
  ylab("Exhuastion Score")
print(rasterise(p,dpi = 300))
dev.off() 


#-----Fig. 4d-----####
library(ggplot2)

Tem_Tex_obj <- readRDS("Tem_Tex_RunPCA_VariableFeatures.rds")

Idents(Tem_Tex_obj) <- "T_main.type"
Tem_Tex_obj <- subset(Tem_Tex_obj,idents=c("CD8_Tex_GZMK","CD8_Tex_CXCL13"))

gene_list <- read.csv("gene_list_plus.txt",stringsAsFactors = F,header = T,sep = "\t")
gene_list <- split(as.matrix(gene_list)[,1], gene_list[,2])
gene_list <- gene_list[c(8)]

Tem_Tex_obj@meta.data$group <- factor(Tem_Tex_obj@meta.data$group,levels = c("Post_treatment_CD8_Tex_CXCL13","Pre_treatment_CD8_Tex_CXCL13",
                                                                             "Post_treatment_CD8_Tex_GZMK","Pre_treatment_CD8_Tex_GZMK"))
DefaultAssay(Tem_Tex_obj) <- "RNA"

Idents(Tem_Tex_obj) <- "group"

pdf("gene_list_exp.pdf",width = 6,height = 3)
DotPlot(object = Tem_Tex_obj, features = gene_list) +
  scale_color_gradient2(low = "#377594", mid = "#ffffff", high = "#9F1626")+
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,color = "black"),
        axis.text.y =element_text(color = "black"),legend.position = "bottom")+
  xlab('')+ylab('')
dev.off()


#-----Fig. 4e-----####

library(Seurat)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(viridis)
library(ggunchull)

CD8_obj <- readRDS("CD8_obj_annotation.rds")

plot_umap_density <- function(data, type, base_size = 12, base_family = "") {
  ggplot(data = data, aes(x = UMAP_1, y = UMAP_2)) +
    theme_black(base_size, base_family) +
    xlim(c(min(data$UMAP_1)-2, max(data$UMAP_1)+2)) +
    ylim(c(min(data$UMAP_2)-2, max(data$UMAP_2)+2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, alpha = 1) +
    geom_point(color = 'grey90', size = 0.05, alpha = 0.4) +
    scale_fill_viridis(option = "magma", alpha = 1) +
    labs(title = paste0("UMAP Density Plot for ", type)) +
    theme(legend.position = "none")
}

theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white",  fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size*0.8, color = "white"),
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
      legend.position = "none",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      panel.background = element_rect(fill = "black", color  =  NA),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0, "lines"),
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size*0.8, color = "white"),
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size*1.2, color = "white"),
      plot.margin = unit(rep(0, 4), "lines")
    )
}

meta <- as.data.frame(CD8_obj@meta.data)
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)

coord <- Embeddings(object =  CD8_obj[["umap"]])

coord = data.frame(ID = rownames(coord), coord)

meta <- merge(meta, coord, by = "ID")

unique_types <- unique(meta$Response)

lapply(unique_types, function(tissue_type) {
  tissue_data <- meta %>% filter(Response == tissue_type)
  p <- plot_umap_density(tissue_data, tissue_type)+theme(aspect.ratio = 1)
  # 保存PDF
  ggsave(filename = paste0("Density/", tissue_type, "_umap_density.pdf"),
         width = 4, height = 4, plot = p)
})


plotData <- as.data.frame(CD8_obj[["umap"]]@cell.embeddings)
plotData$T_main.type <- CD8_obj$T_main.type

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

identity_label_pos <- plotData %>%
  group_by(T_main.type) %>%
  summarise(
    x = median(UMAP_1),
    y = median(UMAP_2)
  )


pdf("CD8_celltype_border.pdf",width = 4,height = 4)
ggplot(plotData, aes(x = UMAP_1, y = UMAP_2, fill = T_main.type, color = T_main.type)) +
  stat_unchull(alpha = 0, size = 0.3, lty = 2,delta=0.1,color = "black") +
  geom_point(size = 0.000001) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(hjust = 0.05, face = "italic")) +
  guides(color = FALSE, x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_manual(values = T_color) +
  scale_color_manual(values = T_color)
dev.off()


#-----Fig. 4j-----####
library(ggplot2)
library(Seurat)
library(cowplot)
library(Matrix)
library(BayesSpace)
library(viridis)

sce.enhanced <- readRDS("BGM_enhanced_expr.rds")
sce <- readRDS("BGM_raw_expr.rds")
markers <- c("CXCL13", "TNFRSF9")
sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                feature_names=markers,
                                nrounds=0)

enhanced.plots <- purrr::map(markers, function(x) featurePlot(sce.enhanced, x))
spot.plots <- purrr::map(markers, function(x) featurePlot(sce, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=2)

pdf("bayes_expr.pdf",width = 5.5,height = 5)
colorscale <- magma(8)
breaks <- c(0.00, 0.50,1.00,1.50,2.00) 
featurePlot(sce.enhanced, "CXCL13",color=NA)+
  scale_fill_gradientn(colors = colorscale, values = scales::rescale(breaks), 
                       limits = c(0,2), oob = scales::squish,
                       na.value = "#FCFDBFFF",
                       name = "Expression",
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black"))
dev.off()


pdf("raw_expr.pdf",width = 5.5,height = 5)
colorscale <- magma(8)
breaks <- c(0.00, 0.50,1.00,1.50,2.00) 
featurePlot(sce, "CXCL13",color=NA)+
  scale_fill_gradientn(colors = colorscale, values = scales::rescale(breaks), 
                       limits = c(0,2), oob = scales::squish,
                       na.value = "#FCFDBFFF",
                       name = "Expression",
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black"))
dev.off()


#-----Extended Data Fig. 5a-----####
library(Seurat)
library(dplyr)

# pre_treatment
pre_CD8_out <- readRDS("pre_CD8_STARTRAC.rds")
pre_CD8_out@pIndex.tran[1:3,]

dat.plot <- as.matrix(subset(pre_CD8_out@pIndex.tran,aid==pre_CD8_out@proj)[,c(-1,-2,-3,-4)])
rownames(dat.plot) <- subset(pre_CD8_out@pIndex.tran,aid==pre_CD8_out@proj)[,4]
dat.plot[is.na(dat.plot)] <- 0

col_fun = colorRamp2(c(0, 0.25, 0.5), c("#F4C9D3", "#F27CA7", "#EA5384"))

pdf("Pre_CD8_pTran.pdf",width = 5.2,height = 5)
Heatmap(dat.plot, 
        name = "pTran",
        col = col_fun,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        show_row_names = TRUE, 
        show_column_names = TRUE,
        rect_gp = gpar(col = "white", lty = 1, lwd = 2)) 
dev.off()


# post_treatment
post_CD8_out <- readRDS("post_CD8_STARTRAC.rds")

post_CD8_out@pIndex.tran[1:3,]

dat.plot <- as.matrix(subset(post_CD8_out@pIndex.tran,aid==post_CD8_out@proj)[,c(-1,-2,-3,-4)])
rownames(dat.plot) <- subset(post_CD8_out@pIndex.tran,aid==post_CD8_out@proj)[,4]
dat.plot[is.na(dat.plot)] <- 0

col_fun = colorRamp2(c(0, 0.25, 0.5), c("#F4C9D3", "#F27CA7", "#EA5384"))

pdf("Post_CD8_pTran.pdf",width = 5.2,height = 5)
Heatmap(dat.plot, 
        name = "pTran",
        col = col_fun,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        show_row_names = TRUE, 
        show_column_names = TRUE,
        rect_gp = gpar(col = "white", lty = 1, lwd = 2)) 
dev.off()


#-----Extended Data Fig. 5c-----####
library(ggthemes)
library(tidydr)

my_cds <- readRDS("my_cds.rds")

pdf("Tem_Tex_Pseudotime.pdf",width = 6,height = 5)
plot_cell_trajectory(my_cds, color_by = "Pseudotime",cell_size = 0.3,show_branch_points=F)+ theme_dr()+
  theme(panel.grid=element_blank(),legend.position = "right",aspect.ratio = 1,
        # axis.text.x =element_text(size=10,color = "black"),
        # axis.text.y =element_text(size=10,color = "black")
  ) + 
  scale_colour_gradient(low = "#FAF5A9", high = "#C52F2D")
dev.off()

Tex_col <- c("CD8_Tem"="#94C54B","CD8_Tex_CXCL13"="#639EBC","CD8_Tex_GZMK"="#E9866D")

pdf("Tem_Tex_main.type.pdf",width = 6,height = 5)
plot_cell_trajectory(my_cds, color_by = "T_main.type",cell_size = 0.3,show_branch_points=F)+ theme_dr()+
  theme(panel.grid=element_blank(),legend.position = "right",aspect.ratio = 1,
        # axis.text.x =element_text(size=10,color = "black"),
        # axis.text.y =element_text(size=10,color = "black")
  ) + 
  scale_color_manual(values=Tex_col)
dev.off()


library(dplyr)
library(tidyverse)

my_cds <- readRDS("Tem_Tex_monocle2_cds.rds")
str(pData(my_cds))

data_df <- t(reducedDimS(my_cds)) %>%
  as.data.frame() %>%
  select_('Component 1' = 1, 'Component 2' = 2) %>%
  rownames_to_column("Cells") %>%
  mutate(pData(my_cds)$State,
         pData(my_cds)$Pseudotime,
         pData(my_cds)$sample,
         pData(my_cds)$T_main.type,
         pData(my_cds)$Treatment)
colnames(data_df) <- c("cells","Component_1","Component_2","State",
                       "Pseudotime","sample","T_main.type","Treatment")

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
x_breaks <- seq(-6, 6, length.out = 11)
y_breaks <- seq(-3, 3, length.out = 11)

pie_data$grid_x <- cut(pie_data$Component_1, breaks = x_breaks, include.lowest = TRUE)
pie_data$grid_y <- cut(pie_data$Component_2, breaks = y_breaks, include.lowest = TRUE)

grid_counts <- table(pie_data$grid_x, pie_data$grid_y, pie_data$T_main.type)
grid_counts <- as.data.frame(grid_counts)
colnames(grid_counts) <- c("grid_x", "grid_y", "T_main.type", "Count")

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
Tex_col <- c("CD8_Tem"="#94C54B","CD8_Tex_CXCL13"="#639EBC","CD8_Tex_GZMK"="#E9866D")

grid_counts$T_main.type <- factor(grid_counts$T_main.type,levels = c("CD8_Tem","CD8_Tex_GZMK","CD8_Tex_CXCL13"))

pdf("Tem_Tex_pie.pdf",width = 6,height = 5)
ggplot(grid_counts, aes(x = "", y = Proportion, fill = T_main.type)) +
  geom_bar(stat = "identity",color = "white", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(grid_y ~ grid_x) +
  scale_fill_manual(values = Tex_col)+
  labs(x = "Grid X",
       y = "Proportion") +
  theme_void() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 1)
dev.off()



#-----Extended Data Fig. 5d-----####
library(ggrepel)
library(data.table)

out <- readRDS("All_T_STARTRAC.rds")
tcr_meta <- tcr_meta[which(tcr_meta$T_main.type %in% c("CD8_Tcm" ,"CD8_Tem",
                                                       "CD8_Tex_GZMK","CD8_Tex_CXCL13","CD8_HSP",
                                                       'CD4_Tfh','CD4_Tn','CD4_Treg_TNFRSF9+',
                                                       'CD4_Treg_TNFRSF9-','CD4_Th17')),]
subT_stat <- tcr_meta %>%
  select(T_main.type,sample,Treatment,Clone_ID, Clone_NUM)

result <- subT_stat %>%
  filter(Clone_NUM >= 2) %>%
  group_by(T_main.type) %>%
  summarize(n_clone = n() ) %>%
  mutate(percentage = n_clone / nrow(subT_stat))


dat.expa <- as.data.table(out@cluster.sig.data)[aid==out@proj,] %>% 
  filter(index == "expa")

T_color <- c("CD8_Tex_CXCL13"="#639EBC","CD8_Tcm"="#D47AAD","CD8_Tem"="#94C54B", "CD8_Tex_GZMK"="#E9866D",
             "CD8_HSP"="#87BEE7",
             "CD4_Treg_TNFRSF9+"="#E9A9AA","CD4_Treg_TNFRSF9-"="#93C1B7","CD4_Tcm"="#C84853","CD4_Th17"="#EFC384","CD4_Tn"="#5B98C7",
             "CD4_Tfh"="#6F983A")

result <- as.data.frame(result)
colnames(result)[1] <- "majorCluster"
data_in <- merge(dat.expa,result,by="majorCluster")


pdf("T_TCR_STARTRAC.pdf",width = 5,height = 5)
ggplot(data_in,aes(x=percentage,y=value,color=majorCluster))+
  geom_point(size=6) +
  labs(x = "Frequency of proliferating cells" , y = "STARTRAC-expa") +
  theme_few() +
  theme(legend.position = "none",
        axis.text.x =element_text(size=10,color = "black"), 
        axis.text.y=element_text(size=10,color = "black"),
        aspect.ratio=1,
        panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid")) +
  scale_color_manual(values = T_color) +
  geom_text_repel(aes(label = majorCluster),
                  color = "black", vjust = -0.3, size = 4)
dev.off()




#-----Extended Data Fig. 5e-----####
library(msigdbr)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)

KEGG <- msigdbr(species = "Homo sapiens", subcategory="CP:KEGG") %>% 
  dplyr::select(gs_name, entrez_gene)
KEGG <- as.data.frame(KEGG)

TCR_pathway <- KEGG[which(KEGG$gs_name=="KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY"),]
TCR_pathway <- bitr(TCR_pathway$entrez_gene,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db)
TCR_pathway$pathway <- "TCR_SIGNALING"
TCR_pathway <- TCR_pathway[,c(3,2)]

TCR_KEGG <- GSEA(geneList,TERM2GENE = TCR_pathway,
                 minGSSize = 3,maxGSSize = 1000,
                 pvalueCutoff = 1)

pdf("CD8_Tex_TCR_signaling_KEGG.pdf",width = 7,height = 5)
gseaNb(object = TCR_KEGG,
       geneSetID = 'TCR_SIGNALING',
       #addGene = mygene,
       addPval = T,
       pvalX = 0.55,pvalY = 0.8,
       pCol = 'black',pvalSize = 5,
       pHjust = 0)+
  theme(
    axis.text.x =element_text(size=10,color = "black"),
    axis.text.y =element_text(size=10,color = "black"))
dev.off()


