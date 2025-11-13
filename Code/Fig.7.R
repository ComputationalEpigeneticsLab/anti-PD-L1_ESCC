# ------------------------------------------------------------------
# Fig. 7
# ------------------------------------------------------------------

#-----Fig. 7a-----####

relative_abundance_matrix <- read.csv("relative_abundance_matrix.txt",row.names = T,col.names = T,sep = '\t',quote = F)

annotation <- read.csv("sample_info.txt",row.names = T,col.names = T,sep = '\t',quote = F)

mycolor_sampletype <- c("N"="#5897A5","P"="#E0B77E","T"="#CD8CAD")
mycolor_treat <- c("Pre_treatment"="#FFB600","Post_treatment"="#77B6BE")
group_color <- c("Pre_PCR"="#A8541E","Pre_nPCR"="#FC8F62","Post_nPCR"="#997DB0","Post_PCR"="#76A3B2")

ann_colors <- list(sample.type = mycolor_sampletype,Treatment = mycolor_treat,Response = group_color)

annotation <- annotation[,-1]
annotation$sample.type <- factor(annotation$sample.type,levels = c("N","P","T"))
annotation$Treatment <- factor(annotation$Treatment,levels = c("Pre_treatment","Post_treatment"))
annotation$Response <- factor(annotation$Response,levels = c("Pre_nPCR","Pre_PCR","Post_nPCR","Post_PCR"))

pdf("cluster.pdf",width = 12,height = 4)
res <- pheatmap(relative_abundance_matrix,
                color = colorRampPalette(c("#053061", "#4393C3", "#92C5DE", "#D1E5F0" ,"#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(50), #颜色的设置
                legend = TRUE, 
                
                #单元格设置
                border_color = NA,
                display_numbers = F, 
                
                #行和列名设置
                fontsize = 10,
                fontsize_row = 10,
                fontsize_col = 10,
                angle_col = 90,
                show_rownames = TRUE,
                show_colnames = T,
                
                #聚类和距离
                scale = "row",
                cluster_cols = TRUE,
                cluster_rows = TRUE,
                cutree_cols = 4,
                treeheight_row = 10,
                treeheight_col = 10,
                clustering_distance_rows = "euclidean",
                clustering_method = "ward.D2",
                
                #分组注释
                annotation_col = annotation,
                annotation_colors = ann_colors,
                annotation_names_col = T,
                annotation_legend = T
)
dev.off()


#-----Fig. 7c-----####
library(pheatmap)
library(circlize)
library(ComplexHeatmap)

gene_cell_exp <- read.csv("expr.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)

gene_use <- c("LAG3","TIGIT","CTLA4","PDCD1")

gene_expr_use <- gene_cell_exp[which(rownames(gene_cell_exp) %in% gene_use),]

post <- gene_expr_use[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29)]

colnames(post) <- c("CD4_Tfh","CD4_Th17","CD4_Tn","CD4_Treg_TNFRSF9-","CD4_Treg_TNFRSF9+","CD8_HSP","CD8_Tcm","CD8_Tem",
                    "CD8_Tex_CXCL13","CD8_Tex_GZMK","MAIT","NK_CD16","NK_XCL1","T_MALAT1","T_Pro")

rownames(post) <- c("CTLA4_post","LAG3_post","PDCD1_post","TIGIT_post")

pre <- gene_expr_use[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)]

colnames(pre) <- c("CD4_Tfh","CD4_Th17","CD4_Tn","CD4_Treg_TNFRSF9-","CD4_Treg_TNFRSF9+","CD8_HSP","CD8_Tcm","CD8_Tem",
                   "CD8_Tex_CXCL13","CD8_Tex_GZMK","MAIT","NK_CD16","NK_XCL1","T_MALAT1","T_Pro")
rownames(pre) <- c("CTLA4_pre","LAG3_pre","PDCD1_pre","TIGIT_pre")

post_pre <- rbind(post,pre)
post_pre <- post_pre[,c('CD4_Tn','CD4_Tfh','CD4_Th17','CD4_Treg_TNFRSF9+','CD4_Treg_TNFRSF9-',
                        'CD8_Tcm','CD8_Tem','CD8_Tex_GZMK','CD8_Tex_CXCL13','CD8_HSP','MAIT',
                        'NK_CD16','NK_XCL1','T_Pro','T_MALAT1')]

post_pre <- post_pre[c("CTLA4_pre","CTLA4_post","LAG3_pre","LAG3_post","PDCD1_pre","PDCD1_post","TIGIT_pre","TIGIT_post"),]

post_pre <- post_pre[c(3,4,7,8,5,6,1,2),]

col_fun = colorRamp2(c(-1, 0, 1), c("#7ACED0", "#FEFCFA", "#F2A15F"))

pdf("immune_checkpoint_expr.pdf",width = 8,height = 4)
pheatmap(as.matrix(post_pre),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames = TRUE,
         legend_breaks = c(-1, 0, 1),
         border_color = "white",
         color = col_fun,
         display_numbers = TRUE, 
         fontsize_number = 8,
         number_format = "%.2f"
)
dev.off()







