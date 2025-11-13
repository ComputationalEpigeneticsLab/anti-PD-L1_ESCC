
library(CellChat)
library(Seurat)
library(patchwork)

seu_obj <- readRDS('ESCC_integrated_all_sub_celltype.rds')

mye_obj <- readRDS('Mye_obj_annotation.rds')
macro_obj<-subset(mye_obj,Mye_main.type %in% c('Macro_APOE','Macro_FOLR2',
                                               'Macro_CSF1R','Macro_HLA-DR+'))

t_obj<-readRDS('T_obj_annotation.rds')

head(seu_obj@meta.data)

head(mye_obj@meta.data)

head(t_obj@meta.data)

t_cell_use<-Cells(t_obj)
t_cell_meta<-t_obj@meta.data[,c('T_main.type','Treatment','Response')]
t_cell_meta$label<-t_cell_meta$T_main.type
t_cell_meta<-t_cell_meta[,c('Treatment','Response','label')]

macro_cell_use<-Cells(macro_obj)
macro_cell_meta<-macro_obj@meta.data[,c('Treatment','Response')]
macro_cell_meta$label<-'Macrophages'

obj_use<-subset(seu_obj,cells=c(t_cell_use,macro_cell_use))

data_use<-as.matrix(GetAssayData(obj_use,slot = 'data',assay='RNA'))
meta_use<-rbind(t_cell_meta,macro_cell_meta)


reticulate::use_condaenv(condaenv = "/software/miniconda3/envs/work", required = TRUE)
cellchat <- createCellChat(object = data_use, meta = meta_use, group.by = "label")
cellchat <- setIdent(cellchat, ident.use = "label")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
mat_count <- cellchat@net$count
saveRDS(cellchat,"cellchat.rds")




