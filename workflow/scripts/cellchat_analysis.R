library(CellChat)
# 1. Create CellChat object from Seurat
seurat_obj <- readRDS("/home/johan/output/skin_pmh_harmony_sctransform2/TN.combined_dim30.rds")
cellchat <- createCellChat(object = seurat_obj, group.by = "SingleR_label")

# 2. Set database to Human Secreted Signaling
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# 3. Compute probabilities
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 4. Save Network Plots (Circle plots for Fibroblasts -> Macrophages)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, title.name = "Number of interactions")