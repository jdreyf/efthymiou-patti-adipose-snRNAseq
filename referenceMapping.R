
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ezRun)
library(qs)

# Load Seurat object
scIntegrated <- qread("snAT.qs")

# Load reference Emont et al dataset
hsAT <- readRDS("../human_all.rds")

DefaultAssay(scIntegrated) <- "integrated"
DefaultAssay(hsAT) <- "integrated"

# Select common variable features
features <- intersect(hsAT$integrated@var.features,scIntegrated$RNA@var.features) # 1629

# Reference mapping via transfer anchors
epiAT.anchors <- FindTransferAnchors(reference = hsAT, query = scIntegrated, normalization.method = "SCT",
                                     dims = 1:30, reference.reduction = "pca", features = features)

# Label transfer
predictions <- TransferData(anchorset = epiAT.anchors, refdata = hsAT$cell_type,
                            dims = 1:30)
scIntegrated <- AddMetaData(scIntegrated, metadata = predictions)

# Query mapping on to UMAP embeddings of reference data set
hsAT <- RunUMAP(hsAT, dims = 1:30, reduction = "pca", return.model = TRUE)
scIntegrated <- MapQuery(anchorset = epiAT.anchors, reference = hsAT, query = scIntegrated,
                       refdata = list(celltype = "cell_type"), reference.reduction = "pca", reduction.model = "umap")

scIntegrated <- AddMetaData(object = scIntegrated, metadata = MappingScore(epiAT.anchors,ndim=30), col.name = "mapping.score")

p1 <- DimPlot(hsAT, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(scIntegrated, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

# Visualize cell type prediction and mapping scores
FeaturePlot(scIntegrated,features = "mapping.score", reduction = "umap")
FeaturePlot(scIntegrated,features = "predicted.celltype.score", reduction = "umap")

VlnPlot(scIntegrated, features = c("mapping.score", "predicted.celltype.score"), pt.size = 0) # Check C13, C16

