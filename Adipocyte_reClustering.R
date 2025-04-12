
library(SingleR)
library(xlsx)
library(AUCell)
library(Seurat)
library(pheatmap)
library(cowplot)
library(ggplot2)
library(harmony)
library(RColorBrewer)
library(clustree)
library(dplyr)
library(ezRun)
library(viridis)
library(SCpubr)
library(enrichR)

# Load Seurat object
scIntegrated <- qread("snAT.qs")

# Select adipocytes
snAd <- subset(snData, seurat_clusters %in% c(0,7))
DimPlot(snAd, split.by = "tissue")
DefaultAssay(snAd) <- "RNA"

# Create list of tissue specific Seurat objects
snSampleList <- SplitObject(snAd, split.by = "tissue")

# Perform SCTransform on each object
for (i in 1:length(snSampleList)) {
  snSampleList[[i]] <- SCTransform(snSampleList[[i]], verbose = FALSE)
}

# Integrate the SCT normalized objects
features <- SelectIntegrationFeatures(object.list = snSampleList)

snSampleList <- PrepSCTIntegration(object.list = snSampleList, anchor.features = features)

int.anchors <- FindIntegrationAnchors(object.list = snSampleList, normalization.method = "SCT",
                                      anchor.features = features)
scIntegrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT")


# Check assay
names(scIntegrated@assays)

scIntegrated@active.assay

# Run the standard workflow for visualization and clustering
scIntegrated <- ScaleData(scIntegrated, verbose = FALSE)
scIntegrated <- RunPCA(scIntegrated, verbose = FALSE)
ElbowPlot(scIntegrated, ndims = 50)

scIntegrated <- RunUMAP(scIntegrated, reduction = "pca", dims = 1:20)
scIntegrated <- RunTSNE(scIntegrated, reduction = "pca", dims = 1:20)
scIntegrated <- FindNeighbors(scIntegrated, reduction = "pca", dims = 1:20)

# Select a range of resolutions
resolution.range <- seq(from = 0.1, to = 1, by = 0.1)

# Find clusters using a range of resolutions
scIntegrated <- Seurat::FindClusters(object = scIntegrated, resolution = resolution.range)

# Clustering overview
clustree(scIntegrated)

scIntegrated <- FindClusters(scIntegrated, resolution = 0.6)

DimPlot(scIntegrated, reduction = "umap", label = TRUE)
DimPlot(scIntegrated, reduction = "tsne", label = TRUE)

DimPlot(scIntegrated, reduction = "umap", label = TRUE, group.by = "tissue")
DimPlot(scIntegrated, reduction = "tsne", label = TRUE, group.by = "tissue")

DimPlot(scIntegrated, reduction = "umap", split.by = "tissue", ncol = 2)
DimPlot(scIntegrated, reduction = "tsne", split.by = "tissue", ncol = 2)

######################################################################################################################################

# Load reference Emont et al dataset
hsAT <- readRDS("/srv/gstore/projects/p3608/hsWAT_Atlas/human_adipocytes.rds")

DefaultAssay(scIntegrated) <- "integrated"
DefaultAssay(hsAT) <- "integrated"

# Select common variable features
features <- intersect(hsAT$integrated@var.features,scIntegrated$integrated@var.features)

# Reference mapping via transfer anchors
epiAT.anchors <- FindTransferAnchors(reference = hsAT, query = scIntegrated, 
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
p2 <- DimPlot(scIntegrated, reduction = "ref.umap", group.by = "ident", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

p1 <- DimPlot(scIntegrated, reduction = "umap",group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p2 <- DimPlot(scIntegrated, reduction = "tsne",group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p1 + p2

# Visualize cell type prediction and mapping scores
FeaturePlot(scIntegrated,features = "mapping.score", reduction = "umap")
FeaturePlot(scIntegrated,features = "predicted.celltype.score", reduction = "umap")

VlnPlot(scIntegrated, features = c("mapping.score", "predicted.celltype.score"), pt.size = 0) # Check C13, C16

######################################################################################################################################

# Perform differential gene expression analysis
degs <- FindMarkers(scIntegrated, ident.1 = "Adipo1", ident.2 = "Adipo2", min.pct = 0.25,
                          logfc.threshold = 0.25, return.thresh = 0.01)

sigDEGs <- subset(degs, abs(avg_log2FC) > 0.5 & p_val_adj < 0.05)
sigDEGs <- rownames_to_column(sigDEGs, "gene")
sigDEGs$pct.diff <- sigDEGs$pct.1 - sigDEGs$pct.2

#####################################################################################################################

# Load reference Lazarescu et al datasets
sAdipo <- readRDS("sAdipo_CZI_NatGenet.rds")
vAdipo <- readRDS("vAdipo_CZI_NatGenet.rds")

# Define Adipo2 signature based on top upregulated genes
adipomat_low_iat <- c("ENSG00000163359", "ENSG00000114494", "ENSG00000115159", "ENSG00000152284", "ENSG00000131018", "ENSG00000142945", "ENSG00000142949",
                      "ENSG00000162600", "ENSG00000135903", "ENSG00000134352", "ENSG00000168078", "ENSG00000143457", "ENSG00000101187", "ENSG00000161960")

adipomat_low_sat <- c("ENSG00000133083", "ENSG00000172260", "ENSG00000163359", "ENSG00000154262", "ENSG00000154258", "ENSG00000005893", "ENSG00000162600",
                       "ENSG00000167755", "ENSG00000145362", "ENSG00000161960", "ENSG00000105647", "ENSG00000142949", "ENSG00000162654")

# Create list of IAT and SAT specific Adipo2 signatures
adipomat_low_sign <- list(adipomat_low2, adipomat_low_iat2, adipomat_low_sat2)

# Create modules in reference data sets based on Adipo2 signatures
sAdipo <- AddModuleScore(sAdipo, features = adipomat_low_sat, name = "Adipo2sig", assay = "RNA")
FeaturePlot(sAdipo, features = c("Adipo2sig1"), reduction = "umap", ncol=3, min.cutoff = 0.5)

vAdipo <- AddModuleScore(vAdipo, features = adipomat_low_iat, name = "Adipo2sig", assay = "RNA")
FeaturePlot(vAdipo, features = c("Adipo2sig1"), reduction = "umap", ncol=3, min.cutoff = 0.5)

DimPlot(sAdipo, group.by = "author_cell_type", label = TRUE, label.size = 4)
DimPlot(vAdipo, group.by = "author_cell_type", label = TRUE, label.size = 4)

######################################################################################################################################

# Generate pseudo bulk expression for Adipocytes
pseudoBulk_obese <- AggregateExpression(scIntegrated, group.by = c("sample", "cellType"), assays = "log_10k_norm")
pseudoBulk_obese <- pseudoBulk_obese$log_10k_norm

# Remove genes with no count across all samples
pseudoBulk <- pseudoBulk_obese[rowSums(pseudoBulk_obese) > 0, ]

# Define sample and conditions
samples <- colnames(pseudoBulk)

metadata <- data.frame(row.names = samples,
						"cellType" = "Adipo1")
						
metadata[grep("Adipo2", samples), "cellType"] <- "Adipo2"

# Create DGEList object
y <- edgeR::DGEList(counts=pseudoBulk, group=metadata$cellType)

# Remove genes with low expression in > 50% samples
keep <- edgeR::filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- edgeR::calcNormFactors(y)

# Extract normalized counts
normCounts <- edgeR::cpm(y, log=TRUE)

# Run PCA
pcDat  <- prcomp(t(normCounts), scale. = TRUE)

# PC1 vs PC2
autoplot(pcDat, data=metadata, label=FALSE,
         colour="cellType", 
         shape=FALSE,
         size=4) + theme_classic()
