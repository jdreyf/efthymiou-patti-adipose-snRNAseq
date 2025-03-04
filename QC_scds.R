args = commandArgs(trailingOnly=TRUE)
# 1nd arg: filtered_feature_bc_matrix.h5 or output_filtered.h5 (output from either from cellranger or cellbender
# 2rd arg: quality control flags on barcodes (doublets, and dead cells)



path_to_h5=args[1]
QC_Results_Dir=args[2]

message("file name is", path_to_h5)
message("QC results dir is ", QC_Results_Dir)



### ----> purpose of the scirpt
# Doublet detection
# Doublets: scds 
### purpose of the scirpt <---|


# Load packages to be used
#library(DoubletFinder)
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(scds)))



# Read cellranger/cellbender filtered_output.h5 file
# Add path to filtered_output.h5 file below
mat=Read10X_h5(filename = path_to_h5)

# As SingleCellExperiment object
sce <- SingleCellExperiment(list(counts=mat))


## Annotate doublet using binary classification based doublet scoring:
sce = bcds(sce, retRes = TRUE, estNdbl=TRUE)

## Annotate doublet using co-expression based doublet scoring:
try({
    sce = cxds(sce, retRes = TRUE, estNdbl=TRUE)
})


### If cxds worked, run hybrid, otherwise use bcds annotations
if ("cxds_score" %in% colnames(colData(sce))) {
        ## Combine both annotations into a hybrid annotation
        sce = cxds_bcds_hybrid(sce, estNdbl=TRUE)
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$hybrid_score, colData(sce)$hybrid_call))
} else {
        print("this pool failed cxds so results are just the bcds calls")
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$bcds_score, colData(sce)$bcds_call))
}


## Doublet scores are now available via colData:
colnames(Doublets) <- c("Barcode","scds_score","scds_DropletType")
Doublets$scds_DropletType <- gsub("FALSE","singlet",Doublets$scds_DropletType) 
Doublets$scds_DropletType <- gsub("TRUE","doublet",Doublets$scds_DropletType)




message("writing output")
write.table(Doublets, paste0(QC_Results_Dir,"/scds_doublets_singlets.tsv"), sep="\t", quote=FALSE, row.names=FALSE)


summary <- as.data.frame(table(Doublets$scds_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
write.table(summary, paste0(QC_Results_Dir,"/scds_doublet_summary.tsv"), sep="\t", quote=TRUE)


