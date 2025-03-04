args = commandArgs(trailingOnly=TRUE)
# 1st arg: filtered_feature_bc_matrix.h5 or output_filtered.h5 (output from either from cellranger or cellbender
# 2nd arg: dir to save the results 
# 3rd arg: mitochondrial transcripts proportion cutoff threshold 

# test >>
path_to_h5="/home/ubuntu/Adipose/Analysis/CellBender_output/RUGBY018_PV_GEX/output_filtered.h5"
QC_Results_Dir="/home/ubuntu/Adipose/Analysis/QC_Process_output/RUGBY018_PV_GEX"
MT_Cutoff=1
## << test 


path_to_h5=args[1]
QC_Results_Dir=args[2]
MT_Cutoff=as.integer(args[3]) # In percent

message("file name is", path_to_h5)
message("QC results dir is ", QC_Results_Dir)


### ----> purpose of the scirpt
# Dead or dying cells: High MT content 
### purpose of the scirpt <---|


# Load packages to be used
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dplyr)))



# Read cellranger/cellbender filtered_output.h5 file
# Add path to filtered_output.h5 file below
mat=Read10X_h5(filename = path_to_h5)
# Create seurat object and perform pre-processing steps.
# No filtering of barcodes based on UMI/gene count is performed at this stage
# No filtering of features is performed at this stage
s.object=CreateSeuratObject(mat)


# Cutoff for MT content is upto the user but you would want to keep it consistent across analysis
df_mt <- PercentageFeatureSet(s.object, pattern = "^MT-")
colnames(df_mt) <- "mt_percentage" 
df_mt["Barcode"] <- rownames(df_mt)
df_mt["dead_cells"] = FALSE
dead_cells_idx = df_mt['mt_percentage'] > (MT_Cutoff) 
dead_cells_idx[which(is.na(dead_cells_idx))] = TRUE  # for 0-count cells, count them as dead as well
df_mt[dead_cells_idx, 'dead_cells' ] = TRUE
df_mt["mt_cutoff"] = MT_Cutoff

mt_df <- df_mt[c("Barcode", "mt_percentage", "dead_cells", "mt_cutoff")]

# Save results
message("writing output")
write.table(mt_df, paste0(QC_Results_Dir,"/mt_dead_alive.tsv"), sep="\t", quote=FALSE, row.names=FALSE)


summary <- as.data.frame(table(mt_df$dead_cells))
colnames(summary) <- c("Classification", "Dead N")
write.table(summary, paste0(QC_Results_Dir,"/mt_summary.tsv"), sep="\t", quote=TRUE)
