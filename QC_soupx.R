args = commandArgs(trailingOnly=TRUE)
# 1st arg: dir to save the results  
# 2nd arg: filtered_feature_bc_matrix.h5 or output_filtered.h5 (output from either from cellranger or cellbender
# 3rd arg: quality control flags on barcodes (doublets, and dead cells)
# 4th arg: broad cell type markers required by SoupX for contamination estimate and removal.


path_to_out_h5=args[1]
path_to_filtered_h5=args[2]
QC_Results_Dir=args[3]
Marker_gene_file=args[4]

flag_file = if ( args[5]=='NA' ) NULL else args[5]
extra_filter_file = if (args[6] == 'NA' ) NULL else args[6]


cat(
	path_to_out_h5, "\n",
	path_to_filtered_h5, 		"\n",
	QC_Results_Dir,		"\n",
	Marker_gene_file,	"\n",

	flag_file, 			"\n",
	extra_filter_file, 	"\n"
)


message("file name is ", path_to_out_h5, "\t", path_to_filtered_h5 )
message("QC results dir is ", QC_Results_Dir)


## Function to load data
load_data = function (path_to_filtered_h5, flag_file=NULL, extra_filter_file=NULL) {
    #### 
	#### extra_filter_file: cellranger's barcode file filtered_feature_bc_matrix/barcodes.tsv.gz
	####
	cmat=Read10X_h5(filename = path_to_filtered_h5)

    extra_filter_df = NULL 
    if (!is.null(extra_filter_file)) {
        extra_filter_df = read.csv(extra_filter_file, header=FALSE)
	}
    if (!is.null(flag_file)) {
        filter_df = read.csv(flag_file, sep=" ")
        cells = filter_df[(filter_df['dead_cells'] == FALSE) & (filter_df['scds_DropletType'] == 'singlet'), ]['Barcode']
		cmat_cells = cmat[, colnames(cmat) %in% cells$Barcode]
		return_value <-  if (!is.null(extra_filter_file))  cmat_cells[, colnames(cmat_cells) %in% extra_filter_df[, 1] ]  else cmat_cells
		return(return_value)
	}
	return_value <- if (!is.null(extra_filter_file))  cmat[, colnames(cmat) %in% extra_filter_df[, 1] ]  else cmat
	return(return_value)

}




### ----> purpose of the scirpt
# Estimate and remove contamination from the background mRNA
# Decontamination: SoupX 
### purpose of the scirpt <---|


# Load packages to be used
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(SoupX)))
library(reticulate)
use_condaenv("scvi")
library(anndata)




### Function to remove non-expressed marker genes in the specification 
### to prevent SoupX from running into error 
update_nonExpressedGeneList = function(sc, nonExpressedGeneList, useToEst) {
  
  # Remove marker genes from the specification when there is no expression in the data 
  # Apadted from SoupX function soupx@calculateContaminationFraction

  if(!is(sc,'SoupChannel')){
    stop("sc must be a SoupChannel object")
  }
  #Check that you've provided the genes in the right format
  if(!is.list(nonExpressedGeneList))
    stop("nonExpressedGeneList must be a list of sets of genes.  e.g. list(HB = c('HBB','HBA2'))")
  #Check we can estimate
  if(sum(useToEst)==0)
    stop("No cells specified as acceptable for estimation.  useToEst must not be all FALSE")
  #Construct the data.frame to perform inferance on


  df = list()
  for(i in seq_along(nonExpressedGeneList)){
    tgts = nonExpressedGeneList[[i]]
    #Get the soup fraction for this set
    sFrac = sum(sc$soupProfile[tgts,'est'])
    w = rownames(useToEst)[useToEst[,i]]
    if(length(w)>0){
      #Get the counts
      cnts = as.matrix(sc$toc[tgts,w,drop=FALSE])
      df[[i]] = data.frame(row.names=NULL,
                           cells=colnames(cnts),
                           geneSet=i,
                           soupFrac = sFrac,
                           counts=colSums(cnts),
                           stringsAsFactors=FALSE)

	  if ( sum(df[[i]]$counts) == 0 |  sum(df[[i]]$soupFrac) == 0 ) { 
		warning(cat("0 counts or 0 `soupFrac` for gene group: ", names(nonExpressedGeneList)[i], ". Removing the group from the specification. \n") )
		df[[i]] = NULL
	  }

    }
  }
	
nonExpressedGeneList_update = list()

for (i in seq_along(df)) { 
	if (!is.null(df[[i]])) {
		group_name = names(nonExpressedGeneList)[i]
		nonExpressedGeneList_update[[group_name]] = nonExpressedGeneList[[group_name]]
	}
}

  return(nonExpressedGeneList_update)
}




# Read cellranger/cellbender filtered_output.h5 file
# Add path to filtered_output.h5 file below
dmat=Read10X_h5(filename = path_to_out_h5)
cmat=load_data(path_to_filtered_h5, flag_file, extra_filter_file)  

## Filter cells that have 0 UMIs
cmat = cmat[, colSums(cmat) > 0]

cxg_cmat_unadjust= t(cmat)
cxg_cmat_unadjust_ad = AnnData(X=cxg_cmat_unadjust)

message("saving unadjusted counts")
write_h5ad(cxg_cmat_unadjust_ad, filename=paste0(QC_Results_Dir,"/raw_counts.h5ad") )


# Create soupx channel
sc = SoupChannel(dmat, cmat, calcSoupProfile = FALSE)
# Calculate soup profile  # no idea why the automatic SoupChannel() function does not work
soupProf = data.frame(row.names = rownames(cmat), est = rowSums(cmat)/sum(cmat), counts = rowSums(cmat))
sc = setSoupProfile(sc, soupProf)

# Estimate contamination using SoupX
## Need to make sure the genes are in the count matrix, otherwise errs
MarkerGene_tb = read.table(Marker_gene_file, sep=",", comment.char="#")
MarkerGeneList = list()
for (i in 1:dim(MarkerGene_tb)[1]) {
        markers = MarkerGene_tb[i, "V2"]
    MarkerGeneList[[MarkerGene_tb[i, "V1"]]] = c(strsplit(markers, " ")[[1]])
}

print("markergenes")
print(MarkerGeneList)


##
whitelist = c()
for (g in unlist(MarkerGeneList) ) {
	if (g %in% rownames(cmat) ) {
		whitelist = c(whitelist, g)
	}
}

##
ExpressedGeneList = list()
for (n in names(MarkerGeneList)) {
	ExpressedGeneList[[n]] = MarkerGeneList[[n]] [MarkerGeneList[[n]] %in% whitelist]
}


# Estimate non-expressing cells based on the given genes
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = ExpressedGeneList, 
                                      clusters = FALSE)

###

## Added to prevent SoupX running into error 
ExpressedGeneList_update = update_nonExpressedGeneList(sc, ExpressedGeneList, useToEst)

###

print(ExpressedGeneList_update)

# Calculate the contamination 
sc = calculateContaminationFraction(sc, ExpressedGeneList_update, useToEst = useToEst, forceAccept=TRUE)




cmat_decont <- adjustCounts(sc)

### 
# from tutorial of soupx: 
##### https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html
# There are two ways to do this: using the automatic autoEstCont method, 
# or manually providing a list of “non expressed genes”. 
# 
# from Sean's script
##### "2022-10-04_Axiom001_SP_SoupX_Script.Rmd"
# sc = setContaminationFraction(sc, 0.1)  
# do you want to specify the contamination percentage beforehand??????
# ^ No.
# when we use specific genes (marker genes) to estimate contamination fraction, we dont have to specifcy manually a fixed percentage





# %% 
## Prepare to save decontaminated count matrix 
cxg_cmat_decont = t(cmat_decont)
cxg_cmat_ad = AnnData(X=cxg_cmat_decont)




message("writing output")
write_h5ad(cxg_cmat_ad, filename=paste0(QC_Results_Dir,"/soupx_adjusted_counts.h5ad") )

message("writing the marker genes that were used for soupX decontamination")
for (i in seq_along(ExpressedGeneList_update)) {
	content = c(names(ExpressedGeneList_update)[i], ",",  ExpressedGeneList_update[[i]])
	write(paste(content, collapse=" "),
		file=paste0(QC_Results_Dir,"/markergenes_used_for_soupx.txt") 
		, append=TRUE, sep=" " )
}


message("writing the global contamination fraction")
write(sprintf("Estimated global contamination fraction of %0.2f%%",100*exp(coef(sc$fit))), 
	file=paste0(QC_Results_Dir,"/soupx_estimated_global_contamination_percentage.txt") 
	, append=FALSE, sep=" " )
