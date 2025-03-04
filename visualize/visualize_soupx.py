###
#### Same code visusalize soupx output
####
#%%
import anndata as ad
import scanpy as sc
import os, sys 
from matplotlib import pyplot as plt, colors 
import numpy as np 
import pandas as pd



de_df_file = sys.argv[1]
soupx_inte_file = sys.argv[2]
figure_dir = sys.argv[3]


# %% 
soupx_adata = ad.read_h5ad(soupx_inte_file)
de_df = pd.read_csv(de_df_file, index_col=0)


soupx_adata.obs['totalUMI'] = soupx_adata.X.sum(axis = 1)
soupx_adata.obs['logUMI'] = sc.pp.log1p(soupx_adata.obs['totalUMI'])



if not os.path.exists(figure_dir):
    os.mkdir(figure_dir)

os.chdir(figure_dir)




# %%

##make red colormap
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colorsComb = np.vstack([colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)
my_cmap = mymap(np.arange(mymap.N))
my_cmap[:,-1] = np.linspace(0, 1, mymap.N)
my_cmap = colors.ListedColormap(my_cmap)


##plot

# %%
# Umap by sample
sc.pl.umap(soupx_adata, color="sample", save="umap_integratedsamples.pdf")

# Umap by UMI
sc.pl.umap(soupx_adata, color='totalUMI', save="umap_integratedsamples_UMI.pdf")
sc.pl.umap(soupx_adata, color='logUMI', save="umap_integratedsamples_logUMI.pdf")


#TODO: take gene list for plot
gene_list = ["PLIN5","PDE3B","ADIPOQ","GHR","PLIN1",
"NR4A1","PDE3B","GHR","PLIN1","LINC00278",
"PDGFRA","DCN","TSHZ2","ABCA10",
"PECAM1","VWF","ADGRL4","ANO2",
"MYH11","RYR2","TAGLN",
"VEGFC","MECOM","PTPRB","FLT1",
"F13A1","CD163",
"MMRN1","PTPRE","NRG3","RELN",
"DCN","FBN1","SEMA3C",
"SKAP1","RIPOR2","CD247",
"CPA3","KIT","MS4A2",
"ADGRB3","COL25A1","STEAP4",
"PKHD1L1","ITLN1","KCTD8","THSD4",
"CNTNAP2","CNTN5","RBFOX1",
"TM4SF19","ALCAM","TPRG1",
"NRXN1","XKR4"]
sc.pl.umap(soupx_adata, color=gene_list, layer = 'log_10k_norm', use_raw=False, save="umap_integratedsamples_genes_scaled.pdf")
sc.pl.umap(soupx_adata, color='leiden_scVI', save="umap_integratedsamples_cluster.pdf")
sc.pl.umap(soupx_adata, color='leiden_scVI', legend_loc='on data', save="umap_integratedsamples_cluster2.pdf")


# %% 
markers = {}
cats = soupx_adata.obs.leiden_scVI.cat.categories
for i, c in enumerate(cats):
    cid = f"{c} vs Rest"
    cell_type_df = de_df.loc[de_df.comparison == cid]

    cell_type_df = cell_type_df[cell_type_df.lfc_mean > 0]

    cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 3]
    cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.1]

    markers[c] = cell_type_df.index.tolist()[:3]


# %%
sc.tl.dendrogram(soupx_adata, groupby="leiden_scVI", use_rep="X_scVI")

sc.pl.dotplot(
    soupx_adata,
    markers,
    groupby="leiden_scVI",
    dendrogram=True,
    color_map="Blues",
    swap_axes=True,
    use_raw=True,
    standard_scale="var",
    save='dendrogram_cluster.pdf'
)


sc.pl.heatmap(
    soupx_adata,
    markers,
    groupby="leiden_scVI",
    layer="log_10k_norm",
    standard_scale="var",
    dendrogram=True,
    save="heatmap_cluster_markers.pdf",
)

