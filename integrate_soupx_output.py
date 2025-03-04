import scvi
import scanpy as sc
import argparse
import os
import pandas as pd
import anndata as ad

from  _utils.load_cellbender_output  import anndata_from_h5  # customized function to deal with scalar datasets in order to read cellbender output



def load_h5(h5_file, flag_file=None, extra_filter_file=None):
    """
       extra_filter_file: cellranger's barcode file filtered_feature_bc_matrix/barcodes.tsv.gz
    """
    if h5_file.endswith('.ad'):
        try:
            adata = sc.read_10x_h5(h5_file, gex_only=True)
        except:
            # In case scanpy can't read cellbender output
            adata = anndata_from_h5(h5_file)
    else:
        adata = sc.read_h5ad(h5_file) # read soupx output saved as .h5ad
    adata.var_names_make_unique()
    extra_filter_df = None
    if extra_filter_file:
        extra_filter_df = pd.read_csv(extra_filter_file, header=None)       
    if flag_file:
        filter_df = pd.read_csv(flag_file, sep=' ', header=0) 
        cells = filter_df[(filter_df['dead_cells'] == False) & (filter_df['scds_DropletType'] == 'singlet')]['Barcode']
        adata_cells = adata[adata.obs_names.isin(cells)]
        return adata_cells[adata_cells.obs_names.isin(extra_filter_df[0])]  if extra_filter_file else adata_cells
    return adata[adata.obs_names.isin(extra_filter_df[0])] if extra_filter_file else adata



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', nargs=1, required=True)
    parser.add_argument('-d', nargs=1,  required=True)
    parser.add_argument('-o', nargs=1, required=True)
    parser.add_argument('-s', nargs='+', required=True)
    parser.add_argument('--filter', nargs=1, default=None)
    parser.add_argument('--xfilter', nargs=2, default=None)
    args = parser.parse_args()
    #args = parser.parse_args('-d /home/ubuntu/Adipose/Analysis/QC_Process_output/  -o SCVI_integration -f soupx_adjusted_counts.h5ad  -s BOM'.split()) #test
    adata_list = {}
    for sample in args.s:
        h5_file=os.path.join(args.d[0], sample, args.f[0])
        flag_file=os.path.join(args.d[0], sample, args.filter[0]) if args.filter else None
        extra_filter_file=os.path.join(args.xfilter[0], sample, args.xfilter[1]) if args.xfilter else None
        adata_list[sample] = load_h5(h5_file, flag_file, extra_filter_file)


    adata_concat = ad.concat(adata_list, join="inner", label="sample")

    sc.pp.filter_genes(adata_concat, min_counts=3)
    adata_concat.raw = adata_concat # keep full dimension safe
    adata_concat.layers['counts'] = adata_concat.X
    sc.pp.highly_variable_genes(
	adata_concat,
	flavor="seurat_v3",
	n_top_genes=2000,
	layer="counts",
	batch_key="sample",
	subset=False)

    scvi.model.SCVI.setup_anndata(adata_concat, layer="counts", batch_key="sample")
    vae = scvi.model.SCVI(adata_concat, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    adata_concat.obsm["X_scVI"] = vae.get_latent_representation()


    print("Save integrated samples into .h5ad file and trained SCVI model into subfolder\n") 
    vae_model_folder=os.path.join(args.o[0], "integrated_scvi_model")
    vae.save(vae_model_folder)
    output_file=os.path.join(args.o[0], "adata_integrated_soupXoutput.h5ad")
    adata_concat.write_h5ad(output_file)
