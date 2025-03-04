import scvi
import scanpy as sc
import argparse
import os
import pandas as pd
import anndata as ad

from  _utils.load_cellbender_output  import anndata_from_h5  # customized function to deal with scalar datasets in order to read cellbender output





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', nargs=1, required=True)
    parser.add_argument('-d', nargs=1,  required=True)
    parser.add_argument('-o', nargs=1, required=True)
    args = parser.parse_args()
    #args = parser.parse_args('-d /home/ubuntu/Adipose/Analysis/SCVI_integration  -o /home/ubuntu/Adipose/Analysis/SCVI_integration -f adata_integrated_soupXoutput.h5ad'.split()) #test

    # load model 
    integrated_soupx_file=os.path.join(args.d[0], args.f[0])
    vae_model_folder=os.path.join(args.d[0], "integrated_scvi_model")

    adata_concat = ad.read_h5ad(integrated_soupx_file)
    try:
        vae = scvi.model.SCVI.load(vae_model_folder, adata=adata_concat, use_gpu=True)
    except:
        vae = scvi.model.SCVI.load(vae_model_folder, adata=adata_concat)


    sc.pp.neighbors(adata_concat, use_rep='X_scVI')
    sc.tl.umap(adata_concat, min_dist=0.3)

    # clustering using scVI latent representation via Leiden
    sc.tl.leiden(adata_concat, key_added="leiden_scVI", resolution=0.5)


    print(adata_concat)

    ##plot
    adata_concat.layers['log_10k_norm'] = sc.pp.normalize_total(adata_concat, target_sum=1e+4, inplace=False)['X']
    sc.pp.log1p(adata_concat, layer='log_10k_norm')

    # Add scvi generated expression that has batch effect removed
    # adata_concat.layers['scvi_denoised_10k'] = vae.get_normalized_expression(adata_concat, library_size=1e4)
    # Not to add bc it's too large of a dense matrix (can't be sparse from generated data)
    
    de_df = vae.differential_expression( groupby="leiden_scVI")

    print("Save DE genes for each cluster in the .csv \n") 
    output_file=os.path.join(args.o[0], "adata_integrated_soupXoutput_with_infer.h5ad")
    output_de_file=os.path.join(args.o[0], "DE_genes_by_cluster.csv")

    de_df.to_csv(output_de_file)
    adata_concat.write_h5ad(output_file)
