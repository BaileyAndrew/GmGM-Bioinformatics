# GmGM

This contains the code used to generate the results in the GmGM paper.  To acquire the GmGM package, run `pip install gmgm`.

## Test dataset

* For the million-cell dataset, the data is available here: https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1
* The COIL dataset (only shown in supplementary material) is available here: https://www.cs.columbia.edu/CAVE/software/softlib/coil-20.php
* If users want additional tests, the `scanpy` package allows easy download of many datasets that will already be in a format accepted by our methodology (AnnData): https://scanpy.readthedocs.io/en/stable/api/datasets.html

To test yourself, independently of our pre-existing scripts, load the dataset you wish to test on as an AnnData object; GmGM natively works on such objects:

```python
import GmGM
# load your dataset as an AnnData object; as an example, here's loading the million-cell dataset if it were stored in data/million-cell/eQTLAutoimmune.h5ad
adata = ad.read_h5ad("data/million-cell/eQTLAutoimmune.h5ad", backed=True, as_sparse=True)

# put your data in the GmGM 'Dataset' format
# this is optional, GmGM will run directly on adata,
# but it will label axes as "obs" and "var" instead of "cell" and "gene"
# as when running in adata it will not know what each axis represents
dataset = GmGM.Dataset(
    dataset={"data": adata.X[:, :]},
    structure={"data": ("cell", "gene")}
)

# Run GmGM - here are some parameters you may play with
# Parameters with '##' aren't strictly needed to run, but
# as in this example we're working with such a large dataset,
# you'll likely want them (especially a low `n_comps`, which
# governs the number of principal components/eigenvectors to keep.
output = GmGM.GmGM(
    dataset,
    threshold_method="bonferroni",
    to_keep={"cell": 0, "gene": 0.05},
    min_edges={"cell": 0, "gene": 0}, ##
    n_comps=50, ##
    verbose=True, ##
    centering_method=None, ##
    use_nonparanormal_skeptic=True, ##
    nonparanormal_evec_backend="COCA", ##
    readonly=False, ##
    calculate_explained_variance=True, ##
    batch_size=1000 ##
)

# This contains your gene network
output.precision_matrices['gene']
```

## lncRNA Experiments

Source data is available on request; to run GmGM on it, put source data in a folder "data/no_harmony" and run `gmgm_pipeline.py`.  
This will generate several files containing information on the data, as well as an `adata` file that can be used for further study.  
Most of the output files (not the `adata` file due to size) are available already in the `results/gmgm_pathway` folder.  
To generate hdWGCNA results, run `hdwgcna.rmd`; again, these are available in `results/gmgm_pathway` already.

The `results/gmgm_pathway` folder contains a module-by-module breakdown of the results of GmGM and hdWGCNA; **this is the source of Tables 9 and 10 in the paper.**
Figures 3, 4, and 8 are all generated in `for_gmgm_paper.ipynb`, as well as Table 1 (in the guise of the last figure of the file).

## Other experiments

Figures 2 and 6 are generated in `synthetic_data.ipynb`

Figure 7 and the data behind Tables 5 and 6 are generated in `coil.ipynb`.

The data behind Table 4 is generated in `million-cell.ipynb`; the `million-cell` folder already contains a module-by-module breakdown of the results of this file.  This requires a lot of memory to run.
