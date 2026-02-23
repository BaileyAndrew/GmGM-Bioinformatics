import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from GmGM import GmGM
from gprofiler import GProfiler
import warnings
import os
import shutil
from utilities import plot_graph
import matplotlib.pyplot as plt
from scipy import sparse
import warnings

results_filepath = "./results/gmgm_pathway"
data_filepath = "./data/no_harmony"

if not os.path.exists(f"{results_filepath}"):
    os.makedirs(results_filepath)


print("Loading All Data")
datas = []
names = []
for path in os.listdir(f"{data_filepath}"):
    name, extension = os.path.splitext(os.path.basename(path))
    if extension != ".csv":
        continue

    print(f"\tLoading {path}...")

    cur_data = pd.read_csv(f"{data_filepath}/{path}", index_col=" ").drop("Unnamed: 0", axis=1)
    names += [name] * cur_data.shape[1]
    datas.append(cur_data)

data = pd.concat(datas, axis=1)
print(f"Data loaded, shape: {data.shape[0]} genes and {data.shape[1]} cells")
print(f"Removing genes occuring in less than ten cells")

appearances = (data != 0).sum(axis=1)
subset_data = data[(appearances > 10)]

adata = ad.AnnData(subset_data.T)
adata.obs["source"] = names
adata.obs["source"] = adata.obs["source"].astype("category")

print(f"Final data shape: {adata.shape[1]} genes and {adata.shape[0]} cells")
print("Annotating genes, needs internet connection to connect to BioMart...")

# We're only annotiating which genes are highly variable with this line of code
# We are NOT removing the non-highly-variable ones
sc.pp.highly_variable_genes(adata, n_top_genes=5000)

mapping = sc.queries.biomart_annotations(
    "hsapiens",
    [
        "ensembl_gene_id", "external_gene_name", "gene_biotype",
        "chromosome_name", "start_position", "end_position"
    ]
)
mapping = mapping.set_index("ensembl_gene_id")
mapping = mapping[mapping.index.isin(adata.var_names)]

# Add to adata
adata.var["external_gene_name"] = adata.var_names
adata.var["gene_biotype"] = "unknown_biotype"
adata.var["raw_chromosome_name"] = "unknown_chromosome"
adata.var["start_position"] = -1000
adata.var["end_position"] = -1000
adata.var.loc[adata.var_names.isin(mapping.index), "external_gene_name"] = mapping["external_gene_name"]
adata.var.loc[
    adata.var["external_gene_name"].isna(),
    "external_gene_name"
] = adata[:, adata.var["external_gene_name"].isna()].var_names
adata.var.loc[adata.var_names.isin(mapping.index), "gene_biotype"] = mapping["gene_biotype"]
adata.var.loc[adata.var_names.isin(mapping.index), "raw_chromosome_name"] = mapping["chromosome_name"]
adata.var.loc[adata.var_names.isin(mapping.index), "start_position"] = mapping["start_position"]
adata.var.loc[adata.var_names.isin(mapping.index), "end_position"] = mapping["end_position"]

# Fix chromosome names to not include the GL000194.2-style unlocalized regions
adata.var["chromosome_name"] = adata.var["raw_chromosome_name"]
adata.var.loc[
    ~adata.var["chromosome_name"].isin([f"{x}" for x in np.arange(1, 24)]+["MT", "X", "Y"]),
    "chromosome_name"
] = "unknown_chromosome"

print("Running GmGM - may take up to a couple minutes")
output = GmGM(
    adata,
    to_keep={"obs": 5, "var": 5},
    random_state=0,
    n_comps=100,
    threshold_method="rowwise-col-weighted",
    verbose=True,
    use_nonparanormal_skeptic=True,
    tol=1e-10
)

print("Creating and analyzing gene clusters - will need internet access to query GProfiler")
bdata = output.T
bdata.obsp["var_gmgm_connectivities"] = output.varp["var_gmgm_connectivities"].toarray()
sc.tl.leiden(bdata, resolution=3, neighbors_key="var_neighbors_gmgm", key_added="leiden_gmgm")
bdata.obs["leiden_gmgm"].unique()

# This cell will download a list of GO terms associated with each cluster.
# It sometimes fails due to server issues, so we can just reload the saved CSV below if that happens.
def create_pathway_df(
    adata: ad.AnnData,
    graph_type: str,
    use_lncs: bool = True
) -> pd.DataFrame:

    leiden_string = f"leiden_{graph_type}"

    gp = GProfiler(return_dataframe=True)
    pathway_df = gp.profile(
        organism="hsapiens",
        query={
            f"Cluster {module}": adata[
                (adata.obs[leiden_string] == module)
                & (use_lncs | (adata.obs["gene_biotype"] == "protein_coding"))
            ].obs_names.values.tolist()
            for module in adata.obs[leiden_string].unique()
        },
        background=adata.obs_names.tolist(),
    )

    return pathway_df

pathway_df_gmgm = create_pathway_df(bdata, "gmgm")
pathway_df_gmgm_no_lnc = create_pathway_df(bdata, "gmgm", use_lncs=False)

print("Analyzing clusters individually...")

# Overwrite old results, if they exist
if os.path.exists(f"{results_filepath}/gmgm"):
    shutil.rmtree(f"{results_filepath}/gmgm")

# Generate new results
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    for module in bdata.obs["leiden_gmgm"].unique():
        if module == "Singleton":
            continue

        (fig, _), _ = plot_graph(
            bdata,
            module,
            key="var_gmgm",
            pathway_df=pathway_df_gmgm,
            pathway_df_no_lnc=pathway_df_gmgm_no_lnc,
            graph_type="gmgm",
            save=True,
            save_path=f"{results_filepath}",
            vertex_names="external_gene_name",
            top_genes=25,
            top_pathways=20,
            display_centrality="degree",
        )

        # Close figures to prevent notebook slowing down.
        plt.close(fig)

# Save to gmgm directory
pathway_df_gmgm.to_csv(f"{results_filepath}/gmgm/gene_pathways.csv")
pathway_df_gmgm_no_lnc.to_csv(f"{results_filepath}/gmgm/gene_pathways_no_lnc.csv")

print("Mapping gene clusters to cell data...")

# Add gene clusters back to original dataframe
output.var["leiden_gmgm"] = bdata.obs["leiden_gmgm"]

# Annoyingly AnnData (well, the version I'm using) requires sparse "matrices" rather than "arrays"
# I think this has been fixed in more modern versions
output.obsp["obs_gmgm_connectivities"] = sparse.csr_matrix(output.obsp["obs_gmgm_connectivities"])
output.varp["var_gmgm_connectivities"] = sparse.csr_matrix(output.varp["var_gmgm_connectivities"])

# Get lnc percentage for each module, and for each module draw a umap
sc.pp.pca(output, n_comps=50)
sc.pp.neighbors(output, n_neighbors=10)
sc.tl.umap(output)

# Overall chromosome breakdown
fig, ax = plt.subplots(figsize=(6, 6))
overall_chrome_df = output.var["chromosome_name"].value_counts().sort_index()
overall_chrome_df /= overall_chrome_df.sum()
overall_chrome_df.plot(kind='bar', ax=ax)
variance = overall_chrome_df.var()
ax.set_title(f"Whole Dataset Chromosome Breakdown (Variance = {variance:.2f})")
ax.set_xlabel("Chromosome")
ax.set_ylabel("Frequency")
fig.savefig(f"{results_filepath}/chromosome_breakdown.png", bbox_inches="tight")

output.var["lnc_percentage"] = 0
lnc_percentages = {}
morans_i_mean = {}
morans_i_percent = {}
max_chromosome_overexpression = {}
max_chromosome = {}
for module in output.var["leiden_gmgm"].unique():
    with warnings.catch_warnings():
        # scanpy's plotting causes lots of needless warnings here on my environment
        warnings.simplefilter("ignore")
        mask = output.var["leiden_gmgm"] == module
        lnc_percentage = (
            output[:, mask].var["gene_biotype"].value_counts().get("lncRNA", 0)
            / (mask).sum()
        )
        output.var.loc[mask, "lnc_percentage"] = lnc_percentage
        lnc_percentages[module] = lnc_percentage

        # Get info on this module for the cells
        output.obs[f"m{module}_mean"] = output[:, mask].X.mean(axis=1)
        output.obs[f"m{module}_percent"] = output[:, mask].X.sum(axis=1) / output.X.sum(axis=1)
        fig = sc.pl.umap(output, color=["source", f"m{module}_mean"], return_fig=True)
        morans_i = sc.metrics.morans_i(output, vals=output.obs[f"m{module}_mean"])
        morans_i_mean[module] = morans_i
        fig.suptitle(f"Moran's I: {morans_i:.2f}")
        fig.savefig(f"{results_filepath}/gmgm/m{module}/expression.png")
        fig = sc.pl.umap(output, color=["source", f"m{module}_percent"], return_fig=True)
        morans_i = sc.metrics.morans_i(output, vals=output.obs[f"m{module}_percent"])
        morans_i_percent[module] = morans_i
        fig.suptitle(f"Moran's I: {morans_i:.2f}")
        fig.savefig(f"{results_filepath}/gmgm/m{module}/expression_percent.png")

        # Plot module chromosome breakdown:
        fig, ax = plt.subplots(figsize=(6, 6))
        chrome_df = output[:, mask].var["chromosome_name"].value_counts().sort_index()
        chrome_df /= chrome_df.sum()
        chrome_df.plot(kind='bar', ax=ax)
        ax.set_title(f"Module m{module} Chromosome Breakdown")
        ax.set_xlabel("Chromosome")
        ax.set_ylabel("Frequency")
        fig.savefig(f"{results_filepath}/gmgm/m{module}/chromosome_breakdown.png", bbox_inches="tight")

        fig, ax = plt.subplots(figsize=(6, 6))
        chrome_rel = chrome_df / overall_chrome_df
        chrome_rel.plot(kind='bar', ax=ax)
        max_change = chrome_rel.max()
        max_chromosome_overexpression[module] = max_change
        max_chromosome[module] = chrome_rel.idxmax()
        ax.set_title(f"Module m{module} Chromosome Breakdown Relative to Baseline (Max Overexpression = {max_change:.2f})")
        ax.set_xlabel("Chromosome")
        ax.set_ylabel("Frequency")
        fig.savefig(f"{results_filepath}/gmgm/m{module}/chromosome_breakdown_relative.png", bbox_inches="tight")


print("Saving Moran's I and Chromosome Data...")

# Create CSV for Moran's I
morans_df = pd.DataFrame({
    "morans_i_mean": morans_i_mean,
    "morans_i_percent": morans_i_percent
})

# The module names (keys) will automatically become the index
morans_df.index.name = "module"
morans_df.index = "m" + morans_df.index
morans_df = morans_df.sort_values(by="morans_i_percent", ascending=False)
morans_df.to_csv(f"{results_filepath}/gmgm/all_gene_module_homogeneities.csv")

# Create CSV for Chromosome presence
chrome_overexp = pd.DataFrame({
    "max_chromosome_overexpression": max_chromosome_overexpression,
    "most_overexpressed_chromosome": max_chromosome
})

# The module names (keys) will automatically become the index
chrome_overexp.index.name = "module"
chrome_overexp.index = chrome_overexp.index
chrome_overexp = chrome_overexp.sort_values(by="max_chromosome_overexpression", ascending=False)
chrome_overexp.to_csv(f"{results_filepath}/gmgm/all_gene_module_chromosome_overexpression.csv")

print("Investigating LNCs and heterogeneity...")

# Create CSV for LNC presence
lnc_df = pd.DataFrame({
    "lnc_percentage": lnc_percentages
})
# The module names (keys) will automatically become the index
lnc_df.index.name = "module"
lnc_df.index = "m" + lnc_df.index
lnc_df = lnc_df.sort_values(by="lnc_percentage", ascending=False)
lnc_df.to_csv(f"{results_filepath}/gmgm/lnc_percentage_per_module.csv")

# Draw lnc percentage plot
fig, ax = plt.subplots(figsize=(6, 6))
modules = list(lnc_percentages.keys())
values = list(lnc_percentages.values())
ax.boxplot(values, vert=True)
ax.scatter([1] * len(values), values)
# re-zip to ensure ordering is the same
for lnc_percentage, module in zip(values, modules):
    ax.text(1.02, lnc_percentage, f"m{module}", va='center', fontsize=8)

ax.set_xticks([1])
ax.set_xticklabels(["Modules"])
ax.set_ylabel("lncRNA Percentage")
ax.set_title("lncRNA Percentage per Module")

fig.tight_layout()
fig.savefig(f"{results_filepath}/gmgm/lncs_per_module.png")

# Draw Moran's I percentage plot
fig, ax = plt.subplots(figsize=(6, 6))
modules = list(morans_i_percent.keys())
values = list(morans_i_percent.values())
ax.boxplot(values, vert=True)
ax.scatter([1] * len(values), values)
# re-zip to ensure ordering is the same
for mip, module in zip(values, modules):
    ax.text(1.02, mip, f"{module}", va='center', fontsize=8)

ax.set_xticks([1])
ax.set_xticklabels(["Modules"])
ax.set_ylabel("Moran's I")
ax.set_title("Moran's I per Module")

fig.tight_layout()
fig.savefig(f"{results_filepath}/gmgm/morans_i_percents.png")

# Draw the joint plot
fig, ax = plt.subplots(figsize=(6, 6))
modules = list(morans_i_percent.keys())
mip = np.array(list(morans_i_percent.values()))
lnc = np.array([lnc_percentages[m] for m in modules])
chromosome_unmarkers = [max_chromosome_overexpression[m] < 3 for m in modules]
chromosome_markers = [max_chromosome_overexpression[m] >= 3 for m in modules]
ax.scatter(mip[chromosome_markers], lnc[chromosome_markers], marker='*', label='Has Overexpressed Chromosome', color='yellow')
ax.scatter(mip[chromosome_unmarkers], lnc[chromosome_unmarkers], marker='.', label='No Overexpression', color='blue')
ax.set_xlabel("Moran's I")
ax.set_ylabel("Percent LNC")
ax.set_title("LNC Percentage vs Moran's I by Module")
for m, x, y in zip(modules, mip, lnc):
    ax.text(x*1.02, y, f"{m}", va='center', fontsize=8)
ax.legend()
fig.tight_layout()
fig.savefig(f"{results_filepath}/gmgm/lnc_and_morans_i.png")


print("Saving the AnnData file...")

output.write_h5ad(f"{results_filepath}/anndata.h5ad")

print("Done!")