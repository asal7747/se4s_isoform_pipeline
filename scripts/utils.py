import sys
import os
import scanpy as sc
import pandas as pd
import anndata
from anndata import AnnData
import matplotlib.pyplot as plt

def get_anndata(fname):
    adata = sc.read(fname)
    return adata

# function to load short-readdatasets in /data and make name the same 
# as filename without extension return individual anndata objects in a dictionary
def load_short_read_datasets(data_dir='data'):
    datasets = {}
    # starts with short
    for fname in os.listdir(data_dir):
        if fname.startswith('short') and fname.endswith('.h5ad'):
            name = fname[:-5]  # remove .h5ad extension
            adata = get_anndata(os.path.join(data_dir, fname))
            datasets[name] = adata
    return datasets


def short_read_adata_to_df(adata):
    """Convert an AnnData to a pandas DataFrame using obs_names and var_names.

    This function performs sanity checks on length and optionally
    allows making names unique (commented out — enable if you want automatic
    uniqueness).
    """
    # convert sparse to dense if needed
    data = adata.X.todense() if hasattr(adata.X, 'todense') else adata.X

    n_obs, n_vars = data.shape
    if len(adata.obs_names) != n_obs:
        raise ValueError(f"adata.obs_names length ({len(adata.obs_names)}) != n_obs ({n_obs})")
    if len(adata.var_names) != n_vars:
        raise ValueError(f"adata.var_names length ({len(adata.var_names)}) != n_vars ({n_vars})")

    # optional: force uniqueness if needed
    if not adata.var_names.is_unique:
        adata.var_names_make_unique(join='-')
    if not adata.obs_names.is_unique:
        adata.obs_names_make_unique(join='-')

    return pd.DataFrame(data, index=adata.obs_names, columns=adata.var_names)

# function to change unnamed: 0 column name to cell_id and x in var to transcript_id
# and return anndata with obs_names as cell_id and var_names as transcript_id
def short_read_adata_index(adata):
    adata.obs.rename({'Unnamed: 0': 'cell_id'}, axis=1, inplace=True)
    adata.var.rename({'x': 'transcript_id'}, axis=1, inplace=True)
    adata.obs_names = adata.obs['cell_id'].values
    adata.var_names = adata.var['transcript_id'].values
    return adata


# function to load long read datasets
def load_long_read_datasets(data_dir='data'):
    datasets = {}
    for fname in os.listdir(data_dir):
        if fname.startswith('long') and fname.endswith('.h5ad'):
            name = fname[:-5]  # remove .h5ad extension
            adata = get_anndata(os.path.join(data_dir, fname))
            datasets[name] = adata
    return datasets


# function to append short_read dfs into a single dataframe
def append_short_read_dfs(dfs):
    combined_df = pd.concat(dfs, axis=0)
    return combined_df

# load all short read datasets
short_read_datasets = load_short_read_datasets()
for name, adata in short_read_datasets.items():
    print(f'Dataset: {name}')
    print(adata)
    print('\n')

# update indexes to be cell_id and transcript_id for all datasets using function
for name, adata in short_read_datasets.items():
    short_read_datasets[name] = short_read_adata_index(adata)

# instantiate individual anndata objects
short_shallow = short_read_datasets['short_shallow']
short_shallow_nuc = short_read_datasets['short_shallow_nuc']
short_deep = short_read_datasets['short_deep']
short_deep_nuc = short_read_datasets['short_deep_nuc']
short_shallow_myotube = short_read_datasets['short_shallow_myotube']
short_deep_myotube = short_read_datasets['short_deep_myotube']


# append short read anndata objects into 1
combined_short_read_adata = anndata.concat(list(short_read_datasets.values()), axis=0)
print(f'Combined short read anndata shape: {combined_short_read_adata.shape}')

# print var_names of combined anndata
print(f'Combined short read anndata var names (transcript IDs): {combined_short_read_adata.var_names[:5]}')

# check uniqueness of obs_names and var_names in combined anndata
print(f'Are obs names unique in combined anndata? {combined_short_read_adata.obs_names.is_unique}')
print(f'Are var names unique in combined anndata? {combined_short_read_adata.var_names.is_unique}')


# convert combined anndata to dataframe (after making names unique)
combined_short_read_df = short_read_adata_to_df(combined_short_read_adata)
print(f'Combined short read dataframe shape: {combined_short_read_df.shape}')

# check for uniqueness of cell_id and transcript_id in combined dataframe
are_cell_ids_unique = combined_short_read_df.index.is_unique
are_transcript_ids_unique = combined_short_read_df.columns.is_unique
print(f'Are cell IDs unique in combined dataframe? {are_cell_ids_unique}')
print(f'Are transcript IDs unique in combined dataframe? {are_transcript_ids_unique}')

# load long read anndata
long_read_datasets = load_long_read_datasets()
for name, adata in long_read_datasets.items():
    print(f'Dataset: {name}')
    print(adata)
    print('\n')

# instantiate individual long read anndata objects
long_gene = long_read_datasets['long_gene']
long_transcript = long_read_datasets['long_transcript']
long_myotube_gene = long_read_datasets['long_myotube_gene']
long_myotube_transcript = long_read_datasets['long_myotube_transcript']
long_nuc_gene = long_read_datasets['long_nuc_gene']
long_nuc_transcript = long_read_datasets['long_nuc_transcript']

# datasets ending in _transcript provide transcript-level quantification
# datasets ending in _gene provide gene-level quantification
# we want to examine both separately, so we will make 2 different combined anndata objects

# combine long read anndata objects ending in "transcript" or "gene" into 1
combined_long_read_adata_transcript = anndata.concat([long_transcript, long_myotube_transcript, long_nuc_transcript], axis=0)
combined_long_read_adata_gene = anndata.concat([long_gene, long_myotube_gene, long_nuc_gene], axis=0)
print(f'Combined long read transcript anndata shape: {combined_long_read_adata_transcript.shape}')
print(f'Combined long read gene anndata shape: {combined_long_read_adata_gene.shape}')

# check uniqueness of obs_names and var_names in combined long read anndata
print(f'Are obs names unique in combined long transcript read anndata? {combined_long_read_adata_transcript.obs_names.is_unique}')
print(f'Are var names unique in combined long read transcript anndata? {combined_long_read_adata_transcript.var_names.is_unique}')
# check uniquess in gene combined long read anndata
print(f'Are obs names unique in combined long read gene anndata? {combined_long_read_adata_gene.obs_names.is_unique}')
print(f'Are var names unique in combined long read gene anndata? {combined_long_read_adata_gene.var_names.is_unique}')

# do these all have matching short reads? This is the 464 cell library
# check if barcode IDs in long read anndata are present in combined short read anndata
combined_long_read_adata_gene_names = set(combined_long_read_adata_gene.obs_names)
combined_short_read_barcodes = set(combined_short_read_adata.obs_names)
matching_barcodes = combined_long_read_adata_gene_names.intersection(combined_short_read_barcodes)
print(f'Number of barcodes in long read gene: {len(combined_long_read_adata_gene_names)}')
print(f'Number of barcodes in combined_short_read_adata: {len(combined_short_read_barcodes)}')
print(f'Number of matching barcodes: {len(matching_barcodes)}')

# do all the below for the combined long read data, make function to avoid repetition

# print highly expressed genes in combined anndata
print(sc.pl.highest_expr_genes(combined_short_read_adata, n_top=20))

# perform some basic filtering
# this returns a tuple of (cell_qc_dataframe, gene_qc_dataframe)
qc = sc.pp.calculate_qc_metrics(combined_short_read_adata)
cell_qc_dataframe = qc[0]
gene_qc_dataframe = qc[1]

print('This is the cell quality control dataframe:')
print(cell_qc_dataframe.head(5))

print('\n\n\n\nThis is the gene quality control dataframe:')
print(gene_qc_dataframe.head(5))

# Cell QC
# cells with few reads should be removed
# plot a histogram showing read counts per cell
plt.hist(cell_qc_dataframe['total_counts'], bins=1000)
plt.xlabel('Total counts')
plt.ylabel('N cells')
plt.xlim(0,50000)
plt.ylim(0,2000)
plt.title('Library size (total counts) per cell')
plt.show()

# zoom in to the cells with fewer reads
plt.hist(cell_qc_dataframe['total_counts'], bins=1000)
plt.xlabel('Total counts')
plt.ylabel('N cells')
plt.xlim(0,20000)
plt.ylim(0,2000)
plt.show()

# add a red line showing an approximate cutoff 
read_count_cutoff = 1000
plt.hist(cell_qc_dataframe['total_counts'], bins=1000)
plt.xlabel('Total counts')
plt.ylabel('N cells')
plt.xlim(0,20000)
plt.ylim(0,1250)
# e.g. chosen cutoff: 1000 reads
plt.axvline(read_count_cutoff, color='red')
plt.title('Library size (total counts) per cell with cutoff')
plt.show()

# remove cells with fewer than N number of reads
print('Started with: \n', combined_short_read_adata)
# stringent cut-off
sc.pp.filter_cells(combined_short_read_adata, min_counts = 1000)
print('Finished with: \n', combined_short_read_adata)

# Detected Genes QC
# want to make sure that reads are distributed across the transcriptome
# e.g. count the total number of unique genes detected in each sample
plt.hist(cell_qc_dataframe['n_genes_by_counts'], bins=100)
plt.xlabel('N genes')
plt.ylabel('N cells')
plt.show()

# zoom in to pick a good cutoff
gene_counts_cutoff = 750
plt.hist(cell_qc_dataframe['n_genes_by_counts'], bins=100)
plt.xlabel('N genes')
plt.ylabel('N cells')
plt.xlim(0,10000) 
plt.ylim(0,2000)
plt.axvline(gene_counts_cutoff, color='red')
plt.show()

print('Started with: \n', combined_short_read_adata)
sc.pp.filter_cells(combined_short_read_adata, min_genes = 750)
print('Finished with: \n', combined_short_read_adata)

# next, remove genes whose expression level is considered "undetectable"
# we define a gene as detectable if at least two cells contain more than 5 reads from the gene
# but the threshold strongly depends on the sequencing depth
# filter genes
# e.g. set the threshold to be at least 2 cells with at least 5 reads

print('Started with: \n', combined_short_read_adata)
sc.pp.filter_genes(combined_short_read_adata, min_cells = 2)
sc.pp.filter_genes(combined_short_read_adata, min_counts = 5)
print('Finished with: \n', combined_short_read_adata)

# find transcript with MT after in transcript ID
mt_transcripts = combined_short_read_adata.var_names.str.contains('mt-')
print(f'Number of MT transcripts: {mt_transcripts.sum()}')
print('MT transcripts:')
print(combined_short_read_adata.var_names[mt_transcripts])

combined_short_read_adata.var['mt'] = mt_transcripts
# calculate qc metrics for mitochondrial genes
sc.pp.calculate_qc_metrics(combined_short_read_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# make a violin plot of some of the computed quality measures:

# the number of genes expressed in the count matrix
# the total counts per cell
# the percentage of counts in mitochondrial genes

sc.pl.violin(combined_short_read_adata, 
             ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# or, visualize mitochondrial counts as scatterplots 
sc.pl.scatter(combined_short_read_adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(combined_short_read_adata, x='total_counts', y='n_genes_by_counts')

# remove cells that have too many mitochondrial genes expressed 
# and cells that have too many total counts
print('Started with: \n', combined_short_read_adata)
combined_short_read_adata = combined_short_read_adata[combined_short_read_adata.obs.n_genes_by_counts < 200000, :]
print('After removing cells with too many total counts: \n', combined_short_read_adata)

combined_short_read_adata = combined_short_read_adata[combined_short_read_adata.obs.pct_counts_mt < 20, :]
print('After removing cells with too many mitochondrial genes: \n', combined_short_read_adata)

# print the final dimensions of the QC'd dataset
print(combined_short_read_adata)

# save the filtered anndata object
combined_short_read_adata.write('qc/combined_short_read_qc.h5ad')

combined_short_read_adata = sc.read('qc/combined_short_read_qc.h5ad')
# check values of count matrix
print(combined_short_read_adata.X[:20,:20])

# normalize (i.e. library-size correct) the data matrix
# so that counts become comparable among cells
sc.pp.normalize_total(combined_short_read_adata, target_sum=1e4)
print(combined_short_read_adata.X[:20,:20])

# log transform the data
sc.pp.log1p(combined_short_read_adata)
print(combined_short_read_adata.X[:20,:20])

# identify highly variable genes
sc.pp.highly_variable_genes(combined_short_read_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(combined_short_read_adata)

# this simply freezes the state of the AnnData object
combined_short_read_adata.raw = combined_short_read_adata
print(combined_short_read_adata)

# save normalized and transformed anndata object
combined_short_read_adata.write('qc/combined_short_read_norm_log_transformed.h5ad')

combined_short_read_adata = sc.read('qc/combined_short_read_norm_log_transformed.h5ad')
# reduce the dimensionality of the data by PCA
sc.tl.pca(combined_short_read_adata, svd_solver='arpack')
# plot first 2 principal components
sc.pl.pca(combined_short_read_adata, annotate_var_explained=True, color='SampleType')

# inspect the contribution of single PCs to the total variance in the data
sc.pl.pca_variance_ratio(combined_short_read_adata, log=True, n_pcs=50)

# make t-SNE plot
sc.tl.tsne(combined_short_read_adata, perplexity=20, learning_rate=1000,
           random_state=42, n_pcs=40)
sc.pl.tsne(combined_short_read_adata, color='SampleType')
sc.pl.tsne(combined_short_read_adata, color='CellType')

# try UMAP plot
sc.pp.neighbors(combined_short_read_adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(combined_short_read_adata)
sc.pl.umap(combined_short_read_adata, color="SampleType")

### long read data exploration and QC ###
# print highly expressed genes in combined anndata
# remove merged_bc from var_names if present
print(sc.pl.highest_expr_genes(combined_long_read_adata_gene, n_top=20))

# perform some basic filtering
# this returns a tuple of (cell_qc_dataframe, gene_qc_dataframe)
qc = sc.pp.calculate_qc_metrics(combined_long_read_adata_gene)
cell_qc_dataframe = qc[0]
gene_qc_dataframe = qc[1]

print('This is the cell quality control dataframe:')
print(cell_qc_dataframe.head(5))

print('\n\n\n\nThis is the gene quality control dataframe:')
print(gene_qc_dataframe.head(5))
# find max counts in a cell
print(gene_qc_dataframe['n_cells_by_counts'].max())


# Cell QC
# cells with few reads should be removed
# plot a histogram showing read counts per cell
# auto adjust x and y limits for better visualization by setting xlim and ylim
lim_lower = 0
xlim_upper = cell_qc_dataframe['n_genes_by_counts'].max()
ylim_upper = cell_qc_dataframe.shape[0] * 0.1
plt.hist(cell_qc_dataframe['total_counts'], bins=1000)
plt.xlabel('Total counts')
plt.ylabel('N cells')
plt.xlim(lim_lower, xlim_upper)
plt.ylim(lim_lower,ylim_upper)
plt.title('Library size (total counts) per cell')
plt.show()

# zoom in to the cells with fewer reads
plt.hist(cell_qc_dataframe['total_counts'], bins=1000)
plt.xlabel('Total counts')
plt.ylabel('N cells')
plt.xlim(0,(xlim_upper/2))
plt.ylim(0,ylim_upper)
plt.show()

# add a red line showing an approximate cutoff 
read_count_cutoff = 500
plt.hist(cell_qc_dataframe['total_counts'], bins=1000)
plt.xlabel('Total counts')
plt.ylabel('N cells')
plt.xlim(0,xlim_upper)
plt.ylim(0,ylim_upper)
# e.g. chosen cutoff: 500 reads
plt.axvline(read_count_cutoff, color='red')
plt.title('Library size (total counts) per cell with cutoff')
plt.show()

# remove cells with fewer than N number of reads
print('Started with: \n', combined_short_read_adata)
# stringent cut-off
sc.pp.filter_cells(combined_short_read_adata, min_counts = 1000)
print('Finished with: \n', combined_short_read_adata)

# Detected Genes QC
# want to make sure that reads are distributed across the transcriptome
# e.g. count the total number of unique genes detected in each sample
plt.hist(cell_qc_dataframe['n_genes_by_counts'], bins=100)
plt.xlabel('N genes')
plt.ylabel('N cells')
plt.show()

# zoom in to pick a good cutoff
gene_counts_cutoff = 750
plt.hist(cell_qc_dataframe['n_genes_by_counts'], bins=100)
plt.xlabel('N genes')
plt.ylabel('N cells')
plt.xlim(0,10000) 
plt.ylim(0,2000)
plt.axvline(gene_counts_cutoff, color='red')
plt.show()

print('Started with: \n', combined_short_read_adata)
sc.pp.filter_cells(combined_short_read_adata, min_genes = 750)
print('Finished with: \n', combined_short_read_adata)

# next, remove genes whose expression level is considered "undetectable"
# we define a gene as detectable if at least two cells contain more than 5 reads from the gene
# but the threshold strongly depends on the sequencing depth
# filter genes
# e.g. set the threshold to be at least 2 cells with at least 5 reads

print('Started with: \n', combined_long_read_adata_transcript)
sc.pp.filter_genes(combined_long_read_adata_transcript, min_cells = 2)
sc.pp.filter_genes(combined_long_read_adata_transcript, min_counts = 5)
print('Finished with: \n', combined_long_read_adata_transcript)

# find transcript with MT after in transcript ID
mt_transcripts = combined_long_read_adata_transcript.var_names.str.contains('mt-')
print(f'Number of MT transcripts: {mt_transcripts.sum()}')
print('MT transcripts:')
print(combined_short_read_adata.var_names[mt_transcripts])

combined_long_read_adata_transcript.var['mt'] = mt_transcripts
# calculate qc metrics for mitochondrial genes
sc.pp.calculate_qc_metrics(combined_long_read_adata_transcript, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# make a violin plot of some of the computed quality measures:

# the number of genes expressed in the count matrix
# the total counts per cell
# the percentage of counts in mitochondrial genes

sc.pl.violin(combined_long_read_adata_transcript, 
             ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# or, visualize mitochondrial counts as scatterplots 
sc.pl.scatter(combined_long_read_adata_transcript, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(combined_long_read_adata_transcript, x='total_counts', y='n_genes_by_counts')

# remove cells that have too many mitochondrial genes expressed 
# and cells that have too many total counts
print('Started with: \n', combined_long_read_adata_transcript)
combined_long_read_adata_transcript = combined_long_read_adata_transcript[combined_long_read_adata_transcript.obs.n_genes_by_counts < 200000, :]
print('After removing cells with too many total counts: \n', combined_long_read_adata_transcript)

combined_long_read_adata_transcript = combined_long_read_adata_transcript[combined_long_read_adata_transcript.obs.pct_counts_mt < 20, :]
print('After removing cells with too many mitochondrial genes: \n', combined_long_read_adata_transcript)

# print the final dimensions of the QC'd dataset
print(combined_long_read_adata_transcript)

# save the filtered anndata object
combined_short_read_adata.write('qc/combined_short_read_qc.h5ad')

combined_short_read_adata = sc.read('qc/combined_short_read_qc.h5ad')
# check values of count matrix
print(combined_short_read_adata.X[:20,:20])

# normalize (i.e. library-size correct) the data matrix
# so that counts become comparable among cells
sc.pp.normalize_total(combined_short_read_adata, target_sum=1e4)
print(combined_short_read_adata.X[:20,:20])

# log transform the data
sc.pp.log1p(combined_short_read_adata)
print(combined_short_read_adata.X[:20,:20])

# identify highly variable genes
sc.pp.highly_variable_genes(combined_short_read_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(combined_short_read_adata)

# this simply freezes the state of the AnnData object
combined_short_read_adata.raw = combined_short_read_adata
print(combined_short_read_adata)

# save normalized and transformed anndata object
combined_short_read_adata.write('qc/combined_short_read_norm_log_transformed.h5ad')

combined_short_read_adata = sc.read('qc/combined_short_read_norm_log_transformed.h5ad')
# reduce the dimensionality of the data by PCA
sc.tl.pca(combined_short_read_adata, svd_solver='arpack')
# plot first 2 principal components
sc.pl.pca(combined_short_read_adata, annotate_var_explained=True, color='SampleType')

# inspect the contribution of single PCs to the total variance in the data
sc.pl.pca_variance_ratio(combined_short_read_adata, log=True, n_pcs=50)

# make t-SNE plot
sc.tl.tsne(combined_short_read_adata, perplexity=20, learning_rate=1000,
           random_state=42, n_pcs=40)
sc.pl.tsne(combined_short_read_adata, color='SampleType')
sc.pl.tsne(combined_short_read_adata, color='CellType')

# try UMAP plot
sc.pp.neighbors(combined_short_read_adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(combined_short_read_adata)
sc.pl.umap(combined_short_read_adata, color="SampleType")



# explore Anndata object
# print(short_read_shallow.obs.head())
# print(short_read_shallow.var.head())
# print(short_read_shallow.X.shape)
# print(short_read_shallow.obs_names[:5])
# print(short_read_shallow.obs['Unnamed: 0'][:5])
# print(short_read_shallow.var_names[:5])
# print(short_read_shallow.obs['nCount_RNA'].describe())
# print(short_read_shallow.var.x)
# print(short_read_shallow.obs.index)
# print(short_read_shallow.var.index)

# df = pd.DataFrame(short_read_shallow.X.todense(),
#                     index=short_read_shallow.obs['Unnamed: 0'],
#                     columns=short_read_shallow.var.x)

# print(df.head())
# print(df.shape)
# # make obs names based on values from column 1
# print("Original adata.obs:")
# print(short_read_shallow.obs.head())
# print("\nOriginal adata.obs_names:")
# print(short_read_shallow.obs_names[:5])

# short_read_shallow.obs_names = short_read_shallow.obs['Unnamed: 0'].values
# del short_read_shallow.obs['Unnamed: 0']
# print("\nUpdated adata.obs_names:")
# print(short_read_shallow.obs_names[:5])
# print("\nUpdated adata.obs:")
# print(short_read_shallow.obs.head())

# # get list of obs names
# obs_names_list = short_read_shallow.obs_names.tolist()

# # check for uniqueness
# are_obs_names_unique = short_read_shallow.obs_names.is_unique
# print(f"\nAre obs names unique? {are_obs_names_unique}")


# # make var_names based on values from column 'x'
# short_read_shallow.var_names = short_read_shallow.var['x'].values
# # rename x to transcript_id to avoid confusion
# short_read_shallow.var.rename({'x': 'transcript_id'}, axis=1, inplace=True)
# print(short_read_shallow.var_names[:5])
# print(short_read_shallow.var.head())
# are_var_names_unique = short_read_shallow.var_names.is_unique
# print(f"\nAre var names unique? {are_var_names_unique}")
# print(short_read_shallow.var)

# # convert to dataframe to make obs_names each cell and var_names each gene
# df = pd.DataFrame(short_read_shallow.X.todense(),
#                     index=short_read_shallow.obs_names,
#                     columns=short_read_shallow.var_names)
# print(df.head())
# print(df.shape)

# print highly expressed genes
# print(sc.pl.highest_expr_genes(short_read_shallow, n_top=20))

# perform some basic filtering
# this returns a tuple of (cell_qc_dataframe, gene_qc_dataframe)
# qc = sc.pp.calculate_qc_metrics(short_read_shallow)
                                
# cell_qc_dataframe = qc[0]
# gene_qc_dataframe = qc[1]

# print('This is the cell quality control dataframe:')
# print(cell_qc_dataframe.head(5))

# print('\n\n\n\nThis is the gene quality control dataframe:')
# print(gene_qc_dataframe.head(5))

# Cell QC
# cells with few reads should be removed
# plot a histogram showing read counts per cell
# plt.hist(cell_qc_dataframe['total_counts'], bins=1000)
# plt.xlabel('Total counts')
# plt.ylabel('N cells')
# plt.xlim(0,50000)
# plt.ylim(0,1250)
# plt.title('Library size (total counts) per cell')
# plt.show()

# # zoom in to the cells with fewer reads
# plt.hist(cell_qc_dataframe['total_counts'], bins=1000)
# plt.xlabel('Total counts')
# plt.ylabel('N cells')
# plt.xlim(0,20000)
# plt.ylim(0,1250)
# plt.show()

# # add a red line showing an approximate cutoff 
# read_count_cutoff = 1000
# plt.hist(cell_qc_dataframe['total_counts'], bins=1000)
# plt.xlabel('Total counts')
# plt.ylabel('N cells')
# plt.xlim(0,20000)
# plt.ylim(0,1250)
# # e.g. chosen cutoff: 1000 reads
# plt.axvline(read_count_cutoff, color='red')
# plt.title('Library size (total counts) per cell with cutoff')
# plt.show()

# # remove cells with fewer than N number of reads
# print('Started with: \n', short_read_shallow)
# sc.pp.filter_cells(short_read_shallow, min_counts = 1000)
# print('Finished with: \n', short_read_shallow)

# # Detected Genes QC
# # want to make sure that reads are distributed across the transcriptome
# # e.g. count the total number of unique genes detected in each sample

# plt.hist(cell_qc_dataframe['n_genes_by_counts'], bins=100)
# plt.xlabel('N genes')
# plt.ylabel('N cells')
# plt.show()

# # zoom in to pick a good cutoff
# gene_counts_cutoff = 750
# plt.hist(cell_qc_dataframe['n_genes_by_counts'], bins=100)
# plt.xlabel('N genes')
# plt.ylabel('N cells')
# plt.xlim(0,10000) 
# plt.ylim(0,2000)
# plt.axvline(gene_counts_cutoff, color='red')
# plt.show()

# print('Started with: \n', short_read_shallow)
# sc.pp.filter_cells(short_read_shallow, min_genes = 750)
# print('Finished with: \n', short_read_shallow)

# # next, remove genes whose expression level is considered "undetectable"
# # we define a gene as detectable if at least two cells contain more than 5 reads from the gene
# # but the threshold strongly depends on the sequencing depth
# # filter genes
# # e.g. set the threshold to be at least 2 cells with at least 5 reads

# print('Started with: \n', short_read_shallow)
# sc.pp.filter_genes(short_read_shallow, min_cells = 2)
# sc.pp.filter_genes(short_read_shallow, min_counts = 5)
# print('Finished with: \n', short_read_shallow)

# # find transcript with MT after in transcript ID
# mt_transcripts = short_read_shallow.var_names.str.contains('mt-')
# print(f'Number of MT transcripts: {mt_transcripts.sum()}')
# print('MT transcripts:')
# print(short_read_shallow.var_names[mt_transcripts])

# short_read_shallow.var['mt'] = mt_transcripts
# # calculate qc metrics for mitochondrial genes
# sc.pp.calculate_qc_metrics(short_read_shallow, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# # make a violin plot of some of the computed quality measures:

# # the number of genes expressed in the count matrix
# # the total counts per cell
# # the percentage of counts in mitochondrial genes

# sc.pl.violin(short_read_shallow, 
#              ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#              jitter=0.4, multi_panel=True)

# # or, visualize mitochondrial counts as scatterplots 
# sc.pl.scatter(short_read_shallow, x='total_counts', y='pct_counts_mt')
# sc.pl.scatter(short_read_shallow, x='total_counts', y='n_genes_by_counts')

# # remove cells that have too many mitochondrial genes expressed 
# # and cells that have too many total counts
# print('Started with: \n', short_read_shallow)
# short_read_shallow = short_read_shallow[short_read_shallow.obs.n_genes_by_counts < 200000, :]
# print('After removing cells with too many total counts: \n', short_read_shallow)

# short_read_shallow = short_read_shallow[short_read_shallow.obs.pct_counts_mt < 20, :]
# print('After removing cells with too many mitochondrial genes: \n', short_read_shallow)

# # print the final dimensions of the QC'd dataset
# print(short_read_shallow) 

# # save the filtered anndata object
# short_read_shallow.write('qc/short_shallow_qc.h5ad')

# short_read_shallow = sc.read('qc/short_shallow_qc.h5ad')
# # check values of count matrix
# print(short_read_shallow.X[:20,:20])

# # normalize (i.e. library-size correct) the data matrix
# # so that counts become comparable among cells
# sc.pp.normalize_total(short_read_shallow, target_sum=1e4)
# print(short_read_shallow.X[:20,:20])

# # log transform the data
# sc.pp.log1p(short_read_shallow)
# print(short_read_shallow.X[:20,:20])

# # identify highly variable genes
# sc.pp.highly_variable_genes(short_read_shallow, min_mean=0.0125, max_mean=3, min_disp=0.5)
# sc.pl.highly_variable_genes(short_read_shallow)

# # this simply freezes the state of the AnnData object
# short_read_shallow.raw = short_read_shallow
# print(short_read_shallow)

# # save normalized and transformed anndata object
# short_read_shallow.write('qc/short_shallow_norm_log_transformed.h5ad')

# short_read_shallow = sc.read('qc/short_shallow_norm_log_transformed.h5ad')
# # reduce the dimensionality of the data by PCA
# sc.tl.pca(short_read_shallow, svd_solver='arpack')
# # plot first 2 principal components
# sc.pl.pca(short_read_shallow, annotate_var_explained=True)

# # inspect the contribution of single PCs to the total variance in the data
# sc.pl.pca_variance_ratio(short_read_shallow, log=True, n_pcs=50)

# # make t-SNE plot
# sc.tl.tsne(short_read_shallow, perplexity=20, learning_rate=1000,
#            random_state=42, n_pcs=40)
# sc.pl.tsne(short_read_shallow)

# # try UMAP plot
# sc.pp.neighbors(short_read_shallow, n_neighbors=10, n_pcs=50)
# sc.tl.umap(short_read_shallow)
# sc.pl.umap(short_read_shallow)