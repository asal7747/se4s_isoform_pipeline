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

fname = 'data/short_shallow.h5ad'

short_read_shallow = get_anndata(fname)

# explore Anndata object
print(short_read_shallow.obs.head())
print(short_read_shallow.var.head())
print(short_read_shallow.X.shape)
print(short_read_shallow.obs_names[:5])
print(short_read_shallow.obs['Unnamed: 0'][:5])
print(short_read_shallow.var_names[:5])
print(short_read_shallow.obs['nCount_RNA'].describe())
print(short_read_shallow.var.x)
print(short_read_shallow.obs.index)
print(short_read_shallow.var.index)

df = pd.DataFrame(short_read_shallow.X.todense(),
                    index=short_read_shallow.obs['Unnamed: 0'],
                    columns=short_read_shallow.var.x)

print(df.head())

# make obs names based on values from column 1
print("Original adata.obs:")
print(short_read_shallow.obs.head())
print("\nOriginal adata.obs_names:")
print(short_read_shallow.obs_names[:5])

short_read_shallow.obs_names = short_read_shallow.obs['Unnamed: 0'].values
del short_read_shallow.obs['Unnamed: 0']
print("\nUpdated adata.obs_names:")
print(short_read_shallow.obs_names[:5])
print("\nUpdated adata.obs:")
print(short_read_shallow.obs.head())

# get list of obs names
obs_names_list = short_read_shallow.obs_names.tolist()

# check for uniqueness
are_obs_names_unique = short_read_shallow.obs_names.is_unique
print(f"\nAre obs names unique? {are_obs_names_unique}")

# make var_names based on values from column 'x'
short_read_shallow.var_names = short_read_shallow.var['x'].values
# rename x to transcript_id to avoid confusion
short_read_shallow.var.rename({'x': 'transcript_id'}, axis=1, inplace=True)
print(short_read_shallow.var_names[:5])
print(short_read_shallow.var.head())
are_var_names_unique = short_read_shallow.var_names.is_unique
print(f"\nAre var names unique? {are_var_names_unique}")
print(short_read_shallow.var)

# convert to dataframe to make obs_names each cell and var_names each gene
df = pd.DataFrame(short_read_shallow.X.todense(),
                    index=short_read_shallow.obs_names,
                    columns=short_read_shallow.var_names)
print(df.head())
print(df.shape)

# print highly expressed genes
print(sc.pl.highest_expr_genes(short_read_shallow, n_top=20))

# perform some basic filtering
# this returns a tuple of (cell_qc_dataframe, gene_qc_dataframe)
qc = sc.pp.calculate_qc_metrics(short_read_shallow)
                                
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
plt.ylim(0,1250)
plt.title('Library size (total counts) per cell')
plt.show()

# zoom in to the cells with fewer reads
plt.hist(cell_qc_dataframe['total_counts'], bins=1000)
plt.xlabel('Total counts')
plt.ylabel('N cells')
plt.xlim(0,20000)
plt.ylim(0,1250)
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
print('Started with: \n', short_read_shallow)
sc.pp.filter_cells(short_read_shallow, min_counts = 1000)
print('Finished with: \n', short_read_shallow)

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

print('Started with: \n', short_read_shallow)
sc.pp.filter_cells(short_read_shallow, min_genes = 750)
print('Finished with: \n', short_read_shallow)

# next, remove genes whose expression level is considered "undetectable"
# we define a gene as detectable if at least two cells contain more than 5 reads from the gene
# but the threshold strongly depends on the sequencing depth
# filter genes
# e.g. set the threshold to be at least 2 cells with at least 5 reads

print('Started with: \n', short_read_shallow)
sc.pp.filter_genes(short_read_shallow, min_cells = 2)
sc.pp.filter_genes(short_read_shallow, min_counts = 5)
print('Finished with: \n', short_read_shallow)

# find transcript with MT after in transcript ID
mt_transcripts = short_read_shallow.var_names.str.contains('mt-')
print(f'Number of MT transcripts: {mt_transcripts.sum()}')
print('MT transcripts:')
print(short_read_shallow.var_names[mt_transcripts])

short_read_shallow.var['mt'] = mt_transcripts
# calculate qc metrics for mitochondrial genes
sc.pp.calculate_qc_metrics(short_read_shallow, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# make a violin plot of some of the computed quality measures:

# the number of genes expressed in the count matrix
# the total counts per cell
# the percentage of counts in mitochondrial genes

sc.pl.violin(short_read_shallow, 
             ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# or, visualize mitochondrial counts as scatterplots 
sc.pl.scatter(short_read_shallow, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(short_read_shallow, x='total_counts', y='n_genes_by_counts')

# remove cells that have too many mitochondrial genes expressed 
# and cells that have too many total counts
print('Started with: \n', short_read_shallow)
short_read_shallow = short_read_shallow[short_read_shallow.obs.n_genes_by_counts < 200000, :]
print('After removing cells with too many total counts: \n', short_read_shallow)

short_read_shallow = short_read_shallow[short_read_shallow.obs.pct_counts_mt < 20, :]
print('After removing cells with too many mitochondrial genes: \n', short_read_shallow)

# print the final dimensions of the QC'd dataset
print(short_read_shallow) 

# save the filtered anndata object
short_read_shallow.write('qc/short_shallow_qc.h5ad')

short_read_shallow = sc.read('qc/short_shallow_qc.h5ad')
# check values of count matrix
print(short_read_shallow.X[:20,:20])

# normalize (i.e. library-size correct) the data matrix
# so that counts become comparable among cells
sc.pp.normalize_total(short_read_shallow, target_sum=1e4)
print(short_read_shallow.X[:20,:20])

# log transform the data
sc.pp.log1p(short_read_shallow)
print(short_read_shallow.X[:20,:20])

# identify highly variable genes
sc.pp.highly_variable_genes(short_read_shallow, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(short_read_shallow)

# this simply freezes the state of the AnnData object
short_read_shallow.raw = short_read_shallow
print(short_read_shallow)

# save normalized and transformed anndata object
short_read_shallow.write('qc/short_shallow_norm_log_transformed.h5ad')

short_read_shallow = sc.read('qc/short_shallow_norm_log_transformed.h5ad')
# reduce the dimensionality of the data by PCA
sc.tl.pca(short_read_shallow, svd_solver='arpack')
# plot first 2 principal components
sc.pl.pca(short_read_shallow, annotate_var_explained=True)

# inspect the contribution of single PCs to the total variance in the data
sc.pl.pca_variance_ratio(short_read_shallow, log=True, n_pcs=50)

# make t-SNE plot
sc.tl.tsne(short_read_shallow, perplexity=20, learning_rate=1000,
           random_state=42, n_pcs=40)
sc.pl.tsne(short_read_shallow)

# try UMAP plot
sc.pp.neighbors(short_read_shallow, n_neighbors=10, n_pcs=50)
sc.tl.umap(short_read_shallow)
sc.pl.umap(short_read_shallow)

# check what authors did with the short)read_shallow dataset