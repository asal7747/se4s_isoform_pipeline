# short-read shallow Split-seq single-cell myoblast 9,000 cell libs
wget -O short_shallow.h5ad https://www.encodeproject.org/files/ENCFF914DEE/@@download/ENCFF914DEE.h5ad

# short-read deep Split-seq single-cell myoblast 1,000 cell lib
wget -O short_deep.h5ad https://www.encodeproject.org/files/ENCFF418TAK/@@download/ENCFF418TAK.h5ad

# long-read Split-seq single-cell myoblast data
# file 1 is a sparse transcript count matrix
wget -O long_transcript.h5ad https://www.encodeproject.org/files/ENCFF513MSS/@@download/ENCFF513MSS.h5ad
# file 2 is a sparse gene count matrix
wget -O long_gene.h5ad https://www.encodeproject.org/files/ENCFF508RNZ/@@download/ENCFF508RNZ.h5ad