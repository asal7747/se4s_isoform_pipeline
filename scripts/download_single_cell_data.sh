# These commands download 3 sets of data from an ENCODE project:
# 1) short-read shallow, short-read deep, and long-read Split-seq single-cell data
# 2) short-read shallow, short-read deep, and long-read Split-seq single-nucleus data
# 3) short-read shallow, short-read deep, and long-read Split-seq single-cell myotube data

### 1 & 2) myoblast data ###

# short-read shallow Split-seq single-cell myoblast 9,000 cell libs
wget -O short_shallow.h5ad https://www.encodeproject.org/files/ENCFF914DEE/@@download/ENCFF914DEE.h5ad

# short-read shallow Split-seq single-nucleus myoblast 9,000 cell libs
wget -O short_shallow_nuc.h5ad https://www.encodeproject.org/files/ENCFF129LKB/@@download/ENCFF129LKB.h5ad

# short-read deep Split-seq single-cell myoblast 1,000 cell lib
wget -O short_deep.h5ad https://www.encodeproject.org/files/ENCFF418TAK/@@download/ENCFF418TAK.h5ad

# short-read deep Split-seq single-nucleus myoblast 1,000 cell lib
wget -O short_deep_nuc.h5ad https://www.encodeproject.org/files/ENCFF755KGW/@@download/ENCFF755KGW.h5ad

# long-read Split-seq single-cell myoblast data
# file 1 is a sparse transcript count matrix
wget -O long_transcript.h5ad https://www.encodeproject.org/files/ENCFF513MSS/@@download/ENCFF513MSS.h5ad
# file 2 is a sparse gene count matrix
wget -O long_gene.h5ad https://www.encodeproject.org/files/ENCFF508RNZ/@@download/ENCFF508RNZ.h5ad

# long-read Split-seq single-nucleus myoblast data
# file 1 is a sparse transcript count matrix
wget -O long_nuc_transcript.h5ad https://www.encodeproject.org/files/ENCFF033NVO/@@download/ENCFF033NVO.h5ad
# file 2 is a sparse gene count matrix
wget -O long_nuc_gene.h5ad https://www.encodeproject.org/files/ENCFF301DZH/@@download/ENCFF301DZH.h5ad

### myotube data ###
# short-read shallow Split-seq single-cell myotube 9,000 cell libs
wget -O short_shallow_myotube.h5ad https://www.encodeproject.org/files/ENCFF545EUV/@@download/ENCFF545EUV.h5ad
# short-read deep Split-seq single-cell myotube 1,000 cell lib
wget -O short_deep_myotube.h5ad https://www.encodeproject.org/files/ENCFF588HWS/@@download/ENCFF588HWS.h5ad
# long-read Split-seq single-cell myotube data
# file 1 is a sparse transcript count matrix
wget -O long_myotube_transcript.h5ad https://www.encodeproject.org/files/ENCFF362NGZ/@@download/ENCFF362NGZ.h5ad
# file 2 is a sparse gene count matrix
wget -O long_myotube_gene.h5ad https://www.encodeproject.org/files/ENCFF049LST/@@download/ENCFF049LST.h5ad