# exprcompareR
exprcompareR allows for the comparison of single cell RNAseq datasets with microarray and bulk RNAseq datasets.  Since there exist datasets of the latter two methods for highly purified populations of cells, these can be used to aid in the labeling of scRNAseq clusters as a particular cell type.

exprcompareR requires a scRNA-seq data in a Seurat object and a GEO or ArrayExpress accession number.  The Seurat object needs to have been scaled and had clusters identified.  exprcompareR works by scaling the expression matrices of both the reference dataset and the average expression of each scRNAseq cluster and performing canonical correlation analysis to identify gene covariance signatures.
