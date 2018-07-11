#' GEOcompare
#' 
#' Compares the averaged expression of each cluster within a Seurat object
#' to that of an entry in the Gene Expression Omnibus database.
#' 
#' Performs a canonical correlation analysis on the two datasets and then
#' assesses the correlation using the CCA data.
#'
#' @param seuratObj 
#' @param GEOaccession A GEO accession number (i.e. 'GSE24759')
#' @param do.plot Plot the correlation matrix
#' @param group.by Identifier or meta.data column by which to group the data. (default: 'ident')
#' @param add.ident An additional identifier to use as a grouping variable.
#' @param cor.function.use = cor
#' @param method.use = 'pearson' 
#'
#' @importFrom GEOquery getGEO
#' @importFrom dplyr as.data.frame rownames_to_column
#' @importFrom Seurat AverageExpression SetAllIdent
#' @importFrom plyr mapvalues ddply numcolwise
#' @importFrom irlba irlba
#' @importFrom heatmaply heatmaply
#' 
#' @return A matrix of correlation values
#' @export
#'
#' @examples
GEOcompare <- function(seuratObj, 
                       GEOaccession, 
                       do.plot = FALSE, 
                       group.by = NULL, 
                       add.ident = NULL,
                       function.use = cor){
  if(!is.null(group.by)){
    seuratObj <- SetAllIdent(seuratObj, group.by)
  }
  seurat.avg <- AverageExpression(seuratObj, add.ident = add.ident)
  
  # access data from a GEO entry
  ref <- getGEO(GEOaccession)
  ref.exprs <- exprs(ref[[1]])
  ref.exprs <- ref.exprs %>% as.data.frame() %>% rownames_to_column(var = "probe_id")
  translate <- data.frame('probe_id' = ref[[1]]@featureData@data$ID, 
                          'gene_name' = ref[[1]]@featureData@data$`Gene Symbol`)
  ref.exprs$gene_name <- mapvalues(x = ref.exprs$probe_id, from = translate$probe_id, to = as.character(translate$gene_name))
  # Speeds things up and prevents certain errors.  No sense in keeping any data we are not going to use downstream.
  ref.exprs <- ref.exprs %>% filter(gene_name %in% rownames(seurat.avg))
  ref.exprs <- ddply(.data = ref.exprs, .variables = "gene_name", .fun = numcolwise(sum)) %>% column_to_rownames(var = 'gene_name')

  # All subsequent steps rely on an overlapping number of variables, so find the genes in both datasets and then
  # subset the matrices
  genes.use <- intersect(rownames(ref.exprs), rownames(seurat.avg))
  ref.exprs <- ref.exprs[genes.use,]
  seurat.avg <- seurat.avg[genes.use,]
  
  # for CCA:
  # given that reference.mat is a gene x sample matrix
  scaled.ref.exprs <- scale(ref.exprs)
  # seurat.avg is the matrix output from Seurat::AverageExpression()
  scaled.seurat.avg <- scale(seurat.avg)
  sample.reference.dot.product <- t(scaled.ref.exprs) %*% scaled.seurat.avg

  cca.results <- irlba(A = sample.reference.dot.product, 
                       nv = min(dim(sample.reference.dot.product))-1) 
  # note: nv *must* be smaller than either dimension of A
  cca.data <- rbind(cca.results$u, cca.results$v)
  colnames(cca.data) <- paste0("CC",seq(1:(min(dim(sample.reference.dot.product))-1)))
  rownames(cca.data) <- c(as.character(ref[[1]]@phenoData@data$title), colnames(seurat.avg))
  
  # originally thought I needed to align the data, but that isn't true.  What we are looking for is a signature, 
  # and we can find that by examining the correlation between the two matrices.
  reference.cca <- cca.data[1:dim(ref.exprs)[[2]], ]
  seurat.cca <- cca.data[dim(ref.exprs)[[2]]+1:dim(seurat.avg)[[2]], ]
  cor.mat <- cor.function.use(t(reference.cca), t(seurat.cca), method = method.use)
  
  if(do.plot){
    heatmaply(cor.mat)
  }
  return(cor.mat)
}
