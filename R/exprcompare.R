#' ExprRefCompare
#'
#' Compares the averaged expression of each cluster within a Seurat object
#' to that of an entry in the Gene Expression Omnibus database or the
#' ArrayExpress Archive of Functional Genomics Data
#'
#' Performs a canonical correlation analysis on the two datasets and then
#' assesses the correlation using the CCA data.
#'
#' @param seuratObj Either a processed Seurat object (with ScaleData run) or the output from AverageExpression.
#' If a Seurat object is passed, AverageExpression will be run.
#' @param accession A GEO or ArrayExpress accession number (i.e. 'GSE24759')
#' @param GEOentry.index Indicate which ExpressionSet to use in a GEO entry containing multiple ExpressionSets
#' @param rename.samples Replace the sample name with that provided in the ExpressionSet's phenoData@data slot
#' @param do.plot Plot the correlation matrix. (default: FALSE)
#' @param group.by Identifier or meta.data column by which to group the data. (default: 'ident')
#' @param add.ident An additional identifier to use as a grouping variable. (default: NULL)
#' @param cor.function.use Function to use to calculate correlation. (default: stats::cor)
#' @param ... Additional parameters to pass to the correlation function.
#'
#' @import stringr
#' @importFrom Seurat AverageExpression SetAllIdent
#' @importFrom heatmaply heatmaply
#' @importFrom stats cor
#'
#' @return A matrix of correlation values
#' @export
#'
#' @examples
ExprRefCompare <- function(seuratObj,
                           accession,
                           GEOentry.index = 1,
                           rename.samples = TRUE,
                           do.plot = FALSE,
                           group.by = NULL,
                           add.ident = NULL,
                           cor.function.use = cor,
                           ...){
  if(class(seuratObj) == 'seurat'){
    if(!is.null(group.by)){
      seuratObj <- SetAllIdent(seuratObj, group.by)
    }
    seurat.avg <- AverageExpression(seuratObj, add.ident = add.ident, use.scale = TRUE)
  } else {
    if(is.data.frame(seuratObj)){
      seurat.avg <- seuratObj
    }
  }

  rownames(seurat.avg) <- toupper(rownames(seurat.avg))

  if (str_sub(accession, 1, 1) == "G"){
    ref.mat <- GEOprep(GEOaccession = accession,
                       var.list = rownames(seurat.avg),
                       index = GEOentry.index,
                       rename.samples = rename.samples)
  } else if (str_sub(accession, 1, 1) == "E"){
    ref.mat <- ArrayExpressionPrep(AEaccession = accession,
                                   var.list = rownames(seurat.avg),
                                   rename.samples = rename.samples)
  }

  cor.result <- CCAcompare(ref.mat = ref.mat,
                           expr.mat = seurat.avg,
                           cor.function.use = cor.function.use)
  if (do.plot){
    return(heatmaply(cor.result))
  } else {
    return(cor.result)
  }
}


#' CCAcompare
#'
#' Compares the an experimental matrix of expression values against a reference matrix of expression values
#' to that of an entry in the Gene Expression Omnibus database.
#'
#' Performs a canonical correlation analysis on the two datasets and then
#' assesses the correlation using the CCA data.
#'
#' @param ref.mat A variable x observation (i.e. genes as rows, samples as columns) matrix of reference values
#' @param expr.mat A variable x observation matrix of experimental values
#' @param do.plot Plot the correlation matrix. (default: FALSE)
#' @param cor.function.use Correlation function to be used. (default: stats::cor)
#' @param ... Additional parameters to pass to the correlation function
#'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom plyr mapvalues ddply numcolwise
#' @importFrom irlba irlba
#' @importFrom stats cor na.omit

#'
#' @return A matrix of correlation values
#' @export
#'
#' @examples
CCAcompare <- function(ref.mat,
                       expr.mat,
                       do.plot = FALSE,
                       cor.function.use = cor,
                       ...){

  # All subsequent steps rely on an overlapping number of variables, so find the genes in both datasets and then
  # subset the matrices
  ref.mat <- na.omit(ref.mat)
  expr.mat <- na.omit(expr.mat)
  common.variables <- intersect(rownames(ref.mat), rownames(expr.mat))
  ref.mat <- ref.mat[common.variables,]
  expr.mat <- expr.mat[common.variables,]

  # for CCA:
  # given that ref.mat and expr.mat are gene x sample matrices
  scaled.ref.mat <- scale(ref.mat)
  scaled.expr.mat <- scale(expr.mat)
  expr.ref.dot.product <- t(scaled.ref.mat) %*% scaled.expr.mat

  cca.results <- irlba(A = expr.ref.dot.product,
                       nv = min(dim(expr.ref.dot.product))-1)
  # note: nv *must* be smaller than either dimension of A
  cca.data <- rbind(cca.results$u, cca.results$v)
  colnames(cca.data) <- paste0("CC",seq(1:(min(dim(expr.ref.dot.product))-1)))
  rownames(cca.data) <- c(colnames(ref.mat), colnames(expr.mat))

  ref.cca <- cca.data[1:dim(ref.mat)[[2]], ]
  expr.cca <- cca.data[dim(ref.mat)[[2]]+1:dim(expr.mat)[[2]], ]
  cor.mat <- cor.function.use(t(ref.cca), t(expr.cca), ...)

  return(cor.mat)
}


#' GEOprep
#'
#' Processes an entry in the Gene Expression Omnibus database for use in GEOcompare
#'
#' @param GEOaccession A GEO accession number (i.e. 'GSE24759')
#' @param var.list List of variables (i.e. genes) to include.
#' If provided, data for other variables is discarded. (default: NULL)
#' @param index Indicate which ExpressionSet to use in a GEO entry containing
#' multiple ExpressionSets. (default: 1)
#' @param rename.samples Replace the GEO sample name with the title
#' provided in the ExpressionSet's phenoData@data slot. (default: TRUE)
#'
#' @importFrom GEOquery getGEO
#' @importFrom Biobase exprs
#' @importFrom magrittr '%>%'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom plyr mapvalues ddply numcolwise
#'
#' @return A matrix of correlation values
#' @export
#'
#' @examples
GEOprep <- function(GEOaccession,
                    var.list = NULL,
                    index = 1,
                    rename.samples= TRUE){

  # access data from a GEO entry
  ref <- getGEO(GEOaccession)
  exprs.mat <- exprs(ref[[index]])
  exprs.mat <- exprs.mat %>% as.data.frame() %>% rownames_to_column(var = "probe_id")
  translate <- data.frame('probe_id' = ref[[index]]@featureData@data$ID,
                          'gene_name' = ref[[index]]@featureData@data$`Gene Symbol`)
  # use toupper() as a quick hack that will let us compare human and mouse genes.
  # sure, going through homologene or the like would be better, but that takes work.
  exprs.mat$gene_name <- mapvalues(x = exprs.mat$probe_id,
                                   from = translate$probe_id,
                                   to = toupper(as.character(translate$gene_name)))

  # Speeds things up and prevents certain errors.  No sense in keeping any data we are not going to use downstream.
  if (!is.null(var.list)){
    exprs.mat <- exprs.mat %>% filter(gene_name %in% var.list)
  }

  exprs.mat <- ddply(.data = exprs.mat, .variables = "gene_name", .fun = numcolwise(sum)) %>% column_to_rownames(var = 'gene_name')

  if(isTRUE(rename.samples)){
    colnames(exprs.mat) <- ref[[index]]@phenoData@data$title
  }

  # TODO need a function that can merge sample replicates.

  return(exprs.mat)
}


#' ArrayExpressionPrep
#'
#' Processes an entry in the Gene Expression Omnibus database for use in GEOcompare
#'
#' @param AEaccession An ArrayExpress accession number (i.e. 'E-MEXP-2360')
#' @param var.list List of variables (i.e. genes) to include.
#' If provided, data for other variables is discarded. (default: NULL)
#' @param rename.samples Replace the GEO sample name with the title
#' provided in the ExpressionSet's phenoData@data$characteristics.celltype. slot. (default: TRUE)
#'
#' @import dplyr
#' @import AnnotationDbi
#' @import stringr
#' @importFrom rlang quo
#' @importFrom ArrayExpress ArrayExpress
#' @importFrom Biobase exprs
#' @importFrom magrittr '%>%'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom plyr mapvalues ddply numcolwise
#' @importFrom BiocInstaller biocLite
#' @importFrom oligo rma
#'
#' @return A matrix of correlation values
#' @export
#'
#' @examples
ArrayExpressionPrep <- function(AEaccession,
                                var.list = NULL,
                                rename.samples= TRUE){

  ae <- ArrayExpress(AEaccession)
  names(ae@phenoData@data) <- tolower(names(ae@phenoData@data))
  ae <- oligo::rma(ae)

  exprs.mat <- exprs(ae)
  exprs.mat <- exprs.mat %>% as.data.frame() %>% rownames_to_column(var = "probe_id")

  platform <- ae@annotation %>%
    str_replace(pattern = "pd\\.", replacement = "") %>%
    str_replace_all(pattern = "\\.", replacement = "") %>%
    str_c(".db") %>%
    get()

  if (!require(platform)){
    should.install <- readline(prompt = paste0("The package ", quo(platform), " is needed but is not installed.  Would you like to attempt to install it?"))
    if (should.install){
      biocLite(platform)
      did.load <- require(platform)
      if (!did.load){
        stop("Sorry, installation failed.  This analysis cannot proceed.")
      }
    } else {
      stop("This analysis requires that package to translate probeIds to gene names.  It will not work without it.")
    }
  }

  translate <- AnnotationDbi::select(platform, keys = keys(platform, keytype="PROBEID"), columns="SYMBOL", keytype="PROBEID")
  exprs.mat$gene_name <- mapvalues(x = exprs.mat$probe_id,
                                   from = translate$PROBEID,
                                   to = toupper(as.character(translate$SYMBOL)))
  exprs.mat <- exprs.mat %>% filter(!is.na(gene_name))

  # Speeds things up and prevents certain errors.  No sense in keeping any data we are not going to use downstream.
  if (!is.null(var.list)){
    exprs.mat <- exprs.mat %>% filter(gene_name %in% var.list)
  }

  exprs.mat <- ddply(.data = exprs.mat, .variables = "gene_name", .fun = numcolwise(sum)) %>% column_to_rownames(var = 'gene_name')

  if(isTRUE(rename.samples)){
    names(ae@phenoData@data) <- names(ae@phenoData@data) %>%
      tolower() %>%
      str_replace_all(pattern = "cell\\.type", replacement = "celltype")
    colnames(exprs.mat) <- ae@phenoData@data$characteristics.celltype.
  }

  # TODO need a function that can merge sample replicates.

  return(exprs.mat)
}
