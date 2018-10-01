#' ExprRefCompare
#'
#' Compares the averaged expression of each cluster within a Seurat object
#' to that of an entry in the Gene Expression Omnibus database or the
#' ArrayExpress Archive of Functional Genomics Data
#'
#' Performs a canonical correlation analysis on the two datasets and then
#' assesses the correlation using the CCA data.
#'
#' @param seuratObj Either a processed Seurat object (with ScaleData run) or
#'   the output from AverageExpression. If a Seurat object is passed,
#'   AverageExpression will be run.
#' @param accession A GEO or ArrayExpress accession number (i.e. 'GSE24759')
#' @param GEOentry.index Indicate which ExpressionSet to use in a GEO entry
#'   containing multiple ExpressionSets
#' @param rename.samples Replace the sample name with that provided in the
#'   ExpressionSet's phenoData@data slot
#' @param do.plot Plot the correlation matrix. Default: FALSE
#' @param group.by Identifier or meta.data column by which to group the data.
#'   Default: 'ident'
#' @param add.ident An additional identifier to use as a grouping variable.
#'   Default: NULL
#' @param cor.function.use Function to use to calculate correlation.
#'   Default: stats::cor
#' @param platform For ArrayExpress datasets, the Annotation Database
#'   corresponding to the platform used in the assay (i.e. hgu133plus2.db).
#'   If not supplied, will attempt to detect the correct database and retrieve
#'   the names from Biomart (which can be SLOW) Default: NULL
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
                           platform = NULL,
                           ...) {
  if (class(x = seuratObj) == "seurat") {
    if (!is.null(x = group.by)) {
      seuratObj <- SetAllIdent(seuratObj, group.by)
    }
    seurat.avg <- AverageExpression(object = seuratObj,
                                    add.ident = add.ident,
                                    use.scale = TRUE)
  } else {
    if (is.data.frame(seuratObj)) {
      seurat.avg <- seuratObj
    }
  }

  rownames(seurat.avg) <- toupper(rownames(seurat.avg))

  if (str_sub(accession, 1, 1) == "G") {
    ref.mat <- GEOprep(
      GEOaccession = accession,
      var.list = rownames(seurat.avg),
      index = GEOentry.index,
      rename.samples = rename.samples
    )
  } else if (str_sub(accession, 1, 1) == "E") {
    ref.mat <- ArrayExpressPrep(
      AEaccession = accession,
      var.list = rownames(seurat.avg),
      rename.samples = rename.samples,
      platform = platform
    )
  }

  cor.result <- CCAcompare(
    ref.mat = ref.mat,
    expr.mat = seurat.avg,
    cor.function.use = cor.function.use
  )
  if (do.plot) {
    return(heatmaply(cor.result))
  } else {
    return(cor.result)
  }
}
