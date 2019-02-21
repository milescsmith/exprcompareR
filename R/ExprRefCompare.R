#' ExprRefCompare
#'
#' Compares the averaged expression of each cluster within a Seurat object
#' to that of an entry in the Gene Expression Omnibus database or the
#' ArrayExpress Archive of Functional Genomics Data
#'
#' Performs a canonical correlation analysis on the two datasets and then
#' assesses the correlation using the CCA data.
#'
#' @param object Either a processed Seurat object (with ScaleData run) or
#'   the output from AverageExpression. If a Seurat object is passed,
#'   AverageExpression will be run.
#' @param accession A GEO or ArrayExpress accession number (i.e. 'GSE24759')
#' @param GEOentry_index Indicate which ExpressionSet to use in a GEO entry
#'   containing multiple ExpressionSets
#' @param rename_samples Replace the sample name with that provided in the
#'   ExpressionSet's phenoData@data slot
#' @param do_plot Plot the correlation matrix. Default: FALSE
#' @param group_by Identifier or meta.data column by which to group the data.
#'   Default: 'ident'
#' @param add_ident An additional identifier to use as a grouping variable.
#'   Default: NULL
#' @param cor_function_use Function to use to calculate correlation.
#'   Default: stats::cor
#' @param annotation_db For entries that do not contain a probeID-to-gene mapping,
#'   it is necessary to supply the name of the annotation package that
#'   corresponds to the assay used in the GEO/AE entry (i.e.
#'   hugene10sttranscriptcluster.db for the [HuGene-1_0-st] Affymetrix Human
#'   Gene 1.0 ST Array [transcript (gene) version].  Default: NULL.
#' @param probe_set If it is difficult/impossible to find the correct annotation
#'   package for the assay, you can use biomaRt as a backup by supplying the
#'   probe_set name and species.  For available probesets, run
#'   biomaRt::listAttributes(mart = useMart(mart = "ENSEMBL_MART_ENSEMBL",
#'   dataset = "{species}_gene_ensembl")), where {species} is the same as below.
#'   Default = NULL
#' @param species = Dataset species, in the form of genus initial, species with
#'   no space between.  Default: "hsapiens"
#' @param platform For ArrayExpress datasets, the Annotation Database
#'   corresponding to the platform used in the assay (i.e. hgu133plus2.db).
#'   If not supplied, will attempt to detect the correct database and retrieve
#'   the names from Biomart (which can be SLOW) Default: NULL
#' @param ... Additional parameters to pass to the correlation function.
#'
#' @import stringr
#' @import Seurat
#' @importFrom heatmaply heatmaply
#' @importFrom stats cor
#'
#' @return A matrix of correlation values
#' @export
#'
#' @examples
ExprRefCompare <- function(object, ...){
  UseMethod("ExprRefCompare")
}

#' @rdname ExprRefCompare
#' @method ExprRefCompare seurat
#' @export
#' @return
ExprRefCompare.seurat <- function(seuratObj,
                           accession,
                           GEOentry_index = 1,
                           rename_samples = TRUE,
                           do_plot = FALSE,
                           group_by = NULL,
                           add_ident = NULL,
                           cor_function_use = cor,
                           platform = NULL,
                           probe_set = NULL,
                           annotation_db = NULL,
                           species = "hsapiens",
                           ...) {
  if (class(x = seuratObj) == "seurat") {
    if (!is.null(x = group_by)) {
      seuratObj <- SetAllIdent(seuratObj, group_by)
    }
    seurat_avg <- AverageExpression(object = seuratObj,
                                    add.ident = add_ident,
                                    use.scale = TRUE)
  } else {
    if (is.data.frame(seuratObj)) {
      seurat_avg <- seuratObj
    }
  }

  rownames(seurat_avg) <- toupper(rownames(seurat_avg))

  if (str_sub(accession, 1, 1) == "G") {
    ref_mat <- GEOprep(
      GEOaccession = accession,
      var_list = rownames(seurat_avg),
      index = GEOentry_index,
      probe_set = probe_set,
      annotation_db = annotation_db,
      species = species,
      rename_samples = rename_samples
    )
  } else if (str_sub(accession, 1, 1) == "E") {
    ref_mat <- ArrayExpressPrep(
      AEaccession = accession,
      var_list = rownames(seurat_avg),
      rename_samples = rename_samples,
      platform = platform
    )
  }

  cor_result <- CCAcompare(
    ref_mat = ref_mat,
    expr_mat = seurat_avg,
    cor_function_use = cor_function_use
  )
  if (do_plot) {
    return(heatmaply(cor_result))
  } else {
    return(cor_result)
  }
}

#' @rdname ExprRefCompare
#' @method ExprRefCompare Seurat
#'
#' @export
#' @return
ExprRefCompare.Seurat <- function(object,
                                  accession,
                                  GEOentry_index = 1,
                                  rename_samples = TRUE,
                                  do_plot = FALSE,
                                  group_by = NULL,
                                  assay = "RNA",
                                  add_ident = NULL,
                                  cor_function_use = cor,
                                  platform = NULL,
                                  probe_set = NULL,
                                  annotation_db = NULL,
                                  species = "hsapiens",
                                  ...) {
  if (class(x = object) == "seurat") {
    if (!is.null(x = group_by)) {
      Idents(object) <- object$group_by
    }
    seurat_avg <- AverageExpression(object = object,
                                    assay = assay,
                                    use.scale = TRUE)
  } else {
    if (is.data.frame(object)) {
      seurat_avg <- object
    }
  }

  rownames(seurat_avg) <- toupper(rownames(seurat_avg))

  if (str_sub(accession, 1, 1) == "G") {
    ref_mat <- GEOprep(
      GEOaccession = accession,
      var_list = rownames(seurat_avg),
      index = GEOentry_index,
      probe_set = probe_set,
      annotation_db = annotation_db,
      species = species,
      rename_samples = rename_samples
    )
  } else if (str_sub(accession, 1, 1) == "E") {
    ref_mat <- ArrayExpressPrep(
      AEaccession = accession,
      var_list = rownames(seurat_avg),
      rename_samples = rename_samples,
      platform = platform
    )
  }

  cor_result <- CCAcompare(
    ref_mat = ref_mat,
    expr_mat = seurat_avg,
    cor_function_use = cor_function_use
  )
  if (do_plot) {
    return(heatmaply(cor_result))
  } else {
    return(cor_result)
  }
}
