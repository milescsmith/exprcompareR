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


#' CCAcompare
#'
#' Compares the an experimental matrix of expression values against a reference
#' matrix of expression values to that of an entry in the Gene Expression
#' Omnibus database.
#'
#' Performs a canonical correlation analysis on the two datasets and then
#' assesses the correlation using the CCA data.
#'
#' @param ref.mat A variable x observation (i.e. genes as rows, samples as
#'   columns) matrix of reference values
#' @param expr.mat A variable x observation matrix of experimental values
#' @param do.plot Plot the correlation matrix. Default: FALSE
#' @param cor.function.use Correlation function to be used.
#'   Default: stats::cor
#' @param ... Additional parameters to pass to the correlation function
#'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom plyr mapvalues ddply numcolwise
#' @importFrom irlba irlba
#' @importFrom stats cor na.omit
#' @importFrom glue glue
#'
#' @return A matrix of correlation values
#' @export
#'
#' @examples
CCAcompare <- function(ref.mat,
                       expr.mat,
                       do.plot = FALSE,
                       cor.function.use = cor,
                       ...) {

  # All subsequent steps rely on an overlapping number of variables, so find
  # the genes in both datasets and then subset the matrices
  ref.mat <- na.omit(ref.mat)
  expr.mat <- na.omit(expr.mat)
  common.variables <- intersect(rownames(ref.mat), rownames(expr.mat))
  ref.mat <- ref.mat[common.variables, ]
  expr.mat <- expr.mat[common.variables, ]

  # for CCA:
  # given that ref.mat and expr.mat are gene x sample matrices
  scaled.ref.mat <- scale(ref.mat)
  scaled.expr.mat <- scale(expr.mat)
  expr.ref.dot.product <- t(scaled.ref.mat) %*% scaled.expr.mat

  cca.results <- irlba(
    A = expr.ref.dot.product,
    nv = min(dim(expr.ref.dot.product)) - 1,
    work = min(dim(expr.ref.dot.product)) - 1 + 35,
    maxit = 50000
  )
  # note: nv *must* be smaller than either dimension of A
  cca.data <- rbind(cca.results$u, cca.results$v)
  colnames(cca.data) <- glue("CC{seq(1:(min(dim(expr.ref.dot.product)) - 1))}")
  rownames(cca.data) <- c(colnames(ref.mat), colnames(expr.mat))

  ref.cca <- cca.data[1:dim(ref.mat)[[2]], ]
  expr.cca <- cca.data[dim(ref.mat)[[2]] + 1:dim(expr.mat)[[2]], ]
  cor.mat <- cor.function.use(t(ref.cca), t(expr.cca), ...)

  return(cor.mat)
}


#' GEOprep
#'
#' Processes an entry in the Gene Expression Omnibus database for use in
#' GEOcompare
#'
#' @param GEOaccession A GEO accession number (i.e. 'GSE24759')
#' @param var.list List of variables (i.e. genes) to include.
#'   If provided, data for other variables is discarded. Default: NULL
#' @param index Indicate which ExpressionSet to use in a GEO entry containing
#'   multiple ExpressionSets. Default: 1
#' @param rename.samples Replace the GEO sample name with the title
#'   provided in the ExpressionSet's phenoData@data slot. Default: TRUE
#'
#' @importFrom GEOquery getGEO
#' @importFrom Biobase exprs
#' @importFrom magrittr '%>%' '%<>%'
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
                    rename.samples = TRUE) {

  # access data from a GEO entry
  ref <- getGEO(GEOaccession)
  exprs.mat <- exprs(ref[[index]])
  exprs.mat %<>% as.data.frame() %>% rownames_to_column(var = "probe_id")
  translate <- data.frame(
    "probe_id" = ref[[index]]@featureData@data$ID,
    "gene_name" = ref[[index]]@featureData@data$`Gene Symbol`
  )
  # use toupper() as a quick hack that will let us compare human and mouse
  # genes.  Sure, going through homologene or the like would be better, but
  #that takes work.
  exprs.mat$gene_name <- mapvalues(
    x = exprs.mat$probe_id,
    from = translate$probe_id,
    to = toupper(as.character(translate$gene_name))
  )

  # Speeds things up and prevents certain errors.  No sense in keeping any data
  # we are not going to use downstream.
  if (!is.null(var.list)) {
    exprs.mat %<>% filter(gene_name %in% var.list)
  }

  exprs.mat <- ddply(.data = exprs.mat,
                     .variables = "gene_name",
                     .fun = numcolwise(sum)) %>%
    column_to_rownames(var = "gene_name")

  if (isTRUE(rename.samples)) {
    colnames(exprs.mat) <- ref[[index]]@phenoData@data$title
  }

  # TODO need a function that can merge sample replicates.

  return(exprs.mat)
}


#' ArrayExpressnPrep
#'
#' Processes an entry in the Gene Expression Omnibus database for use
#' in GEOcompare
#'
#' @param AEaccession An ArrayExpress accession number (i.e. 'E-MEXP-2360')
#' @param var.list List of variables (i.e. genes) to include.
#'   If provided, data for other variables is discarded. Default: NULL
#' @param rename.samples Replace the GEO sample name with the title
#' @param platform The Annotation Database corresponding to the platform used
#'   in the assay (i.e. hgu133plus2.db).
#'   If not supplied, will attempt to detect the correct database and retrieve
#'   the names from Biomart (which can be SLOW). Default: NULL
#'   provided in the ExpressionSet's phenoData@data$characteristics.celltype.
#'   slot. Default: TRUE
#'
#' @import dplyr
#' @import AnnotationDbi
#' @import stringr
#' @import biomaRt
#' @importFrom ArrayExpress ArrayExpress
#' @importFrom Biobase exprs
#' @importFrom magrittr '%>%' '%<>%'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom plyr mapvalues ddply numcolwise
#' @importFrom oligo rma
#'
#' @return A matrix of correlation values
#' @export
#'
#' @examples
ArrayExpressPrep <- function(AEaccession,
                             var.list = NULL,
                             rename.samples = TRUE,
                             platform = NULL) {
  ae <- ArrayExpress(AEaccession)
  names(ae@phenoData@data) <- tolower(names(ae@phenoData@data))
  ae <- oligo::rma(ae)

  exprs.mat <- exprs(ae)
  exprs.mat <- exprs.mat %>%
    as.data.frame() %>%
    rownames_to_column(var = "probe_id")

  if (is.null(platform)) {
    platform <- ae@annotation %>%
      str_replace(pattern = "pd\\.", replacement = "affy_") %>%
      str_replace_all(pattern = "\\.", replacement = "_")
    if (str_detect(string = ae@phenoData@data$characteristics.organism.,
                   pattern = "sapiens")) {
      species_dataset <- "hsapiens_gene_ensembl"
    } else if (str_detect(string = ae@phenoData@data$characteristics.organism.,
                          pattern = "musculus")) {
      species_dataset <- "mmusculus_gene_ensembl"
    }
    ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                       dataset = species_dataset)
    translate <- getBM(
      attributes = c(platform, "hgnc_symbol"),
      filters = platform,
      values = exprs.mat$probe_id,
      mart = ensembl
    )
  } else {
    translate <- AnnotationDbi::select(get(platform),
                                       keys = keys(get(platform),
                                                   keytype = "PROBEID"),
                                       columns = "SYMBOL",
                                       keytype = "PROBEID")
  }


  exprs.mat$gene_name <- mapvalues(
    x = exprs.mat$probe_id,
    from = translate[, 1],
    to = toupper(as.character(translate[, 2]))
  )
  exprs.mat %<>%
    filter(!is.na(gene_name)) %>%
    filter(!gene_name == "") %>%
    filter(!gene_name == probe_id)

  # Speeds things up and prevents certain errors.  No sense in keeping any
  # data we are not going to use downstream.
  if (!is.null(var.list)) {
    exprs.mat %<>% filter(gene_name %in% var.list)
  }

  exprs.mat <- ddply(.data = exprs.mat,
                     .variables = "gene_name",
                     .fun = numcolwise(sum)) %>%
    column_to_rownames(var = "gene_name")

  if (isTRUE(rename.samples)) {
    names(ae@phenoData@data) %<>%
      tolower() %>%
      str_replace_all(pattern = "cell\\.type", replacement = "celltype")
    colnames(exprs.mat) <- ae@phenoData@data$characteristics.celltype.
  }

  # TODO need a function that can merge sample replicates.

  return(exprs.mat)
}

#' SeuratCompareR
#'
#' Uses Canonical Correlation Analysis to compare the expression profile of the
#' identity classes of a Seurat object to those of a second Seurat object
#'
#' @param first.obj Either a processed Seurat object (with ScaleData run) or
#'   the output from AverageExpression.  If a Seurat object is passed,
#'   AverageExpression will be run.
#' @param second.obj A second (or the same) processed Seurat object or the
#'   output from AverageExpression for such an object.
#' @param group.by Identifier or meta.data column by which to group the data.
#'   Default: 'ident'
#' @param add.ident An additional identifier to use as a grouping variable.
#'   Default: NULL
#' @param cor.function.use Function to use to calculate correlation.
#'   Default: stats::cor
#' @param do.plot Plot the correlation matrix. Default: FALSE)
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
SeuratCompareR <- function(first.obj,
                           second.obj,
                           group.by = NULL,
                           add.ident = NULL,
                           cor.function.use = cor,
                           do.plot = FALSE,
                           ...) {
  if (class(x = first.obj) == "seurat") {
    if (!is.null(x = group.by)) {
      first.obj <- SetAllIdent(
        object = first.obj,
        id = group.by
      )
    }
    first.avg <- AverageExpression(
      object = first.obj,
      add.ident = add.ident,
      use.scale = TRUE
    )
  } else {
    if (is.data.frame(x = first.obj)) {
      first.avg <- first.obj
    }
  }

  if (class(x = second.obj) == "seurat") {
    if (!is.null(x = group.by)) {
      second.obj <- SetAllIdent(
        object = second.obj,
        id = group.by
      )
    }
    seurat.avg <- AverageExpression(
      object = second.obj,
      add.ident = add.ident,
      use.scale = TRUE
    )
  } else {
    if (is.data.frame(x = second.obj)) {
      second.avg <- second.obj
    }
  }

  rownames(first.avg) <- toupper(rownames(first.avg))
  rownames(second.avg) <- toupper(rownames(second.avg))

  cor.result <- CCAcompare(
    ref.mat = first.avg,
    expr.mat = second.avg,
    cor.function.use = cor.function.use
  )
  if (do.plot) {
    return(heatmaply(cor.result))
  } else {
    return(cor.result)
  }
}
