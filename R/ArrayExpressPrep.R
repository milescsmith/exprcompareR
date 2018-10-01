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
