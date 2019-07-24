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
#'   no space between.  Default: "hsapiens" .
#' @param rename_samples Replace the GEO sample name with the title
#'   provided in the ExpressionSet's phenoData@data slot. Default: TRUE
#'
#' @import biomaRt
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
                    var_list = NULL,
                    index = 1,
                    probe_set = NULL,
                    annotation_db = NULL,
                    species = "hsapiens",
                    rename_samples = TRUE) {

  # access data from a GEO entry
  ref <- getGEO(GEOaccession)
  exprs_mat <- exprs(ref[[index]])
  exprs_mat %<>% as.data.frame() %>% rownames_to_column(var = "probe_id")

  if (!is.null(annotation_db)) {
    exprs_mat$gene_name <- mapIds(
      annotation_db,
      keys = exprs_mat$probe_id,
      column = "SYMBOL",
      keytype = "PROBEID")
    exprs_mat <- exprs_mat %>% filter(!is.na(gene_name))
  } else if (!is.null(probe_set)) {
    mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = glue("{species}_gene_ensembl"),
                   host = "useast.ensembl.org")
    translation_table = getBM(attributes = c(probe_set, "external_gene_name"),
                              filters = probe_set,
                              values = exprs_mat$probe_id,
                              mart = mart)
  } else {
    translate <- data.frame(
      "probe_id" = ref[[index]]@featureData@data$ID,
      "gene_name" = ref[[index]]@featureData@data[[grep(pattern="[Ss]ymbol", x = names(ref[[2]]@featureData@data), value = TRUE)]]
    )
    # use toupper() as a quick hack that will let us compare human and mouse
    # genes.  Sure, going through homologene or the like would be better, but
    #that takes work.
    exprs_mat$gene_name <- mapvalues(
      x = exprs_mat$probe_id,
      from = translate$probe_id,
      to = toupper(as.character(translate$gene_name))
    )
  }

  # Speeds things up and prevents certain errors.  No sense in keeping any data
  # we are not going to use downstream.
  if (!is.null(var_list)) {
    exprs_mat %<>% filter(gene_name %in% var_list)
  }

  exprs_mat <- ddply(.data = exprs_mat,
                     .variables = "gene_name",
                     .fun = numcolwise(sum)) %>%
    column_to_rownames(var = "gene_name")

  if (isTRUE(rename_samples)) {
    colnames(exprs_mat) <- ref[[index]]@phenoData@data$title
  }

  # TODO need a function that can merge sample replicates.

  return(exprs_mat)
}
