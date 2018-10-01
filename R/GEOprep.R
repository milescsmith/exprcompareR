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
