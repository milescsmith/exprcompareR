#' CCAcompare
#'
#' Compares the an experimental matrix of expression values against a reference
#' matrix of expression values to that of an entry in the Gene Expression
#' Omnibus database.
#'
#' Performs a canonical correlation analysis on the two datasets and then
#' assesses the correlation using the CC scores.
#'
#' @param ref_mat A variable x observation (i.e. genes as rows, samples as
#'   columns) matrix of reference values
#' @param expr_mat A variable x observation matrix of experimental values
#' @param do_plot Plot the correlation matrix. Default: FALSE
#' @param cor_function_use Correlation function to be used.
#'   Default: stats::cor
#' @param ... Additional parameters to pass to the correlation function
#'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom plyr mapvalues ddply numcolwise
#' @importFrom rsvd rsvd
#' @importFrom Matrix tcrossprod
#' @importFrom stats cor na.omit
#' @importFrom glue glue
#'
#' @return A matrix of correlation values
#' @export
#'
#' @examples
CCAcompare <- function(ref_mat,
                       expr_mat,
                       do_plot = FALSE,
                       cor_function_use = cor,
                       ...) {

  # All subsequent steps rely on an overlapping number of variables, so find
  # the genes in both datasets and then subset the matrices
  ref_mat <- na.omit(ref_mat)
  expr_mat <- na.omit(expr_mat)
  common_variables <- intersect(rownames(ref_mat), rownames(expr_mat))
  ref_mat <- ref_mat[common_variables, ]
  expr_mat <- expr_mat[common_variables, ]

  # for CCA:
  # given that ref.mat and expr.mat are gene x sample matrices
  scaled_ref_mat <- scale(ref_mat)
  scaled_expr_mat <- scale(expr_mat)
  expr_ref_dot_product <- tcrossprod(scaled_ref_mat,scaled_expr_mat)

  cca_results <- rsvd(
    A = expr_ref_dot_product,
    k = min(dim(expr_ref_dot_product)) - 1,
  )

  # cca_results <- irlba(
  #   A = expr_ref_dot_product,
  #   nv = min(dim(expr_ref_dot_product)) - 1,
  #   work = min(dim(expr_ref_dot_product)) - 1 + 35,
  #   maxit = 50000
  # )
  # note: nv *must* be smaller than either dimension of A
  cca_data <- rbind(cca_results$u, cca_results$v)
  colnames(cca_data) <- glue("CC{seq(1:(min(dim(expr_ref_dot_product)) - 1))}")
  rownames(cca_data) <- c(colnames(ref_mat), colnames(expr_mat))

  ref_cca <- cca_data[1:dim(ref_mat)[[2]], ]
  expr_cca <- cca_data[dim(ref_mat)[[2]] + 1:dim(expr_mat)[[2]], ]
  cor_mat <- cor_function_use(t(ref_cca), t(expr_cca), ...)

  return(cor_mat)
}
