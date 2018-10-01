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
