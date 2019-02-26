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
#' @import Seurat
#' @importFrom heatmaply heatmaply
#' @importFrom stats cor
#'
#' @return A matrix of correlation values
#' @export
#'
#' @examples
SeuratCompareR <- function(first.obj, ...){
  UseMethod('SeuratCompareR')
}

#' @rdname SeuratCompareR
#' @method SeuratCompareR seurat
#' @export
#' @return
SeuratCompareR.seurat <- function(first.obj,
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
    second.avg <- AverageExpression(
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

#' @rdname SeuratCompareR
#' @method SeuratCompareR Seurat
#' @export
#' @return
SeuratCompareR.Seurat <- function(first.obj,
                                  second.obj,
                                  group.by = NULL,
                                  add.ident = NULL,
                                  cor.function.use = cor,
                                  do.plot = FALSE,
                                  ...) {
  if (class(x = first.obj) == "Seurat") {
    if (!is.null(x = group.by)) {
      Idents(first.obj) <- id = group.by
    }
    first.avg <- AverageExpression(object = first.obj,
                                   add.ident = add.ident,
                                   use.scale = TRUE)
  } else {
    if (is.data.frame(x = first.obj)) {
      first.avg <- first.obj
    }
  }

  if (class(x = second.obj) == "Seurat") {
    if (!is.null(x = group.by)) {
      Idents(second.obj) <- group.by
    }
    second.avg <- AverageExpression(object = second.obj,
                                    add.ident = add.ident,
                                    use.scale = TRUE)
  } else {
    if (is.data.frame(x = second.obj)) {
      second.avg <- second.obj
    }
  }

  rownames(first.avg) <- toupper(rownames(first.avg))
  rownames(second.avg) <- toupper(rownames(second.avg))

  cor.result <- CCAcompare(ref.mat = first.avg,
                           expr.mat = second.avg,
                           cor.function.use = cor.function.use)
  if (do.plot) {
    return(heatmaply(cor.result))
  } else {
    return(cor.result)
  }
}
