#' @title Find threshold for reconstruction error
#' @description Determines threshold for reconstruction errors computed from
#' multiple resamples
#' @author Lieke Michielsen and Jose Alquicira-Hernandez
#' @param object Seurat object used as reference
#' @param labels Response variable
#' @param nfolds Number of resamples
#' @param fn_perc TBD
#' @param seed Numeric seed to create resamples
#' @param verbose Display Seurat messages?
#' @return A numeric value corresponding to the threshold for reconstruction
#' error
#' @importFrom Seurat DefaultAssay Loadings GetAssayData
#' @importFrom caret createFolds
#' @export
#' @examples
#'
#' library(scPred)
#' data(pbmc_1)
#' res <- find_threshold(pbmc_1, pbmc_1$cell_type)
#'

find_threshold <- function(object, labels, nfolds = 5, fn_perc = 0.01,
                           seed = 66, verbose = FALSE){

  default_assay <- DefaultAssay(object)

  # Add a 5fold CV loop to determine the threshold for the CV
  set.seed(seed)
  folds <-  createFolds(labels, k = nfolds)
  build_re <- function(fold, object, fn_perc, verbose){

    train <-  object[,-fold]
    test <-  object[,fold]

    # PCA on train, get feature loadings
    train <-  do_pca(train, verbose = verbose)

    # Get loadings and scaling info
    loadings <- Loadings(train, "pca")
    features <- rownames(loadings)
    train <- GetAssayData(train, "data", assay = default_assay)[features,]
    means <-  Matrix::rowMeans(train)

    rowVar <- function(x, ...) {
      sqrt(Matrix::rowSums((x - means)^2, ...)/(ncol(x) - 1))
    }

    stdevs  <- rowVar(train)

    # Transform test object and calculate RE
    test <- GetAssayData(test, "data", assay = default_assay)[features,]
    test <- scale(Matrix::t(test), means, stdevs)

    test_transformed <- test %*% loadings
    test_inverse <- test_transformed %*% t(loadings)

    RE <- test_inverse - test
    RE <- RE ** 2
    RE <- rowSums(RE)
    RE <- sqrt(RE)

    re_threshold <- quantile(RE, 1 - fn_perc)
    re_threshold
  }


  res <- lapply(folds, build_re, object, fn_perc = fn_perc, verbose)
  res <- median(unlist(res))

  res

}
