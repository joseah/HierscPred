#' @title Create hierarchical tree
#' @description Creates a tree object representing a cell-type hierarchy
#' @author Lieke Michielsen and Jose Alquicira-Hernandez
#' @param branches A character vector containing \code{n} number of branches.
#' Each branch corresponds to a cell type. Each level is delimited by  "/"
#' character. For example, the following branch is composed by three levels
#' and the _terminal_ cell type is \code{pDC}: "Myeloid/DC/pDC"
#' @return A data.tree object containing the cell type hierarchy
#' @importFrom data.tree as.Node
#' @export
#' @examples
#'
#' branches <- c("Myeloid/DC/pDC", "Myeloid/DC/cDC")
#' h <- create_hierarchy(branches)
#'

create_hierarchy <- function(branches){

  # Validate input
  if(!is.character(branches))
    stop("`branches` parameter value must be character")

  # Create dataframe
  h <- data.frame(pathString = branches)

  # Convert dataframe to a tree
  h <- as.Node(h)

  h
}


#' @title Create hierarchical tree
#' @description Creates a tree object representing a cell-type hierarchy
#' @author Lieke Michielsen and Jose Alquicira-Hernandez
#' @param tree Hierarchy created via \code{create_hierarchy}
#' @param data Seurat object containing cells used to train hierarchical models
#' @param pvar Prediction variable. Column name in Seurat object metadata
#' containing the cell labels of the terminal nodes of the hierarchy tree
#' @param reduction Name of reduction in Seurat objet to be used to determine
#' the feature space. Default: "pca"
#' @param scaledata Whether to select variable features and 
#' scale the data of the reference object. Default: FALSE
#' @param model Classification model supported via caret package.
#' A list of all models can be found here in https://topepo.github.io/caret/available-models.html
#' @param allowParallel Allow parallel processing for resampling?
#' @param reconstruction_error Use reconstruction error?
#' @param fn_perc TBD
#' @return A data.tree object containing the cell type hierarchy
#' @export
#' @examples
#'
#' branches <- c("Myeloid/DC/pDC", "Myeloid/DC/cDC")
#' h <- create_hierarchy(branches)
#' train_tree(data, h)
#'

train_tree <- function(data,
                       tree,
                       pvar = 'cell_type',
                       reduction = 'pca', scaledata = FALSE,
                       model = 'svmRadial', allowParallel = FALSE,
                       reconstruction_error = TRUE,
                       fn_perc = 0.01,
                       verbose = FALSE){

  meta_vars <- colnames(data[[]])
  reduction_list <- Reductions(data)

  if(!pvar %in% meta_vars)
    stop("Prediction variable is not present in metadata")

  if(!reduction %in% reduction_list)
    stop(reduction, " is not a reduction present in Seurat object")



  #labels <- data[[pVar]]

  tree <-  train_node(tree, data, pvar, reduction, scaledata, model, allowParallel,
                      reconstruction_error, fn_perc, verbose)
  tree
}

#' @title Create hierarchical tree
#' @description Creates a tree object representing a cell-type hierarchy
#' @author Lieke Michielsen and Jose Alquicira-Hernandez
#' @param tree Hierarchy created via \code{create_hierarchy}
#' @param data Seurat object containing cells used to train hierarchical models
#' @param pvar Prediction variable. Column name in Seurat object metadata
#' containing the cell labels of the terminal nodes of the hierarchy tree
#' @param reduction Name of reduction in Seurat objet to be used to determine
#' the feature space. Default: "pca"
#' #' @param scaledata Whether to select variable features and 
#' scale the data of the reference object. Default: FALSE
#' @param model Classification model supported via caret package.
#' A list of all models can be found here in https://topepo.github.io/caret/available-models.html
#' @param allowParallel Allow parallel processing for resampling?
#' @param reconstruction_error Use reconstruction error?
#' @param fn_perc TBD
#' @return A data.tree object containing the cell type hierarchy
#' @importFrom scPred getFeatureSpace trainModel get_scpred
#' @importFrom Seurat Cells
#' @importFrom data.tree isLeaf
#' @export
#' @examples
#'
#' branches <- c("Myeloid/DC/pDC", "Myeloid/DC/cDC")
#' h <- create_hierarchy(branches)
#' train_tree(data, h)
#'


train_node <- function(tree, data, pvar, reduction, scaledata, model, allowParallel,
                       reconstruction_error, fn_perc, verbose){

  cat("Training parent node: ", tree$name, "\n", sep = "")

  labels <- data[[pvar, drop = TRUE]]
  new_labels <- labels

  if(reconstruction_error){
      message("Determining reconstruction error threshold...")
    tree$RE <- find_threshold(data, labels, fn_perc = fn_perc, scaledata = scaledata, verbose = verbose)
  } else {
    tree$RE <- FALSE
  }

  if(verbose) message("Running PCA...")
  data <- do_pca(data, scaledata = scaledata)

  ## First rewrite the labels
  for(c in tree$children){

    names_children <- as.vector(c$Get('name'))
    idx_children <- labels %in% names_children

    new_labels[idx_children] <- c$name

  }

  data$response <- new_labels

  message("Cell types in this layer of the tree:\n- ",
      paste(unique(new_labels), collapse = "\n- "))
  ## get informative PCs and train classifier
  data <- getFeatureSpace(data, pvar = "response", reduction = reduction)
  data <- trainModel(data, model = model, allowParallel = allowParallel)

  message("Model trained")
  tree$model <- get_scpred(data)

  ## Continue to children if they are not a leaf
  for(c in tree$children){

    if(!isLeaf(c)){
      names_children <- as.vector(c$Get('name'))
      idx_children <- labels %in% names_children

      # Get subset of the data
      data_subset <- subset(data, cells = Cells(data)[idx_children])

      # Do PCA on this node and continue with children
      train_node(c, data_subset, pvar, reduction, scaledata, model, allowParallel,
                reconstruction_error, fn_perc, verbose)

    }

  }

  tree

}
