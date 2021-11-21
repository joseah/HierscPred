#' @title Predict labels of a dataset
#' @description Predict the labels of a dataset using a tree object
#' @author Lieke Michielsen and Jose Alquicira-Hernandez
#' @param tree data.tree object storing the hierarchy and trained classifiers per layer
#' @param newData Seurat object containing cells to annotate
#' @param threshold Threshold used for probabilities to classify cells into classes. All cells below
#' this threshold value will be labels as "unassigned".
#' @param max.iter.harmony Maximum number of clustering iterations by harmony
#' @return A Seurat object with additional metadata columns with prediction probabilities associated to each
#' class, a \code{prediction} column, indicating the classification based on the provided threshold and
#' a \code{generic_class} column without "unassigned" labels.
#' @export
#' @examples
#'
#' branches <- c("Myeloid/DC/pDC", "Myeloid/DC/cDC")
#' h <- create_hierarchy(branches)
#' train_tree(dataTrain, h)
#' dataTest <- predictTree(h, dataTest)


predictTree <- function(tree, newData, threshold = 0, max.iter.harmony = 10){

  newData <- predictNode(tree, newData, threshold, max.iter.harmony = max.iter.harmony)

  newData

}


#' @title Predict labels of a dataset
#' @description Recursive function to predict the labels of a dataset using a tree object
#' @author Lieke Michielsen and Jose Alquicira-Hernandez
#' @param tree data.tree object storing the hierarchy and trained classifiers per layer
#' @param newData Seurat object containing cells to annotate
#' @param threshold Threshold used for probabilities to classify cells into classes. All cells below
#' this threshold value will be labels as "unassigned".
#' @param max.iter.harmony Maximum number of clustering iterations by harmony
#' @return A Seurat object with additional metadata columns with prediction probabilities associated to each
#' class, a \code{prediction} column, indicating the classification based on the provided threshold and
#' a \code{generic_class} column without "unassigned" labels.
#' @importFrom scPred scPredict
#' @export
#' @examples
#'




# Iterate over the nodes recursively
predictNode <- function(tree, newData, threshold, max.iter.harmony){


  newData <- scPredict(newData, tree$model, threshold = threshold, max.iter.harmony = max.iter.harmony)


  # Assign cells to parent node if they are unassigned
  newData$scpred_prediction <- ifelse(newData$scpred_prediction == "unassigned", tree$name, newData$scpred_prediction)

  for(c in tree$children){

    # If c is not a leaf node, continue with predictions
    if(!isLeaf(c)){

      # Get subset of the data
      namesChildren <- as.vector(c$Get('name'))
      idxChildren <- newData$scpred_prediction %in% namesChildren
	  if(sum(idxChildren > 0)){

		  dataSubset <- subset(newData, cells = Cells(newData)[idxChildren])

		  # Predict labels
		  dataSubset <- predictNode(c, dataSubset, threshold, max.iter.harmony)

		  # Add labels to original object
		  newData$scpred_prediction[idxChildren] <- dataSubset$scpred_prediction
      }
    }
  }

  newData

}
