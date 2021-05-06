#' @title Predict labels of a dataset
#' @description Predict the labels of a dataset using a tree object
#' @author Lieke Michielsen and Jose Alquicira-Hernandez
#' @param tree data.tree object storing the hierarchy and trained classifiers per layer
#' @param newData Seurat object containing cells to annotate
#' @param threshold Threshold used for probabilities to classify cells into classes. All cells below 
#' this threshold value will be labels as "unassigned". 
#' @param recompute_alignment Recompute alignment?
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


predictTree <- function(tree, newData, threshold = 0.75, recompute_alignment = TRUE){
  
  newData <- predictNode(tree, newData, threshold, recompute_alignment)
  
  newData
  
}


#' @title Predict labels of a dataset
#' @description Recursive function to predict the labels of a dataset using a tree object
#' @author Lieke Michielsen and Jose Alquicira-Hernandez
#' @param tree data.tree object storing the hierarchy and trained classifiers per layer
#' @param newData Seurat object containing cells to annotate
#' @param threshold Threshold used for probabilities to classify cells into classes. All cells below 
#' this threshold value will be labels as "unassigned". 
#' @param recompute_alignment Recompute alignment?
#' @return A Seurat object with additional metadata columns with prediction probabilities associated to each 
#' class, a \code{prediction} column, indicating the classification based on the provided threshold and 
#' a \code{generic_class} column without "unassigned" labels.
#' @importFrom scPred scPredict
#' @export
#' @examples
#'




# Iterate over the nodes recursively
predictNode <- function(tree, newData, threshold, recompute_alignment){
  
  # Check if reconstruction error was used during training
  if(tree$RE){
    
    print('Reconstruction Error')
    RE_rejected = reconstructionError(newData, tree$model, tree$RE)
    
    RE_rejected_idx = which(RE_rejected)
    print(length(RE_rejected_idx))
    RE_notrejected_idx = which(!RE_rejected)
    print(length(RE_notrejected_idx))
    
    # Make predictions (only for cells that are not rejected!)
    newData$scpred_prediction = ''
    if(length(RE_notrejected_idx) > 0){
      
      newData1 <- newData[,RE_notrejected_idx]
      newData1 <- scPredict(newData1, tree$model, threshold = threshold, recompute_alignment = recompute_alignment)
      
      newData$scpred_prediction[RE_notrejected_idx] = newData1$scpred_prediction
      
    }
    
    newData$scpred_prediction[RE_rejected_idx] = 'unassigned'
    
  } else{
    
    newData <- scPredict(newData, tree$model, threshold = threshold, recompute_alignment = recompute_alignment)
    
  }
  
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
		  dataSubset <- predictNode(c, dataSubset, threshold, recompute_alignment)
		  
		  # Add labels to original object
		  newData$scpred_prediction[idxChildren] <- dataSubset$scpred_prediction
      }
    }
  }
  
  newData    
  
}

#' @title Reject cells based on reconstruction error
#' @description Recursive function to predict the labels of a dataset using a tree object
#' @author Lieke Michielsen and Jose Alquicira-Hernandez
#' @param newData Seurat object containing cells to annotate
#' @param spmodel scPred model
#' @param RE_threshold Threshold for the reconstruction error. If the reconstruction error of a cell is above this threshold
#' it is rejected. This threshold is determined during the training step.
#' @return A boolean vector indicating which cells are rejected.
#' @export
#' @examples
#'



reconstructionError <- function(newData, spmodel, RE_threshold){
  
  ref_loadings <- spmodel@feature_loadings
  new_features <- rownames(newData)
  reference_features <- rownames(ref_loadings)
  
  shared_features <- intersect(reference_features, new_features)
  
  ref_loadings <- ref_loadings[shared_features,]
  
  new_data <- GetAssayData(newData, "data")[shared_features,]
  shared_scaling <- spmodel@scaling[shared_features,]
  means <- shared_scaling$means
  stdevs  <- shared_scaling$stdevs
  new_data <- Matrix::t(new_data)
  
  scaled_data <- scale(new_data, means, stdevs)
  new_embeddings <- scaled_data %*% ref_loadings
  
  new_inverse = new_embeddings %*% t(ref_loadings)
  
  new_inverse = new_inverse[, order(colnames(new_inverse))]
  
  new_original = as.matrix(scaled_data[, order(colnames(scaled_data))])
  
  RE = new_inverse - new_original
  RE = RE ^ 2
  RE = rowSums(RE)
  RE = sqrt(RE)
  
  return(RE > RE_threshold)
  
}

