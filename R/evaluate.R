#' @title Evaluate performance of hierarchical classifier
#' @description Calculated hierarchical F1-score, median F1-score, percentage of rejected cells, 
#' percentage of internally labeled cells, confusion matrix, and population size
#' @author Lieke Michelsen and Jose Alquicira-Hernandez
#' @param y_true vector with true labels
#' @param y_pred vector with predicted labels
#' @param tree data.tree object storing the hierarchy used during classification
#' @return A list containing the different evaluation metrics (median HF1-score, HF1-score per cell population,
#' median F1-score, F1-score per cell population, % rejected cells, % internally labeled cells, confusion matrix
#' and population size)
#' @export
#' @examples
#'
#' library(data.tree)
#' 
#' res <- evaluate(y_true, y_pred, tree)
#'

evaluate <- function(y_true, y_pred, tree){
  
  # %rejected
  root = tree$name
  rejected = 100*sum(y_pred == root)/length(y_pred)
  
  # % internal
  leafs = data.frame(tree$leaves)
  internal = 100*(length(y_pred) - sum(y_pred %in% leafs) - sum(y_pred == root))/length(y_pred)
  
  # Hierarchical F1-score per cell type
  HF1_all = c()
  HF1_names <- unique(y_true)
  
  for(ct in unique(y_true)){
    
    sum_pred = 0
    sum_true = 0
    sum_overlap = 0
    
    set_true = findSet(ct, tree)
    
    if(is.null(set_true)){
      HF1_names <- HF1_names[HF1_names != ct]
      next
    }
    
    #iterate of all the cells of this cell type
    idx = which(y_true == ct)
    
    for(i in idx){
      
      pred = y_pred[i]
      
      if(pred == tree$name){
        set_pred = c()
      } else{
        set_pred = findSet(pred, tree)
      }
      
      overlap = length(intersect(set_true, set_pred))
      
      if((length(set_pred) > length(set_true)) & (length(set_true == overlap))){
        sum_pred = sum_pred + overlap
      } else{
        sum_pred = sum_pred + length(set_pred)
      }
      
      sum_true = sum_true + length(set_true)
      sum_overlap = sum_overlap + overlap
      
    }
    
    HP <- sum_overlap / sum_pred
    HR <- sum_overlap / sum_true
    
    HF1 <- (2*HP*HR)/(HP+HR)
    HF1_all <- append(HF1_all, HF1)
    
  }
  
  names(HF1_all) <- HF1_names
  HF1 <- median(HF1_all)
  
  ## F1-score per cell type
  unique_true <- unique(y_true)
  unique_pred <- unique(y_pred)
  
  unique_all <- unique(c(unique_true,unique_pred))
  conf <- table(y_true, y_pred)
  pop_size <- rowSums(conf)
  
  
  # Here we need to get a vector with the name of the root node and all internal nodes
  nodes_all <- Traverse(tree)
  remove_mf1 <- c()
  for(no in nodes_all){
    if(!no$isLeaf){
      remove_mf1 <- append(remove_mf1, no$name)
    }
  }
  conf_F1 <- table(y_true,y_pred,exclude = remove_mf1)
  
  F1 <- vector()
  
  for (i in c(1:length(unique_true))){
    findLabel = colnames(conf_F1) == row.names(conf_F1)[i]
    if(sum(findLabel)){
      prec <- conf_F1[i,findLabel] / colSums(conf_F1)[findLabel]
      rec <- conf_F1[i,findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0){
        F1[i] = 0
      } else{
        F1[i] <- (2*prec*rec) / (prec + rec)
      }
    } else {
      F1[i] = 0
    }
  }
  
  pop_size <- pop_size[pop_size > 0]
  
  names(F1) <- names(pop_size)
  
  MF1 <- median(F1)
  
  return(list(HF1 = HF1, 
              HF1_all = HF1_all,
              rejected = rejected,
              internal = internal,
              MF1 = MF1, 
              F1_all = F1,
              conf = conf,
              pop_size = pop_size))
  
}


#' @title Find ancestor of a node in the tree
#' @description Return a set containing the node itself and its ancestors
#' @author Lieke Michelsen and Jose Alquicira-Hernandez
#' @param name name of the node
#' @param tree data.tree object storing the hierarchy used during classification
#' @return A set containing the name of the node and its ancestors
#' @export
#' @examples
#'
#' library(data.tree)
#' 
#' nodeSet <- findSet('B cells', tree)
#'

findSet <- function(name, tree){
  
  nodes_set <- c()
  nodes_all <- Traverse(tree)
  
  for(no in nodes_all){
    
    # Correct node found, now look at ancestors
    if(no$name == name){
      nodes_set <- append(nodes_set, name)
      while(no$parent$name != tree$name){
        nodes_set = append(nodes_set, no$parent$name)
        no <- no$parent
      }
    }
  }
  
  nodes_set
  
}
