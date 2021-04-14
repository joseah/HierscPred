do_pca <- function(data, verbose = FALSE, scaledata = FALSE){

  if(ncol(data) < 50){
    npcs <- ncol(data) - 2
  } else {
    npcs <- 50
  }
  
  if (scaledata){
    data <- FindVariableFeatures(data, verbose = verbose)
    data <- ScaleData(data, verbose = verbose)
    
  }

  data <- RunPCA(data, npcs = npcs, verbose = verbose, seed.use = 42)
  data
}
