do_pca <- function(data, verbose = FALSE, scaledata = FALSE, selection.method = 'vst'){

  if(ncol(data) < 50){
    npcs <- ncol(data) - 2
  } else {
    npcs <- 50
  }

  if (scaledata){
    data <- FindVariableFeatures(data, selection.method = selection.method, verbose = verbose)
    data <- ScaleData(data, verbose = verbose)

  }

  data <- RunPCA(data, npcs = npcs, verbose = verbose, seed.use = 42)
  data
}
