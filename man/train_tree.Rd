% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/training.R
\name{train_tree}
\alias{train_tree}
\title{Create hierarchical tree}
\usage{
train_tree(
  data,
  tree,
  pvar = "cell_type",
  reduction = "pca",
  scaledata = FALSE,
  model = "svmRadial",
  allowParallel = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{data}{Seurat object containing cells used to train hierarchical models}

\item{tree}{Hierarchy created via \code{create_hierarchy}}

\item{pvar}{Prediction variable. Column name in Seurat object metadata
containing the cell labels of the terminal nodes of the hierarchy tree}

\item{reduction}{Name of reduction in Seurat objet to be used to determine
the feature space. Default: "pca"}

\item{scaledata}{Whether to select variable features and
scale the data of the reference object. Default: FALSE}

\item{model}{Classification model supported via caret package.
A list of all models can be found here in https://topepo.github.io/caret/available-models.html}

\item{allowParallel}{Allow parallel processing for resampling?}

\item{selection.method}{\code{vst}, \code{mean.var.plot}, or \code{dispersion}}

\item{reconstruction_error}{Use reconstruction error?}

\item{fn_perc}{TBD}
}
\value{
A data.tree object containing the cell type hierarchy
}
\description{
Creates a tree object representing a cell-type hierarchy
}
\examples{

branches <- c("Myeloid/DC/pDC", "Myeloid/DC/cDC")
h <- create_hierarchy(branches)
train_tree(data, h)

}
\author{
Lieke Michielsen and Jose Alquicira-Hernandez
}
