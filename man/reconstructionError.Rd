% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{reconstructionError}
\alias{reconstructionError}
\title{Reject cells based on reconstruction error}
\usage{
reconstructionError(newData, spmodel, RE_threshold)
}
\arguments{
\item{newData}{Seurat object containing cells to annotate}

\item{spmodel}{scPred model}

\item{RE_threshold}{Threshold for the reconstruction error. If the reconstruction error of a cell is above this threshold
it is rejected. This threshold is determined during the training step.}
}
\value{
A boolean vector indicating which cells are rejected.
}
\description{
Recursive function to predict the labels of a dataset using a tree object
}
\examples{

}
\author{
Lieke Michielsen and Jose Alquicira-Hernandez
}
