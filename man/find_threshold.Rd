% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconstruction_error.R
\name{find_threshold}
\alias{find_threshold}
\title{Find threshold for reconstruction error}
\usage{
find_threshold(
  object,
  labels,
  nfolds = 5,
  fn_perc = 0.01,
  seed = 66,
  scaledata = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{object}{Seurat object used as reference}

\item{labels}{Response variable}

\item{nfolds}{Number of resamples}

\item{fn_perc}{TBD}

\item{seed}{Numeric seed to create resamples}

\item{verbose}{Display Seurat messages?}
}
\value{
A numeric value corresponding to the threshold for reconstruction
error
}
\description{
Determines threshold for reconstruction errors computed from
multiple resamples
}
\examples{

library(scPred)
data(pbmc_1)
res <- find_threshold(pbmc_1, pbmc_1$cell_type)

}
\author{
Lieke Michielsen and Jose Alquicira-Hernandez
}
