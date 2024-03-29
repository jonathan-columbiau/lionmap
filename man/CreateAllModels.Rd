% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateAllModels.R
\name{CreateAllModels}
\alias{CreateAllModels}
\title{A helper function that takes in a ge reference dataset that has been transformed
using PCA, and a set of celltype labels in factor format, and produces a set of
models used to differentiate celltypes.}
\usage{
CreateAllModels(reference_dataset, celltype_labels, models_to_include = NULL)
}
\arguments{
\item{reference_dataset}{PC-transformed dataset}

\item{celltype_labels}{Vector of celltype labels}

\item{models_to_include}{Optional vector giving model names to include. If include
this parameter, include a subset of the following (make sure the names match or it won't work):
"linear_svm", "polynomial_svm", "naive_bayes", "ridge", "lasso", "elastic_net",
"linear_da", "knn", "rf", "quadratic_da"}
}
\value{
List of Models
}
\description{
A helper function that takes in a ge reference dataset that has been transformed
using PCA, and a set of celltype labels in factor format, and produces a set of
models used to differentiate celltypes.
}
\examples{
CreateAllModels(reference_dataset, celltype_labels, models_to_include = NULL)
}
\keyword{internal}
