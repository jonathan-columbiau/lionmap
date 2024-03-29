% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binomial_lasso.R
\name{binomial_lasso}
\alias{binomial_lasso}
\title{Binomial lasso helper function}
\usage{
binomial_lasso(reference_dataset, celltype_labels)
}
\arguments{
\item{reference_dataset}{PC-transformed dataset in matrix form}

\item{celltype_labels}{Celltype labels in vector}
}
\value{
Binomial lasso model trained on evenly split dataset (upsamples if classes aren't evenly split)
}
\description{
Binomial lasso helper function
}
\examples{
data("iris")
dataset = iris[iris$Species \%in\% c("setosa","versicolor"),]
labels = as.factor(dataset$Species)
labels = droplevels(labels)
dataset = dataset[,1:4]
ex_model = binomial_lasso(dataset, labels)
}
\keyword{internal}
