% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/naive_bayes.R
\name{naive_bayes}
\alias{naive_bayes}
\title{naive_bayes  helper function}
\usage{
naive_bayes(reference_dataset, celltype_labels)
}
\arguments{
\item{reference_dataset}{PC-transformed dataset in matrix form}

\item{celltype_labels}{Celltype labels in vector}
}
\value{
naive_bayes model trained on evenly split dataset (upsamples if classes aren't evenly split)
}
\description{
naive_bayes  helper function
}
\examples{
data("iris")
dataset = iris[iris$Species \%in\% c("setosa","versicolor"),]
labels = as.factor(dataset$Species)
labels = droplevels(labels)
dataset = dataset[,1:4]
ex_model = naive_bayes(dataset, labels)
}
\keyword{internal}
