% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.bpr.R
\name{predict.bpr}
\alias{predict.bpr}
\title{Predict method for class "bpr"}
\usage{
\method{predict}{bpr}(object, ...)
}
\arguments{
\item{object}{object of class "bpr"}

\item{...}{ignored}
}
\value{
list with components
\itemize{
     \item fitted.vals: posterior mean fitted values for each subject
     \item fitted.distn: posterior distribution of fitted values for each subject  
}
}
\description{
Predict method for class "bpr"
}
