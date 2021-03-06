% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simProfilesResponse.R
\name{simProfilesResponse}
\alias{simProfilesResponse}
\title{Function to simulate response data as a piece-wise constant function of predictor data (fixed profiles)}
\usage{
simProfilesResponse(X, W)
}
\arguments{
\item{X}{predictor data}

\item{W}{covariate data}
}
\value{
list with components
\itemize{ 
       \item Y: response data
       \item h: evaluation of exposure-response function for each individual 
       \item active: main effects selected to be in the exposure-response function
       \item active.ints: interactions selected to be in the exposure-response function
}
}
\description{
Function to simulate response data as a piece-wise constant function of predictor data (fixed profiles)
}
