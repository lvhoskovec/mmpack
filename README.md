# mmpack
This folder contains code for the R package mmpack.

The main purpose of this package is to fit methods for estimationg the association between multipollutant mixtures and health outcomes using the methods nonparamtric Bayes shrinkage, Bayesian kernel machine regression, and Bayesian profile regression. This package relies on the following R packages. Install these packages in the R environment by using the install.packages("") command.

matrixcalc

mvtnorm

bkmr

PReMiuM

units

Then you can install mmpack by running the following lines in the R console:

library(devtools)

install_github( "lvhoskovec/mmpack", build_vignettes = TRUE)

library(mmpack)

vignette(mmpackTutorial")