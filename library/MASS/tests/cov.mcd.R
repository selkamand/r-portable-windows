library(MASS)
set.seed(123)
N <- 9522
rM <- matrix(rnorm(N*4), N)
cvR <- try(cov.rob(rM, method = "mcd", nsamp = "exact"))
cvR <- try(cov.rob(rM, method = "mcd", nsamp = 1e10))

## Things which failed with a gcc12 mis-compile.

## Extracted from package fpc
xg <- structure(c(4, 5, 6, 7, 8, 8, 7, 6, 5, 8), dim = c(5L, 2L),
                dimnames = list(NULL, c("x", "y")))
cov.mcd(xg)
