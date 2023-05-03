set.seed(123)
N <- 9522
rM <- matrix(rnorm(N*4), N)
cvR <- try(MASS::cov.rob(rM, method = "mcd", nsamp = "exact"))
cvR <- try(MASS::cov.rob(rM, method = "mcd", nsamp = 1e10))

