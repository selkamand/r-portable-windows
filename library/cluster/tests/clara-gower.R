quit("no") ## retracted metric = "gower":  bad memory leakage  -- see  ./clara-gower_valgrind.Rout  <<<<<

## Originally inspired by  Kasper Fischer-Rasmussen 's  clara_gower.html  [html from Rmd]

library(cluster)
packageDescription("cluster")

## carefully getting  150 + 200 + 150 = 500 obs. from the 3  xclara clusters :
str(dd <- xclara[c(1:150, 1001:1200, 2101:2250), ])
dim(dd) # 500 2

set.seed(47)
cl_manhat <- clara(dd, 3, metric = "manhattan", rngR=TRUE, pamLike=TRUE, samples = 500)
cl_gower  <- clara(dd, 3, metric = "gower",     rngR=TRUE, pamLike=TRUE, samples = 500)

table(cl_manhat$cluster,
      cl_gower $cluster)

stopifnot(exprs = {
    ## Apart from [188], they are the same
    ##    usually even *including* [188], but not always ???? {FIXME ??? Random? even we use rngR?}
    cl_manhat$cluster[-188] == cl_gower $cluster[-188]
    identical(rle(unname(cl_manhat$cluster)),
              structure(class = "rle",
                        list(lengths = c(29L, 1L, 120L, 80L, 1L, 119L, 150L),
                             values  = c( 1L, 2L,   1L,  2L, 1L,   2L,   3L))))
})
## ==> no distinction between the clusters wrt Manhattan vs. Gower's distance.


## Using {cluster}'s built in tools to compute Gower's distance.

cl_gower_full <- clara(dd, k = 3, metric = "gower", rngR = TRUE, pamLike = TRUE, samples = 500, sampsize = nrow(dd))
dist_cl_full <- as.matrix(cl_gower_full$diss)
i_full <- rownames(dist_cl_full)
d_full <- data.frame(CLARA = as.vector(cl_gower_full$diss),
                     DAISY = as.vector(daisy(dd[i_full, ], metric = "gower")))

## MM: instead of all this, just
all.equal(d_full$CLARA,
          d_full$DAISY, tol=0) # "Mean relative difference: 2.17e-16"
## ... but sometimes *VERY* different (relative diff.   0.5xxx)
if(FALSE)
stopifnot( all.equal(d_full$CLARA,
                     d_full$DAISY, tol = 1e-15) ) ## equal up to  15 digits!

## We can see that the distance measurements are exactly identical when the
## whole data is used in the clustering. This is because the Gower distance
## scales the distances measurements with the range of each feature. Due to
## the subsampling, approximate ranges are calculated based on each
## subsample explaining the deviations.


## MM: compare -- with pam():
dGow <- daisy(dd, metric="gower")
cl_full <- clara(dd, k = 3, metric = "gower", rngR = TRUE, pamLike = TRUE, samples = 1, sampsize = nrow(dd))

all.equal(c(dGow) , c(cl_full$diss), tol=0) # "Mean relative difference: 2.171402e-16"

pam_3 <- pam(dGow, k = 3, variant = "faster")
## FIXME !! -- bug !?
all.equal(pam_3  $ clustering, # we would want *identical* -- bug ??
          cl_full$ clustering)
all.equal(c(dGow) , c(cl_full$diss), tol = 1e-15)
if(FALSE) ## FIXME
stopifnot(exprs = {
    identical(pam_3  $ clustering,
              cl_full$ clustering)
    all.equal(c(dGow) , c(cl_full$diss), tol = 1e-15)
})

