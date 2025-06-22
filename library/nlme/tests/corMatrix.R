## test for PR#16110
library(nlme)
Orth <- subset(Orthodont, Subject %in% c("M01","F01"))
cs1CompSymm <- corCompSymm(value = 0.3, form = ~ 1 | Subject)
cs1CompSymm <- Initialize(cs1CompSymm, data = Orth)

## corFactor should return corMatrix(, corr = FALSE) as a vector
Linvt <- corMatrix(cs1CompSymm, corr = FALSE)
stopifnot(all.equal(unlist(Linvt, use.names = FALSE),
                    as.vector(corFactor(cs1CompSymm))))
## failed in 3.1-145 because the corFactor.corCompSymm method was
## misnamed corFactor.compSymm (a non-existent class), such that the
## general corFactor.corStruct method was called instead, which returned
## a different solution for the (transpose inverse) square-root factor


## test corMatrix() with spatial correlation *and* groups
coords <- cbind(x = (0:5)/10, y = (0:5)/10, z = (0:5)/10)
dists <- dist(coords, method = "manhattan")
spatDatGrouped <- data.frame(ID = factor(rep(c("A", "B"), each = nrow(coords))),
                             coords)
cs1Exp <- corExp(1, form = ~ x + y + z | ID, metric = "manhattan")
cs1Exp <- Initialize(cs1Exp, spatDatGrouped)
stopifnot(identical(getCovariate(cs1Exp),
                    list(A = c(dists), B = c(dists))))
cM <- corMatrix(cs1Exp)
stopifnot(exprs = {
    identical(cM[[1]], cM[[2]])
    cM[[1]][lower.tri(cM[[1]])] == exp(-dists)
})
