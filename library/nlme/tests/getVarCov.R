## tests for PR#16744 (failed in nlme <= 3.1-149)
library("nlme")
fm3 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare),
           weights = ~as.integer(as.character(Mare)))  # fixed variance weights
stopifnot(
    identical(
        getVarCov(fm3, individual = 3),
        getVarCov(fm3, individual = levels(Ovary$Mare)[3])
    ),
    all.equal(
        vapply(as.character(1:11),
               function (id) getVarCov(fm3, individual = id)[1,1],
               0, USE.NAMES = FALSE),
        fm3$sigma^2 * (1:11)
    )
)

## lme method had a similar issue for data not ordered by levels(group)
## is.unsorted(Orthodont$Subject)  # TRUE
fm4 <- lme(distance ~ age, Orthodont, weights = ~as.integer(Subject))
covmats <- getVarCov(fm4, individuals = levels(Orthodont$Subject),
                     type = "conditional")
stopifnot(
    all.equal(
        vapply(covmats, "[", 0, 1, 1),
        fm4$sigma^2 * seq_len(nlevels(Orthodont$Subject)),
        check.attributes = FALSE
    )
)
