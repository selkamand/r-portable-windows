## PR#16744: ordering of variance weights (failed in nlme <= 3.1-149)
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


## PR#16806: 1-observation groups in corSpatial fits (failed in nlme <= 3.1-164)
data("Phenobarb")
pheno <- subset(Phenobarb, !is.na(conc))
stopifnot(sum(getGroups(pheno) == "28") == 1) # Subject 28 has only 1 obs.
lme1 <- lme(conc ~ time, data = pheno, random = ~1 | Subject,
            correlation = corExp(form = ~time | Subject))
condVarCov <- getVarCov(lme1, type = "conditional", individuals = "28")
## gave Error in dimnames(cond.var) <- list(1:nrow(cond.var), 1:ncol(cond.var)) :
##   length of 'dimnames' [2] not equal to array extent
stopifnot(all.equal(unname(diag(condVarCov[[1]])), lme1$sigma^2))
## same for gls():
gls1 <- gls(conc ~ time, data = pheno,
            correlation = corExp(form = ~time | Subject))
margVarCov <- getVarCov(gls1, individual = "28")
## gave Error in rep(1, nrow(S)) : invalid 'times' argument
stopifnot(all.equal(diag(margVarCov), gls1$sigma^2))
