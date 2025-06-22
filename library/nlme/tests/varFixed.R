library(nlme)

## PR#17712: nlme() with a fixed variance function
set.seed(42)
d <- data.frame(
    obs = rnorm(137), # prime, so recycling would fail
    group = factor(rep_len(c("A", "B"), 137))
)
fit_wt <- function (weights = NULL)
    nlme(obs ~ b, fixed = b ~ 1, random = b ~ 1 | group, start = c(b = 0),
         data = d, weights = weights,
         control = nlmeControl(returnObject = TRUE)) # possibly small PNLS steps

m_unweighted <- fit_wt(NULL)
stopifnot(all.equal(fixef(m_unweighted), c(b = mean(d$obs)), tolerance = 0.05))
## equivalent fit with a constant variance function:
m_varIdent <- fit_wt(varIdent())
## nlme <= 3.1-164 failed with Error in recalc.varFunc(object[[i]], conLin) :
##   dims [product 12] do not match the length of object [13]
stopifnot(all.equal(logLik(m_unweighted), logLik(m_varIdent)))

fit0 <- fit_wt(varExp(form = ~obs))
## equivalent fit with the coefficient **fixed** at its estimate:
vf0 <- fit0$modelStruct$varStruct
fit1 <- fit_wt(varExp(form = ~obs, fixed = unname(coef(vf0))))
## fitting failed in nlme <= 3.1-164 (as above)
stopifnot(all.equal(logLik(fit0), logLik(fit1),
                    check.attributes = FALSE)) # df differs by 1
## equivalent fit using the same variance weights via **varFixed**:
##d$wt <- varWeights(vf0)  # beware of different (internal) ordering
d$wt <- fit0$sigma / attr(fit0$residuals, "std")
fit2 <- fit_wt(varFixed(~1/wt^2))
## fitting failed in nlme <= 3.1-164 (as above)
stopifnot(all.equal(logLik(fit0), logLik(fit2),
                    check.attributes = FALSE)) # df differs by 1
