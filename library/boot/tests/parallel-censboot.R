## Reported by Duncan Murdoch (from StackOverflow) 2021-02-04

library(survival)
library(boot)
statMeanSurv <- function(data, var) {
    surv <- survfit(Surv(time, cens) ~ 1, data = data)
    mean(surv$surv) + var
}
## both have datasets aml
res <- censboot(boot::aml, statMeanSurv, R = 5L, var = 1,
                parallel = "multicore", ncpus = 2L)
res$t
