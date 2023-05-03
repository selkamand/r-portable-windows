times <- c(0, 1, 3, 7, 14, 28, 56, 90, 120)
k_in  <- c(0.1, 0.15, 0.16, 0.18, 0.2)
lrc_in <- mean(log(k_in))
lrc_sd_in <- sd(log(k_in))
n_sample <- length(times) * 2
const <- c(s1 = 1,    s2 = 1)
prop  <- c(s1 = 0.07, s2 = 0.03)
const_in <- rep(const, times = c(3, 2) * n_sample)
prop_in  <- rep(prop,  times = c(3, 2) * n_sample)

pred <- as.numeric(sapply(k_in, function(k) rep(100 * exp(- k * times), each = 2)))

set.seed(123456L)
d_syn <- data.frame(
    time  = rep(times, 5, each = 2),
    ds    = rep(paste0("d", 1:5), each = n_sample),
    study = rep(c("s1", "s2"), times = c(3 * n_sample, 2 * n_sample)),
    value = rnorm(length(pred), pred, sqrt(const_in^2 + pred^2 * prop_in^2))
)

library(nlme)
f_nlsList <- nlsList(value ~ SSasymp(time, 0, 100, lrc) | ds,
                     data = d_syn,
                     start = list(lrc = -3))

(fm_tc_study <-
    ## suppressWarnings( ## as the fit seems to be overparameterised
    nlme(f_nlsList,
         weights = varConstProp(form = ~ fitted(.) | study),
         control = list(sigma = 1))
    ## )
)

(ints <- intervals(fm_tc_study))

## Check if intervals include input used for data generation
stopifnot(exprs = {
  ints$fixed["lrc", "lower"] < lrc_in
  ints$fixed["lrc", "upper"] > lrc_in
  all.equal(c(
      ints$fixed["lrc", "est."],
      ints$reStruct$ds["sd(lrc)", "est."]),
    c(lrc_in, lrc_sd_in), tol = 1e-2) # diff. 0.00797 seen
  ints$varStruct[1:2, "lower"] < const
  ints$varStruct[3:4, "lower"] < prop
  ints$varStruct[1:2, "upper"] > const
  ints$varStruct[3:4, "upper"] > prop
  all.equal(
    as.numeric(ints$varStruct[c("prop.s1", "prop.s2"), "est."]),
    as.numeric(prop), tol =  0.15) # diff. 0.062  seen
})

## We do not get warnings if we fix the constant part of the error model
fm_tc_study_CF <-
    nlme(f_nlsList,
         weights = varConstProp(form = ~ fitted(.) | study,
                                fixed = list(const = c(s1 = 1, s2 = 1))))
summary(fm_tc_study_CF)
(ints_cf <- intervals(fm_tc_study_CF))
stopifnot(
  all.equal(
    as.numeric(ints_cf$varStruct[, "est."]),
    as.numeric(prop), tol =  0.15) # diff.  0.1168  seen
)
