### example of glmmPQL with offset.

### glmmPQL age structured sole egg data example
### from "GAMs:An Intro with R" by Simon Wood

library(MASS)
if(!requireNamespace("nlme", quietly = TRUE)) q("no")

options(warn = 2L)
## Chapter 2 stuff...

if(FALSE) { ## script to create solr
#require(gamair) # for 'data(sole)'
data(sole)
sole$off <- log(sole$a.1 - sole$a.0) # model offset term
sole$a <-(sole$a.1 + sole$a.0)/2     # mean stage age
solr <- sole                         # make copy for rescaling
solr$t <- solr$t-mean(sole$t)
solr$t <- solr$t/var(sole$t)^0.5
solr$la <- solr$la - mean(sole$la)
solr$lo <- solr$lo - mean(sole$lo)
save(solr, file = "solr.rda", version = 2)
} else load("solr.rda")

solr$station <- factor(with(solr, paste0(-la, -lo, -t)))
b <- glmmPQL(eggs ~ offset(off) + lo + la + t + I(lo*la) + I(lo^2) +
            I(la^2) + I(t^2) + I(lo*t) + I(la*t) + I(lo^3) + I(la^3) +
            I(t^3) + I(lo*la*t) + I(lo^2*la) + I(lo*la^2) + I(lo^2*t) +
            I(la^2*t) + I(la*t^2) + I(lo*t^2) + a + I(a*t) + I(t^2*a),
            random = list(station = ~1),
            family = quasi(link = log, variance = "mu"),
            data = solr)
summary(b)

b1 <- update(b, ~ . - I(lo*la*t))
summary(b1)

b2 <- update(b1, ~ . - I(lo*t))
summary(b2)

b3 <- update(b2, ~ . - I(lo^2*t))
summary(b3)

b4 <- update(b3, ~ . - I(la*t^2))
summary(b4)

pdf("glmmPQL.pdf")
fv <- exp(fitted(b4))    # fitted values include offset from 3.1-158
resid <- solr$eggs - fv  # raw residuals
sqrt_fv <- sqrt(fv)
plot(sqrt_fv, sqrt(solr$eggs))
abline(0, 1, lwd = 2)
plot(sqrt_fv, resid/sqrt_fv)
plot(sqrt_fv, resid)
fl <- sort(sqrt_fv)
## add 1 s.d. and 2 s.d. reference lines
lines(fl, fl); lines(fl, -fl); lines(fl, 2 * fl, lty = 2)
lines(fl, -2 * fl, lty = 2)

nlme::intervals(b4, which = "var-cov")

head(p1 <- predict(b))
head(p2 <- predict(b, solr))
all.equal(p1, p2, check.attributes = FALSE)
