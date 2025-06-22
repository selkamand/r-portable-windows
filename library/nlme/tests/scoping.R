library("nlme")

## PR#18157
mygnls <- function (mydata)
    gnls(weight ~ SSlogis(Time, Asym, xmid, scal), data = mydata)
fm1 <- mygnls(Soybean) # failed in 3.1-153 with
## Error in stats::nls(formula = weight ~ SSlogis(Time, Asym, xmid, scal),  : 
##   object 'mydata' not found

## similarly, each of the following calls of
## nlme.formula(), nlsList.selfStart(), nlme.nlsList()
## using a self-starting model with local data would fail in 3.1-153 with
## Error in is.data.frame(data) : object 'mydata' not found
local({
    mydata <- subset(Loblolly, Seed < "307")
    fm2 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
                data = mydata, random = Asym ~ 1)
    fml <- nlsList(SSasymp, data = mydata)
    fm3 <- nlme(fml, random = Asym ~ 1)
})


## look for data in the parent frame, not in nlme's namespace
groupedData <- Orthodont
m3 <- lme(distance ~ age, data = groupedData, random = ~1 | Subject)
augPred(m3, length.out = 2)
## gave Error: data in 'm3' call must evaluate to a data frame
simulate(m3, m2 = list(random = ~ age | Subject), seed = 42, method = "ML")
## gave Error: 'data' must be a data.frame, environment, or list
rm(groupedData)


## PR#15892: formula.gls and formula.lme evaluated the call (in bad scope)
## same for predict.lme
invisible(lapply(list(gls, lme), function (FUN) {
    form <- follicles ~ 1
    model <- FUN(form, Ovary)
    stopifnot(identical(formula(model), form))
    stopifnot(all.equal(predict(model, newdata = Ovary[1,]),
                        fitted(model)[1], check.attributes = FALSE))
}))
## first gave Error in eval(x$call$model) : object 'form' not found
## second gave Error in eval(mCall$fixed) : object 'form' not found


## Subject: [Bug 18559] New: nlme -> anova doesn't handle symbolic formulas nicely
## Date: Sat, 08 Jul 2023 --- by Ben Bolker
formula <- height ~ a*exp(-b*age)
fm1 <- nlme(formula, data = Loblolly,
            fixed  = a + b ~ 1,
            random = a ~ 1,   start = c(a = 100, b = 1))
## "same" for gnls:
fmGnl <- gnls(formula, data = Loblolly, start = c(a = 100, b = 1))
stopifnot(exprs = {
    identical(formula(fm1), formula) ## was equal to stats::formula (!)
    identical(getResponseFormula(fm1), ~height)
    identical(formula(fmGnl), formula) # was stats::formula
    identical(getResponseFormula(fmGnl), ~height)
})

## similarly for self-starting models:
formula <- height ~ SSasymp(age, Asym, R0, lrc)
fm2 <- nlme(formula, data = Loblolly, fixed = Asym + R0 + lrc ~ 1, random = Asym ~ 1)
## nlme <= 3.1-164 failed with Error in x$formula :
##   object of type 'symbol' is not subsettable
stopifnot(identical(formula(fm2), formula))
