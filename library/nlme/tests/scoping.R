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
invisible(lapply(list(gls, lme), function (FUN) {
    form <- follicles ~ 1
    stopifnot(identical(formula(FUN(form, Ovary)), form))
}))
## gave Error in eval(x$call$model) : object 'form' not found
