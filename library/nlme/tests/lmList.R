library(nlme)

## PR#13788
qm <-  lmList(height ~ age | Subject, data = Oxboys)
nd <- with(Oxboys,
           expand.grid(age = seq(min(age),max(age),length=50),
                       Subject = levels(Subject))
           )

## failed in 3.1-92
res <- predict(qm, nd, se=TRUE)
stopifnot(is.data.frame(res), dim(res) == c(1300, 3),
          identical(names(res), c("Subject", "fit", "se.fit")))


## plots of ranef() and intervals.lmList() with new arguments 'xlab', 'ylab'
req <- ranef(qm)
(p.req <- plot(req, xlab = "R.Eff.", ylab = "Subj"))
# Fails (p.re2 <- plot(req, age ~ fitted(.)))
iqm <- intervals(qm)
stopifnot(is.array(iqm), dim(iqm) == c(26,3,2))
p.iq <- plot(iqm, ylab = "Subject [factor]")
## Error: formal argument "ylab" matched by multiple .. in 3.1.137
stopifnot(inherits(p.iq, "trellis"),
          inherits(p.req, "trellis"),
          identical( unclass(p.req)[c("xlab","ylab")],
                    list(xlab = "R.Eff.", ylab = "Subj")),
          formula(p.iq) == (group ~ intervals | what))
p.iq


## PR#16542: summary.lmList() with NA coefs in subgroup
lms <- lmList(pixel ~ day + I(day^2) | Dog, data = Pixel)
coefs <- coef(lms)
stopifnot(is.na(coefs["9", "I(day^2)"]),
          identical(dim(coefs), c(10L, 3L)))
summary(lms)  # failed in nlme <= 3.1-155 with
## Error in `[<-`(`*tmp*`, use, use, ii, value = lst[[ii]]) : 
##   subscript out of bounds

## same bug: unused factor levels in subgroup
lms2 <- lmList(pixel ~ Dog + day | Side, data = Pixel,
               subset = !(Side == "R" & Dog == "1") & Dog %in% 1:3)
coef(lms2) # failed in nlme <= 3.1-155 with
## Error in `[<-`(`*tmp*`, i, names(coefs[[i]]), value = if (is.null(coefs[[i]])) { : 
##   subscript out of bounds
