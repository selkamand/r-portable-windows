## from Peter Dalgaard 2020-11-24
## '(without the truncate=FALSE, curves pretty much go through
##   the penultimate point, whatever reasonable bandwith is chosen.)'

library(KernSmooth)
if(require("carData")) {
plot(prestige ~ income, data = Prestige)
with(Prestige, lines(locpoly(income, prestige, bandwidth = 5000)))
with(Prestige, lines(locpoly(income, prestige,
                             bandwidth = 5000, truncate = FALSE)))
}
