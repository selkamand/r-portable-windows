library("nlme")

## PR#18777: corARMA() with p or q exceeding "maxLag"
dat3 <- data.frame(ID = gl(5, 3), visit = 0:2) # => max time diff is 2
csARMA22 <- corARMA(form = ~ visit | ID, p = 2, q = 2)
Initialize(csARMA22, dat3) # ok
csAR3 <- corARMA(form = ~ visit | ID, p = 3)
csMA3 <- corARMA(form = ~ visit | ID, q = 3)
tools::assertError(Initialize(csAR3, dat3), verbose = TRUE)
tools::assertError(Initialize(csMA3, dat3), verbose = TRUE)
## both crashed with "free(): invalid next size (fast)" in nlme <= 3.1-165
