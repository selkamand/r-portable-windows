## IGNORE_RDIFF_BEGIN
.libPaths() # show full library tree {also as check of R CMD check!}
## IGNORE_RDIFF_END
library(cluster)

####---------- Tests for FANNY  i.e., fanny() --------------------------
####
### -- thanks to ../.Rbuildignore , the output of this is
### -- only compared to saved values for the maintainer

###--- An extension of  example(fanny) : -------------------
set.seed(21)
## generate 10+15 objects in two clusters, plus 3 objects lying
## between those clusters.
x <- rbind(cbind(rnorm(10, 0, 0.5), rnorm(10, 0, 0.5)),
           cbind(rnorm(15, 5, 0.5), rnorm(15, 5, 0.5)),
           cbind(rnorm( 3,3.2,0.5), rnorm( 3,3.2,0.5)))

.proctime00 <- proc.time()

(fannyx <- fanny(x, 2))
summary(fannyx)
str(fannyx)
## Different platforms differ (even gcc 3.0.1 vs 3.2 on same platform)!
## {70 or 71 iterations}
## ==> No "fanny-ex.Rout.save" is distributed !
## --------------------------------------------
## IGNORE_RDIFF_BEGIN
summary(fanny(x,3))# one extra cluster
## IGNORE_RDIFF_END

## CRAN-relevant M1 mac: aarch64-apple-darwin24.1.0 / Apple clang version 16.0.0 (clang-1600.0.26.6) / macOS Sequoia 15.1.1
## IGNORE_RDIFF_BEGIN
(fanny(x,2, memb.exp = 1.5))
## IGNORE_RDIFF_END
(fanny(x,2, memb.exp = 1.2))
(fanny(x,2, memb.exp = 1.1))
(fanny(x,2, memb.exp = 3))
## for subsetting, when compairing:
not.conv <- setdiff(names(fannyx),               c("convergence", "call"))
notMconv <- setdiff(names(fannyx), c("membership", "convergence", "call"))

data(ruspini) # < to run under R 1.9.1
summary(fanny(ruspini, 3), digits = 9)
summary(fanny(ruspini, 4), digits = 9)# 'correct' #{clusters}
summary(fanny(ruspini, 5), digits = 9)

cat('Time elapsed: ', proc.time() - .proctime00,'\n')
data(chorSub)
p4cl <- pam(chorSub, k = 4, cluster.only = TRUE)
## The first two are "completely fuzzy" -- and now give a warnings
## IGNORE_RDIFF_BEGIN
f4.20  <- fanny(chorSub, k = 4, trace.lev = 1)
f4.20$coef
## IGNORE_RDIFF_END
stopifnot(exprs = {
    all.equal(f4.20$coef, c(dunn_coeff = 0.25, normalized = 3.330669e-15))
    all.equal(f4.20$objective[["objective"]], 2665.982, tol = 8e-7)
    all.equal(f4.20$silinfo$avg.width,        0.250643, tol = 2e-6)
})

f4.18  <- fanny(chorSub, k = 4,   memb.exp = 1.8) # same problem
f4.18. <- fanny(chorSub, k = 4,   memb.exp = 1.8,
                iniMem.p = f4.20$membership) # very quick convergence
stopifnot(all.equal(f4.18[not.conv], f4.18.[not.conv], tol = 5e-7))

f4.16  <- fanny(chorSub, k = 4, memb.exp = 1.6) # now gives 4 crisp clusters
## IGNORE_RDIFF_BEGIN
f4.16. <- fanny(chorSub, k = 4, memb.exp = 1.6,
                iniMem.p = f4.18$membership, trace.lev = 2)# wrongly "converged" immediately; no longer!
f4.16.2<- fanny(chorSub, k = 4, memb.exp = 1.6,
                iniMem.p = cluster:::as.membership(p4cl), tol = 1e-10, trace.lev = 2)
all.equal((m1 <- f4.16  $membership),
          (m2 <- f4.16.2$membership))
## IGNORE_RDIFF_END
stopifnot(identical(dimnames(m1), dimnames(m2)),
          0 < m1,m1 < 1,  0 < m2,m2 < 1,
          ## the memberships are quite close but have only converged to precision 0.000228
          all.equal(m1, m2, tol = 0.001))
stopifnot(exprs = {
    f4.16$clustering == f4.16.2$clustering
    all.equal(f4.16[notMconv],  f4.16.2[notMconv],  tol = 1e-7)
})

f4.14 <- fanny(chorSub, k = 4,                memb.exp = 1.4)
f4.12 <- fanny(chorSub, k = 4,                memb.exp = 1.2)

table(f4.12$clustering, f4.14$clustering)# close but different
table(f4.16$clustering, f4.14$clustering)# ditto
table(f4.12$clustering, f4.16$clustering)# hence differ even more

symnum(cbind(f4.16$membership, 1, f4.12$membership),
       cutpoints= c(0., 0.2, 0.6, 0.8, 0.9, 0.95, 1 -1e-7, 1 +1e-7),
       symbols   = c(" ", ".", ",", "+",  "*", "B","1"))


## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
