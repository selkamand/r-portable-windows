#### Matrix Factorizations  --- of all kinds

## for R_DEFAULT_PACKAGES=NULL :
library(stats)
library(utils)

library(Matrix)
source(system.file("test-tools.R", package = "Matrix"))# identical3() etc
options(warn = 0)
is64bit <- .Machine$sizeof.pointer == 8
cat("doExtras:", doExtras,";  is64bit:", is64bit, "\n")

### "sparseQR" : Check consistency of methods
##   --------
data(KNex, package = "Matrix")
mm <- KNex$mm
y  <- KNex$y
stopifnot(is((Y <- Matrix(y)), "dgeMatrix"))
md <- as(mm, "matrix")                  # dense

(cS <- system.time(Sq <- qr(mm))) # 0.009
(cD <- system.time(Dq <- qr(md))) # 0.499 (lynne, 2014 f); 1.04 lynne 2019 ?????
cD[1] / cS[1] # dense is  much ( ~ 100--170 times) slower

## chkQR() in ../inst/test-tools-1.R ;

if(doExtras) { ## ~ 20 sec {"large" example}   + 2x qr.R() warnings
    cat("chkQR( <KNex> ) .. takes time .. ")
    system.time(chkQR(mm, y=y, a.qr = Sq, verbose=TRUE))
    system.time(chkQR(md, y=y, a.qr = Dq, verbose=TRUE))
    cat(" done: [Ok]\n")
}

## consistency of results dense and sparse
##	chk.qr.D.S() and  checkQR.DS.both() >>> ../inst/test-tools-Matrix.R
chk.qr.D.S(Dq, Sq, y, Y)

## Another small example with pivoting (and column name "mess"):
suppressWarnings(RNGversion("3.5.0")); set.seed(1)
X <- rsparsematrix(9,5, 1/4, dimnames=list(paste0("r", 1:9), LETTERS[1:5]))
qX <- qr(X); qd <- qr(as(X, "matrix"))
## are the same (now, *including* names):
assert.EQ(print(qr.coef(qX, 1:9)), qr.coef(qd, 1:9), tol=1e-14)
chk.qr.D.S(d. = qd, s. = qX, y = 1:9)



## rank deficient QR cases: ---------------

## From Simon (15 Jul 2009) + dimnames (11 May 2015)
set.seed(10)
a <- matrix(round(10 * runif(90)), 10,9, dimnames =
            list(LETTERS[1:10], paste0("c", 1:9)))
a[a < 7.5] <- 0
(A <- Matrix(a))# first column = all zeros
qD <- chkQR(a, giveRE=TRUE) ## using base qr
qS <- chkQR(A, giveRE=TRUE) ## using Matrix "sparse qr" -- "structurally rank deficient!
validObject(qS)# with the validity now (2012-11-18) -- ok, also for "bad" case
## Here, have illegal access  Up[-1] in ../src/cs.c
try( ## After patch (2016-10-04 - *NOT* committed), this fails
    ## definitely "fails" (with good singularity message) after c3194 (cs.c):
chk.qr.D.S(qD, qS, y = 10 + 1:nrow(A), force=TRUE)# 6 warnings: "structurally rank deficient"
)
try( ## NOTE: *Both* checks  currently fail here:
    chkQR(A, Qinv.chk=TRUE, QtQ.chk=TRUE)
)


## Larger Scale random testing
oo <- options(Matrix.quiet.qr.R = TRUE, Matrix.verbose = TRUE, nwarnings = 1e4)
set.seed(101)

quiet <- doExtras
for(N in 1:(if(doExtras) 1008 else 24)) {
    A <- rsparsematrix(8,5, nnz = rpois(1, lambda=16))
    cat(sprintf(if(quiet) "%d " else "%4d -", N)); if(quiet && N %% 50 == 0) cat("\n")
    checkQR.DS.both(A, Qinv.chk= NA, QtQ.chk=NA, quiet=quiet,
    ##                          --- => FALSE if struct. rank deficient
		    giveRE = FALSE, tol = 1e-12)
    ## with doExtras = TRUE, 64bit (F34, R 4.3.0-dev. 2022-05): seen 8.188e-13
}
summary(warnings())

## Look at single "hard" cases: --------------------------------------

## This is *REALLY* nice and small :
A0 <- new("dgCMatrix", Dim = 4:3, i = c(0:3, 3L), p = c(0L, 3:5), x = rep(1,5))
A0
checkQR.DS.both(A0, Qinv.chk = FALSE, QtQ.chk=FALSE)
##                                           ----- *both* still needed :
try( checkQR.DS.both(A0,  TRUE, FALSE) )
try( checkQR.DS.both(A0, FALSE,  TRUE) )

## and the same when dropping the first row  { --> 3 x 3 }:
A1 <- A0[-1 ,]
checkQR.DS.both(A1, Qinv.chk = FALSE, QtQ.chk=FALSE)
##                                           ----- *both* still needed :
try( checkQR.DS.both(A1,  TRUE, FALSE) )
try( checkQR.DS.both(A1, FALSE,  TRUE) )


qa <- qr(as(A0,"matrix"))
qA <- qr(A0) # -> message: ".. Matrix structurally rank deficient"

drop0(crossprod( Qd <- qr.Q(qa) ), 1e-15) # perfect = diag( 3 )
drop0(crossprod( Qs <- qr.Q(qA) ), 1e-15) # R[3,3] == 0 -- OOPS!
## OTOH, qr.R() is fine, as checked in the checkQR.DS.both(A0, *) above


## zero-row *and* zero-column :
(A2 <- new("dgCMatrix", i = c(0L, 1L, 4L, 7L, 5L, 2L, 4L)
           , p = c(0L, 3L, 4L, 4L, 5L, 7L)
           , Dim = c(8L, 5L)
           , x = c(0.92, 1.06, -1.74, 0.74, 0.19, -0.63, 0.68)))
checkQR.DS.both(A2, Qinv.chk = FALSE, QtQ.chk=FALSE)
##                                           ----- *both* still needed :
try( checkQR.DS.both(A2,  TRUE, FALSE) )
try( checkQR.DS.both(A2, FALSE,  TRUE) )


## Case of *NO* zero-row or zero-column:
(A3 <- new("dgCMatrix", Dim = 6:5
           , i = c(0L, 2L, 4L, 0L, 1L, 5L, 1L, 3L, 0L)
           , p = c(0L, 1L, 3L, 6L, 8L, 9L)
           , x = c(40, -54, -157, -28, 75, 166, 134, 3, -152)))
checkQR.DS.both(A3, Qinv.chk = FALSE, QtQ.chk=FALSE)
##                                           ----- *both* still needed :
try( checkQR.DS.both(A3,  TRUE, FALSE) )
try( checkQR.DS.both(A3, FALSE,  TRUE) )



(A4 <- new("dgCMatrix", Dim = c(7L, 5L)
           , i = c(1:2, 4L, 6L, 1L, 5L, 0:3, 0L, 2:4)
           , p = c(0L, 4L, 6L, 10L, 10L, 14L)
           , x = c(9, -8, 1, -9, 1, 10, -1, -2, 6, 14, 10, 2, 12, -9)))
checkQR.DS.both(A4, Qinv.chk = FALSE, QtQ.chk=FALSE)
##                                           ----- *both* still needed :
try( checkQR.DS.both(A4,  TRUE, FALSE) )
try( checkQR.DS.both(A4, FALSE,  TRUE) )

(A5 <- new("dgCMatrix", Dim = c(4L, 4L)
           , i = c(2L, 2L, 0:1, 0L, 2:3), p = c(0:2, 4L, 7L)
           , x = c(48, 242, 88, 18, -167, -179, 18)))
checkQR.DS.both(A5, Qinv.chk = FALSE, QtQ.chk=FALSE)
##                                           ----- *both* still needed :
try( checkQR.DS.both(A5,  TRUE, FALSE) )
try( checkQR.DS.both(A5, FALSE,  TRUE) )


quiet <- doExtras
for(N in 1:(if(doExtras) 2^12 else 128)) {
    A <- round(100*rsparsematrix(5,3, nnz = min(15,rpois(1, lambda=10))))
    if(any(apply(A, 2, function(x) all(x == 0)))) ## "column of all 0"
        next
    cat(sprintf(if(quiet) "%d " else "%4d -", N)); if(quiet && N %% 50 == 0) cat("\n")
    checkQR.DS.both(A, Qinv.chk=NA, giveRE=FALSE, tol = 1e-12, quiet = quiet)
    ##                         --- => FALSE if struct. rank deficient
}

summary(warnings())


options(oo)



### "denseLU"

## Testing expansions of factorizations {was ./expand.R, then in simple.R }
## new: [m x n]  where m and n  may differ
x. <- c(2^(0:5),9:1,-3:8, round(sqrt(0:16)))
set.seed(1)
for(nnn in 1:100) {
    y <- sample(x., replace=TRUE)
    m <- sample(2:6, 1)
    n <- sample(2:7, 1)
    x <- matrix(seq_len(m*n), m,n)
    lux <- lu(x)# occasionally a warning about exact singularity
    xx <- with(expand(lux), (P %*% L %*% U))
    print(dim(xx))
    assert.EQ.mat(xx, x, tol = 16*.Machine$double.eps)
}

### "sparseLU"
por1 <- readMM(system.file("external/pores_1.mtx", package = "Matrix"))
lu1 <- lu(por1)
pm <- as(por1, "CsparseMatrix")
(pmLU <- lu(pm)) # -> show(<MatrixFactorization>)
xp <- expand(pmLU)
## permute rows and columns of original matrix
ppm <- pm[pmLU@p + 1:1, pmLU@q + 1:1]
Ppm <- pmLU@L %*% pmLU@U
## identical only as long as we don't keep the original class info:
stopifnot(identical3(lu1, pmLU, pm@factors$sparseLU),# TODO === por1@factors$LU
	  identical(ppm, with(xp, P %*% pm %*% t(Q))),
	  sapply(xp, is, class2="Matrix"))

Ipm <- solve(pm, sparse=FALSE)
Spm <- solve(pm, sparse=TRUE)  # is not sparse at all, here
assert.EQ.Mat(Ipm, Spm, giveRE=TRUE, tol = 1e-13)# seen 7.36e-15 only on 32-bit
stopifnot(abs(as.vector(solve(Diagonal(30, x=10) %*% pm) / Ipm) - 1/10) < 1e-7,
	  abs(as.vector(solve(rep.int(4, 30)	  *  pm) / Ipm) - 1/ 4) < 1e-7)

## these two should be the same, and `are' in some ways:
assert.EQ.mat(ppm, as(Ppm, "matrix"), tol = 1e-14, giveRE=TRUE)
## *however*
length(ppm@x)# 180
length(Ppm@x)# 317 !
table(Ppm@x == 0)# (194, 123) - has 123 "zero" and 14 ``almost zero" entries

##-- determinant() and det() --- working via LU ---
m <- matrix(c(0, NA, 0, NA, NA, 0, 0, 0, 1), 3,3)
m0 <- rbind(0,cbind(0,m))
M <- as(m,"Matrix"); M ## "dsCMatrix" ...
M0 <- rbind(0, cbind(0, M))
dM  <- as(M, "denseMatrix")
dM0 <- as(M0,"denseMatrix")
try( lum  <- lu(M) )# Err: "near-singular A"
(lum  <- lu(M,  errSing=FALSE))# NA --- *BUT* it is not stored in @factors
(lum0 <- lu(M0, errSing=FALSE))# NA --- and it is stored in M0@factors[["LU"]]
## "FIXME" - TODO: Consider
replNA <- function(x, value) { x[is.na(x)] <- value ; x }
(EL.1 <- expand(lu.1 <- lu(M.1 <- replNA(M, -10))))
## so it's quite clear how  lu() of the *singular* matrix  M	should work
## but it's not supported by the C code in ../src/cs.c which errors out
stopifnot(all.equal(M.1,  with(EL.1, t(P) %*% L %*% U %*% Q)),
	  is.na(det(M)), is.na(det(dM)),
	  is.na(det(M0)), is.na(det(dM0)) )

###________ Cholesky() ________

##--------  LDL' ---- small exact examples

set.seed(1)
for(n in c(5:12)) {
    cat("\nn = ",n,"\n-------\n")
    rr <- mkLDL(n)
    ##    -------- from 'test-tools.R'
    stopifnot(all(with(rr, A ==
		       as(L %*% D %*% t(L), "symmetricMatrix"))),
	      all(with(rr, A == tcrossprod(L %*% sqrt(D)))))
    d <- rr$d.half
    A <- rr$A
    .A <- as(A, "TsparseMatrix") # 'factors' slot is retained => do chol() _after_ coercion
    R <- chol(A)
    assert.EQ.Mat(R, chol(.A)) # gave infinite recursion
    print(d. <- diag(R))
    D. <- Diagonal(x= d.^2)
    L. <- t(R) %*% Diagonal(x = 1/d.)
    stopifnot(all.equal(as.matrix(D.), as.matrix(rr$ D)),
              all.equal(as.matrix(L.), as.matrix(rr$ L)))
    ##
    CAp <- Cholesky(A)# perm=TRUE --> Permutation:
    validObject(CAp)
    p <- CAp@perm + 1L
    P <- as(p, "pMatrix")
    ## the inverse permutation:
    invP <- solve(P)@perm
    lDet <- sum(2* log(d))# the "true" value
    ldetp  <-         .diag.dsC(Chx = CAp, res.kind = "sumLog")
    ldetp. <- sum(log(.diag.dsC(Chx = CAp, res.kind = "diag") ))
    ##
    CA	<- Cholesky(A,perm=FALSE)
    validObject(CA)
    ldet <- .diag.dsC(Chx = CA, res.kind = "sumLog")
    ## not printing CAp : ends up non-integer for n >= 11
    mCAp <- as(CAp, "CsparseMatrix")
    print(mCA  <- drop0(as(CA, "CsparseMatrix")))
    stopifnot(identical(A[p,p], as(P %*% A %*% t(P),
				   "symmetricMatrix")),
	      relErr(d.^2, .diag.dsC(Chx= CA, res.kind="diag")) < 1e-14,
	      relErr(A[p,p], tcrossprod(mCAp)) < 1e-14)
    if(FALSE)
        rbind(lDet,ldet, ldetp, ldetp.)
    ##  ==> Empirically, I see  lDet = ldet  !=  ldetp == ldetp.
    ## if(rr$rcond.A < ...) warning("condition number of A ..." ## <- TODO
    cat(1,""); assert.EQ.(lDet,  ldet,   tol = 1e-14)
    cat(2,""); assert.EQ.(ldetp, ldetp., tol = 1e-14)
    cat(3,""); assert.EQ.(lDet,  ldetp, tol = n^2* 1e-7)# extreme: have seen 0.0011045 !!
}## for()

mkCholhash <- function(r.all) {
    ## r.all %*% (2^(2:0)), but only those that do not have NA / "?" :
    stopifnot(is.character(rn <- rownames(r.all)),
              is.matrix(r.all), is.logical(r.all))
    c.rn <- vapply(rn, function(ch) strsplit(ch, " ")[[1]], character(3))
    ## Now
    h1 <- function(i) {
        ok <- rep.int(TRUE, 3L)
        if(c.rn[3L, i] == "?")
            ok[2:3] <- FALSE # no supernodal LDL' factorization !!
        r.all[i, ok] %*% 2^((2:0)[ok])
    }
    vapply(seq_len(nrow(r.all)), h1, numeric(1))
}

set.seed(17)
(rr <- mkLDL(4))
(CA <- Cholesky(rr$A))
validObject(CA)
stopifnot(all.equal(determinant(rr$A) -> detA,
                    determinant(as(rr$A, "matrix"))),
          is.all.equal3(c(detA$modulus), log(det(rr$D)), sum(log(rr$D@x))))
A12 <- mkLDL(12, 1/10)
(r12 <- allCholesky(A12$A))[-1]
aCh.hash <- mkCholhash(r12$r.all)
if(requireNamespace("sfsmisc"))
split(rownames(r12$r.all), sfsmisc::Duplicated(aCh.hash))

## TODO: find cases for both choices when we leave it to CHOLMOD to choose
for(n in 1:50) { ## used to seg.fault at n = 10 !
    mkA <- mkLDL(1+rpois(1, 30), 1/10, rcond = FALSE, condest = FALSE)
    cat(sprintf("n = %3d, LDL-dim = %d x %d ", n, nrow(mkA$A), ncol(mkA$A)))
    r <- allCholesky(mkA$A, silentTry=TRUE)
    ## Compare .. apart from the NAs that happen from (perm=FALSE, super=TRUE)
    iNA <- apply(is.na(r$r.all), 1, any)
    cat(sprintf(" -> %3s NAs\n", if(any(iNA)) format(sum(iNA)) else "no"))
    stopifnot(aCh.hash[!iNA] == mkCholhash(r$r.all[!iNA,]))
##     cat("--------\n")
}


## This is a relatively small "critical example" :
A. <-
    new("dsCMatrix", Dim = c(25L, 25L), uplo = "U"
	, i = as.integer(
          c(0, 1, 2, 3, 4, 2, 5, 6, 0, 8, 8, 9, 3, 4, 10, 11, 6, 12, 13, 4,
            10, 14, 15, 1, 2, 5, 16, 17, 0, 7, 8, 18, 9, 19, 10, 11, 16, 20,
            0, 6, 7, 16, 17, 18, 20, 21, 6, 9, 12, 14, 19, 21, 22, 9, 11, 19,
            20, 22, 23, 1, 16, 24))
	##
	, p = c(0:6, 8:10, 12L, 15:16, 18:19, 22:23, 27:28, 32L, 34L, 38L, 46L, 53L, 59L, 62L)
	##
	, x = c(1, 1, 1, 1, 2, 100, 2, 40, 1, 2, 100, 6700, 100, 100, 13200,
	  1, 50, 4100, 1, 5, 400, 20, 1, 40, 100, 5600, 9100, 5000, 5,
	  100, 100, 5900, 100, 6200, 30, 20, 9, 2800, 1, 100, 8, 10, 8000,
	  100, 600, 23900, 30, 100, 2800, 50, 5000, 3100, 15100, 100, 10,
	  5600, 800, 4500, 5500, 7, 600, 18200))
validObject(A.)
## A1: the same pattern as  A.   just simply filled with '1's :
A1 <- A.; A1@x[] <- 1; A1@factors <- list()
A1.8 <- A1; diag(A1.8) <- 8
##
nT. <- as(AT <- as(A., "TsparseMatrix"),"nMatrix")
stopifnot(all(nT.@i <= nT.@j),
	  identical(qr(A1.8), qr(as(A1.8, "generalMatrix"))))

CA <- Cholesky(A. + Diagonal(x = rowSums(abs(A.)) + 1))
validObject(CA)
stopifnotValid(CAinv <- solve(CA), "dsCMatrix")
MA <- as(CA, "CsparseMatrix") # with a confusing warning -- FIXME!
stopifnotValid(MAinv <- solve(MA), "dtCMatrix")
## comparing MAinv with some solve(CA, system="...") .. *not* trivial? - TODO
##
CAinv2 <- solve(CA, Diagonal(nrow(A.)))
CAinv2 <- as(CAinv2, "symmetricMatrix")
stopifnot(identical(CAinv, CAinv2))

## FINALLY fix "TODO": (not implemented *symbolic* factorization of nMatrix)
try(    tc <- Cholesky(nT.)  )

for(p in c(FALSE,TRUE))
    for(L in c(FALSE,TRUE))
        for(s in c(FALSE,TRUE, NA)) {
            cat(sprintf("p,L,S = (%2d,%2d,%2d): ", p,L,s))
            r <- tryCatch(Cholesky(A., perm=p, LDL=L, super=s),
                          error = function(e)e)
            cat(if(inherits(r, "error")) " *** E ***" else
                sprintf("%3d", r@type),"\n", sep="")
        }
str(A., max.level=3) ## look at the 'factors'

facs <- A.@factors
names(facs) <- sub("Cholesky$", "", names(facs))
facs <- facs[order(names(facs))]

sapply(facs, class)
str(lapply(facs, slot, "type"))
## super = TRUE  currently always entails  LDL=FALSE :
## hence isLDL is TRUE for ("D" and not "S"):
sapply(facs, isLDL)

chkCholesky <- function(chmf, A) {
    stopifnot(is(chmf, "CHMfactor"),
              validObject(chmf),
              is(A, "Matrix"), isSymmetric(A))
    if(!is(A, "dsCMatrix"))
        A <- as(as(as(A, "CsparseMatrix"), "symmetricMatrix", "dMatrix"))
    L <- drop0(zapsmall(L. <- as(chmf, "CsparseMatrix")))
    cat("no. nonzeros in L {before / after drop0(zapsmall(.))}: ",
        c(nnzero(L.), nnzero(L)), "\n") ## 112, 95
    ecc <- expand(chmf)
    A... <- with(ecc, crossprod(crossprod(L,P)))
    stopifnot(all.equal(L., ecc$L, tolerance = 1e-14),
              all.equal(A,  A...,  tolerance = 1e-14))
    invisible(ecc)
}

c1.8 <- try(Cholesky(A1.8, super = TRUE))# works "always", interestingly ...
chkCholesky(c1.8, A1.8)



## --- now a "large" (712 x 712) real data example ---------------------------

data(KNex, package = "Matrix")
mtm <- with(KNex, crossprod(mm))
ld.3 <- determinant(Cholesky(mtm, perm = TRUE), sqrt = FALSE)
stopifnot(identical(names(mtm@factors),
                    "sPDCholesky"))
ld.4 <- determinant(Cholesky(mtm, perm = FALSE), sqrt = FALSE)
stopifnot(identical(names(mtm@factors),
                    c("sPDCholesky", "spDCholesky")))
c2 <- Cholesky(mtm, super = TRUE)
validObject(c2)
stopifnot(identical(names(mtm@factors),
                    c("sPDCholesky", "spDCholesky", "SPdCholesky")))

r <- allCholesky(mtm)
r[-1]

## is now taken from cache
c1 <- Cholesky(mtm)

bv <- 1:nrow(mtm) # even integer
b <- matrix(bv)
## solve(c2, b) by default solves  Ax = b, where A = c2'c2 !
x <- solve(c2,b)
stopifnot(identical3(drop(x), solve(c2, bv), drop(solve(c2, b, system = "A"))),
          all.equal(x, solve(mtm, b)))
for(sys in c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt")) {
    x <- solve(c2, b,  system = sys)
    cat(sys,":\n"); print(head(x))
    stopifnot(dim(x) == c(712, 1),
              identical(drop(x), solve(c2, bv, system = sys)))
}

## log(|LL'|) - check if super = TRUE and simplicial give same determinant
(ld.1 <- determinant(mtm))
if(FALSE) {
## MJ: CHMfactor_ldetL2 is unused outside of these tests, so we no longer
##     have it in the namespace { old definition is in ../src/CHMfactor.c }
ld1 <- .Call("CHMfactor_ldetL2", c1)
ld2 <- .Call("CHMfactor_ldetL2", c2)
stopifnot(all.equal(ld1, ld2),
	  all.equal(ld1, as.vector(ld.1$modulus), tolerance = 1e-14),
          all.equal(ld1, as.vector(ld.3$modulus), tolerance = 1e-14),
          all.equal(ld1, as.vector(ld.4$modulus), tolerance = 1e-14))
} else {
stopifnot(all.equal(as.vector(ld.1$modulus), as.vector(ld.3$modulus),
                    tolerance = 1e-14),
          all.equal(as.vector(ld.1$modulus), as.vector(ld.4$modulus),
                    tolerance = 1e-14))
}

## MJ: ldet[123].dsC() are unused outside of these tests, so we no longer
##     have them in the namespace { old definitions are in ../R/determinant.R }
if(FALSE) {
## Some timing measurements
mtm <- with(KNex, crossprod(mm))
I <- .symDiagonal(n=nrow(mtm))
set.seed(101); r <- runif(100)

system.time(D1 <- sapply(r, function(rho) Matrix:::ldet1.dsC(mtm + (1/rho) * I)))
## 0.842 on fast cmath-5
system.time(D2 <- sapply(r, function(rho) Matrix:::ldet2.dsC(mtm + (1/rho) * I)))
## 0.819
system.time(D3 <- sapply(r, function(rho) Matrix:::ldet3.dsC(mtm + (1/rho) * I)))
## 0.810
stopifnot(is.all.equal3(D1,D2,D3, tol = 1e-13))
}

## Updating LL'  should remain LL' and not become  LDL' :
cholCheck <- function(Ut, tol = 1e-12, super = FALSE, LDL = !super) {
    L <- Cholesky(UtU <- tcrossprod(Ut), super=super, LDL=LDL, Imult = 1)
    L1 <- update(L, UtU, mult = 1)
    L2 <- update(L, Ut,  mult = 1)
    stopifnot(is.all.equal3(L, L1, L2, tol = tol),
              all.equal(update(L, UtU, mult = pi),
                        update(L, Ut,  mult = pi), tolerance = tol)
              )
}

## Inspired by
## data(Dyestuff, package = "lme4")
## Zt <- as(Dyestuff$Batch, "sparseMatrix")
Zt <- new("dgCMatrix", Dim = c(6L, 30L), x = 2*1:30,
          i = rep(0:5, each=5),
          p = 0:30, Dimnames = list(LETTERS[1:6], NULL))
cholCheck(0.78 * Zt, tol=1e-14)

oo <- options(Matrix.quiet.qr.R = TRUE, warn = 2)# no warnings allowed
qrZ <- qr(t(Zt))
Rz <- qr.R(qrZ)
stopifnot(exprs = {
    inherits(qrZ, "sparseQR")
    inherits(Rz, "sparseMatrix")
    isTriangular(Rz)
    isDiagonal(Rz) # even though formally a "dtCMatrix"
    qr2rankMatrix(qrZ, do.warn=FALSE) == 6
})
options(oo)

## problematic rank deficient rankMatrix() case -- only seen in large cases ??
## MJ: NA in diag(<sparseQR>@R) not seen with Apple Clang 14.0.3
Z. <- readRDS(system.file("external", "Z_NA_rnk.rds", package="Matrix"))
(rnkZ. <- rankMatrix(Z., method = "qr")) # gave errors; now warns typically, but not on aarm64 (M1)
qrZ. <- qr(Z.)
options(warn=1)
rnk2 <- qr2rankMatrix(qrZ.) # warning ".. only 684 out of 822 finite diag(R) entries"
oo <- options(warn=2)# no warnings allowed from here
di.NA <- anyNA(diag(qrZ.@R))
stopifnot(is(qrZ, "sparseQR"),
          identical(is.na(rnkZ.), di.NA),
          identical(is.na(rnk2), di.NA))

## The above bug fix was partly wrongly extended to  dense matrices for "qr.R":
x <- cbind(1, rep(0:9, 18))
qr.R(qr(x))              # one negative diagonal
qr.R(qr(x, LAPACK=TRUE)) # two negative diagonals
chkRnk <- function(x, rnk) {
    stopifnot(exprs = {
        rankMatrix(x) == rnk
        rankMatrix(x, method="maybeGrad") == rnk ## but "useGrad" is not !
        rankMatrix(x, method="qrLINPACK") == rnk
        rankMatrix(x, method="qr.R"     ) == rnk
    })# the last gave '0' and a warning in Matrix 1.3-0
}
chkRnk(   x,    2)
chkRnk(diag(1), 1) # had "empty stopifnot" (-> Error in MM's experimental setup) +  warning 'min(<empty>)'
(m3 <- cbind(2, rbind(diag(pi, 2), 8)))
chkRnk(m3, 3)
chkRnk(matrix(0, 4,3), 0)
chkRnk(matrix(1, 5,5), 1) # had failed for "maybeGrad"
chkRnk(matrix(1, 5,2), 1)


showSys.time(
for(i in 1:120) {
    set.seed(i)
    M <- rspMat(n=rpois(1,50), m=rpois(1,20), density = 1/(4*rpois(1, 4)))
    cat(sprintf("%3d: dim(M) = %2dx%2d, rank=%2d, k=%9.4g; ",
		i, nrow(M), ncol(M), rankMatrix(M), kappa(M)))
    for(super in c(FALSE,TRUE)) {
        cat("super=",super,"M: ")
        ## 2018-01-04, Avi Adler: needed 1.2e-12 in Windows 64 (for i=55, l.1):
        cholCheck( M  , tol=2e-12, super=super); cat(" M': ")
        cholCheck(t(M), tol=2e-12, super=super)
    }
    cat(" [Ok]\n")
})

.updateCHMfactor
## TODO: (--> ../TODO "Cholesky"):
## ----
## allow Cholesky(A,..) when A is not symmetric *AND*
## we really want to factorize  AA' ( + beta * I)


## Schur() ----------------------
checkSchur <- function(A, SchurA = Schur(A), tol = 1e-14) {
    stopifnot(is(SchurA, "Schur"),
              isOrthogonal(Q <- SchurA@Q),
              all.equal(as.mat(A),
                        as.mat(Q %*% SchurA@T %*% t(Q)), tolerance = tol))
}

SH <- Schur(H5 <- Hilbert(5))
checkSchur(H5, SH)
checkSchur(Diagonal(x = 9:3))

p <- 4L
uTp <- new("dtpMatrix", x=c(2, 3, -1, 4:6, -2:1), Dim = c(p,p))
(uT <- as(uTp, "unpackedMatrix"))
## Schur ( <general> )  <--> Schur( <triangular> )
Su <- Schur(uT) ;   checkSchur(uT, Su)
gT <- as(uT,"generalMatrix")
Sg  <- Schur(gT) ;  checkSchur(gT, Sg)
Stg <- Schur(t(gT));checkSchur(t(gT), Stg)
Stu <- Schur(t(uT));checkSchur(t(uT), Stu)

stopifnot(exprs = {
    identical3(Sg@T, uT, Su@T)
    identical(Sg@Q, as(diag(p), "generalMatrix"))
    ## LaPck 3.12.0: these must be more careful (Q is *different* permutation):
    is.integer(print(ip <- invPerm(pp <- as(Stg@Q, "pMatrix")@perm)))
    identical(Stg@T, as(t(gT[,ip])[,ip], "triangularMatrix"))
    identical(Stg@Q, as(   diag(p)[,ip], "generalMatrix"))
    ## Stu still has p:1 permutation, but should not rely on it
    is.integer(print(i2 <- invPerm(as(Stu@Q, "pMatrix")@perm)))
    identical(Stu@T, as(t(uT[,i2])[,i2], "triangularMatrix"))
    identical(Stu@Q, as(   diag(p)[,i2], "pMatrix")) # Schur(<triangular>) ==> 'Q' is pMatrix
})


## the pedigreemm example where solve(.) failed:
p <- new("dtCMatrix", i = c(2L, 3L, 2L, 5L, 4L, 4:5), p = c(0L, 2L, 4:7, 7L),
	 Dim = c(6L, 6L), Dimnames = list(as.character(1:6), NULL),
	 x = rep.int(-0.5, 7), uplo = "L", diag = "U")
Sp <- Schur(p)
Sp. <- Schur(as(p,"generalMatrix"))
Sp.p <- Schur(crossprod(p))
## the last two failed
ip <- solve(p)
assert.EQ.mat(solve(ip), as(p,"matrix"))


## chol2inv() for a traditional matrix
assert.EQ.mat(     crossprod(chol2inv(chol(Diagonal(x = 5:1)))),
              C <- crossprod(chol2inv(chol(    diag(x = 5:1)))))
stopifnot(all.equal(C, diag((5:1)^-2)))
## failed in some versions because of a "wrong" implicit generic

U <- cbind(1:0, 2*(1:2))
(sU <- as(U, "CsparseMatrix"))
validObject(sS <- crossprod(sU))
C. <- chol(sS)
stopifnot(all.equal(C., sU, tolerance=1e-15))
## chol(<triangular sparse which is diagonal>)
tC7 <- .trDiagonal(7, 7:1)
stopifnotValid(tC7, "dtCMatrix")
ch7  <- chol(tC7) ## this (and the next 2) failed: 'no slot .. "factors" ..."dtCMatrix"'
chT7 <- chol(tT7 <- as(tC7, "TsparseMatrix"))
chR7 <- chol(tR7 <- as(tC7, "RsparseMatrix"))
stopifnot(expr = {
    isDiagonal(ch7)
    identical(chT7, ch7) # "ddiMatrix" all of them
    identical(chR7, ch7) # "ddiMatrix" all of them
    all.equal(sqrt(7:1), diag(ch7 ))
})



## From [Bug 14834] New: chol2inv *** caught segfault ***
n <- 1e6 # was 595362
A <- chol( D <- Diagonal(n) )
stopifnot(identical(A,D)) # A remains (unit)diagonal
is(tA <- as(A,"triangularMatrix"))# currently a dtTMatrix
stopifnotValid(tA, "dsparseMatrix")
CA <- as(tA, "CsparseMatrix")

selectMethod(solve, c("dtCMatrix","missing"))
##--> .Call(dtCMatrix_sparse_solve, a, .trDiagonal(n))  in ../src/dtCMatrix.c
sA  <- solve(CA)## -- R_CheckStack() segfault in Matrix <= 1.0-4
nca <- diagU2N(CA)
stopifnot(identical(sA, nca))
## same check with non-unit-diagonal D :
A <- chol(D <- Diagonal(n, x = 0.5))
ia <- chol2inv(A)
stopifnot(is(ia, "diagonalMatrix"),
	  all.equal(ia@x, rep(2,n), tolerance = 1e-15))

##------- Factor caches must be cleaned - even after scalar-Ops such as "2 *"
set.seed(7)
d <- 5
S <- 10*Diagonal(d) + rsparsematrix(d,d, 1/4)
class(M <- as(S, "denseMatrix")) # dgeMatrix
m <- as.matrix(M)
(dS <- determinant(S))
stopifnot(exprs = {
    all.equal(determinant(m), dS, tolerance=1e-15)
    all.equal(dS, determinant(M), tolerance=1e-15)
    ## These had failed, as the "LU" factor cache was kept unchanged in 2*M :
    all.equal(determinant(2*S), determinant(2*M) -> d2M)
    all.equal(determinant(S^2), determinant(M^2) -> dM2)
    all.equal(determinant(m^2), dM2)
    all.equal(d*log(2), c(d2M$modulus - dS$modulus))
})

## misc. bugs found in Matrix 1.4-1
L. <- new("dtCMatrix", Dim = c(1L, 1L), uplo = "L",
          p = c(0L, 1L), i = 0L, x = 1)
S. <- forceSymmetric(L.)
lu(S.)
stopifnot(validObject(lu(L.)), # was invalid
          identical(names(S.@factors), "sparseLU")) # was "lu"

## chol() should give matrix with 'Dimnames',
## even if 'Dimnames' are not cached
D. <- as(diag(3), "CsparseMatrix")
D.@Dimnames <- dn <- list(zzz = letters[1:3], ZZZ = LETTERS[1:3])
cd1 <- chol(D.) # "fresh"
stopifnot(identical(cd1@Dimnames, rep(dn[2L], 2L)))
cd2 <- chol(D.) # from cache
stopifnot(identical(cd1, cd2))

## lu(<m-by-0>), lu(<0-by-n>), BunchKaufman(<0-by-0>), chol(<0-by-0>)
stopifnot(identical(lu(new("dgeMatrix", Dim = c(2L, 0L))),
                    new("denseLU", Dim = c(2L, 0L))),
          identical(lu(new("dgeMatrix", Dim = c(0L, 2L))),
                    new("denseLU", Dim = c(0L, 2L))),
          identical(BunchKaufman(new("dsyMatrix", uplo = "U")),
                    new("BunchKaufman", uplo = "U")),
          identical(BunchKaufman(new("dspMatrix", uplo = "L")),
                    new("pBunchKaufman", uplo = "L")),
          identical(Cholesky(new("dpoMatrix", uplo = "U")),
                    new("Cholesky", uplo = "U")),
          identical(Cholesky(new("dppMatrix", uplo = "L")),
                    new("pCholesky", uplo = "L")))

## determinant(<ds[yp]Matrix>) going via Bunch-Kaufman
set.seed(15742)
n <- 10L
syU <- syL <- new("dsyMatrix", Dim = c(n, n), x = rnorm(n * n))
spU <- spL <- new("dspMatrix", Dim = c(n, n), x = rnorm((n * (n + 1L)) %/% 2L))
syL@uplo <- spL@uplo <- "L"
for(m in list(syU, syL, spU, spL))
    for(givelog in c(FALSE, TRUE))
        stopifnot(all.equal(determinant(   m,            givelog),
                            determinant(as(m, "matrix"), givelog)))

## was an error at least in Matrix 1.5-4 ...
BunchKaufman(as.matrix(1))


## 'expand2': product of listed factors should reproduce factorized matrix
## FIXME: many of our %*% methods still mangle dimnames or names(dimnames) ...
##        hence for now we coerce the factors to matrix before multiplying
chkMF <- function(X, Y, FUN, ...) {
    ## t(x)@factors may preserve factorizations with x@uplo
    X@factors <- list()

    mf <- FUN(X, ...)
    e2.mf <- expand2(mf)
    e1.mf <- sapply(names(e2.mf), expand1, x = mf, simplify = FALSE)

    m.e2.mf <- lapply(e2.mf, as, "matrix")
    m.e1.mf <- lapply(e1.mf, as, "matrix")

    identical(m.e1.mf, lapply(m.e2.mf, unname)) &&
        isTRUE(all.equal(Reduce(`%*%`, m.e2.mf), Y))
}
set.seed(24831)
n <- 16L
mS <- tcrossprod(matrix(rnorm(n * n), n, n,
                        dimnames = list(A = paste0("s", seq_len(n)), NULL)))
sS <- as(pS <- as(S <- as(mS, "dpoMatrix"), "packedMatrix"), "CsparseMatrix")
stopifnot(exprs = {
    chkMF(   S , mS,    Schur)
    chkMF(  pS , mS,    Schur)
    chkMF(   S , mS,       lu)
    chkMF(  pS , mS,       lu)
    chkMF(  sS , mS,       lu)
    chkMF(  sS , mS,       qr)
    chkMF(   S , mS, BunchKaufman)
    chkMF(  pS , mS, BunchKaufman)
    chkMF(t( S), mS, BunchKaufman)
    chkMF(t(pS), mS, BunchKaufman)
    chkMF(   S , mS, Cholesky)
    chkMF(  pS , mS, Cholesky)
    chkMF(t( S), mS, Cholesky)
    chkMF(t(pS), mS, Cholesky)
    chkMF(  sS , mS, Cholesky, super = FALSE, LDL =  TRUE)
    chkMF(  sS , mS, Cholesky, super = FALSE, LDL = FALSE)
    chkMF(  sS , mS, Cholesky, super =  TRUE, LDL = FALSE)
})


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
