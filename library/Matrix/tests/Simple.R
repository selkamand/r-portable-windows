#### Currently a collection of simple tests
##	(since 'Matrix' takes long to load, rather have fewer source files!)

## for R_DEFAULT_PACKAGES=NULL :
library(methods)
library(stats)
library(utils)

##-------- *BEFORE* attaching Matrix: --------------------------------
str(Matrix::Matrix)# -> load the namespace
T <- new("ngTMatrix", i=0L, j=2L, Dim = c(2L,6L))
T
as(T, "CsparseMatrix")
## gave Error in asMethod(object) : could not find function ".M.classEnv"
## from  0.999375-23 to *-25

## another even shorter version of this:
n <- new("dgCMatrix")
n
## this:
m <- Matrix::Matrix(cbind(1,0,diag(x=2:4)))
m
mt <- m + table(gl(3,5), gl(5,3))# failed in Matrix <= 1.2.9
mt
stopifnot(is(mt, "sparseMatrix"))

##--------------------------------------------------------------------

library(Matrix)
source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

if(interactive()) {
    options(error = recover, Matrix.verbose = TRUE, warn = 1)
} else options(              Matrix.verbose = TRUE, warn = 1)
                                        #   ^^^^^^ to show Matrix.msg()s

(Mv <- Matrix.Version())
stopifnot(identical(packageVersion("Matrix"), Mv[["package"]]))

### Matrix() ''smartness''
(d40 <- Matrix( diag(4)))
(z4 <- Matrix(0*diag(4)))
(o4 <- Matrix(1+diag(4)))
(tr <- Matrix(cbind(1,0:1)))
(M4 <- Matrix(m4 <- cbind(0,rbind(6*diag(3),0))))
dM4 <- Matrix(M4, sparse = FALSE)
d4. <- diag(4); dimnames(d4.) <- dns <- rep(list(LETTERS[1:4]), 2)
d4a <- diag(4); dimnames(d4a) <- dna <- list(LETTERS[1:4], letters[1:4])# "a"symmetric
m1a <- matrix(0, dimnames=list("A","b"))# "a"symmetric
d4di<- as(d4., "diagonalMatrix")
d4da<- as(d4a, "diagonalMatrix")
d4d <- as(d4., "denseMatrix")
d4aS <- Matrix(d4a, sparse=TRUE, doDiag=FALSE)
d1aS <- Matrix(m1a, sparse=TRUE, doDiag=FALSE)
stopifnot(exprs = {
    identical(d4di@x, numeric()) # was "named" unnecessarily
    identical(dimnames(d4 <- Matrix(d4.)), dns)
    identical4(d40, Matrix(diag(4)), unname(d4), unname(d4da))
    identical3(d4, as(d4., "Matrix"), as(d4., "diagonalMatrix"))
    is(d4aS, "dtCMatrix") # not "dsC*", as asymmetric dimnames
    is(d4d, "denseMatrix")
    identical(dimnames(d4T <- as(d4., "TsparseMatrix")), dns) # failed till 2019-09-xx
    ## identical(d4T, as(d4., "dsTMatrix")) # deprecated
})
showProc.time()

class(mN <-  Matrix(NA, 3,4)) # NA *is* logical
validObject(Matrix(NA))
bd4 <- bdiag(M4,dM4,M4)
stopifnotValid(o4, "dsyMatrix")
stopifnotValid(M4, "dtCMatrix")
stopifnot(validObject(dM4), validObject(mN),
          identical(bdiag(M4), bdiag(dM4)),
          identical(bd4@p, c(0L,0:3,3:6,6:9)),
          identical(bd4@i, c(0:2, 4:6, 8:10)), bd4@x == 6
          )
assert.EQ.mat(dM4, m4)
assert.EQ.mat(M4^M4, m4^m4)
assert.EQ.mat(mN, matrix(NA, 3,4))
assert.EQ.mat(bdiag(diag(4)), diag(4))
sL <- Matrix(, 3,4, sparse=TRUE)# -> "lgC"
trS <- Matrix(tr, sparse=TRUE)# failed in 0.9975-11
stopifnotValid(d4, "diagonalMatrix");   stopifnotValid(z4,  "diagonalMatrix")
stopifnotValid(tr, "triangularMatrix"); stopifnotValid(trS, "triangularMatrix")
stopifnot(all(is.na(sL@x)), ## not yet:  all(is.na(sL)),
          !any(sL, na.rm=TRUE), all(!sL, na.rm=TRUE),
          validObject(Matrix(c(NA,0), 4, 3, byrow = TRUE)),
          validObject(Matrix(c(NA,0), 4, 4)))
stopifnotValid(Matrix(c(NA,0,0,0), 4, 4), "sparseMatrix")
I <- i1 <- I1 <- Diagonal(1)
## TODO? stopifnot(identical(I, Matrix(1, sparse=TRUE))) # doDiag=TRUE default
I1[1,1] <- i1[1, ] <- I [ ,1] <- NA
stopifnot(identical3(I,i1,I1))
image(d4) # gave infinite recursion

## Steve Walker, Mar 12, 2014:
n <- 7
(M <- triu(Matrix(seq_len(n^2), n, sparse=TRUE)))
im <- image(M) # was not an n-by-n image plot
stopifnot(n == diff(sort(im$y.limits)))
## ylimits were too small (by 1 on each side)

## Image of *empty*  sparseMatrix (has failed in spite of claim in r3282 | 2018-08-20)
(Z <- Matrix(0, 11, 22)) # empty "dgCMatrix"
Z0 <- Z; Z0[1,1] <- 1; Z0@x <- 0 ## is also all 0, but not "empty"
stopifnot(all(Z == Z0))
image(Z)
## gave Error in seq.default(zrng[1], zrng[2], length.out = cuts + 2) : 'from' must be a finite number
image(Z, useAbs=FALSE) # gave *different* Error in seq.int... : 'length.out' must be ...
image(Z0,useAbs=FALSE) #  (ditto)
image(Z0) # had worked previously already
showProc.time()

assertError( Matrix(factor(letters)) )
n.lsec <- length(.leap.seconds)# 27 (2017-07)
mlp <- matrix(.leap.seconds)## 27 x 1 numeric matrix
Mlp <- Matrix(.leap.seconds)
stopifnot(identical(dim(Mlp), c(n.lsec, 1L)))
assert.EQ.mat(Mlp, mlp)
lt.leap.seconds <- .LS <- as.POSIXlt(.leap.seconds, tz = "GMT") # GMT => sparse HMS
.LS <- unclass(.LS); .LS <- .LS[setdiff(names(.LS), "zone")]
                     # "zone" is character (not there for GMT/UTC in R <= 4.2.x)
(matLS <- data.matrix(data.frame(.LS)))
stopifnot(inherits(MLp <- as(matLS, "Matrix"), "sparseMatrix"),
          is.EQ.mat(MLp, matLS))
printSpMatrix(MLp, col.names = TRUE) # nice sparse  dgC* w/ col.names
(mLp <- matrix(lt.leap.seconds))## prints fine as 27 x 1 matrix of dates (internally is list+dim)
##
E <- rep(c(TRUE,NA,TRUE), length.out=8)
F <- new("nsparseVector", length = 8L, i = c(2L, 5L, 8L))
e <- as(E, "sparseVector"); f <- as(F,"lsparseVector")
Fv <- as.vector(F, "any") # failed in Matrix <= 1.2.0, and base::as.vector(.) failed too:
stopifnot(E | as.vector(F), identical(E | F, F | E),
	  all(e | f), all(E | F), # <- failed  Ops.spv.spv
	  identical(Fv, base::as.vector(F)),
	  is.logical(Fv), which(Fv) == c(2,5,8))
F[-8:-1] # was partly "illegal" (length = 0 is ok; nnz '3' is not)
F.ineg <- lapply(1:8, function(k) F[-8:-k])
F.pos <-  lapply(1:8, function(k) F[seq_len(k-1)])
F.vec <-  lapply(1:8, function(k) Fv[seq_len(k-1)])
str(F.vec, vec.len=8)
(nT <- vapply(F.vec, sum, 1L)) # == 0 0  1 1 1  2 2 2
str(whichT <- lapply(F.vec, which))
i.lengths <- function(L) vapply(L, function(.) length(.@i), 1L)
stopifnot(exprs = {
    identical(lengths(F.vec), 0:7)
    identical(lengths(F.ineg), 0:7)
    identical(lengths(F.pos),  0:7)
    identical(i.lengths(F.pos),  nT)
    identical(i.lengths(F.ineg), nT) # failed before 2018-03-19
    identical(lapply(F.pos, slot, "i"), whichT)
    identical(lapply(F.ineg, slot, "i"), whichT) # failed before 2018-03-19
})
## Here, sparseVector '[' is really wrong:
SV <- new("nsparseVector", length = 30L,
          i = c(1L, 8L, 9L, 12L, 13L, 18L, 21L, 22L))
NI <- -c(1:5, 7:10, 12:16, 18:27, 29,30)
sv <- SV[1:14]; ni <- -(1:14)[-c(6,11)] # smaller example
selectMethod("[", c("nsparseVector","index","missing","missing"))# the culprit
if(FALSE)
  trace("[", browser, signature=c("nsparseVector","index","missing","missing"))
stopifnot(
    SV[NI] == as.vector(SV)[NI] ## badly failed before 2018-03-19
    ,
    sv[ni] == as.vector(sv)[ni] ## ditto
)
if(FALSE)
    untrace("[", signature=c("nsparseVector","index","missing","missing"))

dT <- new("dgTMatrix",
	  i = c(1:2,1:2), j=rep(1:2, each=2), Dim = c(4L, 4L), x = c(1, 1, NA, 2))
dt <- new("dtTMatrix", i = 0:3, j = 0:3, Dim = c(4L, 4L), x = c(1,0,0,0),
	  uplo = "U", diag = "N")
c1 <- as(dT, "CsparseMatrix")
c2 <- as(dt, "CsparseMatrix")
isValid(lc <- c1 > c2,"lgCMatrix")
isValid(lt <- dT > dt,"lgCMatrix")
stopifnot(identical(lc,lt))

## Versions of Diagonal()
dD <- Diagonal(x = 5:1)
sD <- .symDiagonal(5, 5:1)
tD <- .trDiagonal (5, 5:1) # x=<integer> had failed for both
stopifnot(all(sD == dD)
        , all(tD == dD)
        , identical(sD, .symDiagonal(x = 5:1))
        , identical(tD, .sparseDiagonal(x=5:1))
        , identical(tD, .trDiagonal (x = 5:1)))# 'n' now has default

M <- Diagonal(4); M[1,2] <- 2 ; M
cM <- crossprod(M) # >> as_cholmod_l_triplet(): could not reallocate for internal diagU2N()
stopifnot(identical(cM, tcrossprod(t(M))))

S.na <- spMatrix(3, 4, c(1,2,3), c(2,3,3), c(NA,1,0))
(S.na <- S.na - 2 * S.na)
(L <- S.na != 0)

(Ln0 <- S.na != rep(0, prod(dim(L))))
.Lm0 <- S.na != Matrix(0, 3, 4)
stopifnot(Q.eq(L, Ln0), identical(Ln0, .Lm0));  rm(Ln0, .Lm0)
showProc.time()


### Unit-diagonal and unitriangular  {methods need diagU2N() or similar}
I <- Diagonal(3)
(T <- as(I,"TsparseMatrix")) # unitriangular
(C <- as(I,"CsparseMatrix")) #   (ditto)
lT <- as(T,"lMatrix")
lC <- as(C,"lMatrix")
stopifnot(
    identical((n0 <- I != 0), Diagonal(3, TRUE)), I@diag == "U",
    identical(n0, I & TRUE), identical(n0, I | FALSE),
    identical(n0, TRUE & I), identical(n0, FALSE | I),
    all(n0 == !(I == 0)), all(I == n0), identical(n0 == I, I == n0)
    ,
    identical4(lT, as(Diagonal(3, x=TRUE),"TsparseMatrix"), T & TRUE, TRUE & T),
    identical4(lC, as(Diagonal(3, x=TRUE),"CsparseMatrix"), C & TRUE, TRUE & C),
    identical3(lT, T | FALSE, FALSE | T),
    identical3(lC, C | FALSE, FALSE | C),
    TRUE)
I[,1] <- NA; I[2,2] <- NA ; I[3,] <- NaN
stopifnotValid(I, "sparseMatrix")
I # gave error in printSpMatrix() - because of R bug in format.info()

L <- spMatrix(9, 30, i = rep(1:9, 3), 1:27, (1:27) %% 4 != 1)
M <- drop0(crossprod(L))
diag(M) <- diag(M) + 5 # to make it pos.def.
M. <- M[1:12,1:12] # small ex
N3 <- as(Matrix(upper.tri(diag(3))), "nMatrix")
stopifnotValid(bdN <- bdiag(N3, N3),"nsparseMatrix")

stopifnot(identical(L, L == TRUE), ## used to give infinite recursion
          all(drop0((0 - L) != 0) == drop0(L)))
L[sample(length(L), 10)] <- NA
ll <- as(L,"logical")
stopifnot(all.equal(mean(L,  na.rm=TRUE),
		    mean(ll, na.rm=TRUE), tolerance=1e-14),
	  all.equal(mean(L,  na.rm=TRUE, trim=1/4),# <- with a warning
		    mean(ll, na.rm=TRUE, trim=1/4), tolerance=1e-14))


## Examples where  is.na(.) was wrong:
validObject(sc <- new("dsCMatrix", i=as.integer(c(0,0:1,1:2,0:1,3)), Dim=c(4L,4L),
                      p = c(0L,1L,3L,5L,8L), x = c(0,NA,NA,0:1,0,NA,1)))
validObject(gc <- as(sc, "generalMatrix"))

stopifnot(isSymmetric(M), isSymmetric(M.),
	  is(bdiag(M., M.),"symmetricMatrix"),
	  is(bdN, "triangularMatrix"),
          all(sc == gc | (is.na(sc) & is.na(gc))),
	  all.equal(N3,N3),
	  identical(all.equal(N3, t(N3)), all.equal.raw(1:6, 2:7)), # == "6 element mismatches"
	  all((bdN != t(bdN)) == (bdN + t(bdN))), # <nsparse> != <nsparse>  failed to work...
	  !any((0+bdN) > bdN), # <dsparse> o <nsparse>
	  !any(bdN != (0+bdN)), # <nsparse> o <dsparse>
	  length(grep("Length", all.equal(M., (vM <- as.vector(M.))))) > 0,
	  identical(M., (M2 <- Matrix(vM, 12,12))),
	  all.equal(M., M2, tolerance =0)
	  )

Filter(function(.) inherits(get(.), "symmetricMatrix"), ls())
## [1] "M"   "M."  "M2"  "cM"  "d4T" "d4d" "o4"  "sD"  "sc"
cc <- kronecker(cM, Diagonal(x = c(10,1)))
dimnames(cc) <- list(NULL, cn <- letters[1:ncol(cc)])
stopifnotValid(cc, "dsCMatrix")
(tt <- as(cc, "TsparseMatrix")) # shows *symmetric* dimnames
stopifnot(identical3(  cc @Dimnames,   tt @Dimnames, list(NULL, cn)),
          ## t() does not reverse 'Dimnames' slot for symmetricMatrix
	  identical3(t(cc)@Dimnames, t(tt)@Dimnames, list(NULL, cn)),
	  identical3(dimnames(cc), dimnames(tt), list(cn, cn))) # now symmetric!

stopifnot(identical3(dimnames(cc),
                     dimnames(as(cc, "generalMatrix")), ## should fixup dimnames to *symmetric*
                     dimnames(as(tt, "generalMatrix"))))
## -->  .Call(Csparse_symmetric_to_general, from)

mat <- as(cc, "matrix") ## --> should fixup dimnames to *symmetric*
mat # should print *symmetric* dimnames
stopifnot(identical3(dimnames(cc), dimnames(mat), dimnames(as(tt, "matrix"))))

selectMethod(coerce, c("dsCMatrix", "denseMatrix"))
dmat <- as(cc, "denseMatrix") ## --> gave Error (!!) in Matrix 1.1-5
stopifnot(identical3(tt@Dimnames, dmat@Dimnames, list(NULL, cn)))
dmat # should print *symmetric* dimnames (not modifying dmat as it did intermittently)
stopifnot(identical(dmat@Dimnames, list(NULL, cn)))
ttdm <- as(tt, "denseMatrix")
stopifnot(identical(dmat, ttdm),
          identical(dimnames(cc), dimnames(dmat)),
          ## coercing back should give original :
	  identical(cc,            as(dmat,  "sparseMatrix")),
	  identical(asUniqueT(tt), as(ttdm, "TsparseMatrix")))

## MM: now *if* cc is "truly symmetric", these dimnames should be, too:
d5 <- cn[1:5]; dnm5 <- list(d5,d5)
stopifnot(identical(dimnames( cc  [1:5, 1:5]), dnm5),
          identical(dimnames(t(cc)[1:5, 1:5]), dnm5))


## large sparse ones: these now directly "go sparse":
str(m0 <- Matrix(0,     nrow=100, ncol = 1000))
str(l0 <- Matrix(FALSE, nrow=100, ncol = 200))
stopifnot(all(!l0),
          identical(FALSE, any(l0)))

if(!interactive()) warnings()
## really large {length(<dense equivalent>) is beyond R's limits}:
op <- options(warn = 2) # warnings (e.g. integer overflow!) become errors:
n <- 50000L
stopifnot(n^2 > .Machine$integer.max)
## had integer overflow in index constructions:
x <- 1:n
D <- Diagonal(n, x=x[n:1])
summary(D)# special method
summary(D != 0)
stopifnot(identical(x*D, (Dx <- D*x)),
	  identical(D != 0, as(D, "lMatrix")),
	  identical(Dx, local({d <- D; d@x <- d@x * x; d})))

Lrg <- new("dgTMatrix", Dim = c(n,n))
l0 <- as(as(Lrg, "lMatrix"), "CsparseMatrix") # lgC
d0 <- as(l0, "dMatrix")
showProc.time()

if(FALSE) { #_____________________ FIXME: Should use cholmod_l_*() everywhere (?)____
## problem in  Csparse_to_dense :
dl0 <- as(l0, "denseMatrix")
dd0 <- as(d0, "denseMatrix")
## currently, both give --- Error in asMethod(object) :
##   Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 105
##--> And there it is 'Int_max' ==> ../src/CHOLMOD/Include/cholmod_internal.h
## defines 'Int_max'  and does that depending of "mode", and
## MM thinks we should use the "DLONG" mode now -- for 64-bit integers!
##  ==> Then Int_max := SuiteSparse_long_max := LONG_MAX
## (the latter from ../src/SuiteSparse_config/SuiteSparse_config.h )
## ==> use cholmod_l_<foo> instead of cholmod_<foo> in *many places*
##
## check they are ok
stopifnot(identical(dim(dl0), c(n,n)), identical(dim(dd0), c(n,n)),
          !any(dl0), all(dd0 == 0))
rm(dl0, dd0)# too large to keep in memory and pass to checkMatrix()
}

diag(Lrg[2:9,1:8]) <- 1:8
## ==:  Lrg[2:9,1:8] <- `diag<-`(Lrg[2:9,1:8], 1:8)

options(warn = 1) # Matrix 1.5-0 (Sep 2022): we now warn when lots of memory is used
(memGB <- Sys.memGB("MemFree")) # from test-tools-1.R, only works with /proc/*
system.time( # ~10 sec.                            __vv__
    e1 <- if(doExtras && is.finite(memGB) && memGB > 30) { # need around 18 GB
              try(Lrg == Lrg)
              ## had Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 105
## (error message almost ok)
    }) # now works, taking 42.7 sec on ada-20 w/ 504 GB;
if(is(e1, "Matrix")) object.size(e1) # 10000001176 bytes
system.time( # ~10 sec.                            __vv__
    e2 <- if(doExtras && is.finite(memGB) && memGB > 30) { # need around 18 GB
              try(!Lrg) # now *works* on 64-bit machines with enough RAM
              ## and immediately errors if LONG_VECTORs are not available
          }) ## when it works, see
## Warning in .sparse2dense(x) : sparse->dense coercion: allocating vector of size 9.3 GiB
##  user  system elapsed
## 6.812  10.744  17.612

str(e2) # error, NULL or "worked" (=> 50000 x 50000 lgeMatrix)
ina <- is.na(Lrg)# "all FALSE"
stopifnot(if(inherits(e1, "try-error")) grepl("too large", e1)
          else is.null(e1) ||
               (is(e1, "denseMatrix") && is(e1, "lMatrix")),
	  if(inherits(e2, "try-error")) grep("too large", e2) == 1
	  else is.null(e2) || length(e2@x) == n^2,
          !any(ina))# <- gave warning previously
stopifnot(suppressWarnings(any(Lrg)))# (double -> logical  warning)
rm(e1, e2)# too large...
## Matrix bug #6610, did segfault
system.time({ # ... sec                           __vv__
 ## FIXME: reproducible example (not using 'MatrixModels' which triggers) "Out of memory" etc
})
showProc.time()

RNGversion("3.6.0")# future proof
if(doExtras && is.finite(memGB) && memGB > 49) withAutoprint({
    cat("computing SM .. \n")
    showSys.time(m <- matrix(0, 3e6, 1024))
    ##  user  system elapsed
    ## 2.475  10.688  13.196  (faster in past ??)
    set.seed(1); inot0 <- unique(sort(c(1, length(m), sample(length(m), 20))))
    ai0 <- arrayInd(inot0, .dim=dim(m), useNames=FALSE)
    showSys.time(m[inot0] <- 1:22)
    ##  user  system elapsed
    ## 5.931  11.184  17.162
    showSys.time(SM  <- as(m, "sparseMatrix")) # ~ 8 sec
    ## gave 'Error in asMethod(object) : negative length vectors are not allowed'
    ## now works - via C  matrix_to_Csparse()
    showSys.time(n0.m <- c(m) != 0) # logical (full, base R) matrix, 12 GB
    ##   user  system elapsed
    ## 14.901  10.789  25.776
    try( ## _FIXME_ in R: Error ... long vectors not supported yet
        in0.m <- which(n0.m)
    )
    ## DONE: now very fast! [previously did coerce the whole matrix to dense first !]
    subS <- SM[inot0]
    selectMethod("[", c("dgCMatrix","numeric","missing","missing"))# -> .M.vectorSub(x,i)
    ## Directly via arrayInd(), is *FAST*:
    subSij <- SM[ai0]
    stopifnot(subS == 1:22, identical(subS, subSij))
    cat(" [Ok]\n")
    rm(m)
    str(SM)
    ## checking SM:
    TM <- as(SM, "TsparseMatrix")
    stopifnot(as.matrix(summary(TM)) == cbind(ai0, 1:22))
    ## cleanup:
    rm(SM, TM)
})
showProc.time()

## Constructing *packed* dense symmetric (dsp*) | triangular (dtp*) Matrices:
if(doExtras && is.finite(memGB) && memGB > 35) withAutoprint({
    m <- as.integer(2^16) ## = 65536
    showSys.time(x <- rep(as.numeric(1:100), length.out=m*(m+1)/2))
    ##  user  system elapsed
    ## 6.028   8.964  15.074
    gc()
    object.size(x) # 17'180'131'368 bytes: ~ 17 GB
    mat <- new("dspMatrix", x = x, Dim = c(m, m)) # failed with
    ## long vectors not supported yet: ../../src/include/Rinlinedfuns.h:...
    validObject(mat)
    mat <- new("dtpMatrix", x = x, Dim = c(m, m)) # failed .......
    validObject(mat)
    ## cleanup
    rm(mat)
})
showProc.time()

options(warn = 2)# warnings => errors
## with dimnames:
v <- c(a=1, b=2:3)
m <- as.matrix(v)
(M <- as(v, "denseMatrix") )
stopifnot(identical(dimnames(m), list(c("a", "b1", "b2"), NULL)),
	  inherits(M, "dgeMatrix"),
	  identical(dimnames(M), dimnames(m)))

## dimnames(.) of symmpart() / skewpart() :
ns <- c("symmpart", "skewpart", "forceSymmetric")
symFUNs <- setNames(lapply(ns, get), ns); rm(ns)
chkSS <- function(m) {
  r <- lapply(symFUNs, function(fn) fn(m))
  m0 <- as(m, "matrix")
  r0 <- lapply(symFUNs, function(fn) fn(m0))
  stopifnotValid(fS  <- r [["forceSymmetric"]], "symmetricMatrix")
  stopifnotValid(fS0 <- r0[["forceSymmetric"]], "symmetricMatrix")
  dnms <- dimnames(m)
  d.sy <- dimnames(r[["symmpart"]])
  id <- if(is.null(dnms[[2]]) && !is.null(dnms[[1]])) 1 else 2
  stopifnot(identical(d.sy, dnms[c(id,id)]),
	    identical(d.sy, dimnames(r [["skewpart"]])),
            if(identical(dnms, list(NULL, NULL)))
                is.null(dimnames(r0[["skewpart"]]))
	    else identical(d.sy, dimnames(r0[["skewpart"]])),
	    all(m  == with(r,  symmpart + skewpart)),
	    all(m0 == with(r0, symmpart + skewpart)),
	    identical(dS <- dimnames(fS), dimnames(fS0)),
	    identical(dS[1], dS[2]))
}
cat(sprintf("chkSS() {valid %s} for a list of matrices:\n",
            paste(paste0(names(symFUNs), "()"), collapse=", ")))
L <- list(M1 = Matrix(1:4, 2,2),
          M2 = Matrix(c(0, rep(1:0, 3),0:1), 3,3))
L$M3 <- pack(as(forceSymmetric(L$M2), "denseMatrix"))
stopifnotValid(L$M3, "dspMatrix")
for(m in L) {
    cat("\n---\nm:\n"); show(m)
    chkSS(m)
    dn <- list(row = paste0("r", 1:nrow(m)), col = paste0("var.", 1:ncol(m)))
    dimnames(m) <- dn		; chkSS(m)
    colnames(m) <- NULL		; chkSS(m)
    dimnames(m) <- unname(dn)	; chkSS(m)
}

m. <- matrix(c(0, 0, 2:0), 3, 5)
dimnames(m.) <- list(LETTERS[1:3], letters[1:5])
(m0 <- m <- Matrix(m.))
m@Dimnames[[2]] <- m@Dimnames[[1]]
## not valid anymore:
(val <- validObject(m, test=TRUE)); stopifnot(is.character(val))
dm <- as(m0, "denseMatrix"); rm(m)
stopifnot(all.equal(rcond(dm), rcond(m.), tolerance = 1e-14),
	  ##^^^^^^^ dm and m. are both dense, interestingly small differences
	  ## show in at least one case of optimized BLAS
	  all.equal(rcond(dm), 0.4899474520656),
          ## show(<dgRMatrix>) had revealed a bug in C:
          identical(capture.output(show(as(m0, "RsparseMatrix")))[-(1:2)],
                    gsub("0", ".",  capture.output(show(m.))[-1])))
m.1 <- m.; dimnames(m.1) <- list(row=NULL, col=NULL)
M.1 <- Matrix(m.1, sparse=TRUE)
show(M.1)# had bug in .formatSparseSimple()
showProc.time()

###--  Sparse Triangular :

g5 <- new("dgCMatrix", Dim = c(5L, 5L),
          x = c(10, 1, 3, 10, 1, 10, 1, 10, 10),
          i = c(0L,2L,4L, 1L, 3L,2L,4L, 3L, 4L),
          p = c(0L, 3L, 5L, 7:9))
t5 <- as(g5, "triangularMatrix") # fine
stopifnot(class(t5) == "dtCMatrix",
          identical(t5, tril(g5)))
## This is really a regression test for 'methods::selectMethod()'
## Maybe move to R once 'Matrix' is recommended
sm <- selectMethod(coerce, c("dgCMatrix", "triangularMatrix"), verbose=TRUE)
stopifnot(identical(sm(g5), t5))
dimnames(t5) <- list(row=paste0("r",1:5), col=paste0("C.",1:5))
s5 <- symmpart(t5) # gave an error
showProc.time()

(t1 <- new("dtTMatrix", x= c(3,7), i= 0:1, j=3:2,
           Dim= as.integer(c(4,4))))
## Diagonal  o  Sparse
I4 <- Diagonal(4)
D4 <- Diagonal(4, x=1:4)
validObject(t1)
validObject(t2  <-   t1  + I4)
validObject(tt2 <- t(t1) + I4)
validObject(t1c <- as(t1, "CsparseMatrix"))
validObject(t2c <- as(t2, "CsparseMatrix"))
stopifnotValid(2 * I4, "diagonalMatrix")
stopifnotValid(D4 * 3, "diagonalMatrix")
stopifnotValid(I4 / 5, "diagonalMatrix")
stopifnotValid(D4 / 2, "diagonalMatrix")
stopifnotValid(t1c + I4,"triangularMatrix")
stopifnotValid(t2c + I4,"triangularMatrix")
stopifnot(identical(t1, t(t(t1))),
          identical(t1c, t(t(t1c))),
          c(class(t2), class(t1c), class(t2c), class(tt2)) == "dtCMatrix",
          identical(t(tt2), t2))
assert.EQ.mat(t1, as(t1c, "matrix"))
D4. <- D4 * (A4 <- Matrix(1:4, 4,4))
D4p <- A4 + D4
Lg1 <- D4 > 0 & D4 > 1
nLg <- !Lg1
nnLg <- !nLg
D4m <- D4 * 4:1
assert.EQ.mat(D4., diag(x= (1:4)^2))
assert.EQ.mat(D4p, diag(x= (1:4)) + (1:4))
assert.EQ.mat(D4m, diag(x=c(4,6,6,4)))
assert.EQ.mat(Lg1, diag(x= c(FALSE, rep(TRUE,3))))
stopifnot(is(Lg1, "diagonalMatrix"), is(D4m, "diagonalMatrix"),
	  is(D4., "diagonalMatrix"),
          is(nLg, "generalMatrix"), is(nnLg, "generalMatrix"),
          identical3(Lg1,
                     Matrix(nnLg, forceCheck = TRUE),
                     as(nnLg, "diagonalMatrix")),
          all(Lg1 != (!Lg1)))

## tri[lu](<diagonal>)
td3 <- triu(diag(3)); stopifnot(is(td3, "triangularMatrix"), td3@uplo == "U")
Ld3 <- tril(diag(3)); stopifnot(is(Ld3, "triangularMatrix"), Ld3@uplo == "L")
## the latter did fail till 2014-12-20
D3 <- Diagonal(3)
stopifnot(identical3(D3, tril(D3), triu(D3)))
## methods were missing
showProc.time()

## as(<diag>, <anything>) :
str(cls <- names(getClass("Matrix")@subclasses))# all Matrix classes
## Matrix >= 1.4.2 : (basically) only coerce to virtual classes:
table(isVirt <- vapply(cls, isVirtualClass, NA))
for(cl in cls[isVirt])
    if(canCoerce(I4, cl)) {
	cat(cl,":")
	M  <- as(I4, cl)
	M. <- as(D4, cl)
        stopifnot(diag(4) == as(M,"matrix"),
                  if(is(cl,"dMatrix")) diag(x=1:4) == as(M.,"matrix") else TRUE)
	cat(" [Ok]\n")
    }
s4 <- as(D4,"sparseMatrix")
v <- c(11,2,2,12); s4[2:3,2:3] <- v; validObject(s4)
s4. <- D4; s4.[2:3,2:3] <- v; validObject(s4.)
stopifnot(all(s4 == s4.))
## now assign symmetrically to symmetricMatrix
s4 <- as(as(D4,"sparseMatrix"),"symmetricMatrix")
s4[2:3,2:3] <- v
validObject(s4)
stopifnot(is(s4,"symmetricMatrix"))
assert.EQ.mat(s4, as(s4.,"matrix"),tol=0)

## lower-triangular unit-diagonal
L <- new("dtCMatrix", i = 1L, p = c(0:1, 1L), Dim = c(2L, 2L),
         x = 0.5, uplo = "L", diag = "U")
stopifnot(range(L) == 0:1, all.equal(mean(L), 5/8))

## from  0-diagonal to unit-diagonal triangular {low-level step}:
tu <- t1 ; tu@diag <- "U"
tu
validObject(cu <- as(tu, "CsparseMatrix")) # still unitriangular
validObject(cnu <- diagU2N(cu))# <- testing diagU2N
validObject(tt <- as(cu, "TsparseMatrix"))
tu. <- tt ## validObject(tu. <- as(cu, "dtTMatrix"))
stopifnot(exprs = { ## NOT: identical(tu, tu.), # since T* is not unique!
    identical(cu, as(tt, "CsparseMatrix"))
    length(cnu@i) == length(cu@i) + nrow(cu)
    identical(cu, diagN2U(cnu)) # <- testing diagN2U
    all(cu >= 0, na.rm = TRUE)
    all(cu >= 0)
    any(cu >= 7)
})
validObject(tcu <- t(cu))
validObject(ttu <- t(tu))
validObject(ltu <- as(ttu, "lMatrix"))
validObject(ldtu <- as(ltu, "denseMatrix"))
validObject(Cltu <- as(ltu, "CsparseMatrix"))
stopifnot(identical(asCsp(ttu > 0), asCsp(ltu)),
          all(ltu == as(ttu > 0,"denseMatrix")))
ltu - (ttu > 0) # failed
assert.EQ.mat(cu,  as(tu,"matrix"), tol=0)
assert.EQ.mat(cnu, as(tu,"matrix"), tol=0)
showProc.time()

C <- suppressWarnings(Matrix(c(0,1,0,0), 5,5)) + Diagonal(5)
(tU <- diagN2U(tril(C))) # dtC Unitriangular
ntU <- as(tU, "nMatrix")
nT <- as(ntU, "TsparseMatrix")
R <- as(tU, "RsparseMatrix")
Tt <- diagU2N(R) # used to accidentally drop the diag.
stopifnot(R@x == c(1,1,1), diag(Tt) == 1)

lcu <- new("ltCMatrix", Dim = c(4L, 4L), i = c(0:1, 0L), p = c(0L, 0:3),
           x = c(TRUE, FALSE, FALSE), uplo = "U", diag = "U")
(lTu <- as(lcu,"TsparseMatrix"))# prints wrongly (in Matrix 0.999375-31)
stopifnot(identical3(rowSums(lcu), rowSums(lTu), rowSums(drop0(lcu))))
(ncu <- as(lcu, "nMatrix"))# -- gives the "pattern" of lcu, i.e. FALSE are *there*
ncn <- diagU2N(ncu)
(cncn <- crossprod(ncn))# works -> "nsCMatrix"
stopifnot(identical(ncu, as(lcu,"nsparseMatrix")),
	  identical(rowSums(ncu), c(3:1, 1L)),
	  Q.eq(ncn, ncu),
	  Q.eq(crossprod(drop0(lcu)), crossprod(lcu)),# crossprod works -> "dsCMatrix"
	  identical(crossprod(ncu), cncn),
	  Q.eq(cncn, t(ncu) %&% ncu)) #used to seg.fault

U <- new("dtCMatrix", Dim = c(6L, 6L),
	 i = c(0:1, 0L, 2:3, 1L, 4L),
	 p = c(0L,0L,0L, 2:3, 5L, 7L),
	 x = rep.int(-0.5, 7), diag = "U")
validObject(U)
U. <- solve(iU <- solve(U))#-> gave segmentation fault
stopifnot(validObject(U), ## had a case where solve(U) modified U !
	  validObject(iU),
	  validObject(U.),
	  ## no rounding error, since have iU@x * 8 is integer :
	  identical(U, diagN2U(drop0(U.))))
showProc.time()

## <sparse> o <numeric> (of length > 1):
stopifnotValid(tm <- tu * 1:8, "sparseMatrix")
stopifnot(identical4(tm, cu * 1:8, 1:8 * cu, 1:8 * tu))

cu[1,2] <- tu[1,2] <- NA
mu <- as(tu,"matrix")
stopifnotValid(cu, "CsparseMatrix"); stopifnotValid(cu, "triangularMatrix")
stopifnotValid(tu, "TsparseMatrix"); stopifnotValid(tu, "triangularMatrix")
stopifnot(identical(cu * 1:8, tu * 1:8), # but are no longer triangular
          identical(cu > .1, as(tu > .1, "CsparseMatrix")),
          all(cu >= 0, na.rm=TRUE), !all(cu >= 1), is.na(all(tu >= 0)),
          ## Csparse_drop: preserves triangularity incl diag="U"
          ## ^^^^^^^^^^^^ now using more general R_sparse_drop0
          identical(cu, .Call(Matrix:::R_sparse_drop0, cu, 0.))
          )
assert.EQ.mat(cu * 1:8, mu * 1:8)

ina <- is.na(as(cu,"matrix"))
## These 3 were each different (2008-03) !!
stopifnot(all(ina == is.na(cu)),
	  all(ina == is.na(as(cu,"generalMatrix"))),
	  all(ina == as(is.na(as(cu,"matrix")),"nMatrix")))


set.seed(7)
xx <- rpois(10, 50)
Samp <- function(n,size) sample(n, size, replace=TRUE)
Tn <- sparseMatrix(i=Samp(8, 50), j=Samp(9,50), x=xx, repr = "T")
Tn
stopifnot(xx == Tn@x,
	  max(xx) < max(Tn), 0 == min(Tn),
	  (sT <- sum(Tn)) == sum(colSums(Tn)), sT == sum(Tn@x),
	  range(Tn) == range(as(Tn, "CsparseMatrix")))

## tu. is diag "U", but tu2 not:
tu2 <- as(as(tu., "generalMatrix"), "triangularMatrix")
assert.EQ.mat(cu, mu, tol=0)
stopifnot(identical3(cu[cu > 1],  tu [tu > 1], mu [mu > 1]),
          identical3(cu <= 1, tu <= 1, as(mu <= 1, "lMatrix")),# all lgeMatrix
	  identical3(cu[cu <= 1], tu[tu <= 1], mu[mu <= 1]),
	  identical3(cu , triu(cu ), t(t(cu))),
	  identical3(tu , triu(tu ), t(t(tu))),
	  identical3(tu., triu(tu.), t(t(tu.))),
	  identical(tu2, triu(tu2)),
	  identical(tcu , tril(tcu)),
	  identical(ttu , tril(ttu)),
	  identical(t(tu), tril(t(tu)))
          )
assert.EQ.mat(triu(cu), as(triu(as(cu, "matrix")), "matrix"))
for(k in -1:1)
    assert.EQ.mat(tril(cu,k), as(tril(as(cu, "matrix"), k), "matrix"))

(dtr <- Matrix(local({m <- diag(2); m[1,2] <- 3;m})))
identical(dtr, triu(dtr))
assert.EQ.mat(tril(dtr), diag(2))


(t4 <- new("dgTMatrix", i = 3:0, j = 0:3, x = rep(1,4), Dim = as.integer(c(4,4))))
c4 <- as(t4, "CsparseMatrix")
## the same but "dsT" (symmetric)
suppressWarnings(M <- Matrix(c(0, rep(c(0,0:1),4)), 4,4))# warning:.. length [13] is not ..multiple
tt <- as(M, "TsparseMatrix")
stopifnot(all.equal(triu(t4) + tril(t4), c4),
          all.equal(triu(tt) + tril(tt), c4))
showProc.time()


###-- Numeric Dense: Crossprod & Solve

set.seed(123)
mm. <- mm <- Matrix(rnorm(500 * 150), ncol = 150)
stopifnot(validObject(mm))
xpx <- crossprod(mm)
stopifnot(identical(mm, mm.))# once upon a time, mm was altered by crossprod()
stopifnotValid(xpx, "dpoMatrix")
str(mm) # 'dge*"
str(xpx)# 'dpo*"
xpy <- crossprod(mm, rnorm(500))
res <- solve(xpx, xpy)
str(xpx)# now with Cholesky factor
stopifnot(validObject(xpx),
          validObject(xpy),
          validObject(res))
stopifnot(all.equal(xpx %*% res, xpy, tolerance = 1e-12))
lp <- xpx >= 1
slp <- as(lp, "sparseMatrix")

ltlp  <-  lp[ lower.tri(lp) ]
sltlp <- slp[ lower.tri(slp) ]
dim(ij <- which(lower.tri(lp), arr.ind = TRUE))
ss <- slp[ij] # now fast (!)
stopifnot(identical4(lp[ij], ltlp, sltlp, as(lp, "matrix")[ij]),
          identical(ss, sltlp),
          isValid(lp, "lsyMatrix"), lp@uplo == "U")
showProc.time()

###-- more solve() methods  {was ./solve.R }

## first for "dgeMatrix" and all kinds of RHS :
(m6 <- 1 + as(diag(0:5), "generalMatrix"))
rcond(m6)
I6 <- as(diag(6), "generalMatrix")
stopifnot(all.equal(I6, m6 %*% solve(m6)),
          all.equal(I6, solve(m6) %*% m6) )

(i6 <- solve(m6, Matrix(1:6)))
stopifnot(identical(i6, as(cbind(c(-4, rep(1,5))), "generalMatrix")),
          identical(drop(i6), solve(m6, 1:6)),
          identical(i6, solve(m6, matrix(1:6))),
          identical(i6, solve(m6, matrix(c(1,2,3,4,5,6))))
          )

## solve(<sparse>)
(m <- t1+ t(t1) + Diagonal(4))
i.m <- solve(as.mat(m))
I1 <- m %*% i.m
o4 <- diag(I1)
im <- solve(m)# is now sparse {not yet}
(I2 <- m %*% im)
(ms <- as(m, "symmetricMatrix"))
## solve(<sparse>, <sparse>):
s.mm <-  solve(m,m)
s.mms <- solve(m, ms)
## these now work "fully-sparse"
s.ms2 <- solve(ms, ms)
s.msm <- solve(ms, m)
I4c <- as(Matrix(diag(4),sparse=TRUE), "generalMatrix")
stopifnot(isValid(im, "Matrix"), isValid(I2, "Matrix"), class(I4c) == "dgCMatrix",
          all.equal(I1, as(I2,"denseMatrix"), tolerance = 1e-14),
          all.equal(diag(4), as.mat(I2), tolerance = 1e-12),
          all.equal(s.mm,  I2, tolerance = 1e-14),
          all.equal(s.mms, I2, tolerance = 1e-14),
          all.equal(s.ms2, s.msm, tolerance = 4e-15),
          all.equal(s.ms2, I4c  , tolerance = 4e-15),
          abs(o4 - 1) < 1e-14)

image(T125 <- kronecker(kronecker(t5,t5),t5),
      main = paste("T125:",class(T125)))
dim(T3k <- kronecker(t5,kronecker(T125, t5)))
system.time(IT3 <- solve(T3k))# incredibly fast
I. <- drop0(zapsmall(IT3 %*% T3k))
I.. <- diagN2U(I.)
I <- Diagonal(5^5)
stopifnotValid(IT3, "dtCMatrix")
stopifnot(## something like the equivalent of  all(I. == Diagonal(3125)) :
          identical(as(I., "diagonalMatrix"), I),
          identical(as(I..,"diagonalMatrix"), I)
          )
showProc.time()

## printSpMatrix() ;  "suppressing (columns | rows) .." {and do it correctly!}
IT3
op0 <- options(width = 70, max.print = 1000)
T125[-(1:50),] ## suppression ... is it correctly done?

## Still buggy -- FIXME: see ../TODO --- even if we'd require max.print >= 5 or so
for(mm in 1:21) {
    options(max.print=mm)
    cat("----------\n\nmax.print=",mm,":\n", sep="")
    cat("\n>> U:   ") ; show(U)
    cat("\n>> slp: ") ; show(slp)
}

options(op0)# revert to max.print = 1000

###-- row- and column operations  {was ./rowcolOps.R }

set.seed(321)
(m1 <- round(Matrix(rnorm(25), 5), 2))
m1k <- Matrix(round(rnorm(1000), 2), 50, 20)
m.m <- as(m1k, "matrix")
stopifnot(all.equal(colMeans(m1k), colMeans(m.m)),
          all.equal(colSums (m1k), colSums (m.m)),
          all.equal(rowMeans(m1k), rowMeans(m.m)),
          all.equal(rowSums (m1k), rowSums (m.m)),
          all.equal(colMeans(m1k, na.rm=TRUE), colMeans(m.m, na.rm=TRUE)),
          all.equal(colSums (m1k, na.rm=TRUE), colSums (m.m, na.rm=TRUE)),
          all.equal(rowMeans(m1k, na.rm=TRUE), rowMeans(m.m, na.rm=TRUE)),
          all.equal(rowSums (m1k, na.rm=TRUE), rowSums (m.m, na.rm=TRUE)) )

###-- kronecker for nonsparse uses Matrix(.):
stopifnotValid(kr <- kronecker(m1, m6), "Matrix")
assert.EQ.mat(kr,
              kronecker(as(m1, "matrix"),
                        as(m6, "matrix")), tol = 0)

## sparse:
(kt1 <- kronecker(t1, tu))
kt2 <- kronecker(t1c, cu)
stopifnot(identical(as(kt1, "CsparseMatrix"), kt2))
ktf <- kronecker(as(t1, "matrix"), as(tu, "matrix"))
if(FALSE) # FIXME? our kronecker treats "0 * NA" as "0" for structural-0
assert.EQ.mat(kt2, ktf, tol= 0)
(cs1 <- colSums(kt1))
NA.or.True <- function(x) is.na(x) | x
eq <- (cs1 == colSums(as(kt1, "matrix")))
stopifnot(NA.or.True(eq), identical(is.na(eq), is.na(cs1)))
nt1 <- as(kt1, "nMatrix") # no NA's anymore
(ng1 <- as(as(nt1, "generalMatrix"),"CsparseMatrix")) # ngC
dg1 <- as(ng1, "dMatrix")# dgC
lt1 <- kt1 > 5
nt1 <- as(lt1, "nMatrix")
(colSums(nt1, sparseResult = TRUE))
(colSums(kt1, sparseResult = TRUE)) # dsparse, with NA
(colSums(lt1, sparseResult = TRUE)) # isparse, with NA
(colSums(lt1, sparseResult = TRUE, na.rm = TRUE))
(colSums(nt1, sparseResult = TRUE)) # isparse, no NA
## check correct sparseness of both:
for(M in list(kt1, nt1, ng1, dg1, lt1, nt1)) {
    m <- as(M, "matrix")
    for(na.rm in c(FALSE,TRUE)) {
	cs  <- colSums(M, na.rm = na.rm)
	cs. <- colSums(M, na.rm = na.rm, sparseResult = TRUE)
	rs  <- rowSums(M, na.rm = na.rm)
	rs. <- rowSums(M, na.rm = na.rm, sparseResult = TRUE)
	stopifnotValid(cs., "sparseVector")
	stopifnotValid(rs., "sparseVector")
	stopifnot(identical(cs, as(cs., "vector")),
		  identical(rs, as(rs., "vector")),
		  {eq <- cs == colSums(m, na.rm = na.rm) ; ineq <- is.na(eq)
		   all(ineq | eq) && identical(ineq, is.na(cs)) },
		  {eq <- rs == rowSums(m, na.rm = na.rm) ; ineq <- is.na(eq)
		   all(ineq | eq) && identical(ineq, is.na(rs)) } )
    }
}

i1 <- cs. == 1
cs2 <- cs.
cs2[i1] <- 0 # failed in *-31 !!
## now *index* with a NA-sparseVector :
i2 <- i1 ; i2[3] <- NA ; li2 <- as.logical(i2)
cs3 <- cs. ;	       cs3 [i2] <- 0
v3 <- as(cs.,"vector"); v3[li2] <- 0
cs4 <- cs.	     ; cs4[li2] <- 0
stopifnot(length(i1@x) == 2, identical(li2, as(i2,"vector")),
	  identical(cs3, cs4),
	  cs3 == v3, all(as(v3, "sparseVector") == cs3)
	  ## indexing simple "numeric" with sparseVector:
	  ## see 'R_FIXME' in ../R/sparseVector.R
	  ## , identical(v3[i2], v3[li2])
	  ## TODO:
	  ## sub-assigning into simple "numeric" with sparseVector index:
	  )
showProc.time()


M <- Matrix(c(2:0,1),2); M. <- as(M, "sparseMatrix")
(N <- as(crossprod(kronecker(diag(2), M)) > 0,
         "nMatrix"))
(L. <- as(N,"lMatrix"))
stopifnot(identical(N, as(L.,"nMatrix")),
	  identical(kronecker(	  c(1,0), M),
		    kronecker(cbind(1:0), M)))
assert.EQ.mat(kronecker(M,	         c(1,0,0)),
	      kronecker(as(M, "matrix"), c(1,0,0)))

## coercion from "dpo" or "dsy"
xx <- as(xpx, "dsyMatrix")
stopifnot(isSymmetric(xxS  <- as(xx,  "sparseMatrix")),
          isSymmetric(xpxS <- as(xpx, "sparseMatrix")))

tm <- matrix(0, 8,8)
tm[cbind(c(1,1,2,7,8),
         c(3,6,4,8,8))] <- c(2,-30,15,20,80)
(tM <- Matrix(tm))                ## dtC
(mM <- Matrix(m <- (tm + t(tm)))) ## dsC
mT <- as(mM, "TsparseMatrix")
gC <- as(as(mT, "generalMatrix"), "CsparseMatrix")
lT <- as(Matrix(TRUE, 2, 2), "TsparseMatrix")
## Check that mT, lT, and gC print properly :
pr.mT <- capture.output(mT)
pr.lT <- capture.output(lT)[-(1:2)]
nn <- unlist(strsplit(gsub(" +\\.", "", sub("^....", "", pr.mT[-(1:2)])), " "))
stopifnot(as.numeric(nn[nn != ""]) == m[m != 0],
	  identical(1:2, grep("|", pr.lT, fixed=TRUE)),
	  identical(pr.lT, capture.output(as(lT, "nMatrix"))[-(1:2)]),
          capture.output(gC)[-1] == pr.mT[-1])
assert.EQ.mat(tM, tm, tol=0)
assert.EQ.mat(gC, m,  tol=0)
assert.EQ.mat(mT, m,  tol=0)
stopifnotValid(mM, "dsCMatrix")
stopifnotValid(tM, "dtCMatrix")
stopifnot(identical(mT, as(mM, "TsparseMatrix"))
	  , identical(gC, as(mM, "generalMatrix"))
	  ## coercions	general <-> symmetric
	  , identical(as(mM. <- as(mM, "generalMatrix"), "symmetricMatrix"), mM)
	  , identical(as(as(mM., "TsparseMatrix"), "symmetricMatrix"), mT)
	  , identical(as(as(tM, "generalMatrix"),"triangularMatrix"), tM)
          , identical(tM + Diagonal(8), tMD <- Diagonal(8) + tM)
	  )
stopifnotValid(tMD, "dtCMatrix")
eM <- eigen(mM) # works thanks to base::as.matrix hack in ../R/zzz.R
stopifnot(all.equal(eM$values,
                { v <- c(162.462112512353, 30.0665927567458)
                  c(v, 15, 0, 0, 160-v[1], -15, -v[2])}, tolerance=1e-14))

##--- symmetric -> pos.def. needs valid test:
m5 <- Matrix(diag(5) - 1)
assertError(as(m5, "dpoMatrix"))# not pos.definite!
pm5 <- as(m5, "packedMatrix") # packed
assertError(as(pm5, "dppMatrix"))# not pos.definite!
sm <- as(Matrix(diag(5) + 1), "packedMatrix")
pm <- as(sm,"dpoMatrix")## gave infinite recursion (for a day or so)
pp <- as(pm,"dppMatrix")

x <- round(100 * crossprod(Matrix(runif(25),5)))
D <- Diagonal(5, round(1000*runif(5)))
px <- pack(x)
stopifnot(is(x, "dsyMatrix"), is(px, "dspMatrix"), is(D, "ddiMatrix"))

class(x+D)#--> now "dsyMatrix"
stopifnot(is(x+D, "symmetricMatrix"),
	  is(D+px, "dspMatrix"),
	  identical(x+D, D+x), identical(px+D, D+px), identical(pack(x-D), px-D))

tx <- tril(x)
ptx <- pack(tx)
stopifnot(is(tx, "dtrMatrix"), is(ptx, "dtpMatrix"),
          is(t(tx), "dtrMatrix"), is(t(ptx), "dtpMatrix"),
          is(D + tx, "dtrMatrix"), is(tx + D, "dtrMatrix"),
          is(ptx + D, "dtpMatrix"), is(D + ptx, "dtpMatrix"))
showProc.time()


###-- dense nonzero pattern:
class(m <- Matrix(TRUE,2,2)) # lsy
isValid(n <- as(m, "nMatrix"), "nsyMatrix")

## 1)
as(n,"CsparseMatrix") # used to give CHOLMOD error: invalid xtype...
ls2 <- as(m, "CsparseMatrix") # works fine
## and really  'm' and 'n' are interally slot identical (!!!)

as(n,"sparseMatrix")
as(m,"sparseMatrix")

### -- now when starting with nsparse :
nT <- new("ngTMatrix",
          i = as.integer(c(0, 1, 0)),
          j = as.integer(c(0, 0, 1)), Dim = as.integer(c(2,2)))
(nC <- as(nT, "CsparseMatrix"))
str(nC)# of course, no 'x' slot

tt <- as(nT,"denseMatrix") # nge (was lge "wrongly")
stopifnot(is(tt,"ngeMatrix"),
	  identical(as(tt, "lMatrix"),
		    as(as(nT, "lMatrix"), "denseMatrix")))
tt
as(nC,"denseMatrix")
showProc.time()


###-- sparse nonzero pattern : ----------

(nkt <- as(as(as(kt1, "CsparseMatrix"), "generalMatrix"), "nMatrix"))# ok
dkt <- as(nkt, "denseMatrix")
(clt <- crossprod(nkt))
stopifnotValid(nkt, "ngCMatrix")
stopifnotValid(clt, "nsCMatrix")
suppressWarnings(crossprod(clt)) ## warning "crossprod() of symmetric ..."

## a Csparse with *repeated* entry is not valid!
assertError(new("ngCMatrix", p = c(0L,2L), i = c(0L,0L), Dim = 2:1))


### "d" <-> "l"  for (symmetric) sparse : --------------------------------------
data(KNex, package = "Matrix")
mm <- KNex$mm
xpx <- crossprod(mm)
## extract nonzero pattern
nxpx <- as(xpx, "nMatrix")
show(nxpx) ## now ok, since subsetting works
r <- nxpx[1:2,]
lmm <- as(mm, "lMatrix")
nmm <- as(lmm, "nMatrix")
xlx <- crossprod(lmm)
x.x <- crossprod(nmm)
## now A = lxpx and B = xlx should be close, but not quite the same
## since <x,y> = 0 is well possible when x!=0 and y!=0 .
## However,  A[i,j] != 0 ==> B[i,j] != 0:
A <- as(as(nxpx, "lMatrix"), "TsparseMatrix")
B <- as(as(xlx,  "lMatrix"), "TsparseMatrix")
ij <- function(a) a@i + ncol(a) * a@j
stopifnot(all(ij(A) %in% ij(B)))

l3 <- upper.tri(matrix(,3,3))
validObject(c3 <- as(l3, "CsparseMatrix"))
stopifnotValid(c3, "lMatrix")# lgC
(M <- Matrix(l3))
stopifnotValid(M, "ltCMatrix")
stopifnotValid(M2 <- M %x% M, "triangularMatrix") # is "dtT" (why not "dtC" ?)
stopifnot(dim(M2) == c(9,9), identical(M2, kronecker(M,M)))
M3 <- M %x% M2 #ok
(cM3 <- colSums(M3, sparseResult=TRUE))
identical(as.vector(cM3),
          as(rev(rowSums(M3, sparseResult=TRUE)), "vector"))
M. <- M2 %x% M # gave infinite recursion

## diagonal, sparse & interactions
stopifnotValid(as(Diagonal(3), "TsparseMatrix"), "TsparseMatrix")
stopifnotValid(X <- Diagonal(7) + 1.5 * tM[1:7,1:7], "sparseMatrix")
stopifnotValid(X, "triangularMatrix")
stopifnotValid(XX <- X - chol(crossprod(X)), "triangularMatrix")
X
XX
XX <- as(drop0(XX), "symmetricMatrix")
stopifnot(identical(XX, Matrix(0, nrow(X), ncol(X), doDiag=FALSE)))

M <- Matrix(m., sparse = FALSE)
(sM <- Matrix(m.))
class(dlM <- M >= 1)
stopifnot(identical(dlM, !(M < 1)))
stopifnotValid(sM, "sparseMatrix")
stopifnotValid(dlM, "denseMatrix")
(lM  <- as(dlM, "sparseMatrix"))
lM2 <- as(dlM, "CsparseMatrix") #-> now ok
lM0 <- .dense2sparse(dlM, "C")
stopifnot(identical3(lM, lM2, lM0))

selectMethod("coerce",	c("lgeMatrix", "CsparseMatrix"),
	     useInherited = c(from = TRUE, to = FALSE))

ms0 <- Matrix(c(0,1,1,0), 2,2)
ms <- as(ms0, "TsparseMatrix")
cs <- as(ms, "CsparseMatrix")
ll <- as(ms, "lMatrix")
lt <- as(ll, "generalMatrix")
nn <- as(cs, "nsparseMatrix")
l2 <- as(cs, "lsparseMatrix")
nt <- triu(nn)
n3 <- as(nt, "lsparseMatrix")
da <- nt + t(nt)
dm <- nt * t(nt) + da
##
mnt <- as(nt, "matrix")
m <- rbind(v <- 2:3)
assert.EQ.mat(nt %*% v, mnt %*% v)
assert.EQ.mat(v %*% nt, v %*% mnt)
assert.EQ.mat( crossprod(nt, v),  crossprod(mnt,v))
assert.EQ.mat( crossprod(v, nt),  crossprod(v,mnt))
assert.EQ.mat(tcrossprod(v, nt), tcrossprod(v,mnt))
assert.EQ.mat(tcrossprod(nt, m), tcrossprod(mnt, m))
##
stopifnotValid(ms, "dsTMatrix")
stopifnot(as(ms0,"matrix") == as(ll, "matrix"), # coercing num |-> log
	  as(lt, "matrix") == as(ll, "matrix"),
	  identical(ms, as(ll, "dMatrix")),
	  identical4(as(ll, "CsparseMatrix"), as(cs, "lMatrix"),# lsC*
		     as(nn, "lsparseMatrix"), l2),
	  identical3(da, dm, as(cs, "generalMatrix")),		# dgC*
	  identical(as(da, "lMatrix"), as(lt, "CsparseMatrix")) # lgC*
	  )
## Dense *packed* ones:
s4 <- as(D4, "symmetricMatrix")
sp <- as(s4, "packedMatrix")
tp <- triu(sp)
tpL <- tril(sp)
(spL <- t(sp))
stopifnot(sp @uplo=="U", tp @uplo=="U",
	  spL@uplo=="L", tpL@uplo=="L")

## band():
n <- 4 ; m <- 6
r1 <- Matrix(1:24, n,m)
validObject(M1 <- band(r1, 0,0))
(M1 <- as(M1, "sparseMatrix"))
r2 <- Matrix(1:18, 3, 6)
stopifnot(identical(M1, bandSparse(n,m, k=0, diagonals = list(diag(r1)))),
	  identical(band(r2, 0,4),
		    band(r2, 0,3) + band(r2, 4,4)))
s1 <- as(r1, "sparseMatrix") # such that band(s1) is sparse, too
for(k1 in (-n):m)
    for(k2 in k1:m) {
	stopifnotValid(br1 <- band(r1, k1,k2), "ddenseMatrix")
	stopifnotValid(bs1 <- band(s1, k1,k2), "CsparseMatrix")
        stopifnot(all(r1 == s1))
    }
showProc.time()

## large dimensions -- gave integer overflow
## R-forge bug 6743 by Ariel Paulson
## https://r-forge.r-project.org/tracker/?func=detail&atid=294&aid=6743&group_id=61
n <- 47000
stopifnotValid(Mn <- sparseMatrix(i = rep(1:6, 2), dims = c(n,n),
                                  j = c(1L,4L, 6:8, 10:12, 16:19)), "nsparseMatrix")
stopifnotValid(M <- as(Mn, "dMatrix"), "dgCMatrix")
dim(M) # 47000 47000
i <- 46341
stopifnotValid(bM <- band(M, i,   i  ), "dtCMatrix")
## gave Error in  if (sqr && k1 * k2 >= 0) ....
##      In addition: Warning message:
##      In k1 * k2 : NAs produced by integer overflow
x <- 1:999
bM2 <- bandSparse(n, k=i+(0:3), diagonals = list(x,10*x,32*x,5*x), symmetric=TRUE)
stopifnotValid(bM2, "dsCMatrix")
stopifnotValid(bb2 <- band(bM2, k1=i-2, k2=i+5), "dtCMatrix")
stopifnotValid(b0  <- band(bM2, -1000, 1000),    "dsCMatrix")
stopifnotValid(b0a <- band(bM2, -1000, 1001),    "dgCMatrix")
(id <- nrow(M)-i)# 659
colN <- colSums(bM2 != 0)
stopifnot(exprs = {
    identical(bb2, triu(bM2))
    identical(b0 @x, numeric(0))
    identical(b0a@x, numeric(0))
    identical(bM2, band(bM2, -(i+3), i+3))
    assert.EQ(as(bM2, "generalMatrix"),
              band(bM2, -(i+3), i+11), showOnly = TRUE)
    colN == { cN <- c(1:3, rep(4L, id-3)); c(rev(cN), rep(0L, i-id), cN)}
})
showProc.time()
## some of these failed before Matrix 1.4.0 (Oct.7, 2021)


D. <- Diagonal(x= c(-2,3:4)); D.[lower.tri(D.)] <- 1:3 ; D.
D0 <- Diagonal(x= 0:3);       D0[upper.tri(D0)] <- 1:6 ; D0
stopifnot(all.equal(list(modulus = structure(24, logarithm = FALSE), sign = -1L),
                    unclass(determinant(D.,FALSE)), tolerance=1e-15),
	  det(Matrix(0,1)) == 0,
          all.equal(list(modulus = structure(0, logarithm = FALSE), sign = 1L),
                    unclass(determinant(D0,FALSE)), tolerance=0)
          )

### More sparseVector checks: -------------------------------
validObject(new("isparseVector"))
R <- sv <- as(D4, "sparseVector")
## dim(<sparseVector>) <- (n1,n2)  --> sparse Matrix :
dim(R) <- dim(D4)
stopifnotValid(sv,"sparseVector")
stopifnotValid(R, "sparseMatrix")
stopifnot(identical(D4, as(R, "diagonalMatrix")))
iv <- c(rep(0, 5), 3, 0,0,7,0,0,0)
sv <- as(iv, "sparseVector")
sv. <- as(as.integer(iv), "sparseVector")
## Note: Method with signature "numeric#sparseVector" chosen ...
(sv2 <- as(sv, "isparseVector")) ## gave error
as(sv, "zsparseVector")
stopifnot(identical(sv., sv2),
	  identical(  Matrix(sv, 3,4, byrow=TRUE),
		    t(Matrix(sv, 4,3))))
options(warn = 0)# no longer error


## "Large" sparse:
n <- 100000
m <-  50000 ; nnz <- 47
M <- spMatrix(n, m,
              i = sample(n, nnz, replace = TRUE),
              j = sample(m, nnz, replace = TRUE),
              x = round(rnorm(nnz),1))
validObject(Mv <- as(M, "sparseVector"))
validObject(Dv <- as(Diagonal(60000), "sparseVector"))
validObject(LD <- Diagonal(60000, TRUE))
validObject(Lv <- as(LD, "sparseVector"))
Dm <- Dv; dim(Dm) <- c(180000L, 20000L)
stopifnot(!doExtras || isValid(Md <- M * rowSums(M, sparseResult=TRUE), "sparseMatrix"),
	  LD@diag == "U",
          isValid(Dm, "sparseMatrix"),
	  identical(Dv, as(Dm, "sparseVector")))

p. <- new("dtCMatrix", i = c(2:3, 2L), p = c(0L, 2:3, 3L, 3L),
          Dim = c(4L, 4L), x = rep(-0.5, 3), uplo = "L", diag = "U")
assert.EQ.mat(solve(solve(p.)), as(p., "matrix"))
dimnames(p.)[[1]] <- paste(1:4)
ii <- is.na(p.)
stopifnot(all(!ii), !any(as(ii, "denseMatrix")))# used to fail

lst <- ls()
table(istri <- sapply(lst, function(.) is(get(.),"triangularMatrix")))
table(triC <- sapply(lst[istri], function(.) class(get(.))))
table(uniC <- sapply(lst[istri], function(.) get(.)@diag == "U"))
lsUtr <- lst[istri][uniC]
(di <- sapply(lsUtr, function(.) dim(get(.))))
## TODO: use %*%, crossprod(), .. on all those  4 x 4 -- and check "triangular rules"

assertError(new("ltrMatrix", Dim = c(2L,2L), x=TRUE))# gave "illegal" object w/o error
assertError(new("ntrMatrix", Dim = c(2L,2L)))#  dito


showProc.time()# == "stats"

cat("doExtras:",doExtras,"\n")
if(doExtras) {
    cat("checkMatrix() of all: \n---------\n")
    Sys.setlocale("LC_COLLATE", "C")    # to keep ls() reproducible
    for(nm in setdiff(ls(), "d4da")) { # FIXME: checkMatrix(d4da)
        if(is(.m <- get(nm), "Matrix")) {
            cat("\n", rep("-",nchar(nm)),"\n",nm, ":\n", sep='')
            checkMatrix(.m)
        }
    }
    showProc.time()
}
## in any case, test
d4d.2 <- .dense2sparse(!!d4da, "C") ## <<- did wrongly make dimnames symmetric
l4da <- as(d4da, "lMatrix")
assert.EQ.Mat(l4da, as(l4da,"CsparseMatrix"))

dtr <- tr4 <- triu(Matrix(1:16, 4,4))
dtr@x[Matrix:::indTri(4, upper=FALSE, diag=FALSE)] <- 100*(-3:2)
stopifnot(all.equal(dtr, tr4), # because are same *as* simple matrices
	  dtr@x[1:4] == c(1, -(3:1)*100),
	  range(tr4) == c(0,16),
	  range(dtr) == c(0,16)) # <- failed

##  new("nsyMatrix") + new("lgeMatrix") # failed
cln <- sort(outer(c("l","n"), paste0(c("ge","sy"), "Matrix"), paste0))
dim(c.c <- as.matrix(expand.grid(cln, cln, KEEP.OUT.ATTRS=FALSE))) # 16 x 2
## clTry <- function(expr) class(tryCatch(expr, error=identity))[[1]]
## '+' [Arith] failed -- now fixed
cbind(c.c, Res = apply(c.c, 1, function(x) class(new(x[1]) + new(x[2]))))
## '<' [Compare] works fine
cbind(c.c, Res = apply(c.c, 1, function(x) class(new(x[1]) < new(x[2]))))

if(!interactive()) warnings()

## R-forge matrix-Bugs [#6708] (2021-02-25, by David Cortes):
sVec <- sparseVector(c(1,exp(1),pi), c(1,3,7), length=9)
 vec <- c(1, 0, exp(1), 0, 0, 0, pi, 0, 0)
stopifnot(identical(as.matrix(sVec), as.matrix(vec)),
          identical(as.array (sVec), as.array (vec)))

## R-forge matrix-Bugs [#6656] (2020-02-05, by Chun Fung (Jackson) Kwok (kcf.jackson)
## (*is* a bug, but not in kronecker etc, but rather in Arith / Ops)
dC <- sparseMatrix(i=1:4, j=1:4, x=5:2, triangular = TRUE)
(dT <- as(dC, "TsparseMatrix"))
stopifnot(identical(--dC, dC),
          identical(--dT, dT)
          )
## both  - <sparse-triang.> gave : Error .... 'factors' is not a slot in class "dtTMatrix"

## R PR#18250 - by Mikael Jagan
nm2 <- c("a","b")
x <- new("dspMatrix", x = c(3,2,1), Dim = c(2L,2L), Dimnames = list(nm2, NULL))
dn <- list(nm2,nm2)
stopifnotValid(x. <- unpack(x), "dsyMatrix")
validObject( y  <- as(x , "generalMatrix") )
validObject( y. <- as(x., "generalMatrix") )
stopifnotValid( l  <- x > 0,     "lspMatrix")
stopifnotValid( l. <- unpack(l), "lsyMatrix")
stopifnotValid( lg <- as(l, "generalMatrix"), "lgeMatrix")
stopifnotValid( lg2<- as(l.,"generalMatrix"), "lgeMatrix")
stopifnot(exprs = {
    identical(dimnames(x ), dn)
    identical(dimnames(x.), dn)
    identical(dimnames(y ), dn) # was wrong
    identical(dimnames(y.), dn) # was wrong
    identical(dimnames(l ), dn)
    identical(dimnames(l.), dn)
    identical(dimnames(lg), dn) # was wrong
    identical(lg, lg2)
    ## even more cases (?)
})
showProc.time()

dn4 <- list(letters[1:4], LETTERS[1:4])
(D4n <- `dimnames<-`(D4, dn4))
m4 <- as(D4n, "matrix")
stopifnot(identical(dimnames(m4), dn4), Q.eq(D4n, m4, superclasses=NULL))
## as(<ddi>, "matrix")  had lost dimnames before

s24 <- new("dgCMatrix", Dim = c(2L, 4L), p = integer(5L))
triu(s24, k = 4L) # was an error
tril(s24, k = 4L) # was an error

## band(<sparseMatrix>, -k, k) used isSymmetric(tol > 0) to test
## for symmetry of the result, and forcing symmetry lost information
s44 <- new("dgCMatrix", Dim = c(4L, 4L), p = c(0L, 0L, 1L, 1L, 1L),
           i = 0L, x = .Machine$double.xmin)
(bs44 <- band(s44, -1L, 1L))
stopifnot(identical(s44, bs44))

l0.u <- new("ltrMatrix", diag = "U")
!l0.u # was an error

n11 <- new("ngeMatrix", Dim = c(1L, 1L), x = NA)
nn11 <- !n11 # did not respect NA<=>TRUE
stopifnot(is(nn11, "ngeMatrix"), identical(nn11@x, FALSE))

## coercions preserving mathematical equality ought to preserve 'factors' slot;
## though currently only for sparse->sparse and dense->dense
mC <- as(as(Diagonal(x = rlnorm(5)), "CsparseMatrix"), "symmetricMatrix")
stopifnot(identical(mC@factors, list()))
chol(mC)
mR <- as(mC, "RsparseMatrix")
mT <- as(mR, "TsparseMatrix")
stopifnot(!identical(mTf <- mT@factors, list()),
          identical(mTf, mR@factors),
          identical(mTf, mC@factors))

## overallocated l.T should follow usual logic ... NA || TRUE -> TRUE, etc.
lT. <- lT0 <- lT1 <-
    new("lgTMatrix", Dim = c(1L, 1L), i = c(0L, 0L), j = c(0L, 0L),
        x = c(NA, NA))
lT0@x[1L] <- FALSE
lT1@x[1L] <- TRUE
stopifnot(identical(as.vector(lT.), NA),
          identical(as.vector(lT0), NA),
          identical(as.vector(lT1), TRUE),
          identical(as(lT1, "CsparseMatrix")@x, TRUE),
          identical(as(lT1, "dMatrix")@x, 1))

## various is.na(), anyNA(), which() bugs in Matrix <= 1.4-1
.nge <- new("ngeMatrix", Dim = c(2L, 2L), x = rep.int(NA, 4L))
.dtr <- new("dtrMatrix", Dim = c(2L, 2L), x = c(NA, 1, 2, Inf), diag = "U")
.dsy <- new("dsyMatrix", Dim = c(2L, 2L), x = c(1, NA, 2, 3))
stopifnot(!any(is.na(.nge)),
          !any(is.na(.dtr)),
          !any(is.infinite(.dtr)),
          !any(is.na(.dsy)),
          !anyNA(.nge),
          !anyNA(.dtr),
          !anyNA(.dsy),
          identical(which(.nge), 1:4))

## various `dim<-`() bugs in Matrix <= 1.4-1
.dgR <- new("dgRMatrix", Dim = c(2L, 2L), p = integer(3L))
assertError(`dim<-`(.dgR, -x@Dim)) # had yielded an invalid object
stopifnot(is(`dim<-`(.dgR, c(4L, 1L)), "RsparseMatrix")) # was TsparseMatrix

## symmpart(<ldiMatrix>) was not a dMatrix or a symmetricMatrix;
## symmpart(<diagonalMatrix>) did not get symmetrized 'Dimnames'
.ldi.sp <- symmpart(new("ldiMatrix", Dim = c(1L, 1L), Dimnames = list("a", "b"),
                        x = TRUE))
stopifnot(is(.ldi.sp, "dMatrix"),
          is(.ldi.sp, "diagonalMatrix"),
          Matrix:::isSymmetricDN(.ldi.sp@Dimnames))

## as.vector(<ndenseMatrix>), etc. must do NA->TRUE
stopifnot(identical(as.vector(.nge), rep.int(TRUE, 4L)),
          identical(as.logical(.nge), rep.int(TRUE, 4L)),
          identical(as.double(.nge), rep.int(1, 4L)),
          identical(as(.nge, "vector"), rep.int(TRUE, 4L)),
          identical(as(.nge, "matrix"), array(TRUE, .nge@Dim)),
          identical(as(.nge, "dMatrix")@x, rep.int(1, 4L)),
          identical(nnzero(.nge), 4L))

## symmpart(<matrix>) and skewpart(<matrix>) have been documented
## as returning matrix, _not_ Matrix
z0 <- matrix((-2)^(0:3), 2L, 2L)
z1 <- symmpart(z0)
z2 <- skewpart(z0)
stopifnot(!isS4(z1), is.matrix(z1), !isS4(z2), is.matrix(z2))

## various Matrix() bugs in Matrix <= 1.4-1
.d0 <- new("ddiMatrix", Dim = c(1L, 1L), Dimnames = list("A", "B"), diag = "U")
.d1 <- Matrix(`dimnames<-`(diag(1), list("A", "B")), doDiag = TRUE)
.s0 <- new("dsyMatrix", Dim = c(1L, 1L), Dimnames = list("B", "B"), x = 1)
.s1 <- Matrix(.d0, doDiag = FALSE)
stopifnot(identical(.d1, .d0),
          identical(.s1, .s0),
          identical(Matrix(sparseVector(1, 1L, 3L)),
                    new("dgTMatrix", Dim = c(3L, 1L), i = 0L, j = 0L, x = 1)),
          identical(Matrix(new("dgeMatrix"), sparse= TRUE, forceCheck=FALSE),
                    new("dgCMatrix")),
          identical(Matrix(new("dgCMatrix"), sparse=FALSE, forceCheck=FALSE),
                    new("dgeMatrix")),
          identical(Matrix(table(1)),       Matrix(1)),
          identical(Matrix(table(1, 1, 1)), Matrix(1)),
          grepl("too long",
                vapply(alist(Matrix(0, 0.5,    ), Matrix(0,    , 0.5),
                             Matrix(0, 0.0,    ), Matrix(0,    , 0.0)),
                       function(e) conditionMessage(assertError(eval(e))[[1L]]),
                       "")))

## From: Mikael Jagan <jaganmn2@gmail.com> , 16 Aug 2022
validObject(x <- new("dtrMatrix", Dim = c(2L, 2L), diag = "U", x = double(4L)))
x == 0 -> L
x >  2 -> L. # ditto
## gave Error in validObject(.Object) :
##    invalid class "ltrMatrix" object: length of x slot != prod(Dim)
stopifnot(all(L == ((1:4) != 2)),
          all(L. == diag(2)))
## 'Same' with "Logic" instead of "Compare":
x & FALSE -> L.
## gave Error in validObject(.Object) : invalid class "ltrMatrix" object
stopifnot(all(L. == diag(2)))

validObject(y <- new("dgCMatrix", Dim = c(0L, 6L), p = integer(7L)))
y == c(0, 0) -> L2
y == rep(0,6)-> L6
## gave Error in validObject(*): invalid class "lgCMatrix" object: slot p ...
stopifnot(identical(L2, L6), is(L2, "lgCMatrix"), identical(dim(L2), c(0L, 6L)))
## .. "Logic" instead of "Compare":
y & c(FALSE, FALSE) -> L2
## gave Error in validObject(*): invalid class "lgCMatrix" object: slot p ...
stopifnot(identical(L2, L6), is(L2, "lgCMatrix"), identical(dim(L2), c(0L, 6L)))

## Briefly wrong between 1.4-1 and 1.5-0
x.inf <- new("dgCMatrix", Dim = c(2L, 3L),
             p = c(0:2, 2L), i = 0:1, x = c(-Inf, Inf))
stopifnot(identical(is.infinite(x.inf), as(abs(x.inf) == Inf, "nMatrix")))
showProc.time()

## C-level bugs in 1.5-0, detected by full CRAN (incl. ASAN) check
as(new("dgTMatrix"), "CsparseMatrix") # out-of-bounds access
as(seq_len(10000), "pMatrix") # segfault

## <pMatrix> %*% <pMatrix> gave b %*% a prior to Matrix 1.5-2
set.seed(163006)
for(i in 1:6) {
    x <- as(sample.int(10L), "pMatrix")
    y <- as(sample.int(10L), "pMatrix")
    stopifnot(as(x %*% y, "matrix") == as(x, "matrix") %*% as(y, "matrix"))
}

## <indMatrix> %*% <indMatrix> forgot to set 'Dim' briefly prior to 1.5-2
x <- new("indMatrix", Dim = c(5L, 3L),
         perm = sample.int(3L, size = 5L, replace = TRUE))
y <- new("indMatrix", Dim = c(3L, 9L),
         perm = sample.int(9L, size = 3L, replace = TRUE))
validObject(x %*% y)
stopifnot(all(tril(y, -1) == 0)) # was wrong in Matrix 1.5-x

## dimScale(x) (i.e., with 'd1' missing) did not work in 1.5-2;
## same with dimScale(<matrix with NULL dimnames>) ...
set.seed(3054)
V <- matrix(rlnorm(16L), 4L, 4L)
stopifnot(all.equal(as(dimScale(V), "matrix"), cov2cor(V)))

## Diagonal(n, x, names=TRUE) must recycle 'x' _and_ its names
p <- 6L
a0 <- c(a = 0)
stopifnot(identical(unname(nD <- Diagonal(n = p, x = a0, names = TRUE)),
                    Diagonal(n = p, x = a0, names = FALSE)),
          identical(nD, Diagonal(n = p, x = rep(a0, p), names = TRUE)))

## Diagonal(n, names=<character>) should also get 'Dimnames'
stopifnot(identical(Diagonal(1L, names = "a")@Dimnames, list("a", "a")))

## Diagonal(x=<named 0-length>, names = TRUE) should get list(NULL, NULL)
stopifnot(identical(Diagonal(x = a0[0L], names = TRUE)@Dimnames,
                    list(NULL, NULL)))

## names were forgotten prior to 1.5-3
d1 <- Diagonal(1L, names = "b")
stopifnot(identical(colSums(d1), c(b = 1)),
          identical(rowMeans(as(d1, "indMatrix")), c(b = 1)))

## na.rm was ignored prior to 1.5-3
d1 <- Diagonal(x = NaN)
stopifnot(identical(colSums(d1), NaN),
          identical(colSums(d1, na.rm = TRUE), 0),
          identical(rowMeans(d1), NaN),
          identical(rowMeans(d1, na.rm = TRUE), NaN))

## Matrix bug #6810
library(Matrix)
x <- as(matrix(c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE), 4L, 2L),
        "RsparseMatrix")
stopifnot(identical(which(x), seq_along(x)[as.vector(!is.na(x) & x)]),
          identical(which(x, arr.ind = TRUE, useNames = FALSE),
                    arrayInd(which(x), dim(x), dimnames(x), useNames = FALSE)))

## drop0(give.Csparse = FALSE)
R <- new("dtRMatrix", Dim = c(6L, 6L), p = 0:6, j = 0:5, x = 0+0:5)
T <- as(as(R, "TsparseMatrix"), "symmetricMatrix")
for(tol in 0+0:5) {
    i.tol <- R@x > tol
    stopifnot(exprs = {
        identical(drop0(R, tol = tol, give.Csparse = FALSE),
                  new("dtRMatrix", Dim = c(6L, 6L),
                      p = c(0L, cumsum(i.tol)), j = R@j[i.tol], x = R@x[i.tol]))
        identical(drop0(T, tol = tol, give.Csparse = FALSE),
                  new("dsTMatrix", Dim = c(6L, 6L),
                      i = R@j[i.tol], j = R@j[i.tol], x = R@x[i.tol]))
    })
}
## No-op for pattern-like sparse matrices:
i1 <- new("indMatrix", Dim = c(6L, 8L), margin = 1L,
          perm = sample.int(8L, size = 6L, replace = TRUE))
i1.s <- as(i1, "nsparseMatrix")
i1.d <- as(i1,  "ndenseMatrix")
stopifnot(exprs = {
    identical(drop0(i1  , give.Csparse = FALSE), i1)
    identical(drop0(i1.s, give.Csparse = FALSE), i1.s)
    identical(drop0(i1.d, give.Csparse = FALSE), as(i1.d, "CsparseMatrix"))
})

## Setting diagonal elements of non-square RsparseMatrix produced
## an invalid object
x <- new("dgRMatrix", Dim = c(1L, 2L), p = c(0L, 0L))
diag(x) <- 1
validObject(x)

## Subassigning double to logical briefly did not change class
x <- new("lgeMatrix", Dim = c(2L, 3L), x = logical(6L))
y <- new("dgeMatrix", Dim = c(2L, 3L), x = replace(double(6L), 1L, 1))
x[1L, 1L] <- 1
stopifnot(identical(x, y))

## as(<data.frame>, "Matrix") was briefly a error (invalid type "list")
stopifnot(identical(as(data.frame(a = 1:2, b = 3:4), "Matrix"),
                    new("dgeMatrix", x = as.double(1:4),
                        Dim = c(2L, 2L), Dimnames = list(NULL, c("a", "b")))))

## tri[ul](<.t[rp]Matrix>) was often wrong at least in 1.6-1
u <- new("dtrMatrix", Dim = c(8L, 8L), x = as.double(seq_len(64L)))
stopifnot(identical(triu(u, 1L), triu(as(u, "generalMatrix"), 1L)))

## more tril()/triu() woes {introduced after 2021; present till 1.6-4}
for(n in 0:7) {
    cat("n = ", n,"\n----\n")
    ##TODO: for(m in pmax(0, n-3):(n+3)) {
    for(m in pmax(0, n-3):(n+3)) {  #--  n x m  matrix
        cat(" m = ", m,": using k's in  ", (-n),":", m, "\n", sep="")
        symm <- (m == n)
        mn2 <- (mn <- m*n) %/% 2
        ma <- array(seq_len(mn) - mn2 - 1/2, dim = c(n, m))
        if(symm) ma <- crossprod(ma) # crossprod() to be symmetric
        dM <- as(as(ma,   "denseMatrix"), "generalMatrix")
        sM <- as(as(ma, "CsparseMatrix"), "generalMatrix")
        for(k in (-n):m) {
            trum <- triuChk(ma,k);  trlm <- trilChk(ma, k)
            trud <- triuChk(dM,k);  trld <- trilChk(dM, k)
            trus <- triuChk(sM,k);  trls <- trilChk(sM, k)
            if(symm) {
                assert.EQ(    trum,        t(tril(ma,-k)))
                assert.EQ.mat(trud, as.mat(t(tril(dM,-k))))
                assert.EQ.mat(trus, as.mat(t(tril(sM,-k))))
            }
            stopifnot(exprs = {
                ## matrix
                identical(trlm, band(ma, -n, k))
                identical(trum, band(ma,  k, m))
                inherits(trum, "denseMatrix") # "dge" / "dtr" ...
                ## denseMatrix
                identical(trld, band(dM, -n, k))
                identical(trud, band(dM,  k, m))
                inherits(trud, "denseMatrix")
                assert.EQ.Mat(trud, trum, giveRE=TRUE)
                ## sparseMatrix
                identical(trls, band(sM, -n, k))
                identical(trus, band(sM,  k, m))
                inherits(trus, "sparseMatrix")
                assert.EQ.Mat(trus, trum, giveRE=TRUE)
            })
        }
    }
}

## Platform - and other such info -- so we find it in old saved outputs
.libPaths()
SysI <- Sys.info()
structure(Sys.info()[c(4,5,1:3)], class="simple.list")
sessionInfo()
c(Matrix = packageDescription("Matrix")$Built)
if(SysI[["sysname"]] == "Linux" && requireNamespace("sfsmisc")) local({
    nn <- names(.Sc <- sfsmisc::Sys.cpuinfo())
    nn <- names(.Sc <- .Sc[!grepl("^flags", nn)])
    print(.Sc[ grep("\\.[0-9]+$", nn, invert=TRUE) ])
})
showProc.time()
