#### For both 'Extract' ("[") and 'Replace' ("[<-") Method testing
####    aka    subsetting     and  subassignment
####           ~~~~~~~~~~          ~~~~~~~~~~~~~

## for R_DEFAULT_PACKAGES=NULL :
library(stats)
library(utils)

if(interactive()) {
    options(error = recover, warn = 1)
} else if(FALSE) { ## MM / developer  testing *manually* :
    options(error = recover, Matrix.verbose = 2, warn = 1)
} else {
    options(                 Matrix.verbose = 2, warn = 1)
}
## Matrix.verbose = .. (*before* loading 'Matrix' pkg)
## ==> will also show method dispath ambiguity messages: getOption("ambiguousMethodSelection")

#### suppressPackageStartupMessages(...)  as we have an *.Rout.save to Rdiff against
suppressPackageStartupMessages(library(Matrix))

source(system.file("test-tools.R", package = "Matrix"), keep.source = FALSE)
##-> identical3() etc
cat("doExtras:",doExtras,"\n")
if(exists("Sys.setLanguage", mode="function"))
    Sys.setLanguage("en")
englishMsgs <- (lang <- Sys.getenv("LANGUAGE")) %in% c("en", "C")
cat(sprintf("LANGUAGE env.: '%s'; englishMsgs: %s\n",
            lang, englishMsgs))

### Dense Matrices

m <- Matrix(1:28 +0, nrow = 7)
validObject(m)
stopifnot(identical(m, m[]),
	  identical(m[2, 3],  16), # simple number
	  identical(m[2, 3:4], c(16,23)), # simple numeric of length 2
	  identical(m[NA,NA], as(Matrix(NA, 7,4), "dMatrix")))

m[2, 3:4, drop=FALSE] # sub matrix of class 'dgeMatrix'
m[-(4:7), 3:4]        # ditto; the upper right corner of 'm'

## rows or columns only:
m[1,]     # first row, as simple numeric vector
m[,2]     # 2nd column
m[,1:2]   # sub matrix of first two columns
m[-(1:6),, drop=FALSE] # not the first 6 rows, i.e. only the 7th
m[integer(0),] #-> 0 x 4 Matrix
m[2:4, numeric(0)] #-> 3 x 0 Matrix

## logical indexing
stopifnot(identical(m[2,3], m[(1:nrow(m)) == 2, (1:ncol(m)) == 3]),
          identical(m[2,], m[(1:nrow(m)) == 2, ]),
          identical(m[,3:4], m[, (1:4) >= 3]))

## dimnames indexing:
mn <- m
dimnames(mn) <- list(paste("r",letters[1:nrow(mn)],sep=""),
                     LETTERS[1:ncol(mn)])
checkMatrix(mn)
mn["rd", "D"]
msr <- ms <- as(mn,"sparseMatrix")
mnr <- mn
v <- rev(as(ms, "vector"))
mnr[] <- v
msr[] <- v # [<- "sparse" -- not very sensical; did fail w/o a message
z <- msr; z[] <- 0
zz <- as(array(0, dim(z)), "sparseMatrix")
a.m <- as(mnr,"matrix")
stopifnot(identical(mn["rc", "D"], mn[3,4]), mn[3,4] == 24,
	  identical(mn[, "A"], mn[,1]), mn[,1] == 1:7,
	  identical(mn[c("re", "rb"), "B"], mn[c(5,2), 2]),
	  identical(ms["rc", "D"], ms[3,4]), ms[3,4] == 24,
	  identical(ms[, "A"], ms[,1]), ms[,1] == 1:7,
	  identical(ms[ci <- c("re", "rb"), "B"], ms[c(5,2), 2]),
	  identical(rownames(mn[ci, ]), ci),
	  identical(rownames(ms[ci, ]), ci),
	  identical(colnames(mn[,cj <- c("B","D")]), cj),
	  identical(colnames(ms[,cj]), cj),
	  identical(a.m, as(msr,"matrix")),
	  identical(unname(z), zz),
	  identical(a.m, array(v, dim=dim(mn), dimnames=dimnames(mn)))
	  )
showProc.time()

## Bug found thanks to Timothy Mak, Feb 3, 2017:
## sparseMatrix logical indexing with (partial) NA:
a.m <- as(mn,"matrix")
assert.EQ(as(ms,"matrix"), a.m) # incl. dimnames
iN4 <- c(NA, TRUE, FALSE, TRUE)
assert.EQ(as(mn[,iN4],"matrix"), a.m[,iN4]) # (incl. dimnames)
##assert.EQ(as.matrix(ms[,iN4]), a.m[,iN4]) # ms[, <with_NA>]  fails still : _FIXME_
try(ms[,iN4])
try(ms[,iN4] <- 100) ## <- segfaulted in Matrix <= 1.2-8  (!)

## R-forge Matrix bug #2556: Subsetting a sparse matrix did remove  names(dimnames(.)) :
m44 <- matrix(1:16, 4, 4, dimnames=list(row=c('a','b','c','d'), col=c('x','y','z','w')))
## Dense matrix: ------------------------------------------
a <- Matrix(m44)
identical(
    dimnames(m44[,FALSE, drop=FALSE]),
    dimnames(  a[,FALSE, drop=FALSE]))
chk.ndn <- function(a, a0=m44)
    stopifnot(identical(names(dimnames(a)), names(dimnames(a0))))
i <- 1:2
chk.ndn(a[i,]); chk.ndn(a[i, i])
## Sparse matrix: -----------------------------------------
s <- as(a %% 3 == 1, "sparseMatrix")
ts <- as(s,"TsparseMatrix")
b <- sparseMatrix(i=1:3, j=rep(2,3), dims=c(4,4), dimnames=dimnames(s))
tb <- as(b,"TsparseMatrix")
stopifnot(identical5(
    dimnames(a), dimnames(s), dimnames(ts),
    dimnames(b), dimnames(tb)))

chk.ndn(b [i, i]); chk.ndn(b [i, ])
chk.ndn(s [i, i]); chk.ndn(s [i, ])
chk.ndn(tb[i, i]); chk.ndn(tb[i, ])
chk.ndn(ts[i, i]); chk.ndn(ts[i, ])
chk.ndn( b[ , 1, drop=FALSE]); chk.ndn( s[i, 2, drop=FALSE])
chk.ndn(tb[ , 1, drop=FALSE]); chk.ndn(ts[i, 2, drop=FALSE])

L0 <- logical(0)
stopifnot(exprs = {
    identical(dim(b[,L0]), c(4L, 0L))
    identical(dim(b[L0,]), c(0L, 4L)) # failed till 2019-09-x
})

## Printing sparse colnames:
ms[sample(28, 20)] <- 0
ms <- t(rbind2(ms, 3*ms))
cnam1 <- capture.output(show(ms))[2] ; op <- options("sparse.colnames" = "abb3")
cnam2 <- capture.output(show(ms))[2] ; options(op) # revert
stopifnot(## sparse printing
	  grep("^ +$", cnam1) == 1, # cnam1 is empty
	  identical(cnam2,
		    paste(" ", paste(rep(rownames(mn), 2), collapse=" "))))

mo <- m
m[2,3] <- 100
m[1:2, 4] <- 200
m[, 1] <- -1
m[1:3,]

m. <- as(m, "matrix")

## m[ cbind(i,j) ] indexing:
iN <- ij <- cbind(1:6, 2:3)
iN[2:3,] <- iN[5,2] <- NA
stopifnot(identical(m[ij], m.[ij]),
	  identical(m[iN], m.[iN]))

## testing operations on logical Matrices rather more than indexing:
g10 <- m [ m > 10 ]
stopifnot(18 == length(g10))
stopifnot(10 == length(m[ m <= 10 ]))
sel <- (20 <  m) & (m <  150)
sel.<- (20 <  m.)& (m.<  150)
nsel <-(20 >= m) | (m >= 150)
(ssel <- as(sel, "sparseMatrix"))
stopifnot(is(sel, "lMatrix"), is(ssel, "lsparseMatrix"),
	  identical3(as.mat(sel.), as.mat(sel), as.mat(ssel)),
	  identical3(!sel, !ssel, nsel), # !<sparse> is typically dense
	  identical3(m[ sel],  m[ ssel], as(m, "matrix")[as( ssel, "matrix")]),
	  identical3(m[!sel],  m[!ssel], as(m, "matrix")[as(!ssel, "matrix")])
	  )
showProc.time()

## more sparse Matrices --------------------------------------

##' @title Check sparseMatrix sub-assignment   m[i,j] <- v
##' @param ms sparse Matrix
##' @param mm its [traditional matrix]-equivalent
##' @param k (approximate) length of index vectors (i,j)
##' @param n.uniq (approximate) number of unique values in  i,j
##' @param vRNG function(n) for random 'v' generation
##' @param show logical; if TRUE, it will not stop on error
##' @return
##' @author Martin Maechler
chkAssign <- function(ms, mm = as(ms, "matrix"),
                      k = min(20,dim(mm)), n.uniq = k %/% 3,
                      vRNG = { if(is.numeric(mm) || is.complex(mm))
                                   function(n) rpois(n,lambda= 0.75)# <- about 47% zeros
                      else ## logical
                          function(n) runif(n) > 0.8 }, ## 80% zeros
                      showOnly=FALSE)
{
    stopifnot(is(ms,"sparseMatrix"))
    d <- dim(ms)
    s1 <- function(n) sample(n, pmin(n, pmax(1, rpois(1, n.uniq))))
    i <- sample(s1(d[1]), k/2+ rpois(1, k/2), replace = TRUE)
    j <- sample(s1(d[2]), k/2+ rpois(1, k/2), replace = TRUE)
    assert.EQ.mat(ms[i,j], mm[i,j])
    ms2 <- ms. <- ms; mm. <- mm # save
    ## now sub*assign* to these repeated indices, and then compare -----
    v <- vRNG(length(i) * length(j))
    mm[i,j] <- v
    ms[i,j] <- v
    ## useful to see (ii,ij), but confusing R/ESS when additionally debugging:
    ## if(!showOnly && interactive()) { op <- options(error = recover); on.exit(options(op)) }
    assert.EQ.mat(ms, mm, showOnly=showOnly)
    ## vector indexing m[cbind(i,j)] == m[i + N(j-1)] ,  N = nrow(.)
    ii <- seq_len(min(length(i), length(j)))
    i <- i[ii]
    j <- j[ii]
    ij <- cbind(i, j)
    ii <- i + nrow(ms)*(j - 1)
    ord.i <- order(ii)
    iio <- ii[ord.i]
    ui <- unique(iio) # compare these with :
    neg.ii <- - setdiff(seq_len(prod(d)), ii)
    stopifnot(identical(mm[ii], mm[ij]),
              identical(ms.[ui], ms.[neg.ii]),
	      ms.[ij] == mm.[ii], ## M[ cbind(i,j) ] was partly broken; now checking
              ms.[ii] == mm.[ii])
    v <- v[seq_len(length(i))]
    if(is(ms,"nMatrix")) v <- as.logical(v)  # !
    mm.[ij] <- v
    ms.[ii] <- v
    nodup <- (length(ui) == length(ii)) ## <==>  ! anyDuplicated(iio)
    if(nodup) { cat("[-]") # rare, unfortunately
	ms2[neg.ii] <- v[ord.i]
	stopifnot(identical(ms2, ms.))
    }
    assert.EQ.mat(ms., mm., showOnly=showOnly)
} ##{chkAssign}

## Get duplicated index {because these are "hard" (and rare)
getDuplIndex <- function(n, k) {
    repeat {
        i <- sample(n, k, replace=TRUE) # 3 4 6 9 2 9 :  9 is twice
        if(anyDuplicated(i)) break
    }
    i
}

suppressWarnings(RNGversion("3.5.0")); set.seed(101)
m <- 1:800
m[sample(800, 600)] <- 0
m0 <- Matrix(m, nrow = 40)
m1 <- add.simpleDimnames(m0)
for(kind in c("n", "l", "d")) {
 for(m in list(m0,m1)) { ## -- with and without dimnames -------------------------
    kClass <-paste0(kind, "Matrix"  )
    Ckind <- paste0(kind, "gCMatrix")
    Tkind <- paste0(kind, "gTMatrix")
    str(mC <- as(as(as(m, kClass), "CsparseMatrix"), "generalMatrix")) # was as(m, Ckind) deprecated
    #was mT <- as(as(as(m, kClass), "TsparseMatrix"), Tkind))
    str(mT <- as(as(as(m, kClass), "generalMatrix"), "TsparseMatrix"))
    mm <- as(mC, "matrix") # also logical or double
    IDENT <- if(kind == "n") function(x,y) Q.eq2(x,y, tol=0) else identical
    stopifnot(exprs = {
        identical(mT, as(mC, "TsparseMatrix"))
        identical(mC, as(mT, "CsparseMatrix")) # was Ckind; now deprecated
        Qidentical(mC[0,0], new(Ckind)) # M..3 Csp..4
        Qidentical(mT[0,0], new(Tkind)) #  "
        identical(unname(mT[0,]), new(Tkind, Dim = c(0L,ncol(m))))# M.3 T.4 C.4
        identical(unname(mT[,0]), new(Tkind, Dim = c(nrow(m),0L)))# M.3 C.4
        is.null(cat("IDENT():\n"))
        IDENT(mC[0,], as(mT[FALSE,], "CsparseMatrix")) # M.3 C.4 M.3 + Tsp..4 Csp..4 | as(., Ckind) deprecated
        IDENT(mC[,0], as(mT[,FALSE], "CsparseMatrix")) # M.3 C.4 M.3 C.4             | as(., Ckind) deprecated
        is.null(cat("sapply(..):\n"))
        sapply(pmin(min(dim(mC)), c(0:2, 5:10)),
               function(k) {i <- seq_len(k); all(mC[i,i] == mT[i,i])})
    })
    cat("ok\n")
    show(mC[,1])
    show(mC[1:2,])
    show(mC[7,  drop = FALSE])
    assert.EQ.mat(mC[1:2,],   mm[1:2,])
    assert.EQ.mat(mC[0,],     mm[0,])
    assert.EQ.mat(mC[,FALSE], mm[,FALSE])
    ##
    ## *repeated* (aka 'duplicated') indices - did not work at all ...
    i <- pmin(nrow(mC), rep(8:10,2))
    j <- c(2:4, 4:3)
    assert.EQ.mat(mC[i,], mm[i,])
    assert.EQ.mat(mC[,j], mm[,j])
    ## FIXME? assert.EQ.mat(mC[,NA], mm[,NA]) -- mC[,NA] is all 0 "instead" of all NA
    ## MM currently thinks we should  NOT  allow  <sparse>[ <NA> ]
    assert.EQ.mat(mC[i, 2:1], mm[i, 2:1])
    assert.EQ.mat(mC[c(4,1,2:1), j], mm[c(4,1,2:1), j])
    assert.EQ.mat(mC[i,j], mm[i,j])
    ##
    ## set.seed(7)
    op <- options(Matrix.verbose = FALSE)
    cat(" for(): ")
    for(n in 1:(if(doExtras) 50 else 5)) { # (as chkAssign() is random)
	chkAssign(mC, mm)
	chkAssign(mC[-3,-2], mm[-3,-2])
        cat(".")
    }
    options(op)
    cat(sprintf("\n[Ok]%s\n\n", strrep("-", 64)))
 }
 cat(sprintf("\nok( %s )\n== ###%s\n\n", kind, strrep("=", 70)))
}## end{for}---------------------------------------------------------------
showProc.time()

if(doExtras) {### {was ./AAA_index.R, MM-only}
    ## an nsparse-example
    A <- Matrix(c(rep(c(1,0,0),2), rep(c(2,0),7), c(0,0,2), rep(0,4)), 3,9)
    i <- c(3,1:2)
    j <- c(3, 5, 9, 5, 9)
    vv <- logical(length(i)*length(j)); vv[6:9] <- TRUE

    print(An <- as(A,"nMatrix")); an <- as(An, "matrix")
    assert.EQ.mat(An, an)
    An[i, j] <- vv
    an[i, j] <- vv
    assert.EQ.mat(An, an)## error
    if(!all(An == an)) show(drop0(An - an))
    ## all are +1

    options("Matrix.subassign.verbose" = TRUE)# output from C
    An <- as(A,"nMatrix"); An[i, j] <- vv
    ## and compare with this:
    Al <- as(A,"lMatrix"); Al[i, j] <- vv
    options("Matrix.subassign.verbose" = FALSE)

    ##--- An interesting not small not large example  for  M[i,j] <- v ------------
    ##
    M <- Matrix(c(1, rep(0,7), 1:4), 3,4)
    N0 <- kronecker(M,M)
    mkN1 <- function(M) {
        stopifnot(length(d <- dim(M)) == 2)
        isC <- is(M,"CsparseMatrix")
        M[,d[2]] <- c(0,2,0)
        N <- kronecker(diag(x=1:2), M)## remains sparse if 'M' is
        if(isC) N <- as(N, "CsparseMatrix")
        diag(N[-1,]) <- -2
        N[9,]  <- 1:4   # is recycled
        N[,12] <- -7:-9 # ditto
        N
    }

    show(N1 <- t(N <- mkN1(N0)))    # transpose {for display reasons}
    C1 <- t(C <- mkN1(as(N0,"CsparseMatrix")))
    stopifnot(all(C == N))
    assert.EQ.mat(C, mkN1(as(N0, "matrix")))

    C. <- C1
    show(N <- N1) ; n <- as(N, "matrix"); str(N)
    sort(i <- c(6,8,19,11,21,20,10,7,12,9,5,18,17,22,13))## == c(5:13, 17:22))
    sort(j <- c(3,8,6,15,10,4,14,13,16,2,11,17,7,5))## == c(2:8, 10:11, 13:17)
    val <- v.l <- 5*c(0,6,0,7,0,0,8:9, 0,0)
    show(spv <- as(val, "sparseVector")); str(spv)

    n [i,j] <- v.l
    N [i,j] <- val# is recycled, too
    C.[i,j] <- val
    assert.EQ.mat(N,n) ; stopifnot(all(C. == N))
    ## and the same *again*:
    n [i,j] <- v.l
    N [i,j] <- val
    C.[i,j] <- val
    assert.EQ.mat(N,n)
    stopifnot(all(C. == N))

    print(load(system.file("external", "symA.rda", package="Matrix"))) # "As"
    stopifnotValid(As, "dsCMatrix"); stopifnot(identical(As@factors, list()))
    R. <- drop0(chol(As))
    stopifnot(exprs = {
        1:32 == sort(diag(R.))  ## !
        R.@x > 0
        R.@x == as.integer(R.@x)## so it is an integer-valued chol-decomp !
        ## shows that (1) As is *not* singular  (2) the matrix is not random
        all.equal(crossprod(R.), As, tolerance=1e-15)
    })
    print(summary(evA <- eigen(As, only.values=TRUE)$values))
    print(tail(evA)) ## largest three ~= 10^7,  smallest two *negative*
    print(rcond(As)) # 1.722 e-21 == very bad !
    ##-> this *is* a border line case, i.e. very close to singular !
    ## and also determinant(.) is rather random here!
    cc0 <- Cholesky(As)# no problem
    try({
        cc <- Cholesky(As, super=TRUE)
        ## gives --on 32-bit only--
        ## Cholmod error 'matrix not positive definite' at file:../Supernodal/t_cholmod_super_numeric.c, line 613
        ecc <- expand(cc)
        L.P <- with(ecc, crossprod(L,P))  ## == L'P
        ## crossprod(L.P) == (L'P)' L'P == P'LL'P
        stopifnot( all.equal(crossprod(L.P), As) )
    })
    ##---- end{ eigen( As ) -----------

} ## only if(doExtras)


##---- Symmetric indexing of symmetric Matrix ----------
m. <- mC
m.[, c(2, 7:12)] <- 0
stopifnotValid(S <- crossprod(add.simpleDimnames(m.) %% 100), "dsCMatrix")
ss <- as(S, "matrix")
ds <- as(S, "denseMatrix")
## NA-indexing of *dense* Matrices: should work as traditionally
assert.EQ.mat(ds[NA,NA], ss[NA,NA])
assert.EQ.mat(ds[NA,  ], ss[NA,])
assert.EQ.mat(ds[  ,NA], ss[,NA])
T <- as(S, "TsparseMatrix")
stopifnot(identical(ds[2 ,NA], ss[2,NA]),
	  identical(ds[NA, 1], ss[NA, 1]),
	  identical(S, as(T, "CsparseMatrix")) )

## non-repeated indices:
i <- c(7:5, 2:4);assert.EQ.mat(T[i,i], ss[i,i])
## NA in indices  -- check that we get a helpful error message:
i[2] <- NA
er <- tryCatch(T[i,i], error = function(e)e)
if(englishMsgs)
stopifnot(as.logical(grep("indices.*sparse Matrices", er$message)))

N <- nrow(T)
set.seed(11)
for(n in 1:(if(doExtras) 50 else 3)) {
    i <- sample(N, max(2, sample(N,1)), replace = FALSE)
    validObject(Tii <- T[i,i]) ; tTi <- t(T)[i,i]
    stopifnot(is(Tii, "dsTMatrix"), # remained symmetric Tsparse
	      is(tTi, "dsTMatrix"), # may not be identical when *sorted* differently
	      identical(as(t(Tii),"CsparseMatrix"), as(tTi,"CsparseMatrix")))
    assert.EQ.mat(Tii, ss[i,i])
}

b <- diag(1:2)[,c(1,1,2,2)]
cb <- crossprod(b)
cB <- crossprod(Matrix(b, sparse=TRUE))
a <- matrix(0, 6, 6)
a[1:4, 1:4] <- cb
A1 <- A2 <- Matrix(0, 6, 6)#-> ddiMatrix
A1[1:4, 1:4] <- cb
A2[1:4, 1:4] <- cB
assert.EQ.mat(A1, a)# indeed
## "must": symmetric and sparse, i.e., ds*Matrix:
stopifnot(identical(A1, A2), is(A1, "dsCMatrix"))

## repeated ones ``the challenge'' (to do smartly):
j <- c(4, 4, 9, 12, 9, 4, 17, 3, 18, 4, 12, 18, 4, 9)
assert.EQ.mat(T[j,j], ss[j,j])
## and another two sets  (a, A) &  (a., A.) :
a <- matrix(0, 6,6)
a[upper.tri(a)] <- (utr <- c(2, 0,-1, 0,0,5, 7,0,0,0, 0,0,-2,0,8))
ta <- t(a); ta[upper.tri(a)] <- utr; a <- t(ta)
diag(a) <- c(0,3,0,4,6,0)
A <- as(Matrix(a), "TsparseMatrix")
A. <- A
diag(A.) <- 10 * (1:6)
a. <- as(A., "matrix")
## More testing {this was not working for a long time..}
set.seed(1)
for(n in 1:(if(doExtras) 100 else 6)) {
    i <- sample(1:nrow(A), 3+2*rpois(1, lambda=3), replace=TRUE)
    Aii  <- A[i,i]
    A.ii <- A.[i,i]
    stopifnot(class(Aii) == class(A),
              class(A.ii) == class(A.))
    assert.EQ.mat(Aii , a [i,i])
    assert.EQ.mat(A.ii, a.[i,i])
    assert.EQ.mat(T[i,i], ss[i,i])
}
showProc.time()

stopifnot(all.equal(mC[,3], mm[,3]),
	  identical(mC[ij], mC[ij + 0.4]),
	  identical(mC[ij], mm[ij]),
	  identical(mC[iN], mm[iN]))
## out of bound indexing must be detected:
assertError(mC[cbind(ij[,1] - 5, ij[,2])])
assertError(mC[cbind(ij[,1],     ij[,2] + ncol(mC))])

assert.EQ.mat(mC[7, , drop=FALSE], mm[7, , drop=FALSE])
identical    (mC[7,   drop=FALSE], mm[7,   drop=FALSE]) # *vector* indexing

stopifnot(dim(mC[numeric(0), ]) == c(0,20), # used to give warnings
          dim(mC[, integer(0)]) == c(40,0),
          identical(mC[, integer(0)], mC[, FALSE]))
validObject(print(mT[,c(2,4)]))
stopifnot(all.equal(mT[2,], mm[2,]),
	  ## row or column indexing in combination with t() :
	  Q.C.identical(mT[2,], t(mT)[,2]),
	  Q.C.identical(mT[-2,], t(t(mT)[,-2])),
	  Q.C.identical(mT[c(2,5),], t(t(mT)[,c(2,5)])) )
assert.EQ.mat(mT[4,, drop = FALSE], mm[4,, drop = FALSE])
stopifnot(identical3(mm[,1], mC[,1], mT[,1]),
	  identical3(mm[3,], mC[3,], mT[3,]),
	  identical3(mT[2,3], mC[2,3], 0),
	  identical(mT[], mT),
          identical4(       mm[c(3,7), 2:4],  as.mat( m[c(3,7), 2:4]),
                     as.mat(mT[c(3,7), 2:4]), as.mat(mC[c(3,7), 2:4]))
          )

x.x <- crossprod(mC)
stopifnot(class(x.x) == "dsCMatrix",
          class(x.x. <- round(x.x / 10000)) == "dsCMatrix",
          identical(x.x[cbind(2:6, 2:6)],
                    diag(x.x[2:6, 2:6], names=FALSE)))
head(x.x.) # Note the *non*-structural 0's printed as "0"
tail(x.x., -3) # all but the first three lines

lx.x <- as(as(x.x, "lMatrix"), "symmetricMatrix") # FALSE only for "structural" 0
(l10 <- lx.x[1:10, 1:10])# "lsC"
(l3 <-  lx.x[1:3, ])
m.x <- as.mat(x.x) # as.mat() *drops* (NULL,NULL) dimnames
stopifnot(class(l10) == "lsCMatrix", # symmetric indexing -> symmetric !
          identical(as.mat(lx.x), m.x != 0),
          identical(as.logical(lx.x), as.logical(m.x)),
          identical(as.mat(l10), m.x[1:10, 1:10] != 0),
          identical(as.mat(l3 ), m.x[1:3, ] != 0)
          )

##-- Sub*assignment* with repeated / duplicated index:
A <- Matrix(0,4,3) ; A[c(1,2,1), 2] <- 1 ; A
B <- A;              B[c(1,2,1), 2] <- 1:3; B; B. <- B
B.[3,] <- rbind(4:2)
## change the diagonal and the upper and lower subdiagonal :
diag(B.) <- 10 * diag(B.)
diag(B.[,-1]) <- 5* diag(B.[,-1])
diag(B.[-1,]) <- 4* diag(B.[-1,]) ; B.
C <- B.; C[,2] <- C[,2];  C[1,] <- C[1,]; C[2:3,2:1] <- C[2:3,2:1]
stopifnot(identical(unname(as(A, "matrix")),
		    local({a <- matrix(0,4,3); a[c(1,2,1), 2] <-  1 ; a})),
	  identical(unname(as(B, "matrix")),
		    local({a <- matrix(0,4,3); a[c(1,2,1), 2] <- 1:3; a})),
	  identical(C, drop0(B.)))
## <sparse>[<logicalSparse>] <- v  failed in the past
T <- as(C,"TsparseMatrix"); C. <- C
T[T>0] <- 21
C[C>0] <- 21
a. <- local({a <- as(C., "matrix"); a[a>0] <- 21; a})
assert.EQ.mat(C, a.)
stopifnot(identical(C, as(T, "CsparseMatrix")))

## used to fail
n <- 5 ## or much larger
sm <- new("dsTMatrix", i=1L, j=1L, Dim=as.integer(c(n,n)), x = 1)
(cm <- as(sm, "CsparseMatrix"))
sm[2,]
stopifnot(sm[2,] == c(0:1, rep.int(0,ncol(sm)-2)),
	  sm[2,] == cm[2,],
	  sm[,3] == sm[3,],
	  all(sm[,-(1:3)] == t(sm[-(1:3),])), # all(<lge.>)
	  all(sm[,-(1:3)] == 0)
	  )
showProc.time()

##---  "nsparse*" sub-assignment :----------
M <- Matrix(c(1, rep(0,7), 1:4), 3,4)
N0 <- kronecker(M,M)
Nn <- as(N0, "nMatrix"); nn <- as(Nn,"matrix")
(Nn00 <- Nn0 <- Nn); nn00 <- nn0 <- nn

set.seed(1)
Nn0 <- Nn00; nn0 <- nn00
for(i in 1:(if(doExtras) 200 else 25)) {
    Nn <- Nn0
    nn <- nn0
    i. <- getDuplIndex(nrow(N0), 6)
    j. <- getDuplIndex(ncol(N0), 4)
    vv <- sample(c(FALSE,TRUE),
                 length(i.)*length(j.), replace=TRUE)
    cat(",")
    Nn[i., j.] <- vv
    nn[i., j.] <- vv
    assert.EQ.mat(Nn, nn)
    if(!all(Nn == nn)) {
        cat("i=",i,":\n i. <- "); dput(i.)
        cat("j. <- "); dput(j.)
        cat("which(vv): "); dput(which(vv))
        cat("Difference matrix:\n")
        show(drop0(Nn - nn))
    }
    cat("k")
    ## sub-assign double precision to logical sparseMatrices: now *with* warning:
    ##  {earlier: gave *no* warning}:
    assertWarning(Nn[1:2,] <- -pi)
    assertWarning(Nn[, 5] <- -pi)
    assertWarning(Nn[2:4, 5:8] <- -pi)
    stopifnotValid(Nn,"nsparseMatrix")
    ##
    cat(".")
    if(i %% 10 == 0) cat("\n")
    if(i == 100) {
        Nn0 <- as(Nn0, "CsparseMatrix")
        cat("Now: class", class(Nn0)," :\n~~~~~~~~~~~~~~~~~\n")
    }
}
showProc.time()
Nn <- Nn0
## Check that  NA is interpreted as TRUE (with a warning), for "nsparseMatrix":
assertWarning(Nn[ii <-     3  ] <- NA); stopifnot(isValid(Nn,"nsparseMatrix"), Nn[ii])
assertWarning(Nn[ii <-   22:24] <- NA); stopifnot(isValid(Nn,"nsparseMatrix"), Nn[ii])
assertWarning(Nn[ii <- -(1:99)] <- NA); stopifnot(isValid(Nn,"nsparseMatrix"), Nn[ii])
assertWarning(Nn[ii <- 3:4  ] <- c(0,NA))
stopifnot(isValid(Nn,"nsparseMatrix"), Nn[ii] == 0:1)
assertWarning(Nn[ii <- 25:27] <- c(0,1,NA))
stopifnot(isValid(Nn,"nsparseMatrix"), Nn[ii] == c(FALSE,TRUE,TRUE))

m0 <- Diagonal(5)
stopifnot(identical(m0[2,], m0[,2]),
	  identical(m0[,1], c(1,0,0,0,0)))
### Diagonal -- Sparse:
(m1 <- as(m0, "TsparseMatrix")) # dtTMatrix unitriangular
(m2 <- as(m0, "CsparseMatrix")) # dtCMatrix unitriangular
m1g <- as(m1, "generalMatrix")
tr1 <- as(m1, "denseMatrix") # dtrMatrix unitriangular
stopifnotValid(m1g, "dgTMatrix")
diag(tr1) <- 100
stopifnot(diag(tr1) == 100)# failed when 'diag<-' did not recycle
assert.EQ.mat(m2[1:3,],    diag(5)[1:3,])
assert.EQ.mat(m2[,c(4,1)], diag(5)[,c(4,1)])
stopifnot(identical(m2[1:3,], as(m1[1:3,], "CsparseMatrix")),
          identical(asUniqueT(m1[, c(4,2)]),
                    asUniqueT(m2[, c(4,2)]))
          )## failed in 0.9975-11

## 0-dimensional diagonal - subsetting ----------------------------
## before that diagU2N() etc for 0-dim. dtC*:
m0. <- m00 <- matrix(numeric(),0,0)
 tC0.<- new("dtCMatrix")
 tC0 <- new("dtCMatrix", diag="U")
(gC0 <- new("dgCMatrix")) # 0 x 0
D0 <- Diagonal(0)
stopifnot(exprs = {
    identical(m0.,       as(tC0, "matrix")) # failed: Cholmod error 'invalid xtype' ..
    identical(numeric(), as(tC0, "numeric"))#   "
    identical(numeric(), tC0[ 0 ])# --> .M.vectorSub(x, i) failed in as(., "matrix")
    identical(m00[TRUE ], tC0[TRUE ])# (worked already)
    identical(m00[FALSE], tC0[FALSE])#    ditto
    ##
    identical(D0, D0[0,0]) # used to fail --> subCsp_ij  (..)
    identical(gC0, D0[, 0]) #  (ditto)     --> subCsp_cols(..)
    identical(gC0, D0[0, ]) #     "        --> subCsp_rows(..)
    identical(D0, D0[,])          # (worked already)
    identical(m00[ 0 ],  D0[ 0 ] )#      ditto
    identical(m00[TRUE ], D0[TRUE ])#      "
    identical(m00[FALSE], D0[FALSE])#      "
    ##
    identical(tC0, tC0[0,0]) # failed --> subCsp_ij  (..)
    identical(gC0, tC0[ ,0]) #   "    --> subCsp_cols(..)
    identical(gC0, tC0[0, ]) #   "    --> subCsp_rows(..)
    identical(tC0, tC0[,])   # (worked already)
    ## vector indexing
})

expr <- quote({ ## FIXME -- both 'TRUE' and 'FALSE'  should fail "out of bound",etc
  D0[TRUE, TRUE ]
  D0[    , TRUE ]
  D0[TRUE,      ] # worked  but should *NOT*
 tC0[TRUE, TRUE ]
 tC0[    , TRUE ]
 tC0[TRUE,      ] # worked  but should *NOT*
  ##
  D0[FALSE,FALSE] # fails --> subCsp_ij(..) -> intI()
  D0[     ,FALSE] #    ditto ............
  D0[FALSE,     ] #    ditto
 tC0[FALSE,FALSE] #      "
 tC0[FALSE,     ] #      "
 tC0[     ,FALSE] #      "
})
EE <- lapply(expr[-1], function(e)
    list(expr = e,
         r    = tryCatch(eval(e), error = identity)))
exR <- lapply(EE, `[[`, "r")
stopifnot(exprs = {
    vapply(exR, inherits, logical(1), what = "error")
    !englishMsgs ||
    unique( vapply(exR, `[[`, "<msg>", "message")
           ) == "logical subscript too long"
})


(uTr <- new("dtTMatrix", Dim = c(3L,3L), diag="U"))
uTr[1,] <- 0
assert.EQ.mat(uTr, cbind(0, rbind(0,diag(2))))

M <- m0; M[1,] <- 0
Z <- m0; Z[] <- 0; z <- array(0, dim(M))
stopifnot(identical(M, Diagonal(x=c(0, rep(1,4)))),
          all(Z == 0), Qidentical(as(Z, "matrix"), z))
M <- m0; M[,3] <- 3 ; M ; stopifnot(is(M, "sparseMatrix"), M[,3] == 3)
checkMatrix(M)
M <- m0; M[1:3, 3] <- 0 ;M
T <- m0; T[1:3, 3] <- 10
stopifnot(identical(M, Diagonal(x=c(1,1, 0, 1,1))),
          isValid(T, "triangularMatrix"), identical(T[,3], c(10,10,10,0,0)))

M <- m1; M[1,] <- 0 ; M ; assert.EQ.mat(M, diag(c(0,rep(1,4))), tol=0)
M <- m1; M[,3] <- 3 ; stopifnot(is(M,"sparseMatrix"), M[,3] == 3)
Z <- m1; Z[] <- 0
checkMatrix(M)
M <- m1; M[1:3, 3] <- 0 ;M
assert.EQ.mat(M, diag(c(1,1, 0, 1,1)), tol=0)
T <- m1; T[1:3, 3] <- 10; checkMatrix(T)
stopifnot(is(T, "triangularMatrix"), identical(T[,3], c(10,10,10,0,0)),
	  Qidentical(as(Z, "matrix"), z))

M <- m2; M[1,] <- 0 ; M ; assert.EQ.mat(M, diag(c(0,rep(1,4))), tol=0)
M <- m2; M[,3] <- 3 ; stopifnot(is(M,"sparseMatrix"), M[,3] == 3)
checkMatrix(M)
Z <- m2; Z[] <- 0
M <- m2; M[1:3, 3] <- 0 ;M
assert.EQ.mat(M, diag(c(1,1, 0, 1,1)), tol=0)
T <- m2; T[1:3, 3] <- 10; checkMatrix(T)
stopifnot(is(T, "dtCMatrix"), identical(T[,3], c(10,10,10,0,0)),
	  Qidentical(as(Z, "matrix"), z))
showProc.time()


## "Vector indices" -------------------
asLogical <- function(x) {
    stopifnot(is.atomic(x))
    storage.mode(x) <- "logical"
    x
}
.iniDiag.example <- expression({
    D <- Diagonal(6)
    M <- as(D,"generalMatrix") # was "dge", now "dgC"
    d <- as(D,"unpackedMatrix")  # "dtr" (unitri)
    m <- as(D,"matrix")
    s <- as(D,"TsparseMatrix"); N <- as(s,"nMatrix")
    S <- as(s,"CsparseMatrix"); C <- as(S,"nMatrix")
})
eval(.iniDiag.example)
i <- c(3,1,6); v <- c(10,15,20)
## (logical,value) which both are recycled:
L <- c(TRUE, rep(FALSE,8)) ; z <- c(50,99)

## vector subassignment, both with integer & logical
## these now work correctly {though not very efficiently; hence warnings}
m[i] <- v # the role model: only first column is affected
M[i] <- v; assert.EQ.mat(M,m) # dgC
D[i] <- v; assert.EQ.mat(D,m) # ddi -> dtC (new! 2019-07; was dgT)
s[i] <- v; assert.EQ.mat(s,m) # dtT -> dgT
S[i] <- v; assert.EQ.mat(S,m); S # dtC -> dtT -> dgT -> dgC
m.L <- asLogical(m) ; assertWarning(
C[i] <- v, verbose=TRUE) # warning: C is nMatrix, v not T/F
assert.EQ.mat(C,m.L); validObject(C); assertWarning(
N[i] <- v, verbose=TRUE)
assert.EQ.mat(N,m.L); validObject(N)
stopifnot(identical(D, as(as(s, "triangularMatrix"), "CsparseMatrix")))
## logical *vector* indexing
eval(.iniDiag.example)
m[L] <- z; m.L <- asLogical(m)
M[L] <- z; assert.EQ.mat(M,m)
D[L] <- z; assert.EQ.mat(D,m)
s[L] <- z; assert.EQ.mat(s,m)
S[L] <- z; assert.EQ.mat(S,m) ; S ; assertWarning(
C[L] <- z, verbose=TRUE); assert.EQ.mat(C,m.L)  ; assertWarning(
N[L] <- z, verbose=TRUE); assert.EQ.mat(N,m.L)


## indexing [i]  vs  [i,] --- now ok
eval(.iniDiag.example)
stopifnot(identical5(m[i], M[i], D[i], s[i], S[i]), identical3(as.logical(m[i]), C[i], N[i]),
          identical5(m[L], M[L], D[L], s[L], S[L]), identical3(as.logical(m[L]), C[L], N[L]))
## bordercase ' drop = .' *vector* indexing {failed till 2009-04-..)
stopifnot(identical5(m[i,drop=FALSE], M[i,drop=FALSE], D[i,drop=FALSE],
		     s[i,drop=FALSE], S[i,drop=FALSE]),
	  identical3(as.logical(m[i,drop=FALSE]),
		     C[i,drop=FALSE], N[i,drop=FALSE]))
stopifnot(identical5(m[L,drop=FALSE], M[L,drop=FALSE], D[L,drop=FALSE],
		     s[L,drop=FALSE], S[L,drop=FALSE]),
	  identical3(as.logical(m[L,drop=FALSE]),
		     C[L,drop=FALSE], N[L,drop=FALSE]))
## using L for row-indexing should give an error
assertError(m[L,]); assertError(m[L,, drop=FALSE])
## these did not signal an error, upto (including) 0.999375-30:
assertError(s[L,]); assertError(s[L,, drop=FALSE])
assertError(S[L,]); assertError(S[L,, drop=FALSE])
assertError(N[L,]); assertError(N[L,, drop=FALSE])

## row indexing:
assert.EQ.mat(D[i,], m[i,])
assert.EQ.mat(M[i,], m[i,])
assert.EQ.mat(s[i,], m[i,])
assert.EQ.mat(S[i,], m[i,])
assert.EQ.mat(C[i,], asLogical(m[i,]))
assert.EQ.mat(N[i,], asLogical(m[i,]))
## column indexing:
assert.EQ.mat(D[,i], m[,i])
assert.EQ.mat(M[,i], m[,i])
assert.EQ.mat(s[,i], m[,i])
assert.EQ.mat(S[,i], m[,i])
assert.EQ.mat(C[,i], asLogical(m[,i]))
assert.EQ.mat(N[,i], asLogical(m[,i]))


### --- negative indices ----------

## 1) negative *vector* indexing
eval(.iniDiag.example)
i <- -(2:30)
stopifnot(identical5(m[i], M[i], D[i], s[i], S[i]),
          identical3(as.logical(m[i]), C[i], N[i]))
##  negative vector subassignment :
v <- seq_along(m[i])
m[i] <- v; m.L <- asLogical(m)
M[i] <- v; assert.EQ.mat(M,m) # dge
D[i] <- v; assert.EQ.mat(D,m) # ddi -> dtT -> dgT
s[i] <- v; assert.EQ.mat(s,m) # dtT -> dgT
S[i] <- v; assert.EQ.mat(S,m); S ; assertWarning( # dtC -> dtT -> dgT -> dgC
N[i] <- v, verbose=TRUE)
assert.EQ.mat(N,m.L); N ; assertWarning(
C[i] <- v, verbose=TRUE)
assert.EQ.mat(C,m.L); C #

options(warn = 2) #----------------------# NO WARNINGS from here -----------------
					## =====================

## 2) negative [i,j] indices
mc <- mC[1:5, 1:7]
mt <- mT[1:5, 1:7]
## sub matrix
assert.EQ.mat(mC[1:2, 0:3], mm[1:2, 0:3]) # test 0-index
stopifnot(identical(mc[-(3:5), 0:2], mC[1:2, 0:2]),
          identical(mt[-(3:5), 0:2], mT[1:2, 0:2]),
          identical(mC[2:3, 4],      mm[2:3, 4]))
assert.EQ.mat(mC[1:2,], mm[1:2,])
## sub vector
stopifnot(identical4(mc[-(1:4), ], mC[5, 1:7],
                     mt[-(1:4), ], mT[5, 1:7]))
stopifnot(identical4(mc[-(1:4), -(2:4)], mC[5, c(1,5:7)],
                     mt[-(1:4), -(2:4)], mT[5, c(1,5:7)]))

## mixing of negative and positive must give error
assertError(mT[-1:1,])
showProc.time()

## Sub *Assignment* ---- now works (partially):
mt0 <- mt
nt <- as(mt, "nMatrix")
mt[1, 4] <- -99
mt[2:3, 1:6] <- 0
mt
m2 <- mt+mt
m2[1,4] <- -200
m2[c(1,3), c(5:6,2)] <- 1:6
stopifnot(m2[1,4] == -200,
          as.vector(m2[c(1,3), c(5:6,2)]) == 1:6)
mt[,3] <- 30
mt[2:3,] <- 250
mt[1:5 %% 2 == 1, 3] <- 0
mt[3:1, 1:7 > 5] <- 0
mt

tt <- as(mt,"matrix")
ii <- c(0,2,5)
jj <- c(2:3,5)
tt[ii, jj] <- 1:6 # 0 is just "dropped"
mt[ii, jj] <- 1:6
assert.EQ.mat(mt, tt)

mt[1:5, 2:6]
as((mt0 - mt)[1:5,], "dsparseMatrix")# [1,5] and lines 2:3

mt[c(2,4), ] <- 0; stopifnot(as(mt[c(2,4), ],"matrix") == 0)
mt[2:3, 4:7] <- 33
checkMatrix(mt)
mt

mc[1,4] <- -99 ; stopifnot(mc[1,4] == -99)
mc[1,4] <-  00 ; stopifnot(mc[1,4] ==  00)
mc[1,4] <- -99 ; stopifnot(mc[1,4] == -99)
mc[1:2,4:3] <- 4:1; stopifnot(as(mc[1:2,4:3], "matrix") == 4:1)

mc[-1, 3] <- -2:1 # 0 should not be entered; 'value' recycled
mt[-1, 3] <- -2:1
stopifnot(mc@x != 0, mt@x != 0,
	  mc[-1,3] == -2:1, mt[-1,3] == -2:1) ## failed earlier

mc0 <- mc
mt0 <- as(mc0, "TsparseMatrix")
m0  <- as(mc0, "matrix")
set.seed(1); options(Matrix.verbose = FALSE)
for(i in 1:(if(doExtras) 50 else 4)) {
    mc <- mc0; mt <- mt0 ; m <- m0
    ev <- 1:5 %% 2 == round(runif(1))# 0 or 1
    j <- sample(ncol(mc), 1 + round(runif(1)))
    nv <- rpois(sum(ev) * length(j), lambda = 1)
    mc[ev, j] <- nv
     m[ev, j] <- nv
    mt[ev, j] <- nv
    if(i %% 10 == 1) print(mc[ev,j, drop = FALSE])
    stopifnot(as.vector(mc[ev, j]) == nv, ## failed earlier...
              as.vector(mt[ev, j]) == nv)
    validObject(mc) ; assert.EQ.mat(mc, m)
    validObject(mt) ; assert.EQ.mat(mt, m)
}
showProc.time()
options(Matrix.verbose = TRUE)

mc # no longer has non-structural zeros
mc[ii, jj] <- 1:6
mc[c(2,5), c(3,5)] <- 3.2
checkMatrix(mc)
m. <- mc
mc[4,] <- 0
mc

S <- as(Diagonal(5),"TsparseMatrix")
H <- Hilbert(9)
Hc <- as(round(H, 3), "CsparseMatrix")# a sparse matrix with no 0 ...
(trH <- tril(Hc[1:5, 1:5]))
stopifnot(is(trH, "triangularMatrix"), trH@uplo == "L",
          is(S, "triangularMatrix"))

## triangular assignment
## the slick (but inefficient in case of sparse!) way to assign sub-diagonals:
## equivalent to tmp <- `diag<-`(S[,-1], -2:1); S[,-1] <- tmp
## which dispatches to (x="TsparseMatrix", i="missing",j="index", value="replValue")
diag(S[,-1]) <- -2:1 # used to give a wrong warning
S <- as(S,"triangularMatrix")
assert.EQ.mat(S, local({s <- diag(5); diag(s[,-1]) <- -2:1; s}))

trH[c(1:2,4), c(2:3,5)] <- 0 # gave an *error* upto Jan.2008
trH[ lower.tri(trH) ] <- 0   # ditto, because of callNextMethod()

m <- Matrix(0+1:28, nrow = 4)
m[-3,c(2,4:5,7)] <- m[ 3, 1:4] <- m[1:3, 6] <- 0
mT <- as(m, "TsparseMatrix")
stopifnot(identical(mT[lower.tri(mT)],
                    m [lower.tri(m) ]))
lM <- upper.tri(mT, diag=TRUE)
mT[lM] <- 0
 m[lM] <- 0
assert.EQ.mat(mT, as(m,"matrix"))
mT[lM] <- -1:0
 m[lM] <- -1:0
assert.EQ.mat(mT, as(m,"matrix"))
(mT <- drop0(mT))

i <- c(1:2, 4, 6:7); j <- c(2:4,6)
H[i,j] <- 0
(H. <- round(as(H, "sparseMatrix"), 3)[ , 2:7])
Hc. <- Hc
Hc.[i,j] <- 0 ## now "works", but setting "non-structural" 0s
stopifnot(as(Hc.[i,j], "matrix") == 0)
Hc.[, 1:6]

## an example that failed for a long time
sy3 <- new("dsyMatrix", Dim = as.integer(c(2, 2)), x = c(14, -1, 2, -7))
checkMatrix(dm <- kronecker(Diagonal(2), sy3))# now sparse with new kronecker
dm <- Matrix(as(dm, "matrix"))# -> "dsyMatrix"
(s2 <- as(dm, "sparseMatrix"))
checkMatrix(st <- as(s2, "TsparseMatrix"))
stopifnot(is(s2, "symmetricMatrix"),
	  is(st, "symmetricMatrix"))
checkMatrix(s.32  <- st[1:3,1:2]) ## 3 x 2 - and *not* dsTMatrix
checkMatrix(s2.32 <- s2[1:3,1:2])
I <- c(1,4:3)
stopifnot(is(s2.32, "generalMatrix"),
          is(s.32,  "generalMatrix"),
          identical(as.mat(s.32), as.mat(s2.32)),
          identical3(dm[1:3,-1], asD(s2[1:3,-1]), asD(st[1:3,-1])),
          identical4(2, dm[4,3], s2[4,3], st[4,3]),
          identical3(diag(dm), diag(s2), diag(st)),
          is((cI <- s2[I,I]), "dsCMatrix"),
          is((tI <- st[I,I]), "dsTMatrix"),
          identical4(as.mat(dm)[I,I], as.mat(dm[I,I]), as.mat(tI), as.mat(cI))
          )

## now sub-assign  and check for consistency
## symmetric subassign should keep symmetry
st[I,I] <- 0; checkMatrix(st); stopifnot(is(st,"symmetricMatrix"))
s2[I,I] <- 0; checkMatrix(s2); stopifnot(is(s2,"symmetricMatrix"))
##
m <- as.mat(st)
 m[2:1,2:1] <- 4:1
st[2:1,2:1] <- 4:1
s2[2:1,2:1] <- 4:1
stopifnot(identical(m, as.mat(st)),
	  1:4 == as.vector(s2[1:2,1:2]),
	  identical(m, as.mat(s2)))

## now a slightly different situation for 's2' (had bug)
s2 <- as(dm, "sparseMatrix")
s2[I,I] <- 0; diag(s2)[2:3] <- -(1:2)
stopifnot(is(s2,"symmetricMatrix"), diag(s2) == c(0:-2,0))
t2 <- as(s2, "TsparseMatrix")
m <- as.mat(s2)
s2[2:1,2:1] <- 4:1
t2[2:1,2:1] <- 4:1
 m[2:1,2:1] <- 4:1
assert.EQ.mat(t2, m)
assert.EQ.mat(s2, m)
## and the same (for a different s2 !)
s2[2:1,2:1] <- 4:1
t2[2:1,2:1] <- 4:1
assert.EQ.mat(t2, m)# ok
assert.EQ.mat(s2, m)# failed in 0.9975-8
showProc.time()

## sub-assign  RsparseMatrix -- Matrix bug [#6709] by David Cortes
## https://r-forge.r-project.org/tracker/?func=detail&atid=294&aid=6709&group_id=61
## simplified by MM
X <- new("dgCMatrix", i = c(0L,3L), p = c(0L,2L,2L,2L), x = c(100, -20), Dim = c(12L,3L))
R <- as(X, "RsparseMatrix")
T <- as(R, "TsparseMatrix")
T[, 2] <- 22 # works fine
R[, 2] <- 22 # failed, as it called replTmat() giving narg() == 3
## now R is Tsparse (as documented on ../man/RsparseMatrix-class.Rd),
identical(R, T) ## but as this may change, rather R & T should have the same *content*
assert.EQ.Mat(R, T)


## m[cbind(i,j)] <- value: (2-column matrix subassignment): -------------------------
m.[ cbind(3:5, 1:3) ] <- 1:3
stopifnot(m.[3,1] == 1, m.[4,2] == 2)
nt. <- nt ; nt[rbind(2:3, 3:4, c(3,3))] <- FALSE
s. <- m. ; m.[cbind(3,4:6)] <- 0 ## assigning 0 where there *is* 0 ..
stopifnot(identical(nt.,nt),       ## should not have changed
	  identical(s., m.))
x.x[ cbind(2:6, 2:6)] <- 12:16
stopifnot(isValid(x.x, "dsCMatrix"),
	  12:16 == as.mat(x.x)[cbind(2:6, 2:6)])
(ne1 <- (mc - m.) != 0)
stopifnot(identical(ne1, 0 != abs(mc - m.)))
(ge <- m. >= mc) # contains "=" -> result is dense
ne. <- mc != m.  # was wrong (+ warning)
stopifnot(identical(!(m. < mc), m. >= mc),
	  identical(m. < mc, as(!ge, "sparseMatrix")),
	  identical(ne., drop0(ne1)))

d6 <- Diagonal(6)
ii <- c(1:2, 4:5)
d6[cbind(ii,ii)] <- 7*ii
stopifnot(is(d6, "ddiMatrix"), identical(d6, Diagonal(x=c(7*1:2,1,7*4:5,1))))

sclass <- function(obj) as.vector(class(obj)) # as.v*(): drop attr(*,"package")
show2cls <- function(C,D, chr = "")
    cat(sprintf("%s & %s%s: %s %s\n",
                deparse(substitute(C)), deparse(substitute(D)), chr,
                sclass(C), sclass(D)))
for(j in 2:6) { ## even and odd j used to behave differently
    cat("j = ", j, ":\n-------\n")
    M <- Matrix(0, j,j); m <- matrix(0, j,j)
    T  <- as(M, "TsparseMatrix")
    TG <- as(T, "generalMatrix")
    G  <- as(M, "generalMatrix");  show2cls(TG, G)
    stopifnot(is(TG, "TsparseMatrix"),
              is(G,  "CsparseMatrix"))
    id <- cbind(1:j,1:j)
    i2 <- cbind(1:j,j:1)
    m[id] <- 1:j
    M[id] <- 1:j
    T[id] <- 1:j ; show2cls(M, T,' ("diag")')
    stopifnot(is(M, "diagonalMatrix"), # since 2019-07 // FIXME (?!) for j=1
              is(T,"triangularMatrix"), isDiagonal(T)) # was "symmetricMatrix"
    G[id] <- 1:j
    TG[id]<- 1:j
    m[i2] <- 10
    M[i2] <- 10
    T[i2] <- 10 ; show2cls(M, T,' ("symm")')
    G[i2] <- 10
    TG[i2]<- 10
    ##
    assert.EQ.mat(M, m)
    assert.EQ.mat(T, m)
    assert.EQ.mat(G, m)
    assert.EQ.mat(TG,m)
}


## drop, triangular, ...
(M3 <- Matrix(upper.tri(matrix(, 3, 3)))) # ltC; indexing used to fail
T3 <- as(M3, "TsparseMatrix")
stopifnot(identical(drop(M3), M3),
	  identical4(drop(M3[,2, drop = FALSE]), M3[,2, drop = TRUE],
		     drop(T3[,2, drop = FALSE]), T3[,2, drop = TRUE]),
	  is(T3, "triangularMatrix"),
	  !is(T3[,2, drop=FALSE], "triangularMatrix")
	  )

(T6 <- as(as(as(kronecker(Matrix(c(0,0,1,0),2,2), t(T3)),
                "lMatrix"), "triangularMatrix"), "TsparseMatrix"))
T6[1:4, -(1:3)] # failed (trying to coerce back to ltTMatrix)
stopifnot(identical(T6[1:4, -(1:3)][2:3, -3],
		    spMatrix(2,2, i=c(1,2,2), j=c(1,1,2), x=rep(TRUE,3))))

M <- Diagonal(4); M[1,2] <- 2
M. <- as(M, "CsparseMatrix")
(R <- as(M., "RsparseMatrix"))
(Ms <- symmpart(M.))
Rs <- as(Ms, "RsparseMatrix")
stopifnot(isValid(M, "triangularMatrix"),
          isValid(M.,"triangularMatrix"),
          isValid(Ms, "dsCMatrix"),
          isValid(R,  "dtRMatrix"),
          isValid(Rs, "dsRMatrix") )
stopifnot(dim(M[2:3, FALSE]) == c(2,0),
          dim(R[2:3, FALSE]) == c(2,0),
          identical(M [2:3,TRUE], M [2:3,]),
          identical(M.[2:3,TRUE], M.[2:3,]),
          identical(R [2:3,TRUE], R [2:3,]),
          dim(R[FALSE, FALSE]) == c(0,0))

n <- 50000L
Lrg <- new("dgTMatrix", Dim = c(n,n))
diag(Lrg) <- 1:n
dLrg <- as(Lrg, "diagonalMatrix")
stopifnot(identical(Diagonal(x = 1:n), dLrg))
diag(dLrg) <- 1 + diag(dLrg)
Clrg <- as(Lrg,"CsparseMatrix")
Ctrg <- as(Clrg, "triangularMatrix")
diag(Ctrg) <- 1 + diag(Ctrg)
stopifnot(identical(Diagonal(x = 1+ 1:n), dLrg),
          identical(Ctrg, as(dLrg,"CsparseMatrix")))

cc <- capture.output(show(dLrg))# show(<diag>) used to error for large n
showProc.time()

## FIXME: "dspMatrix" (symmetric *packed*) not going via "matrix"


## Large Matrix indexing / subassignment
## ------------------------------------- (from ex. by Imran Rashid)
n <- 7000000
m <-  100000
nnz <- 20000
op <- options(Matrix.verbose = 2, warn = 1)

set.seed(12)
f <- sparseMatrix(i = sample(n, size=nnz, replace=TRUE),
                  j = sample(m, size=nnz, replace=TRUE))
str(f)
dim(f) # 6999863 x 99992
prod(dim(f)) # 699930301096 == 699'930'301'096  (~ 700'000 millions)
str(thisCol <-  f[,5000])# logi [~ 7 mio....]
sv <- as(thisCol, "sparseVector")
str(sv) ## "empty" !
validObject(spCol <- f[,5000, drop=FALSE]) # "empty" [n x 1] ngCmatrix
##
## *not* identical(): as(spCol, "sparseVector")@length is "double"prec:
stopifnot(all.equal(as(spCol, "sparseVector"),
                    as(sv,   "nsparseVector"), tolerance=0))
if(interactive())
    selectMethod("[<-", c("ngCMatrix", "missing","numeric", "logical"))
                                        # -> replCmat() in ../R/Csparse.R
f[,5762] <- thisCol # now "fine" and fast thanks to replCmat() --> replCmat4()

fx <- sparseMatrix(i = sample(n, size=nnz, replace=TRUE),
                   j = sample(m, size=nnz, replace=TRUE),
                   x = round(10*rnorm(nnz)))
class(fx)## dgCMatrix
fx[,6000] <- (tC <- rep(thisCol, length.out=nrow(fx)))# fine
thCol <- fx[,2000]
fx[,5762] <- thCol# fine
stopifnot(is(f, "ngCMatrix"), is(fx, "dgCMatrix"),
	  identical(thisCol, f[,5762]),# perfect
	  identical(as.logical(fx[,6000]), tC),
	  identical(thCol,  fx[,5762]))

showProc.time()
options(op)# revert
##
if(doExtras) {#-----------------------------------------------------------------
cat("checkMatrix() of all: \n---------\n")
Sys.setlocale("LC_COLLATE", "C")# to keep ls() reproducible
for(nm in ls()) if(is(.m <- get(nm), "Matrix")) {
    cat(nm, "\n")
    checkMatrix(.m, verbose = FALSE
		, doDet = nm != "As" ## <- "As" almost singular <=> det() "ill posed"
		)
}
showProc.time()
}#--------------end if(doExtras) -----------------------------------------------

## Bugs found by Peter Ralph
n <- 17
x <- Matrix(0, n,n) # "ddiMatrix" now
## x must have at least three nonzero entries
x[1,1] <- x[2,1:2] <- 1.; class(x) # "dtC"
x0 <- x <- as(x,"generalMatrix")  # if x is dgCMatrix, no error
##
z <- matrix(x) # <== not the "Matrix way":  a (n, 1) matrix
z[1] <- 0

x[1:n, 1:n] <- as(z, "sparseVector")
## gave Error: ... invalid subscript type 'S4'
x2 <- x

dim(zC <- as(z, "CsparseMatrix"))
x <- x0
x[] <- zC # did fail, then gave warning.
x1 <- x
##
x <- x0
x[] <- as(zC, "sparseVector") # did fail, too
x2 <- x
stopifnot(identical(x1,x2))
x <- as(x0, "matrix")
x[] <- z
assert.EQ.mat(x1, x)

i <- 4:7
x1 <- x0; x1[cbind(i, i+10)] <- i^2
x2 <- x0; x2[cbind(i, i+10)] <- as(i^2, "matrix")
## failed: nargs() = 4 ... please report

class(x1) # was "dgT", now "dgC"
stopifnot(isValid(x1, class(x1)), identical(x1, x2))
showProc.time()


## check valid indexing (using *random* indices, often duplicated):
chk_dsp_ind <- function(sv, n=512, negI = FALSE, verbose=FALSE) {
    stopifnot(inherits(sv, "dsparseVector"), n >= 1)
    d <- length(sv)
    ## lambda=2 ==> aiming for short 'i' {easier to work with}
    P <- rpois(n, lambda = if(negI) 5 else 2)
    for(i in seq_len(n)) {
        I <-
            if(negI) { # negative indices: 2 are, 4 neither ... always "valid" !!
                k <- max(4L, d - max(1L, P[i]))
                if(verbose) cat(sprintf("size=k = %2d: ", k))
                - sort(sample.int(d, size=k))# replace=FALSE
            }
            else
                sample.int(d, size=1L+P[i], replace=TRUE)
        ##
        validObject(ss <- sv[I]) # Error if not true
    }
    invisible()
}
s <- as(c(3,5,6), "sparseVector")
set.seed(11); chk_dsp_ind(s)
set.seed(3)
(s2 <- as(rsparsematrix(ncol=1, nrow=37, density=1/4),"sparseVector"))
(s3 <- as(rsparsematrix(ncol=1, nrow=64, density=1/4),"sparseVector"))
set.seed(1)
chk_dsp_ind(s2)
chk_dsp_ind(s3)
##
set.seed(47)
## system.time(e.N2 <- chk_dsp_ind(s2, negI=TRUE, verbose=TRUE))
chk_dsp_ind(s2, negI=TRUE)
chk_dsp_ind(s3, negI=TRUE)

iv <- c(rep(0,100), 3, 0,0,7,0,0,0)
sv0  <- sv  <- as(iv, "sparseVector")
sv.0 <- sv. <- as(as.integer(iv), "sparseVector")
stopifnot(canCoerce("integer", "sparseVector"))
sv2 <- as(sv, "isparseVector")
stopifnot(validObject(sv), validObject(sv2), identical(sv., sv2),
          sv == sv.)
n0 <- sv. != 0 # -> is "lsparseV.."

sv [n0] <- sv [n0]
sv.[n0] <- sv.[n0] # gave error
stopifnot(identical(sv , sv0),
          identical(sv., sv.0))
sv [3:7] <- 0
sv.[3:7] <- 0L
stopifnot(identical(sv , sv0), identical(sv., sv.0))
sv [2:4] <- 2:4
sv.[2:4] <- 2:4
stopifnot(which(sv != 0) == (which(sv. != 0) -> in0),
          in0 == c(2:4, 101L, 104L))
sv [2:6] <- 0L
sv.[2:6] <- 0L
stopifnot(identical(sv , sv0), identical(sv., sv.0))

## the next six *all* gave an error -- but should be no-op's:
for(vv in list(sv, sv.0))
    for(ind in list(0, FALSE, logical(length(vv))))
        vv[ind] <- NA
stopifnot(identical(sv , sv0), identical(sv., sv.0))

## <sparseVector>[i] <- val -- failed to resort @i sometimes:  (R-forge Matrix bug #6659)
y1 <- sparseVector(1:3, 13:15, 16)
y2 <- sparseVector(1:6, c(5, 6, 7, 9, 14, 15), 16)
i <- 1:16*12 # 12 24 36 ... 192
x <- sparseVector(numeric(1), 1, length=200)
x[i] <- y1     ; validObject(x[i]) # TRUE
N <- x[i] + y2 ; validObject( N  ) # TRUE
x[i] <- N ## <== bug was here ..
validObject(x)
## gave 'Error' invalid .."dsparseVector".. 'i' must be sorted strictly increasingly
stopifnot(all.equal(x[i] ,  y1+y2, tolerance=0),
		    x[i] == y1+y2)
showProc.time()

if(!interactive()) warnings()

## [matrix-Bugs][#6720] Subsetting with empty indices does not drop -- 17 Apr 2021, by David Cortes
##  https://r-forge.r-project.org/tracker/?func=detail&atid=294&aid=6720&group_id=61

## extended by MM to all versions of "empty" :
x <- c(1,8)
(m1 <-  rbind(x))
m1[]             # remains matrix
m1[,,drop=FALSE] # ditto
m1[,] # [1] 1 2 -- drops (as default drop=TRUE !)

## Sparse Matrix and  actually *any* Matrix-extending class did not work
(M1 <- as(m1, "denseMatrix")) # "dgeMatrix"
S1 <- as(m1, "CsparseMatrix")
R1 <- as(m1, "RsparseMatrix")
stopifnot(exprs = {
    identical(M1[], M1) # remains
    identical(S1[], S1) # remains
    identical(R1[], R1) # remains
    identical(M1[,,drop=FALSE], M1) # ditto
    identical(S1[,,drop=FALSE], S1) #  "
    identical(R1[,,drop=FALSE], R1) #  "
    ## but drop=TRUE  which is the *default* much be obeyed (also for *empty* (i,j):
    identical(m1[,], x)
    identical(M1[,], x) # should drop, but did not
    identical(S1[,], x) #  "
    identical(R1[,], x) #  "
    identical(m1[,,drop=TRUE], x)
    identical(M1[,,drop=TRUE], x) # should drop, but did not
    identical(S1[,,drop=TRUE], x) #  "
    identical(R1[,,drop=TRUE], x)  # "
})


## [matrix-Bugs][#6721] Assignment to 'dgRMatrix' with missing index takes only first element
## MM: This has been fixed already!
X <- rbind(0, 1:3, 0, c(0,1,0))
Rx <- as(X, "RsparseMatrix")
Cx <- as(X, "CsparseMatrix")
X [2,] <- 0
Cx[2,] <- 0
Rx[2,] <- 0
stopifnot(all(Cx == X),
          all(Rx == X))

## [matrix-Bugs][#6745] show(<large sparseVector>)
## NB: is from a bug in head(*); *only* applies to  *empty* sparseV: length(x@i) == 0
op <- options(max.print=999)
( s0 <- sparseVector(i=integer(), length=2^33)) # show -> head() failed in Matrix <= 1.3-*
(xs0 <- sparseVector(i=integer(), length=2^33, x = numeric()))# ditto
options(op); tail(s0) ; tail(xs0) # (always worked)
## *related* bug in `[` --> needed to fix intIv() for such large sparseVectors
stopifnot(exprs = {
    identical(s0[length(s0) - 3:0], # gave Error in if (any(i < 0L)) { : missing value ....
              new("nsparseVector", i=integer(), length=4L))
    identical(xs0[length(s0) - 3:0], # gave Error ..
              new("dsparseVector", i=integer(), length=4L))
})

## Yielded an invalid object in Matrix <= 1.4-1, instead of throwing error
llc04 <- new("dgCMatrix", Dim = c(4L, 0L))
c40 <- new("dgCMatrix", Dim = c(0L, 4L), p = integer(5L))
assertError(c04[1L, ] <- 1)
assertError(c40[, 1L] <- 1)

## Indexing with nMatrix rather than lMatrix
set.seed(3601)
gC <- rsparsematrix(6, 6, 0.6)
gC@x[sample.int(length(gC@x), 6L)] <- NA
ni <- is.na(gC)
li <- as(ni, "lMatrix")
stopifnot(identical(gC[ ni], gC[ li]),
          identical(gC[!ni], gC[!li]))

## Dispatch thinko in Matrix <= 1.5-4
R0 <- R. <-
new("dgRMatrix",
    Dim = c(12L, 12L),
    p = c(0L, 0L, 4L, 5L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 16L),
    j = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 7L, 8L, 4L, 7L, 9L, 9L, 10L, 11L),
    x = as.double(1:16))
R.[1:2, ] <- R.[1:2, ] # was an error
stopifnot(identical(R0, as(R., "RsparseMatrix")))

## Didn't drop dimensions ...
stopifnot(identical(t(as(1:6,"CsparseMatrix"))[TRUE, ], as.double(1:6)))

## Was briefly wrong prior to Matrix 1.6-0
set.seed(0)
S <- new("dsyMatrix", Dim = c(4L, 4L), x = rnorm(16L))
Sii <- S[4:1, 4:1]
stopifnot(exprs = {
    is(Sii, "dsyMatrix")
    Sii@uplo == "L"
    identical(as(Sii, "matrix"), as(S, "matrix")[4:1, 4:1])
})

## Bug #6839: regression in <.s[CT]Matrix>[<logical>] in Matrix 1.6-z
x <- new("dsCMatrix", Dim = c(4L, 4L),
         p = cumsum(0:4), i = sequence(1:4) - 1L, x = as.double(1:10))
i <- c(TRUE, FALSE)
xi <- as(x, "matrix")[i]
for(cl in paste0(c("C", "R", "T"), "sparseMatrix"))
    stopifnot(identical(as(x, cl)[i], xi))
