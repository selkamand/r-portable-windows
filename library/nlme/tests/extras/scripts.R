## run reproduction scripts from the NLME book chapters
testdir <- system.file("scripts", package = "nlme", mustWork = TRUE)
scripts <- dir(testdir, pattern = "^ch[0-9]*\\.R$")
for(f in scripts) {
    writeLines(c("", strrep("=", nchar(f)), basename(f), strrep("=", nchar(f))))
    set.seed(3)
    options(warn = 1)  # chapters set digits
    source(file.path(testdir, f), echo = TRUE,
           max.deparse.length = Inf, keep.source = TRUE)
}
