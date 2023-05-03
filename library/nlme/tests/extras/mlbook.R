## run reproduction scripts from the installed "mlbook" chapters
testdir <- system.file("mlbook", package = "nlme", mustWork = TRUE)
scripts <- dir(testdir, pattern = "^ch[0-9]*\\.R$")
for(f in scripts) {
    writeLines(c("", strrep("=", nchar(f)), basename(f), strrep("=", nchar(f))))
    set.seed(1)
    options(warn = 1, digits = 5)
    source(file.path(testdir, f), echo = TRUE,
           max.deparse.length = Inf, keep.source = TRUE)
}
