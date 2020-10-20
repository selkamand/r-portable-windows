# README #


### What is this repository for? ###

This repo is simply Rportable (https://sourceforge.net/projects/rportable/) with `r-portable-windows/etc/Rprofile.site` updated to ensure all libraries are installed locally (in r-portable-windows/library)

### How do I get set up? ###

1. Clone into any directory and run R with r-portable-windows/bin/R
2. Install packages (will save to r-portable-windows/library)

### Changing Version of R ###
1. Simply download any version of Rportable
2. Open the file: `etc/Rprofile.site`
3. Add the line: `.libPaths(.Library)`
4. Open bin/R and run `.libPaths()`. Shoud only print:
  * `{path_to_repo}/r-portable-windows/library`

### Current R Version
R 3.6.3


**NOTE**: There is a newer version of R-portable available (R 4.0.0) but it has known issues with Rcpp on Windows


### Who do I talk to? ###

* Repo owner: selkamand@ccia.org.au
