# README #

### What is this repository for? ###

This repo contains a portable version of R (normal R installed with no registry modifications -- see instructions below) packaged with an Rprofile.site file downloaded from: [**Rportable**](https://sourceforge.net/projects/rportable/) updated to ensure all libraries are installed locally (in r-portable-windows/library). Realistically this is a one-line file :) so could just be created yourself based on the instructions below (and is what I'd reccomend doing if changing the version of R)

### How do I get set up? ###

1. Clone into any directory and run R with r-portable-windows/bin/R
2. Install packages (will save to r-portable-windows/library)

### Changing Version of R ###
1. Simply download any version of R
  * install to whatever directory you want (e.g. on an external hdd)
  * When running the installer:
    * **DO NOT** Create a start menu folder
    * **DO NOT** enable storing Version Number In Registry or Associating files with it.
2. Create a file 'Rprofile.site' inside the etc folder of the R install
3. Open the file: `etc/Rprofile.site`
5. Add the line: `.libPaths(.Library)`
6. Open bin/R and run `.libPaths()`. Shoud only print:
  * `{path_to_repo}/r-portable-windows/library`

### Current R Version
R 4.0.3

### How do I install packages?

You can install binary packages you can install using `install.packages('packagename')` as normal.

Packages that require compilation from source, may require Rtools to be installed.
Rtools does not come with the Rportable (its quite large). If you're trying to put together a R-portable with some piece of software preinstalled, my advice is to install rtools on your system, and install all the source packages you'll need. Rtools is only needed for the original compilation, and so you R portable stays portable.

### Who do I talk to? ###
* Repo owner: selkamand@ccia.org.au
