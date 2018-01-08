# *popdemo*
R Package for Demographic Modelling Using Projection Matrices

**This is the DEVELOPMENT page for popdemo**. It includes a development version of the package source code with new files and features currently being worked on. The development package may or may not be installable.

**The development package is folder *x.x-x_dev*, where x.x-x is the most recent version number**.

A list of possible future changes can be found in *devlist.txt*.

To install, run the following lines of R code.  
The first line installs some packages from CRAN: devtools (used to install from GitHub), and some dependencies.  
The second line installs popdemo from GitHub.
```
install.packages(c("devtools", "expm", "MCMCpack", "markovchain"))
devtools::install_github("iainmstott/popdemo/x.x-x_dev/popdemo") #x.x-x is the most recent version number
```

Report bugs and errors to iainmstott@gmail.com
