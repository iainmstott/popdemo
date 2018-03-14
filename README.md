# *popdemo*
R Package for Demographic Modelling Using Projection Matrices.

This GitHub repository contains stable, installable versions of the package. `popdemo` is also available on [CRAN](https://cran.r-project.org/web/packages/popdemo). The latest version on GitHub may be ahead of the latest version on CRAN, so for the newest features you may want to install from GitHub.

To install, run the following lines of R code.  
* The first line installs some packages from CRAN: devtools, plus some dependencies.  
* The second line installs popdemo from GitHub (NOTE remember to replace x.x-x with the desired version number!)
```
install.packages(c("devtools", "expm", "MCMCpack"))
devtools::install_github("iainmstott/popdemo/x.x-x/popdemo") #x.x-x is the desired version number
```
  
  
### *development branch*
The [development](https://github.com/iainmstott/popdemo/tree/development) branch includes a development version of the package source code with new files and features currently being worked on. The development package may or may not be installable.  

The development package is in the *Dev* folder. The package is called *popdemoDev*. A list of possible future changes can be found in *devlist.txt*.  

To install the development package, run the following lines of R code.  
* The first line installs some packages from CRAN: devtools, plus some dependencies.  
* The second line installs popdemoDev from GitHub.
```
install.packages(c("devtools", "expm", "MCMCpack"))
devtools::install_github("iainmstott/popdemo/Dev/popdemoDev")
```

Report bugs and errors to iainmstott@gmail.com
