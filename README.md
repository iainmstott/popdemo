# popdemo
R Package for Demographic Modelling Using Projection Matrices

**This is the MASTER page for popdemo**. It includes stable, installable versions of the package. 

popdemo is also available on CRAN. The latest version released to CRAN is 0.2-3 and can be found [here](https://cran.r-project.org/web/packages/popdemo). The latest version here on Github may be ahead of the latest version on CRAN, so for the newest features you may want to install from Github.

To install: Run these lines of code in R.  
The first line installs some packags from [CRAN](https://cran.r-project.org/): devtools (used to install from github), and some dependencies.  
The second line installs popdemo from github.
```
install.packages(c("devtools", "expm", "MCMCpack", "markovchain"))
devtools::install_github("iainmstott/popdemo/x.x-x/popdemo") #x.x-x is the version number
```

Report bugs and errors to iainmstott@gmail.com
