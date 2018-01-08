# *popdemo*
R Package for Demographic Modelling Using Projection Matrices

**This is the MASTER page for popdemo**. It includes stable, installable versions of the package. 

popdemo is also available on [CRAN](https://cran.r-project.org). The latest version released to CRAN is [0.2-3](https://cran.r-project.org/web/packages/popdemo). The latest version on CRAN may be behind the latest version on GitHub, so for the newest features you may want to install from GitHub.

To install: Run these lines of code in R.  
The first line installs some packags from CRAN: devtools (used to install from GitHub), and some dependencies.  
The second line installs popdemo from GitHub.
```
install.packages(c("devtools", "expm", "MCMCpack", "markovchain"))
devtools::install_github("iainmstott/popdemo/x.x-x/popdemo") #x.x-x is the desired version number
```

Report bugs and errors to iainmstott@gmail.com
