# `popdemo`: Demographic modelling using projection matices. 

`popdemo` provides tools for modelling populations and demography using matrix projection models (MPMs), with deterministic and stochastic model implementations. These tools include population projection, indices of short- and long-term population size and growth, perturbation analysis, convergence to stability or stationarity, and diagnostic and manipulation tools.  

# Installing `popdemo`
`popdemo` is available on [CRAN](https://cran.r-project.org/web/packages/popdemo/index.html) as well as [GitHub](https://github.com/iainmstott/popdemo). The GitHub repository may be ahead of the CRAN version, so for the latest stuff check out GitHub. The GitHub repository also includes a 'development' branch with a development version of the package, which is likely to be unstable but may include new features.

## Installing from [GitHub](https://github.com/iainmstott/popdemo):
```{r eval = FALSE}
# Install dependencies from CRAN:
install.packages(c("devtools", "expm", "MCMCpack", "markovchain"))

# Install stable version from GitHub (recommended):
# NOTE don't forget to change the version number!
devtools::install_github("iainmstott/popdemo/x.x-x/popdemo") #x.x-x is the desired version number

# Install development version 'popdemoDev' (not recommended):
devtools::install_github("iainmstott/popdemo/Dev/popdemoDev", ref = "development")

```
Note that vignettes do not install automatically from GitHub. If you want vignettes to be included, include `build_vignettes = TRUE`.  

## Installing from CRAN:
```{r eval = FALSE}
install.packages("popdemo")
```

Report bugs and errors to iainmstott@gmail.com.