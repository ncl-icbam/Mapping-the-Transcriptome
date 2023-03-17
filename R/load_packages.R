#' A script to load all required packages and functions
#' 

## Load/ install packages ----
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

### Packages to load ----
pkgCRAN <- c("RColorBrewer",  # Data visualisation - colours
             "cols4all",      # Data visualisation - colours
             "tidyverse",     # Data manipulation 
             "dplyr",         # Data manipulation
             "spdep",         # Geocomputation
             "sf",            # Geocomputation
             "GWmodel",       # Geocomputation
             "Seurat",        # Data normalisation
             "flextable",     # Data visualisation - tables
             "scales",        # Data visualisation - scales
             "ggpubr",        # Data visualisation - pub. ready plots
             "tmap",          # Data visualisation - maps
             "cluster",       # Clustering methods
             "grid",          # Data visualisation
             "gridExtra"      # Data visualisation
             )

pkgGit <- c("RachelQueen1/SCFunctionsV3",
            "mtennekes/cols4all",
            "MonashBioinformaticsPlatform/varistran")

### Load or install&load packages ----
pkg.check <- lapply(
  pkgCRAN, 
  function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

pkg.check <- lapply(
  pkgGit, 
  function(x) {
    pkg.name <- sub(".*/", "", x)
    if (!require(pkg.name, character.only = TRUE)) {
      devtools::install_git(paste0("https://github.com/", x),
                            force = TRUE)
      library(pkg.name, character.only = TRUE)
    }
  }
)

source("./R/moran.test.ste.R")
