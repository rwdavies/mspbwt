#!/usr/bin/env Rscript

required_packages <- c("Rcpp", "data.table", "devtools", "RcppArmadillo")
for(package in required_packages) {
    if (!suppressPackageStartupMessages(require(package, character.only = TRUE))) {
        out <- install.packages(package, repos="http://cran.rstudio.com/")
        out <- require(package, character.only = TRUE)
        if (!out) {
            stop(paste0("Failed to install package:", package))
        }
    }
}
if (!suppressPackageStartupMessages(require("rrbgen")))
    install.packages("https://github.com/rwdavies/rrbgen/releases/download/0.0.6/rrbgen_0.0.6.tar.gz", repos=NULL)

if (!suppressPackageStartupMessages(require("STITCH")))
    install.packages("https://github.com/rwdavies/STITCH/releases/download/1.6.10/STITCH_1.6.10.tar.gz", repos=NULL)

