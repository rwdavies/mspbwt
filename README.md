# Multi Symbol Positional Burrows Wheeler Transform 

msPBWT was principally designed to facilitate imputation using QUILT2. QUILT2 is available at [https://github.com/rwdavies/QUILT](https://github.com/rwdavies/QUILT). More details about msPBWT are available in the QUILT2 pre-print available at [https://doi.org/10.1101/2024.07.18.604149](https://doi.org/10.1101/2024.07.18.604149).

A detailed README will follow. 

## High level functions
- ms_encode
- ms_index
- ms_find

## Installation

A release can be installed using code like
```
wget https://github.com/rwdavies/mspbwt/releases/download/0.1.0/mspbwt_0.1.0.tar.gz
R CMD INSTALL mspbwt_0.1.0.tar.gz
```

Alternatively msPBWT can be installed directly from github using code like

``` r
## install.package("pak")
pak::pkg_install("rwdavies/mspbwt/mspbwt")
```

## TODO

- make encoding of `d` an option
- make functions to access `u` easily
