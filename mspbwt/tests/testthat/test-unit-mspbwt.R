if (1 == 0) {

    library("STITCH")
    library("crayon")
    library("testthat")
    library("mspbwt")
    dir <- "~/proj/STITCH/"
    setwd(paste0(dir, "/STITCH/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    dir <- "~/proj/mspbwt/"
    setwd(paste0(dir, "/mspbwt/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)


}


    



## key idea throughout this
## do everything in grid form
## then from this, build equivalent SNP one, and test using previous code as well
    

test_that("multi-version with >2 symbols can work", {

    set.seed(2021)

    K <- 95
    nGrids <- 12
    w <- 6 ## width
    irow <- 1
    icol <- 2

    for(irow in 1:5) {
        for(icol in 1:5) {

            print(irow)
            print(icol)
            print("------")
            
            out <- test_driver_multiple(
                K = K,
                nGrids = 12,
                irow = irow,
                icol = icol,
                w = w
            )

            Xs <- out$Xs
            Zs <- out$Zs
            hapMatcherA <- out$hapMatcherA
            all_symbols <- out$all_symbols
            Z <- out$Z
            
            ## original version
            indices <- BuildIndices_Algorithm5(Xs, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
            top_matches <- MatchZ_Algorithm5(Xs, indices, Zs, verbose = FALSE, do_checks = TRUE)

            check_expected_top_match(top_matches, irow, icol, K, nGrids, w = w, is_grid_check_snps = TRUE)
            exhausive_top_matches_checker(X, top_matches, Z)
                
            if (irow == 3 & icol == 3) {
                ms_indices <- build_and_check_indices(hapMatcherA, all_symbols, check_vs_indices = FALSE)
            } else {
                ms_indices <- Rcpp_ms_BuildIndices_Algorithm5(
                    X1C = hapMatcherA,
                    all_symbols = all_symbols,
                    indices = list()
                )
            }
    
            ms_top_matches <- ms_MatchZ_Algorithm5(
                X = hapMatcherA,
                ms_indices = ms_indices,
                Z = Z,
                ##verbose = TRUE,
                do_checks = FALSE,
                check_vs_indices = FALSE
            )
            ##make_plot = TRUE,
            ##    pdfname = paste0("~/Downloads/tempY.", irow, ".", icol, ".ms.pdf")
            

            check_expected_top_match(ms_top_matches, irow, icol, K, nGrids, w = w, is_grid_check_snps = FALSE )
        }
    }

})


