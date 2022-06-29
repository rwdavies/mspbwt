if (1 == 0) {

    library("STITCH")
    library("crayon")
    library("testthat")
    library("mspbwt")
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

    K <- 95
    nGrids <- 12
    w <- 6 ## width
    irow <- 2
    icol <- 4

    for(irow in 1:5) {
        for(icol in 1:5) {

            set.seed(2021)
            
            out <- test_driver_multiple(
                K = K,
                nGrids = nGrids,
                irow = irow,
                icol = icol,
                w = w
            )

            Xs <- out$Xs
            Zs <- out$Zs
            hapMatcher <- out$hapMatcher
            all_symbols <- out$all_symbols
            Z <- out$Z
            
            ## original version
            indices <- BuildIndices_Algorithm5(Xs, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
            top_matches <- MatchZ_Algorithm5(Xs, indices, Zs, verbose = FALSE, do_checks = TRUE)
            ## check_expected_top_match(top_matches, irow, icol, K, nGrids, w = w, is_grid_check_snps = TRUE)
            
            ## etm <- exhaustive_top_matches_checker(Xs, Zs, top_matches, return_only = TRUE)
            etm <- exhaustive_top_matches_checker(Xs, Zs, top_matches)
            
            if (irow == 3 & icol == 3) {
                ms_indices <- build_and_check_indices(hapMatcher, all_symbols, check_vs_indices = FALSE)
            } else {
                ms_indices <- Rcpp_ms_BuildIndices_Algorithm5(
                    X1C = hapMatcher,
                    all_symbols = all_symbols,
                    indices = list()
                )
            }
    
            ms_top_matches <- ms_MatchZ_Algorithm5(
                X = hapMatcher,
                ms_indices = ms_indices,
                Z = Z,
                ##verbose = TRUE,
                do_checks = FALSE,
                check_vs_indices = FALSE
            )
            ##
            ##make_plot = TRUE,
            ##    pdfname = paste0("~/Downloads/tempY.", irow, ".", icol, ".ms.pdf")
            

            ##etm <- exhaustive_top_matches_checker(hapMatcherA, Z, ms_top_matches, return_only = TRUE)
            etm <- exhaustive_top_matches_checker(hapMatcher, Z, ms_top_matches)
            
        }
    }

})


