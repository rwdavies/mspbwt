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





test_that("multi-symbol with 2 symbols can work anywhere", {

    set.seed(14541)
    irow <- 1
    icol <- 2
    K <- 65
    T <- 35
    w <- 15 ## width
    

            out <- test_driver_simple(irow = irow, icol = icol, K = K, T = T, w = w)
            X <- out$X
            Z <- out$Z

            indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
            top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)

            ## check the result is there
            check_expected_top_match(top_matches, irow, icol, K, T, w) 
            
            out <- make_hapMatcherA(X)
            hapMatcherA <- out[["hapMatcherA"]]
            all_symbols <- out[["all_symbols"]]
            Z1 <- map_Z_to_all_symbols(Z, all_symbols)
            
            ms_indices <- Rcpp_ms_BuildIndices_Algorithm5(
                X1C = hapMatcherA,
                all_symbols = all_symbols,
                indices = list()
            )
            
            ms_top_matches <- ms_MatchZ_Algorithm5(
                X = hapMatcherA,
                ms_indices = ms_indices,
                Z = Z1,
                verbose = TRUE,
                do_checks = FALSE,
                check_vs_indices = FALSE,
                indices = FALSE
            )

            ## argh, can be swapped
            expect_equal(
                top_matches[order(top_matches[, "indexB0"]), -1],
                ms_top_matches[order(ms_top_matches[, "indexB0"]), -1]
            )
            
})





## key idea throughout this
## do everything in grid form
## then from this, build equivalent SNP one, and test using previous code as well
    

test_that("multi-version with >2 symbols can work", {

    set.seed(2021)

    skip("WIP")

    K <- 95
    nGrids <- 12
    w <- 6 ## width
    irow <- 1
    icol <- 2

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
                verbose = TRUE,
                do_checks = FALSE,
                check_vs_indices = FALSE
            )

            check_expected_top_match(ms_top_matches, irow, icol, K, nGrids, w = w, is_grid_check_snps = FALSE )

})


