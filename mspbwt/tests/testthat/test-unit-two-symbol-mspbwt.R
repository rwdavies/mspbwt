##
## this file contains tests of mspbwt vs regular pbwt, only using two symbols in the mspbwt
##


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


test_that("multi-version with 2 symbols can work, when symbols aren't relabelled", {

    out <- test_driver_simple()
    X <- out$X
    Z <- out$Z    

    indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
    top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)

    out <- make_hapMatcherA(X)
    hapMatcherA <- out[["hapMatcherA"]]
    all_symbols <- out[["all_symbols"]]
    Z1 <- map_Z_to_all_symbols(Z, all_symbols)            

    ## checking vs indices only works 
    ms_indices <- build_and_check_indices(hapMatcherA, all_symbols, indices = indices, check_vs_indices = TRUE)
    
    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = hapMatcherA,
        ms_indices = ms_indices,
        Z = Z1,
        verbose = FALSE,
        do_checks = FALSE,
        check_vs_indices = TRUE,
        indices = indices
    )
    expect_equal(top_matches, ms_top_matches)

})


test_that("multi-version with 2 symbols can work, when symbols are relabelled", {

    out <- test_driver_simple()
    X <- out$X
    Z <- out$Z

    ## the small change means that some symbols re-orient
    X <- X[1:150, ]

    indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
    top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)

    out <- make_hapMatcherA(X)
    hapMatcherA <- out[["hapMatcherA"]]
    all_symbols <- out[["all_symbols"]]
    Z1 <- map_Z_to_all_symbols(Z, all_symbols)            

    ## checking vs indices only works 
    ms_indices <- build_and_check_indices(hapMatcherA, all_symbols, check_vs_indices = FALSE)
    
    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = hapMatcherA,
        ms_indices = ms_indices,
        Z = Z1,
        verbose = FALSE,
        do_checks = FALSE,
        check_vs_indices = FALSE,
        indices = indices
    )
    expect_equal(top_matches[, -1], ms_top_matches[, -1])

})


test_that("multi-symbol with 2 symbols can work anywhere", {


    skip("WIP")
    set.seed(14541)
    is_first_run <- TRUE
    
    for(irow in 1:5) {
        for(icol in 1:5) {

            out <- test_driver_simple(irow = irow, icol = icol, K = 65, T = 35)
            X <- out$X
            Z <- out$Z

            indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
            top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)

            ## check the result is there
            check_expected_top_match(top_matches, irow, icol, K, T) 
            
            out <- make_hapMatcherA(X)
            hapMatcherA <- out[["hapMatcherA"]]
            all_symbols <- out[["all_symbols"]]
            Z1 <- map_Z_to_all_symbols(Z, all_symbols)
            
            ms_indices <- build_and_check_indices(hapMatcherA, all_symbols, do_checks = is_first_run, check_vs_indices = is_first_run)
            is_first_run <- FALSE
            
            ms_top_matches <- ms_MatchZ_Algorithm5(
                X = hapMatcherA,
                ms_indices = ms_indices,
                Z = Z1,
                verbose = FALSE,
                do_checks = FALSE,
                check_vs_indices = TRUE,
                indices = indices
            )
            expect_equal(top_matches, ms_top_matches)
            
        }
    }
    
})


