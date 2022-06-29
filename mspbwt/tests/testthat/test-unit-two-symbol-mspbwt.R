##
## this file contains tests of mspbwt vs regular pbwt, only using two symbols in the mspbwt
##


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


test_that("multi-version with 2 symbols can work, when symbols aren't relabelled", {

    out <- test_driver_simple()
    X <- out$X
    Z <- out$Z    

    indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
    top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)

    out <- QUILT::make_rhb_t_equality(
        rhb_t = X,
        nSNPs = ncol(X) * 32,
        nMaxDH = 20,
        ref_error = 0.001
    )
    hapMatcher <- out[["hapMatcher"]]
    all_symbols <- out[["all_symbols"]]
    Z1 <- map_Z_to_all_symbols(Z, all_symbols)            

    ## checking vs indices only works 
    ms_indices <- build_and_check_indices(hapMatcher, all_symbols, indices = indices, check_vs_indices = TRUE)
    
    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = hapMatcher,
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
    check_top_matches(top_matches, X, Z)     

    out <- QUILT::make_rhb_t_equality(
        rhb_t = X,
        nSNPs = ncol(X) * 32,
        nMaxDH = 20,
        ref_error = 0.001
    )
    hapMatcher <- out[["hapMatcher"]]
    all_symbols <- out[["all_symbols"]]
    Z1 <- map_Z_to_all_symbols(Z, all_symbols)            

    ## checking vs indices only works 
    ms_indices <- build_and_check_indices(hapMatcher, all_symbols, check_vs_indices = FALSE)

    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = hapMatcher,
        ms_indices = ms_indices,
        Z = Z1,
        verbose = FALSE,
        do_checks = FALSE,
        check_vs_indices = FALSE,
        indices = indices
    )
    expect_equal(top_matches[, -1], ms_top_matches[, -1])

    check_top_matches(ms_top_matches, X, Z)         
    
})


test_that("multi-symbol with 2 symbols can work anywhere", {

    set.seed(14541)
    irow <- 1
    icol <- 2
    K <- 65
    T <- 35
    w <- 15 ## width
    
    for(irow in 1:5) {
        for(icol in 1:5) {

            out <- test_driver_simple(irow = irow, icol = icol, K = K, T = T, w = w)
            X <- out$X
            Z <- out$Z

            indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
            top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)

            ## check the result is there
            check_expected_top_match(top_matches, irow, icol, K, T, w) 

            out <- QUILT::make_rhb_t_equality(
                rhb_t = X,
                nSNPs = ncol(X) * 32,
                nMaxDH = 20,
                ref_error = 0.001
             )
            hapMatcher <- out[["hapMatcher"]]
            all_symbols <- out[["all_symbols"]]
            Z1 <- map_Z_to_all_symbols(Z, all_symbols)
            
            ms_indices <- Rcpp_ms_BuildIndices_Algorithm5(
                X1C = hapMatcher,
                all_symbols = all_symbols,
                indices = list()
            )
            
            ms_top_matches <- ms_MatchZ_Algorithm5(
                X = hapMatcher,
                ms_indices = ms_indices,
                Z = Z1,
                verbose = FALSE,
                do_checks = FALSE,
                check_vs_indices = FALSE,
                indices = FALSE
            )

            ## argh, can be swapped
            expect_equal(
                top_matches[order(top_matches[, "indexB0"]), -1],
                ms_top_matches[order(ms_top_matches[, "indexB0"]), -1]
            )
            
        }
    }
    
})





test_that("multi-version with 2 symbols can capture clean breaks", {

    X <- array(0L, c(100, 30))
    X[3, 1:10] <- 1L
    X[2, 11:20] <- 1L
    X[1, 21:30] <- 1L
    Z <- rep(1L, 30)
    
    indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
    top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)

    out <- QUILT::make_rhb_t_equality(
        rhb_t = X,
        nSNPs = ncol(X) * 32,
        nMaxDH = 20,
        ref_error = 0.001
    )
    hapMatcher <- out[["hapMatcher"]]
    all_symbols <- out[["all_symbols"]]
    Z1 <- map_Z_to_all_symbols(Z, all_symbols)            

    ## checking vs indices only works 
    ms_indices <- build_and_check_indices(hapMatcher, all_symbols, indices = indices, check_vs_indices = FALSE)
    
    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = hapMatcher,
        ms_indices = ms_indices,
        Z = Z1,
        verbose = FALSE,
        do_checks = FALSE,
        check_vs_indices = TRUE,
        indices = indices
    )
    expect_equal(top_matches, ms_top_matches)

})




