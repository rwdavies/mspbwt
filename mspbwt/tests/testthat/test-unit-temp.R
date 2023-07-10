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
    ## test in R

}



test_that("can avoid use of d", {

    set.seed(2029)

    K <- 90
    nGrids <- 50

    ## s = SNP
    out <- test_driver_multiple(
        K = K,
        nGrids = nGrids
    )
    Xs <- out$Xs
    Zs <- out$Zs
    hapMatcher <- out$hapMatcher
    all_symbols <- out$all_symbols
    Z <- out$Z
    nSNPs <- length(Zs)
    hapMatcherR <- array(as.raw(0), dim(hapMatcher))
    for(i in 1:ncol(hapMatcher)) {
        hapMatcherR[, i] <- as.raw(hapMatcher[, i])
    }

    ## do on only some of them
    which_grids <- seq(1, nGrids, 3)

    ## build indices (just do one)
    ms_indices <- Rcpp_ms_BuildIndices_Algorithm5(
        X1C = hapMatcher[, which_grids],
        all_symbols = all_symbols[which_grids],
        indices = list()
    )

    ## Z here
    Z_local <- map_Z_to_all_symbols(Z[which_grids], ms_indices[["all_symbols"]])

    ## test in R
    R_results <- ms_MatchZ_Algorithm5(
        X = hapMatcher[, which_grids],
        ms_indices = ms_indices,
        Z = Z_local
    )

    ## check
    for(i_k in 1:nrow(R_results)) {
        k1 <- R_results[i_k, "indexB0"] + 1
        start1 <- R_results[i_k, "start1"]
        end1 <- R_results[i_k, "end1"]
        w <- which_grids[start1:end1]
        sum(hapMatcher[k1, w] != Z_local[w])
    }

    ## do test with both ways
    R_results_test <- ms_MatchZ_Algorithm5(
        X = hapMatcher[, which_grids],
        ms_indices = ms_indices,
        Z = Z_local,
        test_d = TRUE
    )

    expect_equal(R_results, R_results_test)

    ## confirm again, now without d
    ms_indices_no_d <- ms_indices
    ms_indices_no_d[["d"]] <- NULL

    ## do test with both ways
    R_results_no_d <- ms_MatchZ_Algorithm5(
        X = hapMatcher[, which_grids],
        ms_indices = ms_indices_no_d,
        Z = Z_local
    )

    expect_equal(R_results, R_results_no_d)

    ## now for Rcpp, do the same thing
    ## test both use of double condition, and simpler

    for(i_test in 1:2) {

        if (i_test == 1) {
            X <- hapMatcher[, which_grids]
            XR <- matrix(as.raw(0), 1, 1)
            cols_to_use0 <- as.integer(1)
            use_cols_to_use0 <- FALSE
            use_XR <- FALSE
        } else {
            X <- matrix(as.integer(1), 1, 1)
            XR <- hapMatcherR
            cols_to_use0 <- as.integer(which_grids - 1)
            use_cols_to_use0 <- TRUE
            use_XR <- TRUE
        }

        ## test one with d
        Rcpp_results_test <- Rcpp_ms_MatchZ_Algorithm5(
            ms_indices = ms_indices,
            X = X,
            XR = XR,
            cols_to_use0 = cols_to_use0,
            use_cols_to_use0 = use_cols_to_use0,
            use_XR = use_XR,            
            Z = Z_local,
            test_d = TRUE,
            have_d = TRUE,
            verbose = FALSE
        )
        
        expect_equal(R_results, Rcpp_results_test)

        ## other test
        ## confirm again, now without d
        ms_indices_tiny_d <- ms_indices
        ms_indices_tiny_d[["d"]] <- matrix(as.integer(1), 1, 1)
        
        Rcpp_results_test2 <- Rcpp_ms_MatchZ_Algorithm5(
            ms_indices = ms_indices_tiny_d,
            X = X,
            XR = XR,
            cols_to_use0 = cols_to_use0,
            use_cols_to_use0 = use_cols_to_use0,
            use_XR = use_XR,
            Z = Z_local,
            test_d = FALSE,
            have_d = FALSE,
            verbose = FALSE
        )
        
        expect_equal(R_results, Rcpp_results_test2)

    }
    
})
