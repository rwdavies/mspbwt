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

            make_plot <- FALSE
            ms_top_matches <- ms_MatchZ_Algorithm5(
                X = hapMatcher,
                ms_indices = ms_indices,
                Z = Z,
                ##verbose = TRUE,
                do_checks = FALSE,
                check_vs_indices = FALSE,
                make_plot = make_plot,
                pdfname = paste0("~/tempY.", irow, ".", icol, ".ms.pdf")
            )
            ##
            ##
            ##


            ##etm <- exhaustive_top_matches_checker(hapMatcherA, Z, ms_top_matches, return_only = TRUE)
            etm <- exhaustive_top_matches_checker(hapMatcher, Z, ms_top_matches)

            Rcpp_ms_top_matches <- Rcpp_ms_MatchZ_Algorithm5(
                X = hapMatcher,
                XR = matrix(raw(0), 1, 1),
                use_XR = FALSE,
                ms_indices = ms_indices,
                Z = Z,
                cols_to_use0 = integer(1),
                do_checks = FALSE,
                check_vs_indices = FALSE
            )
            expect_equal(ms_top_matches, Rcpp_ms_top_matches)

        }
    }

})



## not currently needed
## test_that("can remap rhb_t special value to symbols", {

##     nSNPs <- 32

##     for(nSNPs in c(32, 10)) {

##         symbols <- c(3L, 0L, 15L) ## all 0, two 1's, 4 1's

##         ## 0 =  0000
##         ## 3 =  1100
##         ## 15 = 1111
##         ## now because of order, 3 > 0 > 15 in terms of preference
##         ## e.g. 1 will go to 3

##         test_match     <- c(1L, 7L, 4L, 14L)
##         expected_match <- c(3L, 3L, 0L, 15L)

##         distinctHapsB <- matrix(symbols, nrow = length(symbols), ncol = 1)
##         iGrid1 <- 1

##         expect_equal(
##             sapply(test_match, map_one_binary_value_to_hapMatcher, distinctHapsB, iGrid1, nSNPs),
##             match(expected_match, symbols)
##         )

##     }

## })



##
## here we want to add in a rare symbol
## not covered in hapMatcher
## but contained in Z as well
##
test_that("mspbwt can work with rare symbol not in hapMatcher", {

    set.seed(100110)
    nGrids <- 12
    K <- 40
    X <- array(sample(1L:5L, K * nGrids, replace = TRUE), c(K, nGrids))

    ## for the key point, force it to store this

    for(iKeyGrid in c(6, 1, nGrids)) {

        nMaxDH <- 4
        X[, iKeyGrid] <- rep(1L:4L, K)[1:K]
        X[20:40, iKeyGrid] <- 2000L + 100L:120L

        ## add matches to 3rd and 12th haplotype, which are the same (in hapMatcher), but different underlying
        ## by the current strategy, these will both be captured
        w1 <- max(2, iKeyGrid - 3)
        w2 <- min(nGrids - 1, iKeyGrid + 3)
        w <- w1:w2
        X[12, w] <- X[3, w]
        X[3, iKeyGrid] <- 100L
        X[12, iKeyGrid] <- 1000L
        Z <- X[3, ]
        Z[-w] <- X[1, -w]
        X[12, c(w1 - 1, w2 + 1)] <- 5 - Z[c(w1 - 1, w2 + 1)]
        X[3, c(w1 - 1, w2 + 1)] <- 5 - Z[c(w1 - 1, w2 + 1)]

        ##
        rhb_t <- make_rhb_t_from_rhi_t(X)
        ## QUILT::
        out <- make_rhb_t_equality(
            rhb_t = X,
            nSNPs = 32 * nGrids,
            nMaxDH = nMaxDH,
            ref_error = 0.001
        )
        hapMatcher <- out$hapMatcher
        all_symbols <- out$all_symbols
        ## this now requires re-writing quite a lot about how indices are built
        ## including C++ etc
        ## might need to go back to more fundementals for this!
        ms_indices <- ms_BuildIndices_Algorithm5(
            X1C = hapMatcher,
            all_symbols = all_symbols,
            indices = list(),
            egs = 10,
            n_min_symbols = 3
        )
        ms_indices[["all_usg_check"]] <- NULL

        ## OK - have captured the problem!
        Rcpp_ms_indices <- Rcpp_ms_BuildIndices_Algorithm5(
            X1C = hapMatcher,
            all_symbols = all_symbols,
            indices = list(),
            verbose = FALSE,
            egs = 10,
            n_min_symbols = 3
        )

        expect_equal(ms_indices, Rcpp_ms_indices)

        Z1 <- map_Z_to_all_symbols(Z, all_symbols)
        rbind(X[c(1, 3, 12), ], NA, Z)
        rbind(hapMatcher[c(1, 3, 12), ], NA, Z1)

        ms_top_matches <- ms_MatchZ_Algorithm5(
            X = hapMatcher,
            ms_indices = ms_indices,
            Z = Z1
        )

        ## check that both expected matches are there
        ## we expect 0-based 2 and 11 to match from 1-based SNPs 4 through 8
        a <- ms_top_matches
        expect_equal(sum((a[, 2] == 2) & a[, 3] <= w1 & a[, 4] >=w2), 1)
        expect_equal(sum((a[, 2] == 11) & a[, 3] <= w1 & a[, 4] >= w2), 1)

        Rcpp_ms_top_matches <- Rcpp_ms_MatchZ_Algorithm5(
            X = hapMatcher,
            XR = matrix(raw(0), 1, 1),
            use_XR = FALSE,
            ms_indices = ms_indices,
            Z = Z1,
            cols_to_use0 = integer(1)
        )
        expect_equal(ms_top_matches, Rcpp_ms_top_matches)


    }

})





test_that("can work with different intervals", {

    set.seed(2028)

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

    nWindows <- 3
    ms_indices_multiple <- lapply(1:nWindows, function(iWindow) {
        which_grids <- seq(1, nGrids, nWindows)
        ms_indices <- Rcpp_ms_BuildIndices_Algorithm5(
            X1C = hapMatcher[, which_grids],
            all_symbols = all_symbols[which_grids],
            indices = list()
        )
        ms_indices
    })


    ## check results are the same when run the two ways
    ## 1 = default way
    ## 2 = efficient way
    i_window <- 1

    for(i_window in 1:nWindows) {

        ## default way
        which_grids <- seq(1, nGrids, nWindows)
        Z_local <- map_Z_to_all_symbols(Z[which_grids], ms_indices_multiple[[i_window]][["all_symbols"]])

        ms_top_matches <- Rcpp_ms_MatchZ_Algorithm5(
            X = hapMatcher[, which_grids],
            XR = matrix(raw(0), 1, 1),
            ms_indices = ms_indices_multiple[[i_window]],
            Z = Z_local,
            cols_to_use0 = integer(1),
            verbose = FALSE
        )

        ## efficient way
        ms_top_matches2 <- Rcpp_ms_MatchZ_Algorithm5(
            X = hapMatcher,
            XR = matrix(raw(0), 1, 1),
            ms_indices = ms_indices_multiple[[i_window]],
            Z = Z_local,
            cols_to_use0 = as.integer(which_grids - 1L),
            use_cols_to_use0 = TRUE,
            verbose = FALSE
        )

        expect_equal(
            ms_top_matches,
            ms_top_matches2
        )

    }


})



test_that("can do scan up and down", {

    set.seed(2028)

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
        Z = Z_local,
        do_up_and_down_scan = TRUE,
        mspbwtL = 3,
        mspbwtM = 3
    )

    ## check
    for(i_k in 1:nrow(R_results)) {
        k1 <- R_results[i_k, "index0"] + 1
        end1 <- R_results[i_k, "end1"]
        len1 <- R_results[i_k, "len1"]
        start1 <- end1 - len1 + 1
        w <- which_grids[start1:end1]
        sum(hapMatcher[k1, w] != Z_local[w])
        ##expect_equal(0, sum(hapMatcher[k1, w] != Z_local[w]))
    }

    ## X the simple way
    Rcpp1_results <- Rcpp_ms_MatchZ_Algorithm5(
        X = hapMatcher[, which_grids],
        XR = matrix(as.raw(0), 1, 1),
        ms_indices = ms_indices,
        Z = Z_local,
        do_up_and_down_scan = TRUE,
        mspbwtL = 3,
        mspbwtM = 3,
        cols_to_use0 = integer(1)
    )

    ## XR the simple way
    Rcpp2_results <- Rcpp_ms_MatchZ_Algorithm5(
        X = matrix(0, 1, 1),
        XR = hapMatcherR[, which_grids],
        ms_indices = ms_indices,
        Z = Z_local,
        use_XR = TRUE,
        do_up_and_down_scan = TRUE,
        mspbwtL = 3,
        mspbwtM = 3,
        cols_to_use0 = integer(1)
    )

    ## X the cols0 way
    Rcpp3_results <- Rcpp_ms_MatchZ_Algorithm5(
        X = hapMatcher,
        XR = matrix(as.raw(0), 1, 1),
        ms_indices = ms_indices,
        Z = Z_local,
        do_up_and_down_scan = TRUE,
        mspbwtL = 3,
        mspbwtM = 3,
        cols_to_use0 = as.integer(which_grids - 1),
        use_cols_to_use0 = TRUE,
        verbose = FALSE
    )
    
    ## XR the cols0 way
    Rcpp4_results <- Rcpp_ms_MatchZ_Algorithm5(
        X = matrix(0, 1, 1),
        XR = hapMatcherR,
        ms_indices = ms_indices,
        Z = Z_local,
        use_XR = TRUE,
        do_up_and_down_scan = TRUE,
        mspbwtL = 3,
        mspbwtM = 3,
        cols_to_use0 = as.integer(which_grids - 1),
        use_cols_to_use0 = TRUE
    )

    ## AM HERE
    ## fix these!
    expect_equal(R_results, Rcpp1_results)    
    expect_equal(Rcpp1_results, Rcpp2_results)
    expect_equal(Rcpp1_results, Rcpp3_results)
    expect_equal(Rcpp1_results, Rcpp4_results)    

    ## print("----R-----")
    ## print(R_results)
    ## print("----Rcpp1-----")    
    ## print(Rcpp1_results)
    ## print("----Rcpp2-----")    
    ## print(Rcpp2_results)
    ## print("----Rcpp3-----")        
    ## print(Rcpp3_results)
    ## print("----Rcpp4-----")        
    ## print(Rcpp4_results)        
    
})








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

    ## try removing d
    ms_indices_no_d <- ms_indices
    ms_indices_no_d[["d"]] <- NULL
    
    ## R_results_no_d <- ms_MatchZ_Algorithm5(
    ##     X = hapMatcher[, which_grids],
    ##     ms_indices = ms_indices_no_d,
    ##     Z = Z_local
    ## )

    ## expect_equal(R_results, R_results_no_d)
    
})
