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

            Rcpp_ms_top_matches <- Rcpp_ms_MatchZ_Algorithm5(
                X = hapMatcher,
                ms_indices = ms_indices,
                Z = Z,
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
            ms_indices = ms_indices,
            Z = Z1
        )
        expect_equal(ms_top_matches, Rcpp_ms_top_matches)

        
    }
    
})
