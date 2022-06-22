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


    




    

test_that("can run multi-symbol version with multiple symbols", {

    ## build 32 SNP version of things
    ## test multiple symbol version on this
    K <- 100 ## number of haps
    T <- 32 * 6
    set.seed(2021)
    X <- array(1L, c(K, T))
    Z <- rep(0L, T)
    
    ## fill in randomly with a few symbols per position that differ from each other
    for(j in 1:6) {    
        ## choose 2-5 symbols
        starts <- sort(sample(5:14, 1 + sample(4, 1), replace = FALSE))
        ends <- sort(sample(16:25, length(starts), replace = FALSE))
        for(i in 1:nrow(X)) {
            k <- sample(length(starts), 1)
            X[i, 32 * (j - 1) + starts[k]:ends[k]] <- 0L
        }
    }
    ## now build a clear path to find
    X[2, 1:64] <- 0L
    X[10, 33:96] <- 0L
    X[17, 97:(6 * 32)] <- 0L


    ## original version
    indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
    top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)
    expect_equal(top_matches[, 2] + 1, c(2, 10, 17))
    expect_equal(top_matches[, "start1"], c(1, 33, 97))
    expect_equal(top_matches[, "end1"], c(64, 96, 6 * 32))

    ## check - why not!
    check_Algorithm5(X, Z, top_matches, display = FALSE)

    rhb_t <- make_rhb_t_from_rhi_t(X)    
    out <- make_hapMatcherA(rhb_t)
    hapMatcherA <- out[["hapMatcherA"]]
    all_symbols <- out[["all_symbols"]]
    
    ms_indices <- ms_BuildIndices_Algorithm5(
        X1C = hapMatcherA,
        all_symbols = all_symbols,
        egs = 5,
        n_min_symbols = 5,
        do_checks = TRUE
    )
    ms_indices_with_Rcpp <- ms_BuildIndices_Algorithm5(
        X1C = hapMatcherA,
        all_symbols = all_symbols,
        egs = 5,
        n_min_symbols = 5,
        with_Rcpp = TRUE,
        do_checks = TRUE
    )
    expect_equal(ms_indices, ms_indices_with_Rcpp)

    ms_indices_only_Rcpp <- Rcpp_ms_BuildIndices_Algorithm5(
        X1C = hapMatcherA,
        all_symbols = all_symbols,
        indices = list(),
        egs = 5,
        n_min_symbols = 5
    )
    expect_equal(ms_indices, ms_indices_only_Rcpp)
    
    
    ## OK here as directly taking symbols
    Z <- c(
        hapMatcherA[2, 1:2],
        hapMatcherA[10, 3],
        hapMatcherA[17, 4:6]
    )
        
    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = hapMatcherA,
        ms_indices = ms_indices,
        Z = Z,
        verbose = FALSE,
        do_checks = FALSE,
        check_vs_indices = FALSE
    )

    expect_equal(top_matches[, 2], ms_top_matches[, 2])
    expect_equal(top_matches[, "start1"], 32 * (ms_top_matches[, "start1"] - 1) + 1)
    expect_equal(top_matches[, "end1"], 32 * (ms_top_matches[, "end1"]))
    
})


