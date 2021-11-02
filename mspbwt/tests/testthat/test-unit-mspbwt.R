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


test_that("can run multi-symbol version with 2 symbols", {

    K <- 250 ## number of haps
    T <- 50 ## number of SNPs
    set.seed(2021)
    X <- array(sample(c(0L, 1L), K * T, replace = TRUE), c(K, T))
    Z <- rep(1L, T)
    ## now add things to match to
    X[c(10, 18, 25, 15, 20, 40) , ] <- 0L
    X[c(15, 20), 1:20] <- 1L
    X[c(10, 18, 25), 11:30] <- 1L
    X[c(40), 26:50] <- 1L
    X[170:250, ] <- 0L

    ## change symbols so 0 is always the more common one
    indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
    top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)

    ## now - do multiple symbol version on this
    X1 <- X + 1L
    Z1 <- Z + 1L

    ## include in this repo?
    out <- make_hapMatcherA(X)
    hapMatcherA <- out[["hapMatcherA"]]
    all_symbols <- out[["all_symbols"]]
    ## but for checks - want to keep the same order, so swap if not the same
    ## for(iGrid in 1:ncol(hapMatcherA)) {
    ##     if (sum(abs(hapMatcherA[, iGrid] - X1[, iGrid])) > 0) {
    ##         hapMatcherA[, iGrid] <- 3L - hapMatcherA[, iGrid]
    ##         all_symbols[[iGrid]][1:2, ] <- all_symbols[[iGrid]][2:1, ]
    ##     }
    ## }
    ## stopifnot(sum(X1 != hapMatcherA) == 0)
    
    ## now build indices
    X1C <- hapMatcherA
    
    ms_indices <- ms_BuildIndices_Algorithm5(
        X1C = X1C,
        all_symbols = all_symbols,
        check_vs_indices = TRUE,
        indices = indices
    )

    ## this should still work, 
    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = X1,
        ms_indices = ms_indices,
        Z = Z1,
        verbose = FALSE,
        do_checks = FALSE,
        check_vs_indices = TRUE,
        indices = indices
    )
    
    expect_equal(top_matches, ms_top_matches)
    
    
})


    

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


test_that("can build and test efficient multi-symbol version", {

    skip("not for routine use")
    load("~/Download/rhb_t_small.RData")
    rhb_t <- a[, 1:20]
    ## can I just use hapMatcher, and special lookup?

    ## simple hapMatcher
    out <- make_hapMatcherA(rhb_t)
    hapMatcherA <- out[["hapMatcherA"]]
    all_symbols <- out[["all_symbols"]]

    ##
    ms_indices <- ms_BuildIndices_Algorithm5(
        X1C = hapMatcherA,
        all_symbols = all_symbols,
        egs = 100,
        n_min_symbols = 100
    )
    usge_all <- ms_indices$usge_all

    ## kind of surprised so slow? 
    object.size(usge_all)
    object.size(ms_indices$d)
    
    ## check it can work!
    Z <- c(hapMatcherA[100, 1:20])
    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = hapMatcherA,
        ms_indices = ms_indices,
        Z = Z,
        verbose = FALSE,
        do_checks = FALSE,
        check_vs_indices = FALSE
    )
    
    ## does this make sense?
    
    object.size(ms_indices$a) 
    object.size(usge_all) ## aw yeah
    
    ## maybe
    iGrid <- 100
    usL_encoding <- encode_usL(
        usA[[iGrid]],
        symbol_count_at_grid = all_symbols[[iGrid]],
        egs = 256,
        n_min_symbols = 100
    )
    object.size(usL_encoding)
    object.size(1:nrow(usA[[iGrid]]))
    ## OK, so pretty similar to just what a will be

    ## so how many does this require
    iGrid <- 100
    all_symbols[[iGrid]]
    ## so originally, requires 10000 * many integers
    
    sum(sapply(1:8, function(i) {
        out <- encode_column_of_u(usA[[iGrid]][, i], egs = 1000)
    }))



    
    

    ## 
    diff(usA[[iGrid]][, i])

    ## or should I just use encoding / decoding?
    ## store every 32nd value, then store 31 ups/downs?
    
    ## so requires ~1/10th - good, not amazing?
    
    
    lapply(ms_indices, dim)

    b <- hapMatcherA[ms_indices$a[, 4] + 1, ]
    ## nrow(unique(b)) wow only 222 here already
    
})

