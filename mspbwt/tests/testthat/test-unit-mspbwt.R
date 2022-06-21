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

    ms_indices_with_Rcpp <- ms_BuildIndices_Algorithm5(
        X1C = X1C,
        all_symbols = all_symbols,
        check_vs_indices = TRUE,
        indices = indices,
        with_Rcpp = TRUE
    )
    expect_equal(ms_indices, ms_indices_with_Rcpp)

    ms_indices_only_Rcpp <- Rcpp_ms_BuildIndices_Algorithm5(
        X1C = X1C,
        all_symbols = all_symbols,
        indices = list()
    )
    expect_equal(ms_indices, ms_indices_only_Rcpp)
    
    
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
    

    ## check can stretch without issue
    X <- X * 10 - 5
    out <- make_hapMatcherA(X)
    hapMatcherA <- out[["hapMatcherA"]]
    all_symbols <- out[["all_symbols"]]
    
    ms_indices2 <- ms_BuildIndices_Algorithm5(
        X1C = hapMatcherA,
        all_symbols = all_symbols,
        check_vs_indices = TRUE,
        indices = indices
    )
    
    ms_top_matches2 <- ms_MatchZ_Algorithm5(
        X = X1,
        ms_indices = ms_indices2,
        Z = Z1,
        verbose = FALSE,
        do_checks = FALSE,
        check_vs_indices = TRUE,
        indices = indices
    )
    expect_equal(top_matches, ms_top_matches2)
    
    
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


test_that("can build and test efficient multi-symbol version", {

    skip("not for routine use")
    
    load("~/Downloads/rhb_t_small.RData")
    
    load("~/Downloads/test_nicola/quilt_output/RData/QUILT_prepared_reference.chr20.10000001.12000000.RData")

    K <- nrow(rhb_t)
    to_keep <- c(100, 200, 29963 + 1)
    out <- make_hapMatcherA(rhb_t[to_keep, 1:100])
    hapMatcherA <- out[["hapMatcherA"]]
    all_symbols <- out[["all_symbols"]]
    system.time({    
        ms_indices_only_Rcpp <- Rcpp_ms_BuildIndices_Algorithm5(
            X1C = hapMatcherA,
            all_symbols = all_symbols,
            indices = list(),
            egs = 100,
            n_min_symbols = 100
        )
    })["elapsed"]

    ## can I make the problem smaller
    N <- 15
    ## i1 <- 100
    ## i2 <- 200
    i1 <-1
    i2 <- 2
    Z <- c(hapMatcherA[i1, 1:N], hapMatcherA[i2, (N + 1):ncol(hapMatcherA)])
    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = hapMatcherA,
        ms_indices = ms_indices_only_Rcpp,
        Z = Z,
        verbose = FALSE,
        do_checks = FALSE,
        check_vs_indices = FALSE
    )
    ms_top_matches[ms_top_matches[, "indexB0"] %in% c(i1 - 1, i2 - 1), ]
    ## 

    ## check each one?
    for(i in 1:nrow(ms_top_matches)) {
        s <- ms_top_matches[i, "start1"]
        e <- ms_top_matches[i, "end1"]
        
        ## off on first one?
        Z[s:e]
        hapMatcherA[ms_top_matches[i, "indexB0"] + 1, s:e]
        
        stopifnot(sum(Z[s:e] != hapMatcherA[ms_top_matches[i, "indexB0"] + 1, s:e]) == 0)
    }
    
    ## OK, legit better
    ms_top_matches
    (Z == hapMatcherA[29963 + 1, ])[(N + 1):ncol(hapMatcherA)]
    
    
    cbind(Z, NA, t(hapMatcherA[c(i1, i2), ]))

    

    ##
    ## check speed
    ##
    f <- function(T) {
        K <- nrow(rhb_t)
        out <- make_hapMatcherA(rhb_t[, 1:T])
        hapMatcherA <- out[["hapMatcherA"]]
        all_symbols <- out[["all_symbols"]]
        system.time({    
            ms_indices_only_Rcpp <- Rcpp_ms_BuildIndices_Algorithm5(
                X1C = hapMatcherA,
                all_symbols = all_symbols,
                indices = list(),
                egs = 100,
                n_min_symbols = 100
            )
        })["elapsed"]
    }
    ## mostly linear
    o <- sapply(seq(10, 100, 10), f)

    K <- nrow(rhb_t)
    out <- make_hapMatcherA(rhb_t)
    hapMatcherA <- out[["hapMatcherA"]]
    all_symbols <- out[["all_symbols"]]
    system.time({    
        ms_indices_only_Rcpp <- Rcpp_ms_BuildIndices_Algorithm5(
            X1C = hapMatcherA,
            all_symbols = all_symbols,
            indices = list(),
            egs = 100,
            n_min_symbols = 100
        )
    })["elapsed"]
    
    usge_all <- ms_indices$usge_all

    ## so for whatever reason "u" is now a matrix per-site
    ## d is more straightforward - mostly 0
    ## hence the crazy storage approach
    ## other ones more straightforward? 

    if (1 == 0) {

        ## 
        t <- sapply(ms_indices_only_Rcpp, object.size) / 1e6
        data.frame(cbind(names(t), t))

        object.size(ms_indices_only_Rcpp$a) / 1024 / 1024
        object.size(ms_indices_only_Rcpp$d) / 1024 / 1024        
        
        a <- ms_indices_only_Rcpp$a
        d <- ms_indices_only_Rcpp$d

        ## difference between d's? can I rebuild easily?

        ## what does a look like - some long runs, but not always
        ## store difference somehow?
        a <- table(d, useNA = "always")
        a <- a[order(-a)]
        a <- a[a > 0]
        a <- cbind(as.integer(names(a)), a)
        rownames(a) <- NULL
        colnames(a) <- c("symbol", "count")

        ## OK so d is accessed one column at a time, so want efficient coding of a column
        ## hmm, d really isn't the same thing
        encoded_d <- encode_usg(
            usg = d,
            symbol_count_at_grid = a,
            egs = 100,
            n_min_symbols = 100
        )

    }


    
    ## did this truly work? it didn't find 999?
    
    ## kind of surprised so slow? 
    object.size(usge_all)
    object.size(ms_indices$d)
    
    ## check it can work!
    Z <- c(hapMatcherA[100, ])
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

