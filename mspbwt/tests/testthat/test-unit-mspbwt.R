if (1 == 0) {

    library("STITCH")
    library("crayon")
    library("testthat")
    library("mspbwt")
    dir <- "~/proj/QUILT/"
    setwd(paste0(dir, "/QUILT/R"))
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

    K <- 200 ## number of haps
    T <- 50 ## number of SNPs
    set.seed(2021)
    X <- array(sample(c(0L, 1L), K * T, replace = TRUE), c(K, T))
    Z <- rep(1L, T)
    ## now add things to match to
    X[c(10, 18, 25, 15, 20, 40) , ] <- 0L
    X[c(15, 20), 1:20] <- 1L
    X[c(10, 18, 25), 11:30] <- 1L
    X[c(40), 26:50] <- 1L

    indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
    top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)

    ## now - do multiple symbol version on this
    X1 <- X + 1L
    Z1 <- Z + 1L

    ms_indices <- ms_BuildIndices_Algorithm5(
        X1,
        check_vs_indices = TRUE,
        indices = indices
    )
    
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

    ## now - build symbol version
    rhb_t <- make_rhb_t_from_rhi_t(X)
    ## now build multi-symbol PBWT
    out <- QUILT::make_rhb_t_equality(
        rhb_t = rhb_t,
        nSNPs = 32 * 6,
        ref_error = 1e-3,
        nMaxDH = 100,
        verbose = FALSE
    )
    hapMatcher <- out$hapMatcher

    ms_indices <- ms_BuildIndices_Algorithm5(
        X = hapMatcher
    )

    Z <- c(
        hapMatcher[2, 1:2],
        hapMatcher[10, 3],
        hapMatcher[17, 4:6]
    )
        
    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = hapMatcher,
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



## 

## AM HERE
## IMPLEMENT ALTERNATIVE VERSION THAT IS MORE EFFICIENT
## CONSIDER GETTING LARGER DATA TO TEST ON?
