## w = 10
get_row_col <- function(irow, icol, T, K, w = 10) {
    row <- c(1, 2, round(K / 2), K - 1, K)[irow]
    c1 <- c(1, 2, round(T / 2), T - w + 1 - 1, T - w + 1)[icol]
    c2 <- c1 + w - 1
    c(row = row, c1 = c1, c2 = c2)
}

##   K <- 250 ## number of haps
##   T <- 50 ## number of SNPs
test_driver_simple <- function(
    K = 250,
    T = 50,
    irow = NULL,
    icol = NULL,
    w = 10
) {
    X <- array(sample(c(0L, 1L), K * T, replace = TRUE), c(K, T))
    Z <- rep(1L, T)
    if (is.null(irow)) {
        ## now add things to match to
        X[c(10, 18, 25, 15, 20, 40) , ] <- 0L
        X[c(15, 20), 1:20] <- 1L
        X[c(10, 18, 25), 11:30] <- 1L
        X[c(40), 26:50] <- 1L
        X[170:250, ] <- 0L
    } else {
        ## choose where to put the match
        a <- get_row_col(irow, icol, T, K, w)
        X[a["row"], ] <- 0L ## blank out
        X[a["row"], a["c1"]:a["c2"]] <- 1L ## make the match!
    }
    X1 <- X + 1L
    Z1 <- Z + 1L
    list(
        X = X,
        Z = Z,
        X1 = X1,
        Z1 = Z1
    )
}


check_expected_top_match <- function(top_matches, irow, icol, K, T, w = 10) {
    a <- get_row_col(irow, icol, T, K, w) 
    i <- match(a["row"], top_matches[, "indexB0"] + 1)
    expect_equal(length(i), 1)
    expect_equivalent(top_matches[i, "start1"], a["c1"])
    expect_equivalent(top_matches[i, "end1"], a["c2"])    
}

## at a minimum, these should all match
check_top_matches <- function(top_matches, X, Z, continue = FALSE) {
    for(irow in 1:nrow(top_matches)) {
        i <- top_matches[irow, "indexB0"] + 1
        c1 <- top_matches[irow, "start1"]
        c2 <- top_matches[irow, "end1"]
        s <- sum(X[i, c1:c2] != Z[c1:c2])
        if (s > 0) {
            print("top match inaccurate")
            print(top_matches[irow, ])
            print(paste0(X[i, c1:c2], collapse = ","))
            print(paste0(Z[c1:c2], collapse = ","))
            if (!continue) {
                stop("top match innacurate")
            }
        }
    }
}

build_and_check_indices <- function(
    hapMatcherA,
    all_symbols,
    indices = NULL,
    check_vs_indices = TRUE,
    do_checks = TRUE
) {
    X1C <- hapMatcherA
    ms_indices <- ms_BuildIndices_Algorithm5(
        X1C = X1C,
        all_symbols = all_symbols,
        check_vs_indices = check_vs_indices,
        indices = indices,
        do_checks = do_checks
    )
    ## check simple version with Rcpp
    ms_indices_with_Rcpp <- ms_BuildIndices_Algorithm5(
        X1C = X1C,
        all_symbols = all_symbols,
        check_vs_indices = check_vs_indices,
        indices = indices,
        with_Rcpp = TRUE
    )
    expect_equal(ms_indices, ms_indices_with_Rcpp)
    ## check complete version with Rcpp
    ms_indices_only_Rcpp <- Rcpp_ms_BuildIndices_Algorithm5(
        X1C = X1C,
        all_symbols = all_symbols,
        indices = list()
    )
    expect_equal(ms_indices, ms_indices_only_Rcpp)
    ms_indices
}
