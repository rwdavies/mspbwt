## these tests are specifically designed to help with intuition and understanding the method
## the primary functionality of these tests is covered elsewhere

## work on intuition of the two
test_that("intuition of building", {

    pdf("~/plots/mspbwt.example.indices.pdf", height = 12, width = 12)
    
    ## build example using 0, 1 for here
    set.seed(10)
    nGrids <- 8
    X <- matrix(as.integer(runif(6 * nGrids) > 0.5), 6, nGrids)
    ## re-ordered so that 0s more frequent than 1s
    X <- apply(X, 2, function(x) {
        a <- table(x + 1)
        as.integer(order(-a)[x + 1]- 1)
    })
    ## 
    indices <- BuildIndices_Algorithm5(
        X,
        verbose = TRUE,
        do_checks = FALSE,
        do_var_check = TRUE
    )

    ## here, first just do the same thing
    ## just to show how similar
    ## then do another example
    for(i_example in 1:2) {

        if (i_example == 1) {
            X1 <- X + 1L
        } else {
            ## make a few random points 3's
            X1 <- X + 1L
            ## change one of the entries from a 2 to a 3
            X1[2, 5] <- 3L
        }

        all_symbols <- list(1:nGrids)
        for (iGrid in 1:nGrids) {
            ## re-label
            a <- table(X1[, iGrid], useNA = "always")
            a <- a[order(-a)]
            a <- a[a > 0]
            names_a <- as.integer(names(a))
            a <- cbind(names_a, a)
            rownames(a) <- NULL
            colnames(a) <- c("symbol", "count")
            all_symbols[[iGrid]] <- a
        }

        ## am here!
        ms_indices <- ms_BuildIndices_Algorithm5(
            X1C = X1,
            all_symbols = all_symbols,
            do_checks = TRUE,
            return_all_usg_check = TRUE,
            verbose = TRUE
        )
        if (i_example == 1) {
            ms_indices1 <- ms_indices    
        } else {
            ms_indices2 <- ms_indices    
        }

    }

    save(indices, ms_indices1, ms_indices2, X, X1, file = "~/temp2.RData")




})


test_that("multi-version with >2 symbols can work, for simple version good for plotting", {

    set.seed(2021)
    
    K <- 15
    nGrids <- 20
    w <- 6
    irow <- 3
    icol <- 3
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
            
    ms_indices <- build_and_check_indices(hapMatcher, all_symbols, check_vs_indices = FALSE)

    make_plot <- TRUE
    ms_top_matches <- ms_MatchZ_Algorithm5(
        X = hapMatcher,
        ms_indices = ms_indices,
        Z = Z,
        ##verbose = TRUE,
        do_checks = FALSE,
        check_vs_indices = FALSE,
        make_plot = make_plot,
        pdfname = paste0("~/temp.simple.", irow, ".", icol, ".ms.pdf")
    )

    ##etm <- exhaustive_top_matches_checker(hapMatcherA, Z, ms_top_matches, return_only = TRUE)
    etm <- exhaustive_top_matches_checker(hapMatcher, Z, ms_top_matches)

    Rcpp_ms_top_matches <- Rcpp_ms_MatchZ_Algorithm5(
        X = hapMatcher,
        ms_indices = ms_indices,
        Z = Z,
        cols_to_use0 = integer(1),
        do_checks = FALSE,
        check_vs_indices = FALSE
    )
    expect_equal(ms_top_matches, Rcpp_ms_top_matches)

})

