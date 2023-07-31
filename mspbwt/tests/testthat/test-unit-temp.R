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

    K <- 1000
    nGrids <- 50
    ##
    X1C <- matrix(0, K, nGrids)

    ## simulate
    out <- lapply(1:nGrids, function(iGrid) {
        m <- sample(5:10, 1)
        vals <- as.integer(sample(1:m, K, replace = TRUE, prob = (sample(m) ** 2) / sum((1:m) ** 2)))
        ## make sure at least one
        vals[sample(1:m, replace = FALSE)] <- 1:m
        if (iGrid == 3) {
            vals[sample(1:K, 5)] <- 0L
        }
        a <- table(vals)
        ## put 0s at the end if they exist
        if ("0" %in% names(a)) {
            a <- a[c(2:length(a), 1)]
            names(a)[length(a)] <- names(a)[1]
        }
        names_a <- as.integer(names(a))
        a <- cbind(names_a, a)
        rownames(a) <- NULL
        colnames(a) <- c("symbol", "count")
        return(list(vals, a))
    })
    X1C <- sapply(out, function(x) x[[1]])
    all_symbols <- lapply(out, function(x) x[[2]])
    Z <- c(X1C[10, 1:10], X1C[20, -(1:10)])
    
    ## make something big here
    n_min_symbols <- 50
    egs <- 10
    ms_indices <- ms_BuildIndices_Algorithm5(
        X1C = X1C,
        all_symbols = all_symbols,
        return_all_usg_check = TRUE,
        do_checks = TRUE,
        n_min_symbols = n_min_symbols,
        egs = egs
    )

    ## all of these 0-based
    s <- 2
    v2 <- 12
    g <- 1 ## think about this as 0-based
    usge <- ms_indices$usge_all[[g + 1]]
    C <- ms_indices$all_symbols[[g + 1]]
    U <- ms_indices$all_usg_check[[g + 1]]
    ##
    out <- encode_maximal_column_of_u(U[, s], egs = egs, efficient = FALSE)
    out_mat2 <- out[["out_mat"]]
    out_vec2 <- out[["out_vec"]]
    k1 <- get_k_given_matrix_u(s, v2, U)
    which.max(U[, s] == (v2 + 1)) - 1 - 1

    ## exhaustive check! decently slow, but sure, why not
    for(s in 1:nrow(C)) {
        for(v2 in 0:5) {
            ## (C[s, 2] - 1)) {
            ## print(paste0("s = ", s, " v2 = ", v2))
            k1 <- get_k_given_matrix_u(s, v2, U) ## 0-based
            k2 <- get_k_given_encoded_u(s, v2, usge, C, K, egs, do_checks = TRUE, U = U)
            k3 <- Rcpp_get_k_given_encoded_u(s, v2, usge, C, K, egs, do_checks = FALSE, U = matrix(0, 1, 1))
            print(paste0("k1 = ", k1, ", k2 = ", k2, ", k3 = ", k3))
            expect_equal(k1, k2)
            expect_equal(k1, k3)
        }
    }

    
})
