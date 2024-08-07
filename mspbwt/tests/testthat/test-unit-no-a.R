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



test_that("no a initial evaluation", {

    ## from some simulated data
    ## a good way to gain experience, walk through small data
    ms_indices <- list(a = structure(c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 1L, 2L,
    5L, 3L, 7L, 0L, 4L, 6L, 2L, 5L, 7L, 6L, 1L, 0L, 3L, 4L, 2L, 7L,
    6L, 4L, 1L, 0L, 3L, 5L, 2L, 7L, 6L, 5L, 4L, 0L, 1L, 3L, 7L, 5L,
    4L, 0L, 1L, 3L, 2L, 6L, 7L, 0L, 1L, 3L, 6L, 5L, 4L, 2L), dim = 8:7),
    d = structure(c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 2L, 0L,
    0L, 1L, 0L, 1L, 0L, 1L, 2L, 2L, 0L, 1L, 1L, 2L, 1L, 2L, 2L,
    3L, 3L, 1L, 1L, 2L, 3L, 1L, 2L, 3L, 4L, 4L, 1L, 1L, 3L, 4L,
    3L, 4L, 4L, 5L, 5L, 3L, 4L, 3L, 4L, 4L, 5L, 5L, 6L, 6L, 4L,
    4L, 4L, 5L, 6L, 4L, 5L, 7L), dim = c(9L, 7L)), usge_all = list(
    list(c(2, 3, 6), c(4, 8), c(1, 5), 7), list(c(2, 3, 5,
        8), c(1, 6), 4, 7), list(c(1, 3, 4, 8), c(5, 6, 7), 2),
        list(c(1, 2, 3, 8), c(4, 6), 5, 7), list(c(2, 4, 5, 6,
        7, 8), 1, 3), list(c(1, 4, 5, 6, 8), c(2, 3, 7))), egs = 100,
    n_min_symbols = 100, all_symbols = list(structure(c(5L, 7L,
    8L, 9L, 3L, 2L, 2L, 1L), dim = c(4L, 2L), dimnames = list(
        NULL, c("symbol", "count"))), structure(c(10L, 11L, 13L,
    14L, 4L, 2L, 1L, 1L), dim = c(4L, 2L), dimnames = list(NULL,
        c("symbol", "count"))), structure(c(0L, 4L, 10L, 4L,
    3L, 1L), dim = 3:2, dimnames = list(NULL, c("symbol", "count"
    ))), structure(c(2L, 3L, 9L, 15L, 4L, 2L, 1L, 1L), dim = c(4L,
    2L), dimnames = list(NULL, c("symbol", "count"))), structure(c(3L,
    7L, 13L, 6L, 1L, 1L), dim = 3:2, dimnames = list(NULL, c("symbol",
    "count"))), structure(c(0L, 12L, 5L, 3L), dim = c(2L, 2L), dimnames = list(
        NULL, c("symbol", "count")))), all_usg_check = list(structure(c(0,
    0, 1, 2, 2, 2, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 0, 1,
    1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1), dim = c(9L,
    4L)), structure(c(0, 0, 1, 2, 2, 3, 3, 3, 4, 0, 1, 1, 1,
    1, 1, 2, 2, 2, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 1, 1), dim = c(9L, 4L)), structure(c(0, 1, 1, 2, 3,
    3, 3, 3, 4, 0, 0, 0, 0, 0, 1, 2, 3, 3, 0, 0, 1, 1, 1, 1,
    1, 1, 1), dim = c(9L, 3L)), structure(c(0, 1, 2, 3, 3, 3,
    3, 3, 4, 0, 0, 0, 0, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 1, 1,
    1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1), dim = c(9L, 4L)), structure(c(0,
    0, 1, 1, 2, 3, 4, 5, 6, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 1), dim = c(9L, 3L)), structure(c(0, 1,
    1, 1, 2, 3, 4, 4, 5, 0, 0, 1, 2, 2, 2, 2, 3, 3), dim = c(9L,
                                                             2L))))

    K <- sum(ms_indices$all_symbols[[1]][, 2])
    ##A_first_col <- 0:(K - 1) ## WOW this is easy
    A_last_col <- ms_indices$a[, ncol(ms_indices$a)]
    all_usg_check <- ms_indices$all_usg_check
    all_symbols <- ms_indices$all_symbols
    usge_all <- ms_indices$usge_all
    egs <- ms_indices$egs
    n_min_symbols <- ms_indices$n_min_symbols
    K <- sum(ms_indices$all_symbols[[1]][, 2])

    ## backward
    g_in <- 2
    for(v_in in 0:7) {
        index <- find_index_backward(
            g_in = g_in,
            v_in = v_in,
            all_usg_check = all_usg_check,
            all_symbols = all_symbols,
            use_U = TRUE
        )
        v <- v_in
        expect_equal(ms_indices$a[v + 1, g_in + 1 + 1], index)
    }

    ## test using compressed version, works
    for(v_in in 0:7) {
        index <- find_index_backward(
            g_in = g_in,
            v_in = v_in,
            all_symbols = all_symbols,
            usge_all = usge_all,
            egs = egs,
            use_U = FALSE
        )
        v <- v_in
        expect_equal(ms_indices$a[v + 1, g_in + 1 + 1], index)
    }

    ## forward
    g_in <- 2
    for(k_in in 0:7) {
        index <- find_index_forward(g_in, k_in, all_usg_check, all_symbols, A_last_col)
        k <- k_in
        expect_equal(ms_indices$a[k + 1, g_in + 1 + 1], index)
    }

    ## test using Rcpp version, works
    for(v_in in 0:7) {
        index <- Rcpp_find_index_backward(
            g_in = g_in,
            v_in = v_in,
            all_symbols = all_symbols,
            usge_all = usge_all,
            egs = egs,
            K = K,
            list_of_columns_of_A = list(),
            use_list_of_columns_of_A = FALSE
        )
        v <- v_in
        expect_equal(ms_indices$a[v + 1, g_in + 1 + 1], index)
    }
    
    ## test using compressed version, with columns
    list_of_columns_of_A <- as.list(1:ncol(ms_indices$a))
    list_of_columns_of_A[[3]] <- ms_indices$a[, 3]
    ##
    g_in <- 2    
    for(v_in in 0:7) {
        index <- find_index_backward(
            g_in = g_in,
            v_in = v_in,
            all_symbols = all_symbols,
            usge_all = usge_all,
            egs = egs,
            use_U = FALSE,
            list_of_columns_of_A = list_of_columns_of_A
        )
        v <- v_in
        expect_equal(ms_indices$a[v + 1, g_in + 1 + 1], index)
    }

    ## same as above, but using Rcpp
    g_in <- 2    
    for(v_in in 0:7) {
        index <- Rcpp_find_index_backward(
            g_in = g_in,
            v_in = v_in,
            all_symbols = all_symbols,
            usge_all = usge_all,
            egs = egs,
            K = K,
            list_of_columns_of_A = list_of_columns_of_A,
            use_list_of_columns_of_A = TRUE
        )
        v <- v_in
        expect_equal(ms_indices$a[v + 1, g_in + 1 + 1], index)
    }


})




test_that("no a slightly larger experiments", {

    set.seed(10)
    ## for(seed in 1:100) {
    ##     print(paste0("seed = ", seed))
    ##     set.seed(seed)
        
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
        a <- cbind(names_a, as.integer(a))
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
    print("start exhaustive check")
    for(s in 1:nrow(C)) {
        for(v2 in 0:(C[s, 2] - 1)) {
            k1 <- get_k_given_matrix_u(s, v2, U) ## 0-based
            k2 <- get_k_given_encoded_u(s, v2, usge, C, K, egs, do_checks = TRUE, U = U)
            k3 <- Rcpp_get_k_given_encoded_u(s, v2, usge, C, K, egs)
            expect_equal(k1, k2)
            expect_equal(k1, k3)
            ## print("we gooooood")
        }
    }
    print("end exhaustive check")

    hapMatcherR <- matrix(as.raw(0), nrow(X1C), ncol(X1C))
    for(icol in 1:ncol(hapMatcherR)) {
        hapMatcherR[, icol] <- as.raw(X1C[, icol])
    }

    ## check this 
    f <- get_f_given_Z(Z, ms_indices$all_symbols, ms_indices$usge_all, ms_indices$egs)
    f2 <- Rcpp_get_f_given_Z(Z, ms_indices$all_symbols, ms_indices$usge_all, ms_indices$egs)
    expect_equal(f, f2)
    
    out <-  find_good_matches_without_a(
        Z = Z,
        all_symbols = ms_indices$all_symbols,
        usge_all = ms_indices$usge_all,
        egs = ms_indices$egs,
        pbwtL = 10,
        pbwtM = 2,
        hapMatcherR = hapMatcherR,
        do_checks = TRUE,
        A = ms_indices$a
    )
    mat_out <- out[["mat_out"]]

    ## check it has what we expect given construction
    i <- which(mat_out[, "index"] == 19)
    expect_true(length(i) > 0)
    i <- which(mat_out[, "index"] == 9)
    expect_true(length(i) > 0)

    ## check version that uses some rcpp
    out <-  find_good_matches_without_a(
        Z = Z,
        all_symbols = ms_indices$all_symbols,
        usge_all = ms_indices$usge_all,
        egs = ms_indices$egs,
        pbwtL = 10,
        pbwtM = 2,
        hapMatcherR = hapMatcherR,
        use_rcpp = TRUE
    )
    mat_out2 <- out

    ## check it is the same
    expect_equal(mat_out, mat_out2)

    ## check version that uses only Rcpp
    out <-  Rcpp_find_good_matches_without_a(
        Z = Z,
        all_symbols = ms_indices$all_symbols,
        usge_all = ms_indices$usge_all,
        egs = ms_indices$egs,
        pbwtL = 10L,
        pbwtM = 2L,
        K = as.integer(sum(ms_indices$all_symbols[[1]][, 2])),
        hapMatcherR = hapMatcherR,
        which_snps_in_hapMatcherR = 1L:ncol(hapMatcherR),
        list_of_columns_of_A = list(),
        use_list_of_columns_of_A = FALSE,
        verbose = FALSE
    )
    mat_out3 <- out
    expect_equal(mat_out, mat_out3)

})


test_that("explore", {

    skip("WIP")

    ## so OK quite a bit better?
    ## could drop even further if encode decode was say raw, and egs was 256
    ## could try to make it smaller? less RAM maybe

    load("/data/smew1/rdavies/tempAtempAtempAaA.RData")

    ##for(i in 1:4) {
    ##    ms_indices[[i]]$a <- NULL
   ## }
    ##object.size(ms_indices) / 1024 / 1024

    pbwtM <- 5
    pbwtL <- 10
    
    all_symbols <- ms_indices[[1]]$all_symbols
    usge_all <- ms_indices[[1]]$usge_all
    egs <- ms_indices[[1]]$egs

    ## fake one
    Z <- as.integer(hapMatcherR[100000, seq(1, ncol(hapMatcherR), 4)])
    K <- nrow(hapMatcherR)
    A <- ms_indices[[1]]$a

    which_snps_in_hapMatcherR <- seq(1, ncol(hapMatcherR), length(ms_indices))

    ## keep first few, then rest
    list_of_columns_of_A <- as.list(1:ncol(ms_indices[[1]]$a))
    for(i in seq(1, ncol(ms_indices[[1]]$a), 10)) {
        list_of_columns_of_A[[i]] <- ms_indices[[1]]$a[, i]
    }

    
    out <-  find_good_matches_without_a(
        Z = Z,
        all_symbols = ms_indices[[1]]$all_symbols,
        usge_all = ms_indices[[1]]$usge_all,
        egs = ms_indices[[1]]$egs,
        pbwtL = 20,
        pbwtM = 5,
        hapMatcherR = hapMatcherR,
        which_snps_in_hapMatcherR = which_snps_in_hapMatcherR,
        A = ms_indices[[1]]$a,
        do_checks = TRUE,
        verbose = TRUE,
        list_of_columns_of_A = list_of_columns_of_A,
        use_rcpp = TRUE
    )

    mat_out <- out[["mat_out"]]
    list_of_mats <- out$list_of_mats
    
    ## check that they work, yup, good
    for(i_row in 1:nrow(mat_out)) {
        g <- mat_out[i_row, 1]
        k <- mat_out[i_row, 2]
        len <- mat_out[i_row, 3]
        w <- (g + 1):(g + len)
        ## 
        stopifnot(sum(as.integer(hapMatcherR[k + 1, which_snps_in_hapMatcherR[w]]) != Z[w]) == 0)
    }

    ##
    list_of_mats[[114]][[1]]    
    list_of_mats[[115]][[1]]    

    out <- NULL
    for(g in 80:116) {
        out <- rbind(out, list_of_mats[[g - 1]][[1]][, "v"])
        out <- rbind(out, list_of_mats[[g]][[1]][, "k"])        
    }
    
    
m1 <- sapply(list_of_mats, function(x) x[[1]][, c("v")])
    m2 <- sapply(list_of_mats, function(x) x[[1]][, c("k")])    
    
    ## check it out here


    
    ## compare the best matches?
    ## can probably do pretty small pbwtL honestly
    for(i in 1:5) {
        x <- out[[i]]
        print(head(x[order(-x[, 3]), ], 10))
    }
    
    ## ## seems to work
    
    ## mat_out[mat_out[, 2] == 99999, ]
    
    ## ## check matches? seems good! not sure if this should necessarily be perfect always
    ## t(sapply(1:nrow(mat_out), function(i_row) {
    ##     g <- mat_out[i_row, "g"]
    ##     index <- mat_out[i_row, "index"]
    ##     l <- mat_out[i_row, "len"]
    ##     a <- as.integer(hapMatcherR[index + 1, seq(1, 461, 4)])[g:(g + l) + 1] == Z[g:(g + l) + 1]
    ##     c(sum(a), sum(!a))
    ## }))

    ## Z[81:116]
    

    ## A[f[100] + 1 + -3:3, 100]

    ## i <- 105
    ## A[f[length(f) - i]+-1:3, ncol(A) - i]
    
    ## ## yup looks good
    ## m <- cbind(
    ##     A[f[length(f) - 3]+-10:10, ncol(A) - 3],        
    ##     A[f[length(f) - 2]+-10:10, ncol(A) - 2],
    ##     A[f[length(f) - 1]+-10:10, ncol(A) - 1],
    ##     A[f[length(f)]+-10:10, ncol(A)]
    ## )

    ## m2 <- sapply((length(f) - 50):length(f), function(i) A[f[i]+-10:10, i + 1])

    ## for(i in 51:20) {
    ##     stopifnot(sum(m2[, i] == m2[, i - 1]) == 21)
    ## }
    
    ## apply(hapMatcherR[m2[, 1] + 1, seq(1, 461, 4)][, -10:0 + 116], 1, function(x) paste0(as.character(as.integer(x)), collapse = "-"))
    
    ## paste0(tail(Z, 10), collapse = "-")
    
    ## tail(Z, 10)
    ## ## so all perfect matches
    ## ## so why is it moving, does that make sense
    ## ## ALSO

    
    ## ## OK so given some pbwtL, try scanning up and down
    ## ## record indices that are needed
    ## ## extra, can I also get s, and determine when it ends, to keep / reject based on pbwtM?

    ## ## haha am at figure me out bit

    ## ## OK TRY AGAIN
    ## ## SOMETHING NOT QUITE RIGH
    ## ## PROMISING
    ## ## AM HERE
    
    ## ## breaks at one point
    ## A[cbind(out$trajectory[, 2] + 1, out$trajectory[, 1] + 1)]
    
    ## stuff <- sapply(-3:3, function(v_in) {
    ##     out <- find_index_backward(g_in, v_in, all_symbols, usge_all, egs = egs, K, return_trajectory = TRUE)
    ##     out[[2]][, 2]
    ## })
    
    ## t(stuff)


})
