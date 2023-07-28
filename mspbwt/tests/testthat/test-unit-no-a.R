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



})




test_that("no a slightly larger experiments", {

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
        for(v2 in 0:(C[s, 2] - 1)) {
            ## print(paste0("s = ", s, " v2 = ", v2))
            k1 <- get_k_given_matrix_u(s, v2, U) ## 0-based
            k2 <- get_k_given_encoded_u(s, v2, usge, C, K, egs, do_checks = TRUE, U = U)
            expect_equal(k1, k2)
        }
    }

    out_mat <-  find_good_matches_without_a(
        Z = Z,
        all_symbols,
        usge_all,
        egs,
        pbwtL,
        pbwtM,
        hapMatcherR,
        do_checks = TRUE,
        A = NULL
    ) 
    skip("rest not checked")
    
    ##
    ## now try to find f, see what it looks like
    ##
    usge_all <- ms_indices$usge_all


    f <- get_f_given_Z(Z, all_symbols, usge_all, egs)
    m <- cbind(0:(length(all_symbols) - 1), f)

    ## am here
    ## tomorrow, write algorithm to go backwards from this
    ## determine what I need to scan up and down
    ##
    K <- sum(all_symbols[[1]][, 2])

    
    g_in <- 49
    stuff <- sapply(806 + -1:1, function(v_in) {
        out <- find_index_backward(g_in, v_in, all_symbols, usge_all, egs = egs, K, return_trajectory = TRUE)
        out[[2]][, 2]
    })
    t(stuff)

    ## see if this is enough, as is. hopefully?
    ## implement as is, test as is
    
    ## 805 is the one (above?) whatevs

    rbind(X1C[1 + t(stuff)[, 51], ], NA, Z)

    ## looks OK, can I get some real, or semi-real, stuff in here
    rbind(X1C[20,], Z[])
    
    
    cbind(stuff[[4]]$trajectory, stuff[[3]]$trajectory[, 2], stuff[[2]]$trajectory[, 2], stuff[[1]]$trajectory[, 2])


    load("/data/smew1/rdavies/tempAtempAtempAaA.RData")
    all_symbols <- ms_indices[[1]]$all_symbols
    usge_all <- ms_indices[[1]]$usge_all
    egs <- ms_indices[[1]]$egs

    ## fake one
    Z <- hapMatcherR[100000, ]
    f <- get_f_given_Z(Z, all_symbols, usge_all, egs)
    

})


test_that("explore", {

    skip("WIP")


    load("/data/smew1/rdavies/tempAtempAtempAaA.RData")

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
    

    ## seems to work
    
    mat_out[mat_out[, 2] == 99999, ]
    
    ## check matches? seems good! not sure if this should necessarily be perfect always
    t(sapply(1:nrow(mat_out), function(i_row) {
        g <- mat_out[i_row, "g"]
        index <- mat_out[i_row, "index"]
        l <- mat_out[i_row, "len"]
        a <- as.integer(hapMatcherR[index + 1, seq(1, 461, 4)])[g:(g + l) + 1] == Z[g:(g + l) + 1]
        c(sum(a), sum(!a))
    }))

    Z[81:116]
    

    A[f[100] + 1 + -3:3, 100]

    i <- 105
    A[f[length(f) - i]+-1:3, ncol(A) - i]
    
    ## yup looks good
    m <- cbind(
        A[f[length(f) - 3]+-10:10, ncol(A) - 3],        
        A[f[length(f) - 2]+-10:10, ncol(A) - 2],
        A[f[length(f) - 1]+-10:10, ncol(A) - 1],
        A[f[length(f)]+-10:10, ncol(A)]
    )

    m2 <- sapply((length(f) - 50):length(f), function(i) A[f[i]+-10:10, i + 1])

    for(i in 51:20) {
        stopifnot(sum(m2[, i] == m2[, i - 1]) == 21)
    }
    
    apply(hapMatcherR[m2[, 1] + 1, seq(1, 461, 4)][, -10:0 + 116], 1, function(x) paste0(as.character(as.integer(x)), collapse = "-"))
    
    paste0(tail(Z, 10), collapse = "-")
    
    tail(Z, 10)
    ## so all perfect matches
    ## so why is it moving, does that make sense
    ## ALSO

    
    ## OK so given some pbwtL, try scanning up and down
    ## record indices that are needed
    ## extra, can I also get s, and determine when it ends, to keep / reject based on pbwtM?

    ## haha am at figure me out bit

    ## OK TRY AGAIN
    ## SOMETHING NOT QUITE RIGH
    ## PROMISING
    ## AM HERE
    
    ## breaks at one point
    A[cbind(out$trajectory[, 2] + 1, out$trajectory[, 1] + 1)]
    
    stuff <- sapply(-3:3, function(v_in) {
        out <- find_index_backward(g_in, v_in, all_symbols, usge_all, egs = egs, K, return_trajectory = TRUE)
        out[[2]][, 2]
    })
    
    t(stuff)


})
