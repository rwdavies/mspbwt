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



get_k_given_matrix_u <- function(s, v2, U) {
    k <- which.max(U[, s] == (v2 + 1)) - 1 - 1
    k
}
get_k_given_encoded_u <- function(s, v2, usge, C, K, egs, do_checks = FALSE, U = NULL) {
    if (class(usge[[s]]) == "list") {
        ##
        out_mat <- usge[[s]][[1]]
        out_vec <- usge[[s]][[2]]
        ##
        ## first, need the 0-based first position lower than it in out_mat
        ## call it i_row
        ##
        done <- FALSE
        c <- 0 ## count
        i_row <- floor(v2 / C[s, 2] * K / egs)
        while(!done) {
            ## not entirely sure what to do about last entry, we'll see
            ## here let's make an estimate
            ## check below then above
            ## if neither we're good
            if (v2 < out_mat[i_row + 1, 1]) {
                ## eventually, re-estimate, something more efficient
                ## here, re-estimate based on new value?
                ## inc <- floor(v2 / out_mat[i + 1, 1] * K / egs)
                i_row <- i_row - 1
            } else if (out_mat[i_row + 1 + 1, 1] <= v2) {
                i_row <- i_row + 1
            } else {
                done <- TRUE
            }
            c <- c + 1
            if (c > 10000) {
                stop("problem")
            }
        }
        ##
        ## special cases
        ## if all 0 we can skip
        ## if all 1 (increasing) we might need to consider
        ##
        if ((out_mat[i_row + 1 + 1, 1] - out_mat[i_row + 1, 1]) == egs) {
            stop("figure me out")
        }
        ##
        ##
        ## now second part of decoding
        ## slightly trickier
        ## need to figure out how much further to go
        ##
        ## figure out where to start in vector
        if (i_row == 0) {
            vec_pos <- 0
        } else {
            vec_pos <- out_mat[i_row, "vec_pos"] + 1 ## still 0-based
        }
        ## now second part of decoding (different from original decoding idea)
        ##
        ##
        val <- out_mat[i_row + 1, 1] ## at this point
        remainder <- v2 - val
        ## keep adding until it's reached then go back one
        u <- 0 ## is how much we've added
        steps <- 0 ## count of how many steps forward we've gone in original vector
        is_plus <- TRUE
        done <- FALSE
        while(!done) {
            ## print(paste0("before: vec_pos = ", vec_pos, ", u = ", u, ", steps = ", steps))
            if (is_plus) {
                u <- u + out_vec[vec_pos + 1]
            }
            steps <- steps + out_vec[vec_pos + 1]
            if (do_checks) {
                stopifnot(U[egs * i_row + steps + 1, s] == (val + u))
            }
            ## print(paste0("consider:vec_pos = ", vec_pos, ", u = ", u, ", steps = ", steps))
            if (u > remainder) {
                ## something like this
                k <- egs * i_row + steps - (u - remainder)
                if (do_checks) {
                    stopifnot(U[k + 1, s] == v2)
                    stopifnot(U[k + 1 + 1, s] == (v2 + 1))
                }
                done <- TRUE
            }
            vec_pos <- vec_pos + 1
            is_plus <- !is_plus
            if (!done && (vec_pos) > (out_mat[i_row + 1, 2])) {
                ## must need all steps
                k <- egs * (i_row + 1) - 1
                done <- TRUE
            }
        }
        return(as.integer(k))
    } else {
        return(as.integer(usge[[s]][v2 + 1] - 1))
    }
}



find_index_backward <- function(g_in, v_in, all_symbols, usge = NULL, egs = NULL, K = NULL, all_usg_check = NULL, do_checks = FALSE, A = NULL, use_U = TRUE) {
    g <- g_in
    v <- v_in
    if (do_checks) {
        g <- g_in
        message(paste0("Truth is (0-based) A[v, g + 1] = A[", v, ", ", g + 1, "] = ", A[v + 1, g + 1 + 1]))
        to_out <- matrix(0, g_in + 1, 4)
    }
    ## now this g needs to be one lower
    g <- g_in - 1
    for(g in (g_in - 1):(-1)) {
        C <- all_symbols[[g + 1 + 1]] ## symbols at g+1
        ## find k
        ## first find column (gives s)
        symbol_counts <- cumsum(C[, 2])
        s <- which.max((v + 1) <= symbol_counts) ## 1-based (as usual)s
        ## v2 is 0-based instance of first time this symbol occurs in that column
        if (s == 1) {
            v2 <- v ## 0-based
        } else {
            v2 <- v - symbol_counts[s - 1] ## 0-based instance
        }
        ## find first instance of this
        if (use_U) {
            U <- all_usg_check[[g + 1 + 1]] ## "U" i.e. FM index at g+1
            k <- get_k_given_matrix_u(s, v2, U) ## 0-based
            ## k <- which.max(U[, s] == (v2 + 1)) - 1 - 1
        } else {
            if (do_checks) {
                U <- all_usg_check[[g + 1 + 1]] ## "U" i.e. FM index at g+1
            } else {
                U <- NULL
            }
            k <- get_k_given_encoded_u(s, v2, usge[[g + 1 + 1]], C, K, egs, do_checks = do_checks, U = U)
        }
        if (do_checks) {
            ## print(A[k + 1, g + 1 + 1])
            ## print(A[v + 1, g + 1 + 1 + 1])
            stopifnot(A[k + 1, g + 1 + 1] == A[v + 1, g + 1 + 1 + 1])
            to_out[g_in - g, 1] <- g + 1
            to_out[g_in - g, 2] <- k
            to_out[g_in - g, 3] <- A[k + 1, g + 1 + 1]
        }
        ## update
        v <- k ## 0-based
        ## checks
    }
    index <- k ## first col is 0:(K - 1) so this is OK!
    if (do_checks) {
        message(paste0("Inferred is ", index))
    }
    index
}


find_index_forward <- function(g_in, k_in, all_usg_check, all_symbols, A_last_col, do_checks = FALSE, A = NULL) {
    g <- g_in
    k <- k_in
    nGrids <- length(all_symbols)
    if (do_checks) {
        message(paste0("Truth is (0-based) A[k, g + 1] = A[", k, ", ", g + 1, "] = ", A[k + 1, g + 1 + 1]))
        to_out <- matrix(0, nGrids - g_in + 1, 4)
    }
    ## now this g needs to be one lower
    for(g in (g_in:(nGrids - 2))) {
        U <- all_usg_check[[g + 1 + 1]] ## "U" i.e. FM index at g+1
        C <- all_symbols[[g + 1 + 1]] ## symbols at g+1
        symbol_counts <- c(0, cumsum(C[, 2]))
        ##
        s <- which.max(U[k + 1 + 1, ] - U[k + 1, ]) ## symbol
        ##
        v <- symbol_counts[s] + U[k + 1, s]
        ##
        if (do_checks) {
            stopifnot(A[k + 1, g + 1 + 1] == A[v + 1, g + 1 + 1 + 1])
            to_out[g - g_in + 1, 1] <- g + 1
            to_out[g - g_in + 1, 2] <- k
            to_out[g - g_in + 1, 3] <- A[k + 1, g + 1 + 1]
        }
        ## update
        k <- v
        ## checks
    }
    index <- A_last_col[k + 1]
    if (do_checks) {
        message(paste0("Inferred is ", index))
    }
    index
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
    usge <- ms_indices$usge
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
            usge = usge,
            egs = egs,
            K = K,
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

    ## exhaustive check!
    for(s in 1:nrow(C)) {
        for(v2 in 0:(C[s, 2] - 1)) {
            ## print(paste0("s = ", s, " v2 = ", v2))
            k1 <- get_k_given_matrix_u(s, v2, U) ## 0-based
            k2 <- get_k_given_encoded_u(s, v2, usge, C, K, egs, do_checks = TRUE, U = U)
            expect_equal(k1, k2)
        }
    }

    ## OK this works!

})
