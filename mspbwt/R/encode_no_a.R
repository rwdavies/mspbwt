get_f_given_Z <- function(Z, all_symbols, usge_all, egs, all_usg_check, use_U = FALSE) {
    ##
    ## f update
    ##
    nGrids <- length(all_symbols)
    f <- rep(0L, nGrids)
    x <- c(0, cumsum(all_symbols[[1]][, 2]))
    Z_1 <- Z[1]
    if (Z_1 == 0) {
        Z_1 <- nrow(all_symbols[[1]])
    }
    f[1] <- x[Z_1]
    ##
    ## update
    ##
    for(g in 1:(nGrids - 1)) {
        C <- all_symbols[[g + 1]] ## effectively
        s <- Z[g + 1]
        k <- f[g] + 1
        if (use_U) {
            U <- all_usg_check[[g + 1]]
            u <- U[f[g] + 1, s]
        } else {
            u <- decode_value_of_usge(usge_all[[g + 1]], s = s, v = k, egs = egs)
        }
        c <- 0
        if (s > 1) {
            for(s2 in 2:s) {
                c <- c + C[s2 - 1, 2]
            }
        }
        u <- u + c
        f[g + 1] <- u
    }
    f
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
            ## so it's just whatever remains really
            k <- egs * (i_row) + v2 - out_mat[i_row + 1, 1]
            return(k)
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



## f1 <- function(v)  {
##     symbol_counts <- cumsum(C[, 2])
##     s <- which.max((v + 1) <= symbol_counts) ## 1-based (as usual)s
##     s
##     if (s == 1) {
##         v2 <- v
##     } else {
##         v2 <- v - symbol_counts[s - 1] ## 0-based instance
##     }
##     v2
## }

## f2 <- function(v) {
##     symbol_counts <- cumsum(C[, 2])    
##     s <- 1    
##     count <- C[1, 2]
##     while(v > (count - 1)) {
##         s <- s + 1
##         count <- count + C[s, 2]
##     }
##     s
##     if (s == 1) {
##         v2 <- v
##     } else {
##         v2 <- v - (count - C[s, 2])
##     }
##     v2
## }

## r <- sort(unique(c(0, rep(symbol_counts, each = 7) + -3:3)))
## r <- r[r >= 0 & r <= (sum(C[, 2] - 1))]

## for(v in r) {
##     expect_equivalent(f1(v), f2(v))
## }


## g is 0-based
## v is 0-based
go_backwards_one_step <- function(
    g,
    v,
    C,
    usge,
    all_usg_check = NULL,
    do_checks = FALSE,
    use_U = FALSE,
    egs = NULL,
    U = NULL
) {
    ## get s (1-based)
    s <- 1    
    count <- C[1, 2]
    while(v > (count - 1)) {
        s <- s + 1
        count <- count + C[s, 2]
    }
    s
    if (s == 1) {
        v2 <- v
    } else {
        v2 <- v - (count - C[s, 2])
    }
    v2
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
        K <- sum(C[, 2])
        k <- get_k_given_encoded_u(s, v2, usge, C, K, egs, do_checks = do_checks, U = U)
    }
    k
}


find_index_backward <- function(g_in, v_in, all_symbols, usge_all = NULL, egs = NULL, K = NULL, all_usg_check = NULL, do_checks = FALSE, A = NULL, use_U = FALSE, return_trajectory = FALSE) {
    g <- g_in
    v <- v_in
    if (do_checks) {
        g <- g_in
        message(paste0("Truth is (0-based) A[v, g + 1] = A[", v, ", ", g + 1, "] = ", A[v + 1, g + 1 + 1]))
        to_out <- matrix(0, g_in + 1, 4)
    }
    if (return_trajectory) {
        trajectory <- matrix(0, g_in + 1 + 1, 2)
        trajectory[1, 1] <- g + 1
        trajectory[1, 2] <- v ## 0-based
    }
    ## now this g needs to be one lower
    g <- g_in - 1
    for(g in (g_in - 1):(-1)) {
        usge <- usge_all[[g + 1 + 1]]
        C <- all_symbols[[g + 1 + 1]] ## symbols at g+1
        k <- go_backwards_one_step(
            g = g,
            v = v,
            C = C,
            usge = usge,
            all_usg_check = all_usg_check,
            do_checks = do_checks,
            use_U = use_U,
            egs = egs,
            U = U
        )
        if (do_checks) {
            stopifnot(A[k + 1, g + 1 + 1] == A[v + 1, g + 1 + 1 + 1])
            to_out[g_in - g, 1] <- g + 1
            to_out[g_in - g, 2] <- k
            to_out[g_in - g, 3] <- A[k + 1, g + 1 + 1]
        }
        if (return_trajectory) {
            trajectory[g_in - g + 1, 1] <- g + 1
            trajectory[g_in - g + 1, 2] <- k
        }
        ## update
        v <- k ## 0-based
        ## checks
    }
    index <- k ## first col is 0:(K - 1) so this is OK!
    if (do_checks) {
        message(paste0("Inferred is ", index))
    }
    if (return_trajectory) {
        return(list(index = index, trajectory = trajectory))
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


