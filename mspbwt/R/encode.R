encode_maximal_column_of_u <- function(u, egs, efficient = TRUE, verbose = FALSE) {
    n_rows <- ceiling(length(u) / egs)
    out_mat <- array(NA, c(n_rows, 4))
    colnames(out_mat) <- c("start1", "start0", "value", "vec_pos")
    out_vec <- integer(length(u))
    if (!efficient) {
        names(out_vec) <- 1:length(u)
    }
    vec_pos <- 0
    for(i in 0:(n_rows - 1)) {
        s0 <- egs * i
        e0 <- egs * (i + 1) - 1
        e0c <- e0 ## check value for no change
        if ((e0 + 1) > length(u)) {
            e0 <- length(u) - 1
            e0c <- e0 - 1
        }
        out_mat[i + 1, "start1"] <- s0 + 1 ## 1-based start
        out_mat[i + 1, "start0"] <- s0     ## 0-based start
        out_mat[i + 1, "value"] <- u[s0 + 1]
        do_encoding_this_SNP <- TRUE
        if (i < (n_rows - 1)) {
            if (
            (u[e0c + 1 + 1] - u[s0 + 1]) == 0 |
            (u[e0c + 1 + 1] - u[s0 + 1]) == (e0c + 1 - s0)
            ) {
                out_mat[i + 1, "vec_pos"] <- vec_pos - 1
                do_encoding_this_SNP <- FALSE
            }
        }
        if (verbose) {
            print(paste0("i = ", i, ", do_encoding_this_SNP = ", do_encoding_this_SNP))
        }
        if (do_encoding_this_SNP) {
            ## now build runs from this
            cs <- 0
            rt <- TRUE ## true = +1, FALSE = 0+
            for(j in s0:e0) {
                if (verbose) {
                    print(paste0("rt = ", rt, ", cs = ", cs, " j = ", j))
                }
                if (j == (length(u) - 1)) {
                    ## trigger a storage here
                    d <- 2
                } else {
                    d <- u[j + 1 + 1] - u[j + 1]
                }
                final <- j == e0
                if (rt) {
                    if (d == 1 & !final) {
                        cs <- cs + 1
                    } else {
                        if (verbose) {
                            print(paste0("save 1, j = ", j))
                        }
                        out_vec[vec_pos + 1] <- cs
                        names(out_vec)[vec_pos + 1] <- i                    
                        vec_pos <- vec_pos + 1
                        cs <- 1
                        rt <- FALSE
                    }
                } else {
                    if (d == 0 & !final) {
                        cs <- cs + 1                    
                    } else {
                        if (verbose) {
                            print(paste0("save 0, j = ", j))
                        }
                        out_vec[vec_pos + 1] <- cs
                        names(out_vec)[vec_pos + 1] <- i
                        vec_pos <- vec_pos + 1
                        cs <- 1
                        rt <- TRUE
                    }
                }
            }
            out_mat[i + 1, "vec_pos"] <- vec_pos - 1            
        }
    }
    out_vec <- out_vec[1:(vec_pos)]
    if (efficient) {
        names(out_vec) <- NULL
        out_mat <- out_mat[, c("value", "vec_pos")]
    }
    return(
        list(
            out_mat = out_mat,
            out_vec = out_vec
        )
    )
}



## 0-based v here
decode_maximal_value_of_u <- function(out_mat, out_vec, v, egs, do_checks = FALSE) {
    ## remainder is how many more to go
    ## e.g. 14 means you have to go 14 into the creation of this thing
    remainder <- (v) %% egs
    i_row <- (v - remainder) / egs ## how many multiples you have
    if (do_checks) {
        stopifnot((i_row * egs + remainder) == v)
    }
    if (remainder == 0) {
        return(out_mat[i_row + 1, "value"])
    } else {
        val <- out_mat[i_row + 1, "value"]
        if ((i_row + 1) < nrow(out_mat)) {
            ## these are constant - so don't need to check in between!
            next_val <- out_mat[i_row + 1 + 1, "value"]
            if (next_val == val) {
                return(val)
            } else if (((next_val) - val) == egs) {
                return(val + remainder)
            }
        }
        if (i_row == 0) {
            vec_pos <- 0
        } else {
            vec_pos <- out_mat[i_row, "vec_pos"] + 1 ## 0-based
        }
        if (do_checks) {
            stopifnot(which.max(names(out_vec) == i_row) == (vec_pos + 1))
        }
        u <- 0
        is_plus <- TRUE
        while(u < remainder) {
            next_u <- u + out_vec[vec_pos + 1]
            if (remainder <= next_u) {
                if (is_plus) {
                    val <- val + remainder - u
                }
                return(val)
            } else {
                u <- next_u
                if (is_plus) {
                    val <- val + out_vec[vec_pos + 1]
                    is_plus <- FALSE
                } else {
                    is_plus <- TRUE
                }
                vec_pos <- vec_pos + 1
            }
        }
    }
}

encode_minimal_column_of_u <- function(u) {
    which(diff(u) > 0)
}


encode_usg <- function(
    usg,
    symbol_count_at_grid,
    egs,                                
    n_min_symbols
) {
    lapply(1:nrow(symbol_count_at_grid), function(i_symbol) {
        if (symbol_count_at_grid[i_symbol, 2] > n_min_symbols) {
            return(encode_maximal_column_of_u(usg[, i_symbol], egs = egs))
        } else {
            return(encode_minimal_column_of_u(usg[, i_symbol]))
        }
    })
}



## x is
## v is 0-based on position
decode_minimal_value_of_u <- function(x, v) {
    i <- 0
    r <- x[i + 1]
    while(v >= r) {
        i <- i + 1
        if ((i + 1) > length(x)) {
            return(i)
        }
        r <- x[i + 1]
    }
    i
}


## s is the symbol
## v is the position down we are getting from
## e.g. if usg was the full matrix for a site, we would return
## usg[v + 1, s] for 0-based v and 1-based s
decode_value_of_usge <- function(
    usge,
    symbol_count_at_grid,
    s,
    v,
    egs,
    n_min_symbols
) {
    if (class(usge[[s]]) == "list") {
        decode_maximal_value_of_u(
            out_mat = usge[[s]][[1]],
            out_vec = usge[[s]][[2]],
            v = v,
            egs = egs
        )
    } else {
        decode_minimal_value_of_u(
            x = usge[[s]],
            v = v
        )
    }
}
