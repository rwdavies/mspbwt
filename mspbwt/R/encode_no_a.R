#' @export
find_good_matches_without_a <- function(
                                        Z,
                                        all_symbols,
                                        usge_all,
                                        egs,
                                        pbwtL,
                                        pbwtM,
                                        hapMatcherR,
                                        do_checks = FALSE,
                                        A = NULL,
                                        which_snps_in_hapMatcherR = NULL,
                                        verbose = FALSE,
                                        list_of_columns_of_A = NULL,
                                        use_rcpp = FALSE
                                        ) {
  ##
  stopifnot(class(all_symbols[[1]][1, 1]) == "integer")
  ## check class, should be only ints
  stopifnot(!("numeric" %in% unlist(lapply(usge_all[[1]], function(x) sapply(x, class)))))
  f <- get_f_given_Z(Z, all_symbols, usge_all, egs)
  ## 
  if (is.null(which_snps_in_hapMatcherR)) {
    which_snps_in_hapMatcherR <- 1:ncol(hapMatcherR)
  }
  ## print(paste0("pbwtL = ", pbwtL))
  K <- sum(all_symbols[[1]][, 2])
  ## 
  nam <- c("v", "s", "k", "l", "trueA")
  mat_up <- matrix(-1, pbwtL, length(nam))
  colnames(mat_up) <- nam
  mat_down <- matrix(-1, pbwtL, length(nam))
  colnames(mat_down) <- nam
  a <- -1
  to_out <- as.list(1:length(f))
  g <- length(f) - 1
  mat_up_prev <- mat_up
  mat_down_prev <- mat_down
  if (do_checks) {
    list_of_mats <- as.list(1:length(f))
  }
  K <- as.integer(sum(all_symbols[[1]][, 2]))
  ##
  for(g in (length(f) - 1):0) {
    if (verbose) {
      print(paste0(g, ", ", date()))
    }
    fc <- f[g + 1]
    Zc <- Z[g + 1]
    C <- all_symbols[[g + 1]]
    usge <- usge_all[[g + 1]]
    ## intialize the ones we want to store here
    c_up <- 0
    c_down <- 0
    for(l in 0:(pbwtL - 1)) {
      v_up <- fc - l - 1
      if (v_up >= 0 && v_up <= (K - 1)) {
        out <- go_backwards_one_step(g = g + 1, v = v_up, C = C, usge = usge, egs = egs, K = K)
        ## s <- out[1]
        s <- NA
        k <- out[1]
        if (do_checks) {
          a <- A[v_up + 1, g + 1 + 1]
        }
        mat_up[l + 1, ] <- c(v_up, s, k, 0, a)
        cur_v <- mat_up[l + 1, "v"]
        prev_k <- mat_up_prev[c_up + 1, "k"]
        if (cur_v == prev_k) {
          mat_up[l + 1, "l"] <- mat_up_prev[c_up + 1, "l"] + 1
          c_up <- c_up + 1
        }
      } else {
        mat_up[l + 1, ] <- c(-1, -1, -1, -1, -1)
      }
      v_down <- fc + l
      if (v_down >= 0 && v_down <= (K - 1)) {
        out <- go_backwards_one_step(g = g + 1, v = v_down, C = C, usge = usge, egs = egs, K = K)
        ##s <- out[1]
        s <- NA                
        k <- out[1]
        if (do_checks) {
          a <- A[v_down + 1, g + 1 + 1]
        }
        mat_down[l + 1, ] <- c(v_down, s, k, 0, a)
        cur_v <- mat_down[l + 1, "v"]
        prev_k <- mat_down_prev[c_down + 1, "k"]
        if (cur_v == prev_k) {
          mat_down[l + 1, "l"] <- mat_down_prev[c_down + 1, "l"] + 1
          c_down <- c_down + 1
        }
      } else {
        mat_down[l + 1, ] <- c(-1, -1, -1, -1, -1)
      }
    }
    ##
    ## check if need save condition
    ##
    if ((g < (length(f) - 1) & (c_up < pbwtL | c_down < pbwtL)) | g == 0) {
      if (verbose) {
        print(paste0("to check: up ", pbwtL - c_up, ", down = ", pbwtL - c_down))
      }
      ## force test on everything
      if (g == 0) {
        c_up <- 0
        c_down <- 0
      }
      mat_out <- matrix(0, 2 * pbwtL - c_up - c_down, 3) ## somehow, g, index, l
      c_mat <- 0
      if (c_up < pbwtL) {
        for(ic in (c_up):(pbwtL - 1)) {
          if (mat_up_prev[ic + 1, "v"] >= 0) {
            v = mat_up_prev[ic + 1, "v"]
            k = mat_up_prev[ic + 1, "k"]
            len = mat_up_prev[ic + 1, "l"]
            if (verbose) {
              print(paste0("Find index backward: g = ", g, ", v_in = ", k))
            }
            if (use_rcpp) {
              ## print("yup using Rcpp")
              ##print(paste0("g = ", g, ", k = ", k))
              ##print("saving")
              ##save(g, k, all_symbols, usge_all, egs, K, list_of_columns_of_A, file = "/dev/shm/rwdavies/temp.RData", compress = FALSE)
              ## stop("WER")
              index <- Rcpp_find_index_backward(g_in = g, v_in = k, all_symbols = all_symbols, usge_all = usge_all, egs = egs, K = K, list_of_columns_of_A = list_of_columns_of_A, use_list_of_columns_of_A = !is.null(list_of_columns_of_A))
            } else {
              index <- find_index_backward(g_in = g, v_in = k, all_symbols = all_symbols, usge_all = usge_all, egs = egs, list_of_columns_of_A = list_of_columns_of_A)
            }
            ## so might disagree
            g2 <- g
            done <- FALSE
            while(!done & g2 >= 0) {
              if (as.integer(hapMatcherR[index + 1, which_snps_in_hapMatcherR[g2 + 1]]) == Z[g2 + 1]) {
                len <- len + 1
                g2 <- g2 - 1
              } else {
                done <- TRUE
              }
            }
            ## so might disagree
            if (len > pbwtM) {
              ## store this value
              if (do_checks) {
                stopifnot(index == mat_up_prev[ic + 1, "trueA"])
              }
              mat_out[c_mat + 1, ] <- c(g2 + 1, index, len)
              c_mat <- c_mat + 1
            }
          }
        }
      }
      ## down
      if (c_down < pbwtL) {
        for(ic in (c_down):(pbwtL - 1)) {
          if (mat_down_prev[ic + 1, "v"] >= (0)) {                    
            v = mat_down_prev[ic + 1, "v"]
            k = mat_down_prev[ic + 1, "k"]
            len = mat_down_prev[ic + 1, "l"]
            index <- find_index_backward(g_in = g, v_in = k, all_symbols = all_symbols, usge_all = usge_all, egs = egs)
            ## so might disagree
            g2 <- g
            done <- FALSE                        
            while(!done & g2 >= 0) {
              if (as.integer(hapMatcherR[index + 1, which_snps_in_hapMatcherR[g2 + 1]]) == Z[g2 + 1]) {
                len <- len + 1
                g2 <- g2 - 1
              } else {
                done <- TRUE                                
              }
            }
            ## so might disagree
            if (len > pbwtM) {
              ## store this value
              if (do_checks) {
                stopifnot(index == mat_down_prev[ic + 1, "trueA"])
              }
              mat_out[c_mat + 1, ] <- c(g2 + 1, index, len)
              c_mat <- c_mat + 1
            }
          }
        }
      }
      if (c_mat > 0) {
        to_out[[g + 1]] <- mat_out[1:c_mat, , drop = FALSE]
      }
      ## save these results
    }
    mat_up_prev <- mat_up
    mat_down_prev <- mat_down
    if (do_checks) {
      list_of_mats[[g + 1]] <- list(mat_up, mat_down)
    }
  }
  ##
  ## merge them together
  ## 
  w <- sapply(to_out, length) > 1
  if (sum(w) == 0) {
    n <- 0
  } else {
    n <- sum(sapply(to_out[w], nrow))
  }
  mat_out <- matrix(0, n, 3)
  colnames(mat_out) <- c("g", "index", "len")
  if (sum(w) == 0) {
    return(mat_out)
  }
  c <- 1
  for(m in to_out[w]) {
    mat_out[c + 0:(nrow(m) - 1), ] <- m
    c <- c + nrow(m)
  }
  if (do_checks) {
    return(
      list(
        mat_out = mat_out,
        list_of_mats = list_of_mats
      )
    )
  }
  ## check it out
  mat_out
}


find_start <- function(v, k, len, g, all_symbols, usge_all, egs) {
  return(g2)
}

get_f_given_Z <- function(Z, all_symbols, usge_all, egs, all_usg_check, use_U = FALSE) {
  ##
  ## f update
  ##
  nGrids <- length(all_symbols)
  f <- rep(0L, nGrids)
  Z_1 <- Z[1]
  if (Z_1 == 0) {
    Z_1 <- nrow(all_symbols[[1]])
  }
  C <- all_symbols[[1]]
  ##x <- c(0, cumsum(C[, 2]))[Z_1]
  x <- 0
  if (Z_1 > 1) {
    for(s in 1:(Z_1 - 1)) {
      x <- x + C[s, 2]
    }
  }
  f[1] <- x
  
  ##
  ## update
  ##
  for(g in 1:(nGrids - 1)) {
    C <- all_symbols[[g + 1]] ## effectively
    s <- Z[g + 1]
    k <- f[g] ## keep 0-based
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
  as.integer(f)
}



get_k_given_matrix_u <- function(s, v2, U) {
  k <- which.max(U[, s] == (v2 + 1)) - 1 - 1
  k
}
get_k_given_encoded_u <- function(s, v2, usge, C, K, egs, do_checks = FALSE, U = NULL, verbose = FALSE) {
  if (class(usge[[s]]) == "list") {
    ##
    out_mat <- usge[[s]][[1]]
    out_vec <- usge[[s]][[2]]
    ##
    ## first, need the 0-based first position lower than it in out_mat
    ## call it i_row
    ##
    ## if right at the end, call it
    done <- FALSE
    c <- 0 ## count
    up_i <- 0
    up_val <- 0
    down_i <- nrow(out_mat) - 1
    down_val <- C[s, 2]
    frac <- (v2 - up_val) / (down_val - up_val)
    i_row <- min(nrow(out_mat) - 2, round(frac * K /egs))
    if (verbose) {
      print(paste0("in R, initial frac = ", frac, ", i_row = ", i_row))
    }
    ## special case, at the end
    if (v2 >= out_mat[nrow(out_mat), 1]) {
      done <- TRUE
      i_row <- nrow(out_mat) - 1
    } 
    while(!done) {
      ##message(paste0("c = ", c, ", ", "i_row = ", i_row, ", ", "frac = ", frac, ", ",
      ##                "up_i = ", up_i, ", up_val = ", up_val, ", ",
      ##                "down_i = ", down_i, ", down_val = ", down_val
      ##                ))
      ## not entirely sure what to do about last entry, we'll see
      ## here let's make an estimate
      ## check below then above
      ## if neither we're good
      if ((i_row + 1) > nrow(out_mat)) {
        save(s, v2, usge, C, K, egs, file = "~/temp.RData")
        stop("something went wrong, see ~/temp.RData")
      }
      if (v2 < out_mat[i_row + 1, 1]) {
        ## OK, so we're under
        ## so re-set the bounds
        down_i <- i_row
        down_val <- out_mat[i_row + 1, 1]
        ## re-estimate
        i_row_proposed <- round((v2 - up_val) / (down_val - up_val) * (down_i - up_i)) + up_i
        i_row <- min(i_row - 1, i_row_proposed)
        ## 
      } else if (out_mat[i_row + 1 + 1, 1] <= v2) {
        up_i <- i_row
        up_val <- out_mat[i_row + 1 + 1, 1]
        ## re-estimate
        i_row_proposed <- round((v2 - up_val) / (down_val - up_val) * (down_i - up_i)) + up_i
        i_row <- max(i_row + 1, i_row_proposed)
        ## 
      } else {
        done <- TRUE
      }
      c <- c + 1
      if (c > 100000) {
        stop("problem")
      }
    }
    if (verbose) {
      print(paste0("in R, i_row = ", i_row))
    }
    ## message(paste0("c is ", c))
    ##
    ## special cases
    ## if all 0 we can skip
    ## if all 1 (increasing) we might need to consider
    ##
    if (i_row < (nrow(out_mat) - 1)) {
      if ((out_mat[i_row + 1 + 1, 1] - out_mat[i_row + 1, 1]) == egs) {
        ## so it's just whatever remains really
        k <- egs * (i_row) + v2 - out_mat[i_row + 1, 1]
        return(k)
      }
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
      if (verbose) {
        print(paste0("before: vec_pos = ", vec_pos, ", u = ", u, ", steps = ", steps))
      }
      if (is_plus) {
        u <- u + out_vec[vec_pos + 1]
      }
      steps <- steps + out_vec[vec_pos + 1]
      if (do_checks) {
        stopifnot(U[egs * i_row + steps + 1, s] == (val + u))
      }
      if (verbose) {
        print(paste0("consider:vec_pos = ", vec_pos, ", u = ", u, ", steps = ", steps))
      }
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
    if (verbose) {
      print(paste0("in R, steps = ", steps))
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
                                  egs,
                                  K,
                                  all_usg_check = NULL,
                                  do_checks = FALSE,
                                  use_U = FALSE,
                                  U = NULL,
                                  verbose = FALSE
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
    if (verbose) {
      print(paste0("s = ", s, ", v2 = ", v2))
    }
    k <- get_k_given_encoded_u(s, v2, usge, C, K, egs, do_checks = do_checks, U = U)
  }
  return(k)
  ## c(s, k) ## previou return, I don't think necessary now
}


find_index_backward <- function(
                                g_in,
                                v_in,
                                all_symbols,
                                usge_all = NULL,
                                egs = NULL,
                                all_usg_check = NULL,
                                do_checks = FALSE,
                                A = NULL,
                                use_U = FALSE,
                                return_trajectory = FALSE,
                                list_of_columns_of_A = NULL,
                                verbose = FALSE
                                ) {
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
    if (!is.null(list_of_columns_of_A)) {
      ## print(paste0("g = ", g))
      if (length(list_of_columns_of_A[[g + 1 + 1 + 1]]) > 1) {
        ## print(paste0("g = ", g, ", v = ", v))
        return(list_of_columns_of_A[[g + 1 + 1 + 1]][[v + 1]])
      }
    }
    usge <- usge_all[[g + 1 + 1]]
    C <- all_symbols[[g + 1 + 1]] ## symbols at g+1
    if (verbose) {
      print(paste0("g = ", g, ", v = ", v))
    }
    out <- go_backwards_one_step(
      g = g,
      v = v,
      C = C,
      usge = usge,
      all_usg_check = all_usg_check,
      do_checks = FALSE,
      use_U = use_U,
      egs = egs,
      U = U,
      verbose = FALSE
    )
    s <- NA
    ## s <- out[1]
    k <- out[1]
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


