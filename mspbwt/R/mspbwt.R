##
## where X1C means it is X, but 1-based, and complete
## i.e. every column of X contains entries from 1 to j for some arbitrary j >= 1
## and every entry from 1 to j occurs at least once
##
## NOTE - X1C can specially contain zeros, representing missing values
## these ARE recorded in all_symbols, but placed in the final entry in the per-grid matrix
## these are to be dynamically re-mapped to be the new largest symbol
##
#' @export
ms_BuildIndices_Algorithm5 <- function(
    X1C,
    all_symbols,
    verbose = FALSE,
    do_checks = FALSE,
    check_vs_indices = FALSE,
    indices = NULL,
    egs = 100,
    n_min_symbols = 100,
    with_Rcpp = FALSE,
    return_all_usg_check = FALSE
) {
    if (1 == 0) {
        verbose = FALSE
        do_checks = FALSE
        check_vs_indices = FALSE
        indices = NULL
        egs = 100
        n_min_symbols = 100
        with_Rcpp = FALSE
        return_all_usg_check = FALSE
    }
    if (do_checks | check_vs_indices) {
        return_d <- TRUE
    }
    if (do_checks) {
        ## here check validity of X1X
        ## needs to be 1:maxEntry per column
        ## each one used in decreasing order
        stopifnot(class(X1C[1, 1]) == "integer")
        for(iCol in 1:ncol(X1C)) {
            m1 <- min(X1C[, iCol])
            m2 <- max(X1C[, iCol])
            ## stopifnot(m1 == 1)
            stopifnot(sort(unique(X1C[, iCol])) == (m1:m2))
            ## also, check they are in decreasing order
            ##t <- table(X1C[, iCol])
            ##stopifnot(sort(t, decreasing = TRUE) == t)
        }
    }
    nGrids <- ncol(X1C)
    all_usg_check <- as.list(1:nGrids)
    ## thoughts - what to do if "too" rare - bin into too rare category? keep working with?
    ## inefficient?
    n_symbols_per_grid <- as.integer(sapply(all_symbols, nrow))
    Smax <- max(n_symbols_per_grid)
    ## build arrays including a, d and now u, v, c
    K <- nrow(X1C)
    T <- ncol(X1C)
    a <- array(as.integer(NA), c(K, T + 1)) ## orders
    c <- array(NA, T + 1) ## not sure why
    b <- array(NA, K)
    X1C_col1 <- X1C[, 1]
    X1C_col1[X1C_col1 == 0] <- nrow(all_symbols[[1]])
    a[, 1] <- as.integer(0:(K - 1)) ## by definition for some reason
    a[, 2] <- as.integer(order(X1C_col1) - 1) ## 0-based
    ## related to d
    d <- array(as.integer(NA), c(K + 1, T + 1)) ## distances
    d[, 1:2] <- 0L
    for(k in 1:(K - 1)) {
        if (X1C_col1[a[k, 2] + 1] != X1C_col1[a[k + 1, 2] + 1]) {
            d[k + 1, 2] <- 1L
        }
    }
    ##
    d[1, ] <- as.integer(1:(T + 1))
    d[K + 1, ] <- as.integer(1:(T + 1))
    ##
    ##
    ns_obs <- rep(0L, Smax)
    ##
    usge_all <- list(1:T) ## us, but all of them
    ## Smax_for_usl <- 0
    ## for(t in 1:T) {
    ##     x <- sum(all_symbols[[t]][, 2] > n_min_symbols)
    ##     if (x > Smax_for_usl) {
    ##         Smax_for_usl <- x
    ##     }
    ## }
    ## change it up, make it bigger, who cares, only done once
    Smax_for_usl <- max(sapply(all_symbols, nrow))
    usg <- array(0L, c(K + 1, Smax_for_usl))
    ##
    if (check_vs_indices) {
        t <- 1
        expect_equal(a[, t + 1], indices$a[, t + 1])
        expect_equal(d[, t + 1], indices$d[, t + 1])
    }
    ##
    for(t in 1L:as.integer(T)) {
        ##
        ## re-set
        ##
        St <- as.integer(n_symbols_per_grid[t]) ## number of symbols in this grid
        symbol_count <- as.integer(all_symbols[[t]][, "count"])
        ## if (sum(symbol_count) < K) {
        ##     St <- St + 1
        ##     symbol_count <- c(symbol_count, K - sum(symbol_count))
        ## }
        ##
        if (do_checks) {
            usg_check <- array(0L, c(K + 1, St))
        } else {
            usg_check <- array(0L, c(1, 1))
        }
        if (with_Rcpp) {
            f <- Rcpp_one_move_forward_buildindices
        } else {
            f <- one_move_forward_buildindices
        }
        things_to_check <- c(
            "X1C", "a", "usg",
            "usg_check", "d"
        )
        for(thing in things_to_check) {
            result <- eval(parse(text = paste0("class(", thing, "[1])")))
            if (result != "integer") {
                eval(parse(text = paste0("class(", thing, ")")))
                eval(parse(text = paste0("class(", thing, "[1])")))
                stop(paste0("class of ", thing, " has changed"))
            }
        }
        if (verbose) {
            print(paste0("t = ", t))
            print(paste0("St = ", St))
            print(symbol_count)
        }
        ##op <- options(digits.secs = 6); print("--in --"); print(Sys.time())
        out <- f(
            X1C = X1C,
            a = a,
            d = d,
            usg = usg,
            usg_check = usg_check,
            t = t,
            K = K,
            symbol_count = symbol_count,
            egs = egs,
            St = St,
            n_min_symbols = n_min_symbols,
            do_checks = do_checks
        )
        if (verbose) {
            print("done one move forward")
        }
        ## op <- options(digits.secs = 6); print("--out--"); print(Sys.time())
        if (!with_Rcpp) {
            a <- out$a
            d <- out$d
            usg <- out$usg
            usg_check <- out$usg_check
        }
        if (return_all_usg_check) {
            all_usg_check[[t]] <- usg_check
        }
        usge <- out[["usge"]]
        ##
        usge_all[[t]] <- usge
        if (do_checks) {
            expect_equal(
                usge,
                encode_usg(
                    usg = usg_check,
                    symbol_count_at_grid = all_symbols[[t]],
                    egs = egs,
                    n_min_symbols = n_min_symbols
                )
            )
        }
        ##
        if (t == 1) {
            d[0 + 1, t + 1] <- as.integer(t + 1) ## don't override
        }
        if (check_vs_indices) {
            expect_equal(a[, t + 1], indices$a[, t + 1])
            expect_equal(d[, t + 1], indices$d[, t + 1]) ## 160
        }
        ##
        ## do checks
        ##
        if (do_checks) {
            if (verbose) {
                print("do so checks")
            }
            for(k in 0:(K - 1)) {
                s <- X1C[a[k + 1, t] + 1, t] ## symbol here
                if (s == 0) {
                    s <- St
                }
                c <- c(0, cumsum(all_symbols[[t]][, 2]))[s]
                w <- decode_value_of_usge(
                    usge = usge_all[[t]],
                    s = s,
                    v = k + 1,
                    egs = egs
                ) + c
                ## if (verbose) {
                ##     message(paste0(
                ##         "k=", k, ", ",
                ##         "X=", X1C[a[k + 1, t] + 1, t], ", ",
                ##         "a[w,t+1]=", a[w, t + 1], ", ",
                ##         "a[k+1,t]=", a[k + 1, t]
                ##     ))
                ## }
                stopifnot(a[w, t + 1] == a[k + 1, t])
            }
            if (verbose) {
                print("done those  checks")
            }
        }
    }
    to_return <- list(
        a = a,
        d = d,
        usge_all = usge_all,
        egs = egs,
        n_min_symbols = n_min_symbols,
        all_symbols = all_symbols,
        all_usg_check = all_usg_check
    )
    to_return
}



## implements most of the move-forwardness of the algorithm building
one_move_forward_buildindices <- function(
    X1C,
    a,
    d,
    usg,
    usg_check,
    t,
    K,
    symbol_count,
    egs,
    St,
    n_min_symbols,
    do_checks,
    X1C_can_have_zeroes,
    verbose = FALSE
) {

    ## spave(
    ## X1C,
    ## a,
    ## d,
    ## usg,
    ## usg_check,
    ## t,
    ## K,
    ## symbol_count,
    ## egs,
    ## St,
    ## n_min_symbols,
    ## do_checks,
    ## X1C_can_have_zeroes,
    ## file = paste0("~/temp.", t, ".RData"))

    ## load(paste0("~/temp.", t, ".RData"))
    ## verbose <- TRUE
    ## n_min_symbols <- -1
    ##
    ##
    ## get count of number of each
    usge <- as.list(1:St) ## us for this (g)rid (e)ncoded
    first_usg_minimal_symbol <- 1 ## 1-based
    prev_value <- symbol_count[1]
    for(s in 1:St) {
        if (symbol_count[s] <= n_min_symbols) {
            usge[[s]] <- rep(-1, symbol_count[s])
        }
        ##     & symbol_count[s] <= prev_value) {
        ##     first_usg_minimal_symbol <- first_usg_minimal_symbol + 1
        ## } else {
        ##     usge[[s]] <- rep(-1, symbol_count[s])
        ## }
        ## ## see if this works! do for ALL of them
        ## usge[[s]] <- rep(-1, symbol_count[s])
        ## prev_value <- symbol_count[s]
    }
    start_count <- c(0, cumsum(symbol_count))
    ##
    nso <- rep(0L, St) ## n_symbols_observed
    pqs <- rep(as.integer(t), St) ## pqs - vector analogue to pq
    val <- c()
    usg[] <- 0L
    ##
    for(k in 0:(K - 1)) { ## haps (1-based)
        s <- X1C[a[k + 1, t] + 1, t] ## this symbol to consider
        if (s == 0) {
            s <- St
        }
        if (verbose) {
            message("-----")
            message(paste0("k = ", k, ", s = ", s))
        }
        ## match_start <- prev_d[k + 1]
        match_start <- d[k + 1, t]
        for(i in 0:(St - 1)) {
            if (match_start > pqs[i + 1]) {
                pqs[i + 1] <- match_start
            }
        }
        ## now - where it goes - 0 based
        stopifnot(s <= length(nso))
        val <- c(val, start_count[s] + nso[s] + 1)
        v <- start_count[s] + nso[s]
        a[start_count[s] + nso[s] + 1, t + 1] <- a[k + 1, t]
        d[start_count[s] + nso[s] + 1, t + 1] <- as.integer(pqs[s])
        if (verbose) {
            message(paste0("v = ", v, ", a[v + 1, t + 1] <- a[k + 1, t] = ", a[k + 1, t]))
        }
        ## d_vec[start_count[s] + nso[s] + 1] <- as.integer(pqs[s])
        usg[k + 1 + 1,] <- usg[k + 1, ]
        if (symbol_count[s] > n_min_symbols) {
            usg[k + 1 + 1, s] <- usg[k + 1 + 1, s] + 1L
        } else {
            usge[[s]][nso[s] + 1] <- k + 1
        }
        if (verbose) {
            message(paste0(usg[k + 1, ],collapse = "-"))
            message(paste0(usg[k + 1 + 1, ],collapse = "-"))
        }
        pqs[s] <- 0
        nso[s] <- nso[s] + 1
        if (do_checks) {
            usg_check[k + 1 + 1,] <- usg_check[k + 1, ]
            usg_check[k + 1 + 1, s] <- usg_check[k + 1 + 1, s] + 1
        }
    }

    ## encode the rest of them
    for(s in 1:St) {
        if (symbol_count[s] > n_min_symbols) {
            ## usge[[s]] <- Rcpp_encode_maximal_column_of_u(usg[, s], egs = egs)
            usge[[s]] <- encode_maximal_column_of_u(usg[, s], egs = egs)
        }
    }
    list(
        a = a,
        d = d,
        usg = usg,
        usge = usge,
        usg_check = usg_check
    )
}




wf <- function(k, t, s, usge_all, all_symbols, egs, indices, check_vs_indices = FALSE) {
    if (s == 0) {
        s <- nrow(all_symbols[[t]])
    }
    u <- decode_value_of_usge(
        usge = usge_all[[t]],
        s = s,
        v = k,
        egs = egs
    )
    ## u <- usge[k + 1, s] + c
    c <- 0
    if (s > 1) {
        for(s2 in 2:s) {
            c <- c + all_symbols[[t]][s2 - 1, 2]
        }
    }
    u <- u + c
    if (check_vs_indices) {
        if (s == 1) {
            return_val <- indices$u[k + 1, t]
        } else {
            return_val <- indices$v[k + 1, t] + indices$c[t]
        }
        if (u != return_val) {
            stop(paste0(
                "s = ", s,
                ", u = ", u,
                ", u(original) = ", return_val))
        }
    }
    u
}




#' @export
ms_MatchZ_Algorithm5 <- function(
    X,
    ms_indices,
    Z,
    verbose = FALSE,
    do_checks = FALSE,
    check_vs_indices = FALSE,
    indices = FALSE,
    print_or_message = print,
    pdfname = "~/Downloads/temp.pdf",
    make_plot = FALSE,
    do_up_and_down_scan = FALSE,
    mspbwtL = 3,
    mspbwtM = 3,
    XR = NULL,
    use_XR = FALSE,
    test_d = FALSE,
    height = NA,
    width = NA
) {
    ##
    if (is.na(height)) {
        height <- nrow(X) / 2 * 1.25 / 2
    }
    if (is.na(width)) {
        width <- 8
    }
    if (make_plot) pdf(pdfname, height = height, width = width)
    ##
    if (!use_XR) {
        K <- nrow(X)
        T <- ncol(X)
    } else {
        K <- nrow(XR)
        T <- ncol(XR)
    }
    ## indices
    a <- ms_indices[["a"]]
    if (test_d) {
        d <- ms_indices[["d"]]
        use_d <- TRUE ## shouldn't matter
    } else {
        if ("d" %in% names(ms_indices)) {
            d <- ms_indices[["d"]]
            use_d <- TRUE
        } else {
            d <- NULL
            use_d <- FALSE
        }
    }
    ##d_store <- ms_indices[["d_store"]]
    usge_all <- ms_indices[["usge_all"]]
    ## things needed as well
    egs <- ms_indices[["egs"]]
    n_min_symbols <- ms_indices[["n_min_symbols"]]
    all_symbols <- ms_indices[["all_symbols"]]
    ##
    if (length(Z) != (ncol(a) - 1)) {
        stop("Z not the right size")
    }
    e <- array(NA, T) ## keep this 0-based (like d)
    f <- array(NA, T) ## keep these 0-based
    g <- array(NA, T) ## keep these 0-based
    e[1] <- 0
    x <- c(0, cumsum(all_symbols[[1]][, 2]))
    Z_1 <- Z[1]
    if (Z_1 == 0) {
        Z_1 <- nrow(all_symbols[[1]])
    }
    f[1] <- x[Z_1]
    g[1] <- x[Z_1 + 1]
    if (verbose) {
        print(paste0("init, f[1] = ", f[1], ", g[1] = ", g[1]))
    }
    ## set this up
    if (do_up_and_down_scan) {
        ## uppy downy = ud
        ## a = above, b = below
        ud_up_prev <- integer(mspbwtL)
        ud_up_cur <- integer(mspbwtL)
        ud_up_length_prev <- integer(mspbwtL)
        ud_up_length_cur <- integer(mspbwtL)
        ##
        ud_down_prev <- integer(mspbwtL)
        ud_down_cur <- integer(mspbwtL)
        ud_down_length_prev <- integer(mspbwtL)
        ud_down_length_cur <- integer(mspbwtL)
        ##
        ud_up_prev[] <- -1
        ud_down_prev[] <- -1
        ## initialize up
        ## fg <- floor((f[1] + g[1] - 1) / 2) ## 0-based, include in "up"
        fg <- f[1] ## not ideal but consistent
        ## do up, include first entry fg
        i0 <- 0 ## 0-based
        while(i0 <= (mspbwtL - 1)) {
            if (0 <= (fg - i0)) {
                ud_up_prev[i0 + 1] <- a[fg - i0 + 1, 1 + 1] ## go up, so subtract
                ud_up_length_prev[i0 + 1] <- 0
            } else {
                ud_up_prev[i0 + 1] <- -1
                ud_up_length_prev[i0 + 1] <- -1
            }
            i0 <- i0 + 1
        }
        ## do down, go after first entry fg
        i0 <- 0
        while(i0 <= (mspbwtL - 1)) {
            if ((fg + i0 + 1) <= (K - 1)) {
                ud_down_prev[i0 + 1] <- a[fg + i0 + 1 + 1, 1 + 1] ## go up, so subtract
                ud_down_length_prev[i0 + 1] <- 0
            } else {
                ud_down_prev[i0 + 1] <- -1
                ud_down_length_prev[i0 + 1] <- -1
            }
            i0 <- i0 + 1
        }
        ##
        uppy_downy_reporter <- NULL
    }
    ##
    ## just do easy bit for now
    ##
    fc <- f[1]
    gc <- g[1]
    ec <- e[1]
    e1 <- NA
    top_matches <- NULL
    for(t in 2:T) {
        f1 <- wf(k = fc, t = t, s = Z[t], usge_all = usge_all, all_symbols = all_symbols, egs = egs, indices = indices, check_vs_indices = check_vs_indices)
        g1 <- wf(k = gc, t = t, s = Z[t], usge_all = usge_all, all_symbols = all_symbols, egs = egs, indices = indices, check_vs_indices = check_vs_indices)
        if (verbose) {
            print_or_message(paste0("Start of loop t=", t, ", fc = ", fc, ", gc = ", gc, ", ec = ", ec, ", Z[t] = ", Z[t],", f1=", f1, ", g1=", g1, ", e1 = ", e1))
        }
        if (do_up_and_down_scan) {
            ## fg <- floor((f1 + g1 - 1) / 2) ## 0-based, include in "up"
            fg <- f1
            ## if would point out of bounds, re-set
            if (fg == K) {
                fg <- K - 1
            }
            ##
            ## go "up" i.e. above i.e. up in the matrix
            ##
            i0_cur <- 0 ## 0-based, through local
            i0_prev <- 0
            ## go through previous values
            while((i0_prev <= (mspbwtL - 1)) && (-1 < ud_up_prev[i0_cur + 1]) && (0 <= (fg - i0_cur))) {
                ## focus on going through past list
                prev <- ud_up_prev[i0_prev + 1]
                ## now what is the current, does that work
                ## print("A")
                ## print(paste0("f1 = ", f1, ", g1 = ", g1))
                ## print(paste0("fg = ", fg, ", i0_cur = ", i0_cur))
                ## print(paste0("Z[t] = ", Z[t]))
                ## print(X[, t])
                ## print(fg - i0_cur)
                ## print(dim(a))
                cur <- a[fg - i0_cur + 1, t + 1] ## go up, so subtract
                if (cur == prev) {
                    ## it is a match. save and increment match
                    ud_up_cur[i0_cur + 1] <- cur
                    ## do the same for length. use previous length
                    ud_up_length_cur[i0_cur + 1] <- ud_up_length_prev[i0_prev + 1] + 1
                    ## increment up one
                    i0_cur <- i0_cur + 1
                } else {
                    ## it is not a match, report it
                    ## do not increment cur
                    len <- ud_up_length_prev[i0_prev + 1] ## 0-based
                    ## print(paste0("losing:", prev, ", with len = ", len))
                    if (mspbwtM <= len) {
                        uppy_downy_reporter <- rbind(
                            uppy_downy_reporter,
                            matrix(c(prev, t - 1, len), nrow = 1)
                        )
                    }
                }
                i0_prev <- i0_prev + 1
            }
            ## now fill in what was not set
            while((i0_cur <= (mspbwtL - 1))) {
                if (0 <= (fg - i0_cur)) {
                    cur <- a[fg - i0_cur + 1, t + 1] ## go up, so subtract
                    ud_up_cur[i0_cur + 1] <- cur
                    ## go backward, find start
                    e1 <- t
                    if (!use_XR) {
                        while((1 <= e1) && (X[cur + 1, e1] == Z[e1])) {
                            e1 <- e1 - 1
                        }
                    } else {
                        while((1 <= e1) && (XR[cur + 1, e1] == as.raw(Z[e1]))) {
                            e1 <- e1 - 1
                        }
                    }
                    ## go backwards, sort out start, how many before
                    ud_up_length_cur[i0_cur + 1] <- t - e1
                } else {
                    ud_up_cur[i0_cur + 1] <- -1
                    ud_up_length_cur[i0_cur + 1] <- -1
                }
                i0_cur <- i0_cur + 1
            }
            ## now reset
            ##if (verbose) {
            ## if (t <= 3) {
            ##     print(paste0("---------- t = ", t, " ------ reset up"))
            ##     print(paste0("fg = ", fg))
            ##     print(paste0("ud_up_prev = ", paste0(ud_up_prev, collapse = ", ")))
            ##     print(paste0("ud_up_length_prev = ", paste0(ud_up_length_prev, collapse = ", ")))
            ##     print(paste0("ud_up_cur = ", paste0(ud_up_cur, collapse = ", ")))
            ##     print(paste0("ud_up_length_cur = ", paste0(ud_up_length_cur, collapse = ", ")))
            ## }
            ##}
            ud_up_prev <- ud_up_cur
            ud_up_length_prev <- ud_up_length_cur
            ##
            ## go "down" i.e. below
            ##
            i0_cur <- 0 ## 0-based, through local
            i0_prev <- 0
            ## go through previous values
            while((i0_prev <= (mspbwtL - 1)) && (-1 < ud_down_prev[i0_cur + 1]) && ((fg + i0_cur + 1 + 1) <= K)) {
                ## focus on going through past list
                prev <- ud_down_prev[i0_prev + 1]
                ## now what is the current, does that work
                cur <- a[fg + i0_cur + 1 + 1, t + 1]
                if (cur == prev) {
                    ## it is a match. save and increment match
                    ud_down_cur[i0_cur + 1] <- cur
                    ## do the same for length. use previous length
                    ud_down_length_cur[i0_cur + 1] <- ud_down_length_prev[i0_prev + 1] + 1
                    ## increment up one
                    i0_cur <- i0_cur + 1
                } else {
                    ## it is not a match, report it
                    ## do not increment cur
                    len <- ud_down_length_prev[i0_prev + 1] ## 0-based
                    if (mspbwtM <= len) {
                        uppy_downy_reporter <- rbind(
                            uppy_downy_reporter,
                            matrix(c(prev, t - 1, len), nrow = 1)
                        )
                    }
                }
                i0_prev <- i0_prev + 1
            }
            ## now fill in what was not set
            while((i0_cur <= (mspbwtL - 1))) {
                if ((fg + i0_cur + 1) <= (K - 1)) {
                    cur <- a[fg + i0_cur + 1 + 1, t + 1] ## go up, so subtract
                    ud_down_cur[i0_cur + 1] <- cur
                    ## go backward, find start
                    e1 <- t
                    if (!use_XR) {
                        while((1 <= e1) && (X[cur + 1, e1] == Z[e1])) {
                            e1 <- e1 - 1
                        }
                    } else {
                        while((1 <= e1) && (XR[cur + 1, e1] == as.raw(Z[e1]))) {
                            e1 <- e1 - 1
                        }
                    }
                    ## go backwards, sort out start, how many before
                    ud_down_length_cur[i0_cur + 1] <- t - e1
                } else {
                    ud_down_cur[i0_cur + 1] <- -1
                    ud_down_length_cur[i0_cur + 1] <- -1
                }
                i0_cur <- i0_cur + 1
            }
            ##
            ud_down_prev <- ud_down_cur
            ud_down_length_prev <- ud_down_length_cur
        }
        if (g1 > f1) {
            ## nothing to do
            if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches, use_fc = FALSE)
        } else {
            if (verbose) print("save start")
            ## we have reached a maximum - need to report and update e, f, g
            for(k in fc:(gc - 1)) {
                top_matches <- rbind(
                    top_matches,
                    matrix(c(k, a[k + 1, t], ec + 1, t - 1), nrow = 1)
                )
            }
            ## for visualization efficiency
            out <- find_restart(e1 = e1, f1 = f1, g1 = g1, ec = ec, fc = fc, gc = gc, X = X, a = a, Z = Z, t = t, d = d, top_matches = top_matches, K = K, verbose = verbose, make_plot = make_plot, use_d = use_d, test_d = test_d)
            e1 <- out[["e1"]]
            f1 <- out[["f1"]]
            g1 <- out[["g1"]]
            ec <- e1
            print_or_message(paste0("After re-start f1=", f1, ", g1=", g1))
        }
        fc <- f1
        gc <- g1
        e[t] <- ec
        f[t] <- fc
        g[t] <- gc
    }
    t <- t + 1
    if (fc < gc) {
        if (verbose) {
            print_or_message(paste0("final fc = ", fc))
            print_or_message(paste0("final gc = ", gc))
        }
        for(k in fc:(gc - 1)) {
            if (verbose) {
                print_or_message("save final match")
                print_or_message(paste0("k = ", k, ", i0 = ", a[k + 1, t], ", s1 = ", ec + 1, ", e1 = ", t - 1))
            }
            top_matches <- rbind(
                top_matches,
                matrix(c(k, a[k + 1, t], ec + 1, t - 1), nrow = 1)
            )
        }
    }
    colnames(top_matches) <- c("k0", "indexB0", "start1", "end1")
    ##
    if (make_plot) dev.off()
    if (!do_up_and_down_scan) {
        return(top_matches)
    } else {
        ## report everything
        ## up
        i0_cur <- 0
        while((i0_cur <= (mspbwtL - 1))) {
            ## up
            prev <- ud_up_cur[i0_cur + 1]
            len <- ud_up_length_cur[i0_cur + 1]
            if (mspbwtM <= len) {
                uppy_downy_reporter <- rbind(
                    uppy_downy_reporter,
                    matrix(c(prev, t - 1, len), nrow = 1)
                )
            }
            i0_cur <- i0_cur + 1
        }
        ## down
        i0_cur <- 0
        while((i0_cur <= (mspbwtL - 1))) {
            ## up
            prev <- ud_down_cur[i0_cur + 1]
            len <- ud_down_length_cur[i0_cur + 1]
            if (mspbwtM <= len) {
                uppy_downy_reporter <- rbind(
                    uppy_downy_reporter,
                    matrix(c(prev, t - 1, len), nrow = 1)
                )
            }
            i0_cur <- i0_cur + 1
        }
        colnames(uppy_downy_reporter) <- c("index0", "end1", "len1")
        ##
        uppy_downy_reporter <- cbind(uppy_downy_reporter, start1 = uppy_downy_reporter[, "end1"] - uppy_downy_reporter[, "len1"] + 1)
        uppy_downy_reporter <- uppy_downy_reporter[, c("index0", "start1", "end1", "len1")]
        return(uppy_downy_reporter)
        ##     list(
        ##         uppy_downy_reporter = uppy_downy_reporter,
        ##         top_matches = top_matches
        ##     )
        ## )
    }
}













##
## more of a QUILT style function
## here take a matrix with arbitrary integer symbols
## and build 1-based symbol version
##
##  #' @export
## make_hapMatcherA <- function(
##     rhb_t
## ) {
##     K <- nrow(rhb_t)
##     nGrids <- ncol(rhb_t)
##     ## --- hapMatcherA
##     ## matrix K x nGrids
##     ## i is 1-based index of symbol
##     ## this includes all possible symbols
##     hapMatcherA <- array(0L, c(K, nGrids))
##     all_symbols <- list(1:nGrids)
##     for(iGrid in 1:nGrids) {
##         ## can safely ignore the end, it will be zeros what is not captured
##         a <- table(rhb_t[, iGrid], useNA = "always")
##         a <- a[order(-a)]
##         a <- a[a > 0]
##         a <- cbind(as.integer(names(a)), a)
##         rownames(a) <- NULL
##         colnames(a) <- c("symbol", "count")
##         hapMatcherA[, iGrid] <- as.integer(match(rhb_t[, iGrid], a[, "symbol"]))
##         ## here store both the count and the label
##         all_symbols[[iGrid]] <- a
##     }
##     return(
##         list(
##             hapMatcherA = hapMatcherA,
##             all_symbols = all_symbols
##         )
##     )
## }






## here if we have some X1C (original matrix, 1-based increasing possibly with 0s),
## or a regular X (original matrix, no 0s)
## and some Z with the same encoding
## defined at the SAME GRIDs (i.e. both are for grids and neither for SNPs)
## and we've mapped X to this new encoding, into hapMatcherA, and all_symbols
## we want to do the same with Z
#' @export
map_Z_to_all_symbols <- function(Z, all_symbols) {
    Z1 <- Z
    Z1[] <- 0L
    for(i in 1:length(Z)) {
        a <- all_symbols[[i]]
        Z1[i] <- match(Z[i], a[, "symbol"])
        ## if it is NA, something needs to be done
        if (is.na(Z1[i])) {
            ## if the first and last entries agree, then there is a missing character, which we can use
            ## this is represented by a 0 in X1C
            ## though note, can both be NA, so need to be careful
            the_same <- a[1, 1] == a[nrow(a), 1]
            if (is.na(the_same)) {
                if (is.na(a[1, 1]) & is.na(a[nrow(a), 1])) {
                    the_same <- TRUE
                } else {
                    the_same <- FALSE
                }
            }
            if (the_same) {
                Z1[i] <- 0
            } else {
                ## else, there is no missing, so we match to the closest
                dist <- calc_dist_between_rhb_t_and_hap(
                    matrix(a[, 1], ncol = 1),
                    STITCH::rcpp_int_expand(Z[i], 32),
                    32
                )
                Z1[i] <- which(dist == min(dist))[1]
            }
        }
    }
    Z1 <- as.integer(Z1)
    Z1
}


## where Zs is some SNP level vector of 0-1s (then rounded)
## ## then we want it in symbol form, matching when appropriate
## map_snp_hap_to_hapMatcher_grid_symbols <- function(Zs) {

##     Z1 <- Z
##     Z1[] <- 0L
##     for(i in 1:length(Z)) {
##         Z1[i] <- match(Z[i], all_symbols[[i]][, "symbol"])
##         if (is.na(Z1[i])) {
##             stop("have not figured this out yet!")
##         }
##     }
##     Z1 <- as.integer(Z1)
##     Z1


## }


## hapc is the "haplotype" (value) in binary form
## a is available symbols in that grid in matrix form
## so we need to return 1:(nrow(a)) because those are the available entries
map_one_binary_value_to_hapMatcher <- function(
    hapc,
    distinctHapsB,
    iGrid1,
    nSNPs
) {
    dist <- calc_dist_between_rhb_t_and_hap(
        distinctHapsB[, iGrid1, drop = FALSE],
        STITCH::rcpp_int_expand(hapc, nSNPs),
        nSNPs
    )
    which(dist == min(dist))[1]
}


find_restart <- function(
    e1,
    f1,
    g1,
    ec,
    fc,
    gc,
    X,
    a,
    Z,
    t,
    d,
    top_matches,
    K,
    verbose,
    make_plot,
    use_d,
    test_d
) {
    ## save(e1,
    ## f1,
    ## g1,
    ## ec,
    ## fc,
    ## gc,
    ## X,
    ## a,
    ## Z,
    ## t,
    ## d,
    ## top_matches,
    ## K,
    ## verbose,
    ## make_plot,
    ## use_d,
    ## test_d,
    ## file = "~/temp.RData")
    ## ## stop("WER")
    ## load("~/temp.RData")
    ## load("~/Downloads/tempAAA/temp.RData")
    ## test_d <- TRUE
    fc <- f1
    gc <- g1
    if (use_d) {
        e1 <- d[f1 + 1, t + 1] - 1 ## this is 0-based, probably!
    } else {
        e1 <- t - 1
    }
    if (verbose) print(paste0("e1 = ", e1))
    if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches, use_fc = FALSE)
    ##
    ## do complicated matching here, need to potentially come back more
    ##
    ##save(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches, file = "~/temp.RData")
    ##load("~/temp.RData")
    matches_lower <- FALSE
    matches_upper <- FALSE
    if ((e1 == t) && (f1 == K)) {
        e1 <- t - 1
    }
    while((!matches_lower) && (!matches_upper)) {
        ## here I check if I can do ABOVE i.e. one value up (reduced)
        if (f1 > 0) {
            matches_upper <- Z[e1 + 1] == X[a[f1 - 1 + 1, t + 1] + 1, e1 + 1]
        } else {
            matches_upper <- FALSE
        }
        ## here I check if I can go BELOW i.e. one value down
        if (f1 < K) {
            matches_lower <- Z[e1 + 1] == X[a[f1 + 1, t + 1] + 1, e1 + 1]
        } else {
            matches_lower <- FALSE
        }
        if (!matches_lower & !matches_upper) {
            e1 <- e1 + 1
        }
        if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
    }
    ##
    ## this CAN happen, if there is a symbol mis-match, and have to go forward
    ##
    if (matches_upper) {
        f1 <- f1 - 1
        f1_init <- f1
        index <- a[f1 + 1, t + 1] ## a
        if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
        ## should skip if both
        while (Z[e1 - 1 + 1] == X[index + 1, e1 - 1 + 1]) {
            e1 <- e1 - 1
            if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
        }
        ## so should match from (1-based) e1 + 1
        if (use_d | test_d) {
            f1 <- f1_init
            while (d[f1 + 1, t + 1] <= e1) {
                f1 <- f1 - 1
                if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
            }
            f1_with_d <- f1
        }
        if (!use_d | test_d) {
            f1 <- f1_init
            ##} else {
            ## so here, we need to check everything between e1 and t
            ## only do if we can access the next index
            cond <- TRUE
            while(0 <= (f1 - 1) && cond) {
                ## only
                next_index <- a[f1 - 1 + 1, t + 1]
                cond <- TRUE
                for(tt0 in (e1):(t - 1)) {
                    if (X[index + 1, tt0 + 1] != X[next_index + 1, tt0 + 1]) {
                        cond <- FALSE
                    }
                }
                if (cond) {
                    f1 <- f1 - 1
                    if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
                }
            }
            if (test_d) {
                if (f1_with_d != f1) {
                    stop(paste0("Problem:f1_with_d = ", f1_with_d, " and f1 = ", f1))
                }
            }
        }
    }
    if (matches_lower) {
        g1 <- g1 + 1
        g1_init <- g1
        index <- a[f1 + 1, t + 1] ## a
        if (verbose) {
            print(paste0(
                "e1 = ", e1, ", ",
                "Z[e1 - 1 + 1] = ", Z[e1 - 1 + 1], ", X[index + 1, e1 - 1 + 1] = ", X[index + 1, e1 - 1 + 1]
            ))
        }
        if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
        while (Z[e1 - 1 + 1] == X[index + 1, e1 - 1 + 1]) {
            e1 <- e1 - 1
            if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
        }
        ## so should match from (1-based) e1 + 1
        if (use_d | test_d) {
            g1 <- g1_init
            while ((g1 < K) && (d[g1 + 1, t + 1] <= e1)) { ## d
                g1 <- g1 + 1
                if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
            }
            g1_with_d <- g1
        }
        if (!use_d | test_d) {
            g1 <- g1_init
            ## so here, we need to check everything between e1 and t
            cond <- TRUE
            while((g1 < K) && cond) {
                ## only
                next_index <- a[f1 + 1 + 1, t + 1]
                cond <- TRUE
                for(tt0 in (e1):(t - 1)) {
                    if (X[index + 1, tt0 + 1] != X[next_index + 1, tt0 + 1]) {
                        cond <- FALSE
                    }
                }
                if (cond) {
                    g1 <- g1 + 1
                    if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
                }
            }
            if (test_d) {
                if (g1_with_d != g1) {
                    stop(paste0("Problem:g1_with_d = ", g1_with_d, " and g1 = ", g1))
                }
            }
        }
    }
    if (verbose) print(paste0("e1 = ", e1))
    if (verbose) print("save stop")
    return(
        list(
            e1 = e1,
            f1 = f1,
            g1 = g1
        )
    )
}






