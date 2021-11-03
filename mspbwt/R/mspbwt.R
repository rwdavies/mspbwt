##
## where X1C means it is X, but 1-based, and complete
## i.e. every column of X contains entries from 1 to j for some arbitrary j >= 1
## and every entry from 1 to j occurs at least once
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
    n_min_symbols = 100
) {
    ## thoughts - what to do if "too" rare - bin into too rare category? keep working with?
    ## inefficient?
    n_symbols_per_grid <- as.integer(sapply(all_symbols, nrow))
    Smax <- max(n_symbols_per_grid)
    ## build arrays including a, d and now u, v, c
    K <- nrow(X1C)
    T <- ncol(X1C)
    a <- array(NA, c(K, T + 1)) ## orders
    d <- array(NA, c(K + 1, T + 1)) ## distances
    dtemp <- array(NA, K) ## temp for filler
    c <- array(NA, T + 1) ## not sure why
    b <- array(NA, K)
    a[, 1] <- 0:(K - 1) ## by definition for some reason
    a[, 2] <- order(X1C[, 1]) - 1 ## 0-based
    ## d is a 0, except first entry, and on (ordered) switch
    d[, 1] <- 0
    d[, 2] <- 0
    for(k in 1:(K - 1)) {
        if (X1C[a[k, 2] + 1, 1] != X1C[a[k + 1, 2] + 1, 1]) {
            d[k + 1, 2] <- 1
        }
    }
    ## 
    d[1, ] <- 1:(T + 1)
    d[K + 1, ] <- 1:(T + 1)
    ##
    ns_obs <- rep(0L, Smax)
    ##
    usge_all <- list(1:T) ## us, but all of them    
    Smax_for_usl <- 0
    for(t in 1:T) {
        x <- sum(all_symbols[[t]][, 2] > n_min_symbols)
        if (x > Smax_for_usl) {
            Smax_for_usl <- x
        }
    }
    usg <- array(0L, c(K + 1, Smax_for_usl))
    ## 
    if (check_vs_indices) {
        t <- 1
        expect_equal(a[, t + 1], indices$a[, t + 1])
        expect_equal(d[, t + 1], indices$d[, t + 1])
    }
    for(t in 1:T) {
        ##
        ## OK, whatevs
        St <- n_symbols_per_grid[t] ## number of symbols in this grid
        ## 
        ## get count of number of each
        usge <- list(1:St) ## us for this (g)rid (e)ncoded
        symbol_count <- all_symbols[[t]][, "count"]
        first_usg_minimal_symbol <- 1 ## 1-based
        for(s in 1:St) {
            if (symbol_count[s] > n_min_symbols) {
                first_usg_minimal_symbol <- first_usg_minimal_symbol + 1
            } else {
                usge[[s]] <- rep(-1, symbol_count[s])
            }
        }
        start_count <- c(0, cumsum(symbol_count))
        ##
        ## 
        ##
        nso <- rep(0L, St) ## n_symbols_observed
        pqs <- rep(t, St) ## pqs - vector analogue to pq
        val <- c()
        usg[] <- 0
        if (do_checks) {
            usg_check <- array(0L, c(K + 1, St))
        }
        ## fill in the rest of this one!
        for(k in 0:(K - 1)) { ## haps (1-based)
            s <- X1C[a[k + 1, t] + 1, t] ## this symbol to consider
            match_start <- d[k + 1, t]
            for(i in 0:(St - 1)) {
                if (match_start > pqs[i + 1]) {
                    pqs[i + 1] <- match_start
                }
            }
            ## now - where it goes - 0 based
            val <- c(val, start_count[s] + nso[s] + 1)
            a[start_count[s] + nso[s] + 1, t + 1] <- a[k + 1, t]
            d[start_count[s] + nso[s] + 1, t + 1] <- pqs[s]
            usg[k + 1 + 1,] <- usg[k + 1, ]
            if (s < first_usg_minimal_symbol) {
                usg[k + 1 + 1, s] <- usg[k + 1 + 1, s] + 1
            } else {
                usge[[s]][nso[s] + 1] <- k + 1
            }
            pqs[s] <- 0
            nso[s] <- nso[s] + 1
            if (do_checks) {
                usg_check[k + 1 + 1,] <- usg_check[k + 1, ]
                usg_check[k + 1 + 1, s] <- usg_check[k + 1 + 1, s] + 1
            }
        }
        if ((first_usg_minimal_symbol - 1 - 1) >= 0) {
            for(s in 0:(first_usg_minimal_symbol - 1 - 1)) {
                usge[[s + 1]] <- encode_maximal_column_of_u(usg[, s + 1], egs = egs)
            }
        }
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
        if (t == 1) {
            ## Not sure
            d[0 + 1, t + 1] <- t + 1 ## don't override
        }
        if (check_vs_indices) {
            expect_equal(a[, t + 1], indices$a[, t + 1])
            expect_equal(d[, t + 1], indices$d[, t + 1]) ## 160
        }
        ## 
        ## do checks
        ##
        for(k in 0:(K - 1)) {
            s <- X1C[a[k + 1, t] + 1, t] ## symbol here
            c <- c(0, cumsum(all_symbols[[t]][, 2]))[s]
            w <- decode_value_of_usge(
                usge = usge_all[[t]],
                symbol_count_at_grid = all_symbols[[t]],
                s = s,
                v = k + 1,
                egs = egs,
                n_min_symbols = n_min_symbols
            ) + c
            if (verbose) {
                message(paste0(
                    "k=", k, ", ",
                    "X=", X[a[k + 1, t] + 1, t], ", ",
                    "a[w,t+1]=", a[w, t + 1], ", ",
                    "a[k+1,t]=", a[k + 1, t]
                ))
            }
            if (do_checks) {
                stopifnot(a[w, t + 1] == a[k + 1, t])
            }
        }
    }
    return(
        list(
            a = a,
            d = d,
            usge_all = usge_all,
            egs = egs,
            n_min_symbols = n_min_symbols,
            all_symbols = all_symbols
        )
    )
}




#' @export
ms_MatchZ_Algorithm5 <- function(
    X,
    ms_indices,
    Z,
    verbose = FALSE,
    do_checks = FALSE,
    check_vs_indices = FALSE,
    indices = FALSE
) {
    K <- nrow(X)
    T <- ncol(X)
    a <- ms_indices[["a"]]
    d <- ms_indices[["d"]]
    usge_all <- ms_indices[["usge_all"]]
    egs <- ms_indices[["egs"]]
    n_min_symbols <- ms_indices[["n_min_symbols"]]
    all_symbols <- ms_indices[["all_symbols"]]
    if (length(Z) != (ncol(a) - 1)) {
        stop("Z not the right size")
    }
    e <- array(NA, T) ## keep this 0-based (like d)
    f <- array(NA, T) ## keep these 0-based
    g <- array(NA, T) ## keep these 0-based
    e[1] <- 0
    x <- c(0, cumsum(all_symbols[[1]][, 2]))
    f[1] <- x[Z[1]]
    g[1] <- x[Z[1] + 1]
    ##
    ## just do easy bit for now
    ##
    wf <- function(k, t, s, usge_all, all_symbols) {
        c <- c(0, cumsum(all_symbols[[t]][, 2]))[s]
        ## u <- usge[k + 1, s] + c
        u <- decode_value_of_usge(
            usge = usge_all[[t]],
            symbol_count_at_grid = all_symbols[[t]],
            s = s,
            v = k,
            egs = egs,
            n_min_symbols = n_min_symbols
        ) + c
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
    fc <- f[1]
    gc <- g[1]
    ec <- e[1]
    e1 <- NA
    top_matches <- NULL
    for(t in 2:T) {
        f1 <- wf(fc, t, Z[t], usge_all, all_symbols)
        g1 <- wf(gc, t, Z[t], usge_all, all_symbols)
        if (verbose) {
            message(paste0("Start of loop t=", t, ", fc = ", fc, ", gc = ", gc, ", ec = ", ec, ", Z[t] = ", Z[t],", f1=", f1, ", g1=", g1, ", e1 = ", e1))
        } 
        if (g1 > f1) {
            ## nothing to do
        } else {
            ## we have reached a maximum - need to report and update e, f, g
            for(k in fc:(gc - 1)) {
                top_matches <- rbind(
                    top_matches,
                    matrix(c(k, a[k + 1, t], ec + 1, t - 1), nrow = 1)
                )
            }
            e1 <- d[f1 + 1, t + 1] - 1 ## this is 0-based, probably!
            if ((Z[e1 + 1] == 0 && f1 > 0) || f1 == K) {
                f1 <- g1 - 1
                index <- a[f1 + 1, t + 1]
                while (Z[e1 - 1 + 1] == X[index + 1, e1 - 1 + 1]) {
                    e1 <- e1 - 1
                }
                while (d[f1 + 1, t + 1] <= e1) {
                    f1 <- f1 - 1
                }
            } else if (f1 < K) {
                g1 <- f1 + 1
                index <- a[f1 + 1, t + 1]
                while (Z[e1 - 1 + 1] == X[index + 1, e1 - 1 + 1]) {
                    e1 <- e1 - 1
                }
                while ((g1 < K) && (d[g1 + 1, t + 1] <= e1)) {
                    g1 <- g1 + 1
                }
            }
            ec <- e1
        }
        fc <- f1
        gc <- g1
        e[t] <- ec
        f[t] <- fc
        g[t] <- gc
    }
    t <- t + 1
    if (fc < gc) {
        for(k in fc:(gc - 1)) {
            if (verbose) {
                message("save final match")
            }
            top_matches <- rbind(
                top_matches,
                matrix(c(k, a[k + 1, t], ec + 1, t - 1), nrow = 1)
            )
        }
    }
    colnames(top_matches) <- c("k0", "indexB0", "start1", "end1")
    ##
    ## perform checks if wanted
    ##
    return(top_matches)
}













##
## more of a QUILT style function
## here take a matrix with arbitrary integer symbols
## and build 1-based symbol version
##
#' @export
make_hapMatcherA <- function(
    rhb_t
) {
    K <- nrow(rhb_t)
    nGrids <- ncol(rhb_t)
    ## --- hapMatcherA
    ## matrix K x nGrids
    ## i is 1-based index of symbol
    ## this includes all possible symbols
    hapMatcherA <- array(0L, c(K, nGrids))
    all_symbols <- list(1:nGrids)
    for(iGrid in 1:nGrids) {
        ## can safely ignore the end, it will be zeros what is not captured
        a <- table(rhb_t[, iGrid], useNA = "always")
        a <- a[order(-a)]
        a <- a[a > 0]
        a <- cbind(as.integer(names(a)), a)
        rownames(a) <- NULL
        colnames(a) <- c("symbol", "count")
        hapMatcherA[, iGrid] <- as.integer(match(rhb_t[, iGrid], a[, "symbol"]))
        ## here store both the count and the label
        all_symbols[[iGrid]] <- a
    }
    return(
        list(
            hapMatcherA = hapMatcherA,
            all_symbols = all_symbols
        )
    )
}



