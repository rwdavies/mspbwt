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
    n_min_symbols = 100,
    with_Rcpp = FALSE
) {
    if (1 == 0) {
        egs = 100
        n_min_symbols = 100
        do_checks <- FALSE
        check_vs_indices <- TRUE
        with_Rcpp <- FALSE
        return_d <- TRUE        
    }
    if (do_checks | check_vs_indices) {
        return_d <- TRUE
    }
    if (do_checks) {
        ## here check validity of X1X
        stopifnot(class(X1C[1, 1]) == "integer")
        for(iCol in 1:ncol(X1C)) {
            m1 <- min(X1C[, iCol])            
            m2 <- max(X1C[, iCol])
            stopifnot(m1 == 1)
            stopifnot(sort(unique(X1C[, iCol])) == (m1:m2))
        }
    }
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
    ## 
    for(t in 1L:as.integer(T)) {
        ##
        ## re-set
        ##
        St <- as.integer(n_symbols_per_grid[t]) ## number of symbols in this grid
        symbol_count <- as.integer(all_symbols[[t]][, "count"])
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
        ## op <- options(digits.secs = 6); print("--out--"); print(Sys.time())
        if (!with_Rcpp) {
            a <- out$a
            d <- out$d
            usg <- out$usg
            usg_check <- out$usg_check
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
                stopifnot(a[w, t + 1] == a[k + 1, t])
            }
        }
    }
    to_return <- list(
        a = a,
        d = d,
        usge_all = usge_all,
        egs = egs,
        n_min_symbols = n_min_symbols,
        all_symbols = all_symbols
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
    do_checks
) {
    ##
    ## 
    ## get count of number of each
    usge <- list(1:St) ## us for this (g)rid (e)ncoded
    first_usg_minimal_symbol <- 1 ## 1-based
    prev_value <- symbol_count[1]
    for(s in 1:St) {
        if (symbol_count[s] > n_min_symbols & symbol_count[s] <= prev_value) {
            first_usg_minimal_symbol <- first_usg_minimal_symbol + 1
        } else {
            usge[[s]] <- rep(-1, symbol_count[s])
        }
        prev_value <- symbol_count[s]
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
        ## match_start <- prev_d[k + 1]
        match_start <- d[k + 1, t]
        for(i in 0:(St - 1)) {
            if (match_start > pqs[i + 1]) {
                pqs[i + 1] <- match_start
            }
        }
        ## now - where it goes - 0 based
        val <- c(val, start_count[s] + nso[s] + 1)
        a[start_count[s] + nso[s] + 1, t + 1] <- a[k + 1, t]
        d[start_count[s] + nso[s] + 1, t + 1] <- as.integer(pqs[s])
        ## d_vec[start_count[s] + nso[s] + 1] <- as.integer(pqs[s])
        usg[k + 1 + 1,] <- usg[k + 1, ]
        if (s < first_usg_minimal_symbol) {
            usg[k + 1 + 1, s] <- usg[k + 1 + 1, s] + 1L
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
            usge[[s + 1]] <- Rcpp_encode_maximal_column_of_u(usg[, s + 1], egs = egs)
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
    make_plot = FALSE
) {
    ## 
    if (make_plot) pdf(pdfname, height = nrow(X) / 2 * 1.25 / 2, width = 8)
    ## 
    K <- nrow(X)
    T <- ncol(X)
    ## indices
    a <- ms_indices[["a"]]
    d <- ms_indices[["d"]]
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
    ##
    ## just do easy bit for now
    ##
    wf <- function(k, t, s, usge_all, all_symbols) {
        if (s == 0) {
            s <- nrow(all_symbols[[t]])
        }
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
        if (verbose && t <= 120) {
            print_or_message(paste0("Start of loop t=", t, ", fc = ", fc, ", gc = ", gc, ", ec = ", ec, ", Z[t] = ", Z[t],", f1=", f1, ", g1=", g1, ", e1 = ", e1))
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
            ##d_vec <- decompress_d(d_store, t + 1, K)
            e1 <- d[f1 + 1, t + 1] - 1 ## this is 0-based, probably!
            if (verbose) print(paste0("e1 = ", e1))
            fc <- f1; gc <- g1 ## for visualization efficiency
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
                index <- a[f1 + 1, t + 1] ## a
                if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
                ## should skip if both
                while (Z[e1 - 1 + 1] == X[index + 1, e1 - 1 + 1]) {
                    e1 <- e1 - 1
                    if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)                    
                }
                while (d[f1 + 1, t + 1] <= e1) {
                    f1 <- f1 - 1
                    if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
                }
            }
            if (matches_lower) {
                g1 <- g1 + 1
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
                while ((g1 < K) && (d[g1 + 1, t + 1] <= e1)) { ## d
                    g1 <- g1 + 1
                    if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)                    
                }
            }
            if (verbose) print(paste0("e1 = ", e1))
            if (verbose) print("save stop")
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
    ##
    return(top_matches)
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






## here if we have some X (original matrix with symbols), and some Z with the same encoding
## defined at the SAME GRIDs (i.e. both are for grids and neither for SNPs)
## and we've mapped X to this new encoding, into hapMatcherA, and all_symbols
## we want to do the same with Z
map_Z_to_all_symbols <- function(Z, all_symbols) {
    Z1 <- Z
    Z1[] <- 0L
    for(i in 1:length(Z)) {
        a <- all_symbols[[i]]
        Z1[i] <- match(Z[i], a[, "symbol"])
        if (is.na(Z1[i])) {
            Z1[i] <- 0
        }
        ## now - we are OK with this if this is allowed in hapMatcher
        ## i.e. if hapMatcher has 0's, this is allowed
        ## however, if in all_symbols, the first and last symbol are DIFFERENT, there ARE no hapMatcher 0 values
        ## HOWEVER, if hapMatcher does NOT, we want to map to available values
        if (Z1[i] == 0 & (a[1, 1] != a[nrow(a), 1])) {
            dist <- calc_dist_between_rhb_t_and_hap(
                matrix(a[, 1], ncol = 1),
                STITCH::rcpp_int_expand(Z[i], 32),
                32
            )
            Z1[i] <- which(dist == min(dist))[1]
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
