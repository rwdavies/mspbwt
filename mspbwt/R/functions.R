BuildPrefixArray <- function(X) {
    K <- nrow(X)
    T <- ncol(X)
    a <- array(NA, c(K, T))
    b <- array(NA, K)
    a[, 1] <- order(X[, 1])
    for(t in 2:T) { ## SNPs
        u <- 1
        v <- 1
        for(k in 1:K) { ## haps
            if (X[a[k, t - 1], t] == 0) {
                a[u, t] <- a[k, t - 1]
                u <- u + 1
            } else {
                b[v] <- a[k, t - 1]
                v <- v + 1
            }
        }
        ## add in b, if v > 1
        if (v > 1) {
            a[u - 1 + 1:(v - 1), t] <- b[1:(v - 1)]
        }
    }
    return(a)
}

displayPrefixArray <- function(X, a, t) {
    ## so show - before k, k, after k
    if (t != 1) {
        to_cat <- paste0(apply(X[a[, t], 1:(t - 1), drop = FALSE], 1, paste0, collapse = ""))
    } else {
        to_cat <- rep("", nrow(X))
    }
    ## Middle bit
    to_cat <- paste0(to_cat, X[a[, t], t], sep = "")
    ## end bit
    if (t != ncol(X)) {
        to_cat <- paste0(to_cat, " ", apply(X[a[, t], (t + 1):ncol(X), drop = FALSE], 1, paste0, collapse = ""))
    }
    ## add order
    to_cat <- paste0(formatC(a[, t], width = nchar(nrow(X))), ":", to_cat)
    cat(to_cat, sep = "\n")
}




BuildPrefixAndDivergenceArray <- function(X) {
    K <- nrow(X)
    T <- ncol(X)
    a <- array(NA, c(K, T)) ## positions
    d <- array(NA, c(K, T)) ## distances
    b <- array(NA, K)
    e <- array(NA, K)
    a[, 1] <- order(X[, 1])
    ## d is a 0, except first entry, and on (ordered) switch
    d[, 1] <- 0
    d[unique(c(1, which.max(X[a[, 1], 1]))), 1] <- 1
    for(t in 2:T) { ## SNPs
        u <- 1
        v <- 1
        p <- t
        q <- t
        for(k in 1:K) { ## haps
            if (d[k, t - 1] > p) {
                p <- d[k, t - 1]
            }
            if (d[k, t - 1] > q) {
                q <- d[k, t - 1]
            }
            if (X[a[k, t - 1], t] == 0) {
                a[u, t] <- a[k, t - 1]
                d[u, t] <- p
                u <- u + 1
                p <- 0 ## correct
            } else {
                b[v] <- a[k, t - 1]
                e[v] <- q
                v <- v + 1
                q <- 0
            }
        }
        ## add in b, if v > 1
        if (v > 1) {
            a[u - 1 + 1:(v - 1), t] <- b[1:(v - 1)]
            d[u - 1 + 1:(v - 1), t] <- e[1:(v - 1)]
        }
    }
    return(list(a = a, d = d))
}






displayPrefixAndDivergenceArray <- function(X, a, d, t) {
    ## if current k, don't do anything
    ## else, start at k + 1 (in 1-based)
    ## so show - before t, t and after t
    if (t != 1) {
        to_cat <- paste0(apply(X[a[, t], 1:(t - 1), drop = FALSE], 1, paste0, collapse = ""))
    } else {
        to_cat <- rep("", nrow(X))
    }
    to_cat <- paste0(to_cat, X[a[, t], t], sep = "")
    ## end bit
    if (t != ncol(X)) {
        to_cat <- paste0(to_cat, " ", apply(X[a[, t], (t + 1):ncol(X), drop = FALSE], 1, paste0, collapse = ""))
    }
    ## add order
    to_cat <- paste0(formatC(a[, t], width = nchar(nrow(X))), ":", to_cat)
    offset <- nchar(nrow(X)) + 1
    ## now - I think OK - bold appropriately!
    for(i in 1:length(to_cat)) {
        if (d[i, t] == t) {
            cat(to_cat[i]) ## No bold
        } else {
            cat(
                crayon::black(substr(to_cat[i], 1, offset + d[i, t])),
                crayon::black$bold$underline(substr(to_cat[i], offset + d[i, t] + 1, offset + t)),
                crayon::black(substr(to_cat[i], offset + t + 1, offset + t + 1 + 1e6)),
                sep = ""
            )
        }
        cat(sep = "\n")
    }
}


BuildPrefixAndDivergenceArrayAndReportLongMatches <- function(X, L) {
    K <- nrow(X)
    T <- ncol(X)
    a <- array(NA, c(K, T)) ## positions
    d <- array(NA, c(K, T)) ## distances
    b <- array(NA, K)
    e <- array(NA, K)
    a[, 1] <- order(X[, 1])
    ## d is a 0, except first entry, and on (ordered) switch
    d[, 1] <- 0
    d[unique(c(1, which.max(X[a[, 1], 1]))), 1] <- 1
    ##
    long_matches <- list()
    for(t in 2:T) { ## SNPs
        u <- 1
        v <- 1
        p <- t
        q <- t
        us <- u
        vs <- v
        for(k in 1:K) { ## haps
            ## do the check
            if (d[k, t - 1] > (t - 1 - L)) {
                if (((u - us) > 0) & ((v - vs) > 0)) {
                    ## if ((t - 1) == 5) stop("WER")
                    to_store <- list(list(
                        t, ## where it terminates
                        a[us:(u - 1), t],
                        b[vs:(v - 1)]
                    ))
                    long_matches <- append(long_matches, to_store)
                }
                us <- u
                vs <- v
            }
            if (d[k, t - 1] > p) {
                p <- d[k, t - 1]
            }
            if (d[k, t - 1] > q) {
                q <- d[k, t - 1]
            }
            if (X[a[k, t - 1], t] == 0) {
                a[u, t] <- a[k, t - 1]
                d[u, t] <- p
                u <- u + 1
                p <- 0
            } else {
                b[v] <- a[k, t - 1]
                e[v] <- q
                v <- v + 1
                q <- 0
            }
        }
        ## add in b, if v > 1
        if (v > 1) {
            a[u - 1 + 1:(v - 1), t] <- b[1:(v - 1)]
            d[u - 1 + 1:(v - 1), t] <- e[1:(v - 1)]
        }
        ## last bit
        if (((u - us) > 0) & ((v - vs) > 0)) {
            ## if ((t - 1) == 5) stop("WER")
            to_store <- list(list(
                t, ## where it terminates
                a[us:(u - 1), t],
                b[vs:(v - 1)]
            ))
            long_matches <- append(long_matches, to_store)
        }
    }
    return(
        list(
            a = a,
            d = d,
            long_matches = long_matches
        )
    )
}

displayLongMatches <- function(X, long_matches, L) {
    for(long_match in long_matches) {
        p <- long_match[[1]]
        cat(paste0(
            "match ending at 1-based position ", p,
            " between haps ", paste0(long_match[[2]], collapse = ", "),
            " and ", paste0(long_match[[3]], collapse = ", "),
            ":"
        ))
        cat("\n")
        for(a in long_match[[2]]) {
            cat(formatC(a, width = nchar(nrow(X))), ":", sep = "")
            cat(X[a, p + (-L:0)], sep = "")
            cat("\n")
        }
        for(a in long_match[[3]]) {
            cat(formatC(a, width = nchar(nrow(X))), ":", sep = "")
            cat(X[a, p + (-L:0)], sep = "")
            cat("\n")
        }
    }
}



BuildPrefixAndDivergenceArrayAndReportSetMaximalMatches <- function(X) {
    K <- nrow(X)
    T <- ncol(X)
    a <- array(NA, c(K, T)) ## positions
    d <- array(NA, c(K + 1, T)) ## distances
    b <- array(NA, K)
    e <- array(NA, K)
    a[, 1] <- order(X[, 1])
    ## d is a 0, except first entry, and on (ordered) switch
    d[, 1] <- 0
    d[unique(c(1, which.max(X[a[, 1], 1]))), 1] <- 1
    ## add sentinels
    d[1, ] <- 1:T
    d[K + 1, ] <- 1:T
    ##
    maximal_matches <- NULL
    for(t in 2:(T + 1)) { ## SNPs
        ##
        ## do sweeping bit
        ##
        for(k in 1:K) {
            match_continues <- FALSE
            m <- k - 1
            n <- k + 1
            dl <- d[k, t - 1]
            du <- d[k + 1, t - 1]
            ##
            ## match to previous neighbour is longer, scan down the haplotypes, decreasing indices
            ##
            if (dl <= du) { ## if lower one is at least min
                while(m >= 1 && d[m + 1, t - 1] <= dl) {
                    ## RHS will be out of bound BUT this is caught by first bit
                    if (t < (T + 1) && X[a[m, t - 1], t] == X[a[k, t - 1], t]) {
                        match_continues <- TRUE ## i.e. we won't be saving anything!
                        break ## should just break out of while loop - as
                    }
                    m <- m - 1
                }
            }
            ##
            ## match to next neighbour is longer, scale up the haplotypes, increasing indices
            ##
            if (!match_continues & (dl >= du)) {
                while((n <= K) && d[n, t - 1] <= du) {
                    if (t < (T + 1) && X[a[n, t - 1], t] == X[a[k, t - 1], t]) {
                        match_continues <- TRUE
                        break ## should just break out of while loop - as
                    }
                    n <- n + 1
                }
            }
            ## if match doesn't continue - i.e. it ends here - record
            if (!match_continues) {
                if ((m + 1) <= (k - 1)) {
                    for(j in (m + 1):(k - 1)) {
                        ##if (d[j, t - 1] < t) {
                        to_return <- matrix(c(a[k, t - 1], a[j, t - 1], dl + 1, t - 1), nrow = 1)
                        maximal_matches <- rbind(maximal_matches, to_return)
                        ##}
                    }
                }
                if ((k + 1) <= (n - 1)) {
                    for(j in (k + 1):(n - 1)) {
                        ##if (d[j + 1, t - 1] < t) {
                        to_return <- matrix(c(a[k, t - 1], a[j, t - 1], du + 1, t - 1), nrow = 1)
                        maximal_matches <- rbind(maximal_matches, to_return)
                        ##}
                    }
                }
            }
        }
        ##
        ## do normal prefix and divergence arrays
        ##
        if (t <= T) {
            u <- 1
            v <- 1
            p <- t
            q <- t
            for(k in 1:K) { ## haps
                match_start <- d[k, t - 1]
                if (match_start > p) {
                    p <- match_start
                }
                if (match_start > q) {
                    q <- match_start
                }
                if (X[a[k, t - 1], t] == 0) {
                    a[u, t] <- a[k, t - 1]
                    d[u, t] <- p
                    u <- u + 1
                    p <- 0
                } else {
                    b[v] <- a[k, t - 1]
                    e[v] <- q
                    v <- v + 1
                    q <- 0
                }
            }
            ## add in b, if v > 1
            if (v > 1) {
                a[u - 1 + 1:(v - 1), t] <- b[1:(v - 1)]
                d[u - 1 + 1:(v - 1), t] <- e[1:(v - 1)]
            }
        }
    }
    colnames(maximal_matches) <- c("indexA1", "indexB1", "start1", "end1")
    return(
        list(
            a = a,
            d = d,
            maximal_matches = maximal_matches
        )
    )
}

DisplaySetMaximalMatches <- function(X, maximal_matches, k) {
    cat(
        paste0(
            "set maximal matches for haplotype: ", k
        ), sep = "\n"
    )
    cat(formatC(k, width = nchar(nrow(X))), ":", sep = "")
    cat(X[k, ], sep = "")
    cat("\n")
    ##
    mm <- maximal_matches[which(maximal_matches[, 1] == k), , drop = FALSE]
    ##
    for(i in 1:nrow(mm)) {
        m <- mm[i, ]
        otherk <- m[2]
        from <- m[3]
        to <- m[4]
        ##
        cat(formatC(otherk, width = nchar(nrow(X))), ":", sep = "")
        ## now bold bit that overlaps otherwise not bold
        if (from > 1) {
            ## start bit
            cat(X[otherk, 1:(from - 1)], sep = "")
        }
        ## middle bit
        cat(crayon::black$bold$underline(X[otherk, from:to]), sep = "")
        if (to < ncol(X)) {
            cat(X[otherk, (to + 1):ncol(X)], sep = "")
        }
        cat("\n")
    }
}



BuildIndices_Algorithm5 <- function(
    X,
    verbose = FALSE,
    do_checks = FALSE,
    do_var_check = TRUE
) {
    ## build arrays including a, d and now u, v, c
    temp <- colSums(X)
    if (do_var_check) {
        if ((sum(temp == 0) > 0) | (sum(temp == nrow(X)) > 0)) {
            stop("At least one column has no variation")
        }
    }
    K <- nrow(X)
    T <- ncol(X)
    a <- array(NA, c(K, T + 1)) ## positions
    d <- array(NA, c(K + 1, T + 1)) ## distances
    a[, 1] <- 0:(K - 1)
    dtemp <- array(NA, K) ## temp for filler
    u <- array(NA, c(K + 1, T))
    u[1, ] <- 0
    v <- array(NA, c(K + 1, T))
    v[1, ] <- 0
    w <- array(NA, c(K, T)) ## eventually unnecessary
    c <- array(NA, T + 1) ## not sure why
    b <- array(NA, K)
    ## a[, 1] <- order(X[, 1])
    ## d is a 0, except first entry, and on (ordered) switch
    ## d[, 1] <- 0
    ##d[unique(c(1, which.max(X[a[, 1], 1]))), 1] <- 1
    ## build first column of a, u, v, d
    uu <- 0
    vv <- 0
    d[, 1] <- 0
    d[, 2] <- 0    
    for(k in 0:(K - 1)) {
        if (X[k + 1, 1] == 0) {
            a[uu + 1, 2] <- k
            uu <- uu + 1
        } else {
            b[vv + 1] <- k
            vv <- vv + 1
        }
        u[k + 1 + 1, 1] <- uu
        v[k + 1 + 1, 1] <- vv
    }
    if (vv > 0) {
        a[uu + 0:(vv - 1) + 1, 2] <- b[0:(vv - 1) + 1]
    }
    ## d is 1 in first entry and on first swap i.e. first v
    d[uu + 1, 2] <- 1
    c[1] <- uu
    c[T + 1] <- 0
    ## add sentinels
    d[1, ] <- 1:(T + 1)
    d[K + 1, ] <- 1:(T + 1)
    ##
    for(t in 2:T) {
        ##
        ## do sweeping bit
        ##
        uu <- 0
        vv <- 0
        p <- t ## 0-based SNP index + 1
        q <- t ## 0-based SNP index + 1
        for(k in 0:(K - 1)) { ## haps (1-based)
            match_start <- d[k + 1, t]
            if (match_start > p) {
                p <- match_start
            }
            if (match_start > q) {
                q <- match_start
            }
            if (X[a[k + 1, t] + 1, t] == 0) {
                a[uu + 1, t + 1] <- a[k + 1, t]
                d[uu + 1, t + 1] <- p
                uu <- uu + 1
                p <- 0
            } else {
                b[vv + 1] <- a[k + 1, t]
                dtemp[vv + 1] <- q
                vv <- vv + 1
                q <- 0
            }
            u[k + 1 + 1, t] <- uu
            v[k + 1 + 1, t] <- vv
        }
        ## add in b, if v > 1
        if (vv > 0) {
            a[uu + 0:(vv - 1) + 1, t + 1] <- b[0:(vv - 1) + 1]
            d[uu + 0:(vv - 1) + 1, t + 1] <- dtemp[0:(vv - 1) + 1]
        }
        if (do_checks) {
            dtemp[] <- NA
            stopifnot(sum(is.na(a[, t + 1])) == 0)
            stopifnot(sum(is.na(d[, t + 1])) == 0)            
        }
        ##
        c[t] <- uu
        ## perform the check
        if (verbose) {
            message(paste0("t = ", t))
        }
        ##
        ## do checks
        ##
        for(k in 0:(K - 1)) {
            if (X[a[k + 1, t] + 1, t] == 0) {
                w[k + 1, t] <- u[k + 1 + 1, t]
            } else {
                w[k + 1, t] <- v[k + 1 + 1, t] + c[t]
            }
            if (verbose) {
                message(paste0(
                    "k=", k, ", ",
                    "u=", u[k + 1 + 1, t], ", ",
                    "v+c=", v[k + 1 + 1, t] + c[t], ", ",
                    "X=", X[a[k + 1, t] + 1, t], ", ",
                    "w=", w[k + 1,t], ", ",                    
                    "a[w[k+1,t],t+1]=", a[w[k + 1, t], t + 1], ", ",
                    "a[k+1,t]=", a[k + 1, t]
                ))
            }
            if (do_checks) {
                stopifnot(a[w[k + 1, t], t + 1] == a[k + 1, t])
            }
        }
    }
    return(
        list(
            a = a,
            u = u,
            v = v,
            c = c,
            d = d
        )
    )
}


MatchZ_Algorithm5 <- function(
    X,
    indices,
    Z,
    verbose = FALSE,
    do_checks = FALSE
) {
    K <- nrow(X)
    T <- ncol(X)
    if (sum(c("a", "u", "v", "c", "d") %in% names(indices)) != 5) {
        stop("bad indices list")
    }
    a <- indices[["a"]]
    u <- indices[["u"]]
    v <- indices[["v"]]
    c <- indices[["c"]]
    d <- indices[["d"]]
    if (length(Z) != (ncol(a) - 1)) {
        stop("Z not the right size")
    }
    ##
    ## OK these all look fine
    ## 
    ## u
    ## t(t(v) + c(c[-length(c)]))
    ## a
    ## d
    ##
    e <- array(NA, T) ## keep this 0-based (like d)
    f <- array(NA, T) ## keep these 0-based
    g <- array(NA, T) ## keep these 0-based
    ##
    ## init - get same match (if it exists)
    ##
    if (Z[1] == 0) {
        ## then f[1] is the start of the 0's i.e. 0
        ## and goes 1 over i.e. to first 1
        f[1] <- 0
        g[1] <- c[1]
        e[1] <- 0
    } else {
        ## then f[1] is the start of the 1's
        f[1] <- c[1]
        g[1] <- K
        e[1] <- 0
    }
    ##
    ## just do easy bit for now
    ##
    wf <- function(k, t, z) {
        if (z == 0) {
            return(u[k + 1, t])
        } else {
            return(v[k + 1, t] + c[t])
        }
    }
    fc <- f[1]
    gc <- g[1]
    ec <- e[1]
    e1 <- NA
    top_matches <- NULL
    for(t in 2:T) {
        ##
        ## so this tells it where to move to next time (this time)
        ## remember here fc and gc are indices in a
        ## recall that w is set up so that the following holds
        ## (a[w[k + 1, t], t + 1] == a[k + 1, t])
        ## i.e. if we're at index k in a[k + 1, t]
        ## then next time we go to index w[k + 1, t] in a at t + 1
        ## 
        ## NOW, when does this collapse? what does this mean?
        ## if there is a match, at least one, we'll always continue
        ## 
        ## if there are no matches, we collapse
        ## why? recall either u or v is increasing or constant
        ## i.e. those samples have a 0 or 1
        ## if all those samples (between fc and gc) move to something
        ## then the opposite of that i.e. the opposite Z will give (u or v) that are constant
        ## hence they will collapse to a single value, indicating the long match is broken
        ##
        ## now, where are we after that collapse?
        ## we are somewhere just before, or just after, where Z would have matched to
        ## 
        ## that new index of a represents where the previous matches would have gone, had they continued
        ##
        f1 <- wf(fc, t, Z[t])
        g1 <- wf(gc, t, Z[t])
        if (verbose) {
            message(paste0("Start of loop t=", t, ", fc = ", fc, ", gc = ", gc, ", ec = ", ec, ", Z[t] = ", Z[t],", f1=", f1, ", g1=", g1, ", e1 = ", e1))
        } 
        if (g1 > f1) {
            ## nothing to do
        } else {
            ## we have reached a maximum - need to report and update e, f, g
            for(k in fc:(gc - 1)) {
                ## so here, we know they're now broken now
                ## so previous SNP was the end of the run
                ## we know they've all broken, and have previous fc and gc to save
                ## we loop from fc to (gc - 1) in 0-based
                ##
                ## together this gives
                ## 
                ## 1st k is the 0-based what hap we're looking at
                ## 2nd is a[k + 1, t] is the 0-based index. this is effectively back a SNP
                ##     as a[, t + 1] is the index for 1-based SNP t
                ## 3rd ec is 0-based, so ec + 1 is 1-based where we start
                ## 4th this is now 1-based t, so we stopped on previous SNP, hence t - 1
                if (verbose) {
                    message("save top match(es)")
                }
                top_matches <- rbind(
                    top_matches,
                    matrix(c(k, a[k + 1, t], ec + 1, t - 1), nrow = 1)
                )
            }
            ## so what are f1 and g1 here? where this thing wants to go to?
            e1 <- d[f1 + 1, t + 1] - 1 ## this is 0-based, probably!
            if ((Z[e1 + 1] == 0 && f1 > 0) || f1 == K) {
                f1 <- g1 - 1
                index <- a[f1 + 1, t + 1]
                ## figure out where they stop matching
                if (verbose) {
                    print(paste0("f1 = ", f1))
                    print(paste0("index = ", index))
                    print(paste0("Z[e1 - 1 + 1] = ", Z[e1 - 1 + 1]))
                }
                while (Z[e1 - 1 + 1] == X[index + 1, e1 - 1 + 1]) {
                    e1 <- e1 - 1
                }
                ## keep all of those that agree over the length
                while (d[f1 + 1, t + 1] <= e1) {
                    f1 <- f1 - 1
                }
            } else if (f1 < K) {
                ## am here, work to understand this
                ## in basic form, if Z of the new startpoint (?) is a 0
                ## figure out new e1
                ## then go up for new g1?
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
        ## perform switch over
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


check_Algorithm5 <- function(X, Z, top_matches, display = FALSE) {
    X2 <- rbind(X, Z)
    out <- BuildPrefixAndDivergenceArrayAndReportSetMaximalMatches(X2)
    maximal_matches <- out$maximal_matches
    slow_matches <- maximal_matches[maximal_matches[, 1] == (nrow(X) + 1), , drop = FALSE]
    ## check they are the same
    expect_equivalent(slow_matches[, "indexB1"], top_matches[, "indexB0"] + 1)
    expect_equal(slow_matches[, "start1"], top_matches[, "start1"])
    expect_equal(slow_matches[, "end1"], top_matches[, "end1"])
    ## visualize them!
    if (display) {
        message("The original way")
        DisplaySetMaximalMatches(X2, slow_matches, k = nrow(X) + 1)
        ## other way
        message("The faster way")
        top_matchesL <- top_matches
        top_matchesL[, 1] <- nrow(X) + 1
        top_matchesL[, 2] <- top_matchesL[, 2] + 1
        DisplaySetMaximalMatches(X2, top_matchesL, k = nrow(X) + 1)
    }
}

