## w = 10
get_row_col <- function(irow, icol, T, K, w = 10) {
    row <- c(1, 2, round(K / 2), K - 1, K)[irow]
    c1 <- c(1, 2, round(T / 2), T - w + 1 - 1, T - w + 1)[icol]
    c2 <- c1 + w - 1
    c(row = row, c1 = c1, c2 = c2)
}

##   K <- 250 ## number of haps
##   T <- 50 ## number of SNPs
test_driver_simple <- function(
    K = 250,
    T = 50,
    irow = NULL,
    icol = NULL,
    w = 10
) {
    X <- array(sample(c(0L, 1L), K * T, replace = TRUE), c(K, T))
    Z <- rep(1L, T)
    if (is.null(irow)) {
        ## now add things to match to
        X[c(10, 18, 25, 15, 20, 40) , ] <- 0L
        X[c(15, 20), 1:20] <- 1L
        X[c(10, 18, 25), 11:30] <- 1L
        X[c(40), 26:50] <- 1L
        X[170:250, ] <- 0L
    } else {
        ## choose where to put the match
        a <- get_row_col(irow, icol, T, K, w)
        X[a["row"], ] <- 0L ## blank out
        X[a["row"], a["c1"]:a["c2"]] <- 1L ## make the match!
    }
    X1 <- X + 1L
    Z1 <- Z + 1L
    list(
        X = X,
        Z = Z,
        X1 = X1,
        Z1 = Z1
    )
}


check_expected_top_match <- function(top_matches, irow, icol, K, T, w = 10, is_grid_check_snps = FALSE) {
    a <- get_row_col(irow, icol, T, K, w)
    i <- match(a["row"], top_matches[, "indexB0"] + 1)
    expect_equal(length(i), 1)
    ## expect_true(sum(a["row"] %in% top_matches[, "indexB0"]) == 1)
    if (!is_grid_check_snps) {
        expect_equivalent(top_matches[i, "start1"], a["c1"])
        expect_equivalent(top_matches[i, "end1"], a["c2"])    
    } else {
        a["c1"] <- 32 * (a["c1"] - 1) + 1
        a["c2"] <- 32 * a["c2"]
        expect_true(top_matches[i, "start1"] <= a["c1"])
        expect_true(a["c2"] <= top_matches[i, "end1"]) 
    }
}

## at a minimum, these should all match
check_top_matches <- function(top_matches, X, Z, continue = FALSE) {
    for(irow in 1:nrow(top_matches)) {
        i <- top_matches[irow, "indexB0"] + 1
        c1 <- top_matches[irow, "start1"]
        c2 <- top_matches[irow, "end1"]
        s <- sum(X[i, c1:c2] != Z[c1:c2])
        if (s > 0) {
            print("top match inaccurate")
            print(top_matches[irow, ])
            print(paste0(X[i, c1:c2], collapse = ","))
            print(paste0(Z[c1:c2], collapse = ","))
            if (!continue) {
                stop("top match innacurate")
            }
        }
    }
}

build_and_check_indices <- function(
    hapMatcher,
    all_symbols,
    indices = NULL,
    check_vs_indices = TRUE,
    do_checks = TRUE
) {
    X1C <- hapMatcher
    ms_indices <- ms_BuildIndices_Algorithm5(
        X1C = X1C,
        all_symbols = all_symbols,
        check_vs_indices = check_vs_indices,
        indices = indices,
        do_checks = do_checks
    )
    ## check simple version with Rcpp
    ms_indices_with_Rcpp <- ms_BuildIndices_Algorithm5(
        X1C = X1C,
        all_symbols = all_symbols,
        check_vs_indices = check_vs_indices,
        indices = indices,
        with_Rcpp = TRUE
    )
    expect_equal(ms_indices, ms_indices_with_Rcpp)
    ## check complete version with Rcpp
    ms_indices_only_Rcpp <- Rcpp_ms_BuildIndices_Algorithm5(
        X1C = X1C,
        all_symbols = all_symbols,
        indices = list()
    )
    expect_equal(ms_indices, ms_indices_only_Rcpp)
    ms_indices
}




##  
##  
test_driver_multiple <- function(
    K = 95,
    nGrids = 12,
    irow = 1,
    icol = 1,
    w = 6,
    nMaxDH = 6
) {
    T <- nGrids * 32
    ## build SNP X
    X <- array(1L, c(K, T))
    Z <- rep(1L, T)
    a <- get_row_col(irow, icol, nGrids, K, w) ## in grid notation
    for(j in 1:nGrids) {    
        ## choose 2-5 symbols to be present
        starts <- sort(sample(5:14, 1 + sample(4, 1), replace = FALSE))
        ends <- sort(sample(16:25, length(starts), replace = FALSE))
        for(i in 1:nrow(X)) {
            k <- sample(length(starts), 1)
            X[i, 32 * (j - 1) + starts[k]:ends[k]] <- 0L
        }
        ## do Z as well, always choose "first" one
        k <- 1
        Z[32 * (j - 1) + starts[k]:ends[k]] <- 0L
        ## make sure focal X otherwise is distinct
        k <- 2
        X[a["row"], 32 * (j - 1) + starts[k]:ends[k]] <- 0L
    }
    ## now make the match
    which_snps <- (32 * (a["c1"] - 1) + 1):(32 * a["c2"])
    X[a["row"], which_snps] <- Z[which_snps]
    ## re-name with "S" for SNP to make clear
    Xs <- X
    Zs <- Z 
    ## rest
    rhb_t <- make_rhb_t_from_rhi_t(X)    
    out <- QUILT::make_rhb_t_equality(
        rhb_t = rhb_t,
        nSNPs = T,
        nMaxDH = nMaxDH,
        ref_error = 0.001
    )
    hapMatcher <- out[["hapMatcher"]]
    all_symbols <- out[["all_symbols"]]
    Zg <- make_rhb_t_from_rhi_t(matrix(Z, nrow = 1))
    Z <- map_Z_to_all_symbols(Zg, all_symbols) 
    ##    
    list(
        Xs = Xs,
        Zs = Zs,
        hapMatcher = hapMatcher,
        all_symbols = all_symbols,
        Z = Z
    )
}



## slow but probably works
exhaustive_top_matches_checker <- function(X, Z, top_matches, return_only = FALSE) {
    ## 
    ## get index, start and end of each match
    ##
    Y <- sapply(1:nrow(X), function(i) {
        a <- rle(X[i, ] == Z)
        b <- a$lengths
        c <- a$values
        s <- c(1, 1 + cumsum(b[-length(b)]))
        e <- cumsum(b)
        cbind(index = i, length = b[c], start1 = s[c], end1 = e[c])
    })
    Y <- cbind(
        index = unlist(sapply(Y, function(x) x[, 1])),
        length = unlist(sapply(Y, function(x) x[, 2])),
        start1 = unlist(sapply(Y, function(x) x[, 3])),
        end1 = unlist(sapply(Y, function(x) x[, 4]))
    )
    etm <- NULL ## expected top matches
    while(nrow(Y) > 0) {
        ## 
        w <- which(Y[, 2] == max(Y[, 2]))
        ## here keep only those with first stop and end
        w <- w[which(
            Y[w, "start1"] == Y[w[1], "start1"] &
            Y[w, "end1"] == Y[w[1], "end1"]
        )]
        etm <- rbind(
            etm,
            Y[w, ]
        )
        w <- w[1]
        ## remove anything contained in this
        s <- Y[w, "start1"]
        e <- Y[w, "end1"]
        to_remove <- (Y[, "start1"] >= s) & (Y[, "end1"] <= e)
        Y <- Y[!to_remove, ]
    }
    etm <- etm[order(etm[, "start1"], etm[, "index"]), ]
    top_matches <- top_matches[order(top_matches[, "start1"], top_matches[, "indexB0"]), ]
    if (return_only) {
        return(etm)
    }
    ## 
    expect_equivalent(etm[, "index"], top_matches[, "indexB0"] + 1)
    expect_equivalent(etm[, "start1"], top_matches[, "start1"])
    expect_equivalent(etm[, "end1"], top_matches[, "end1"])
    etm
}





    ## M <- array(0, dim(X))
    ## for(i in 1:nrow(X)) {
    ##     a <- rle(X[i, ] == Z)
    ##     s <- 0
    ##     for(j in 1:length(a$l)) {
    ##         s <- s + 1
    ##         e <- s + a$l[j] - 1
    ##         M[i, s:e] <- as.integer(a$l[j]) * as.integer(a$v[j])
    ##         s <- e
    ##     }
    ## }

    
    ## Y <- lapply(1:ncol(X), function(i) {
    ##     ## OK so this gives me (set) of longest
    ##     w <- which(M[, i] == max(M[, i]))
    ##     starts <- sapply(w, function(ww) {
    ##         e <- i
    ##         while((e > 0) && (M[ww, i] == M[ww, e])) {
    ##             e <- e - 1
    ##         }
    ##         e + 1
    ##     })
    ##     a <- cbind(i, max(M[, i]), w, starts)
    ##     ## now need to know if anything also long not contained

    ##     ## AM HERE
    ##     ## FIX THIS
        
        
    ## })
    ## Y <- cbind(
    ##     grid = unlist(sapply(Y, function(x) x[, 1])),
    ##     length = unlist(sapply(Y, function(x) x[, 2])),
    ##     who = unlist(sapply(Y, function(x) x[, 3])),
    ##     starts = unlist(sapply(Y, function(x) x[, 4]))
    ## )
    ## Y <- Y[order(Y[, 3], Y[, 1], Y[, 2]), ]
    ## built_top_matches <- NULL
    ## sp <- -1
    ## ep <- -1
    ## kp <- -1
    ## for(i in 1:nrow(Y)) {
    ##     k <- Y[i, "who"] - 1 ## make 0-based
    ##     s <- Y[i, "starts"]
    ##     e <- s + Y[i, "length"] - 1
    ##     ## if anything is different, do!
    ##     if (sp != s & ep != e & kp != k) {
    ##         built_top_matches <- rbind(
    ##             built_top_matches,
    ##             c(k, s, e)
    ##         )
    ##     }
    ##     sp <- s
    ##     ep <- e
    ##     kp <- k
    ## }



