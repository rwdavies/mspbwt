if (1 == 0) {

    library("testthat")    
    dir <- "~/proj/mspbwt/"
    setwd(paste0(dir, "/mspbwt/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)

}

test_that("can do simple case of encoding and decoding, just in R", {
    
    ## helps with logic!

    u <- c(0, 0, 0, 1, 1, 2, 2, 3)
    u <- as.integer(u)
    egs <- 4

    ## maximal encoding
    out <- encode_maximal_column_of_u(u, egs, efficient = FALSE)
    out_mat <- out[["out_mat"]]
    out_vec <- out[["out_vec"]]

    expect_equal(out_mat[, "value"], c(0, 1))
    expect_equal(out_mat[, "vec_pos"], c(2, 5))
    expect_equal(as.integer(out_vec), as.integer(c(0, 2, 1, 1, 1, 1)))
    
    
    ## check all values
    recoded <- sapply(0:(length(u) - 1), function(v) {
        as.integer(decode_maximal_value_of_u(out_mat, out_vec, v, egs, do_checks = TRUE))
    })
    
    expect_equal(u, recoded)
    

})


test_that("can more exhaustively encode and decode columns of u, stored maximally", {

    ## these are vectors with entries 0 through some maximum value
    ## here we focus on the complete ones
    ## the sparse ones are easy to encode separately
    u <- integer(102)
    u[c(2:11, 13:14, 16:17, 50:55)] <- 1L
    uori <- cumsum(u)

    egs <- 20 ## encode grid size
    i_egs <- 1
    i_end <- 1
    
    for(i_egs in 1:7) {
        
        for(i_end in 1:2) {

            u <- uori
            ## make sure to test exact, minus 1, plus 1, and also 2's on length
            egs <- c(20, 20, 20, 20, 20, 5, 17)[i_egs]
            if (egs == 20) {
                u <- u[1:(97 + i_egs)]
            }
            ## also check tricky final value behaviour
            if (i_end == 2) {
                u[length(u)] <- u[length(u) - 1] + 1
            }
            
            ## maximal encoding
            out <- encode_maximal_column_of_u(u, egs, efficient = TRUE)
            out_mat <- out[["out_mat"]]
            out_vec <- out[["out_vec"]]
            
            ## check all values
            recoded <- sapply(0:(length(u) - 1), function(v) {
                as.integer(decode_maximal_value_of_u(out_mat, out_vec, v, egs, do_checks = FALSE))
            })
            
            expect_equal(u, recoded)
            
            ## check C++
            out2 <- Rcpp_encode_maximal_column_of_u(u, egs, efficient = TRUE)
            expect_equal(out[["out_mat"]], out2[["out_mat"]])
            expect_equal(out[["out_vec"]], out2[["out_vec"]])

            ## check all values using C++
            recoded2 <- sapply(0:(length(u) - 1), function(v) {
                as.integer(Rcpp_decode_maximal_value_of_u(out_mat, out_vec, v, egs, do_checks = FALSE))
            })
            expect_equal(u, recoded2)            
            
            

        }
    }

})


test_that("can encode and decode full u column, stored maximally, with offsets", {

    ## these are vectors with entries 0 through some maximum value
    ## here we focus on the complete ones
    ## the sparse ones are easy to encode separately
    u <- integer(102)
    u[c(2:11, 13:14, 16:17, 50:55)] <- 1L
    uori <- cumsum(u)
    u <- uori
    egs <- 20 ## encode grid size
            
    ## maximal encoding
    out <- encode_maximal_column_of_u(u, egs, efficient = TRUE)
    out_mat <- out[["out_mat"]]
    out_vec <- out[["out_vec"]]

    ## put into new things
    bigger_out_mat <- matrix(as.integer(NA), nrow = 30, ncol = 2)
    bigger_out_mat[10 + 1:nrow(out_mat), ] <- out_mat
    colnames(bigger_out_mat) <- colnames(out_mat)
    bigger_out_vec <- integer(50)
    bigger_out_vec[] <- as.integer(NA)
    bigger_out_vec[20 + 1:length(out_vec)] <- out_vec
    offsets <- as.integer(c(range(10 + 1:nrow(out_mat)) - 1, range(20 + 1:length(out_vec)) - 1))

    ## check all values
    recoded1 <- sapply(0:(length(u) - 1), function(v) {
        as.integer(decode_maximal_value_of_u(out_mat, out_vec, v, egs, do_checks = FALSE))
    })
    recoded2 <- sapply(0:(length(u) - 1), function(v) {
        as.integer(decode_maximal_value_of_u(
            bigger_out_mat, bigger_out_vec,
            v, egs,
            do_checks = FALSE,
            use_offsets = TRUE,
            offsets = offsets
        ))
    })
    
    expect_equal(u, recoded1)
    expect_equal(u, recoded2)    
            
    ## check C++
    out2 <- Rcpp_encode_maximal_column_of_u(u, egs, efficient = TRUE)
    expect_equal(out[["out_mat"]], out2[["out_mat"]])
    expect_equal(out[["out_vec"]], out2[["out_vec"]])
    
    ## check all values using C++
    recoded2 <- sapply(0:(length(u) - 1), function(v) {
        as.integer(Rcpp_decode_maximal_value_of_u(out_mat, out_vec, v, egs, do_checks = FALSE))
    })
    expect_equal(u, recoded2)            
    

})


test_that("can encode and decode columns of u, stored minimally", {

    for(i in 1:3) {

        if (i == 1) {
            ## so the encoding is something like this
            u <- integer(102)
            u[c(4:11, 13:14, 16:17, 50:55)] <- 1L
        } else if (i == 2) {
            u <- integer(102)
            u[2] <- 1L
        } else if (i == 3) {
            u <- integer(102)
            u[length(u)] <- 1L
        }
        uori <- cumsum(u)
        u <- uori

        out1 <- encode_minimal_column_of_u(u)
        out2 <- as.integer(which(diff(u) > 0))
        out3 <- Rcpp_encode_minimal_column_of_u(u)    
        expect_equal(out1, out2)
        expect_equal(out1, out3)
        
        ## check all values
        recoded1 <- sapply(0:(length(u) - 1), function(v) {
            as.integer(decode_minimal_value_of_u(out1, v))
        })
        recoded2 <- sapply(0:(length(u) - 1), function(v) {
            Rcpp_decode_minimal_value_of_u(out3, v)
        })
        
        expect_equal(u, recoded1)
        expect_equal(u, recoded2)    

    }
})


test_that("can decode columns of U stored with offsets", {

    ## so the encoding is something like this
    u <- integer(102)
    u[c(4:11, 13:14, 16:17, 50:55)] <- 1L
    uori <- cumsum(u)
    u <- uori

    out1 <- encode_minimal_column_of_u(u)
    out1_bigger <- integer(200)
    out1_bigger[100 + 1:length(out1)] <- out1
    ## offsets are 0-based, these are inclusive
    offsets <- as.integer(range(100 + 1:length(out1)) - 1)

    ## check all values
    recoded1 <- sapply(0:(length(u) - 1), function(v) {
        as.integer(decode_minimal_value_of_u(out1, v))
    })
    recoded2 <- sapply(0:(length(u) - 1), function(v) {
        as.integer(decode_minimal_value_of_u(out1_bigger, v, use_offsets = TRUE, offsets))
    })
        
    expect_equal(u, recoded1)
    expect_equal(u, recoded2)
        
})





test_that("can encode and decode a usL", {

    ##
    ## build a usL
    ##
    K <- 10000
    Smax <- 10
    usg <- array(0L, c(K + 1, Smax)) ## us, for a grid
    usg[] <- 0
    prob <- 2 ** (10:1)
    prob <- prob / 10
    sx <- sample(1:10, 10000, prob = prob, replace = TRUE)
    ## want to semi-sort s
    sx <- sx[order(sx + rnorm(10000) / 5)]
    
    a <- table(sx, useNA = "always")
    a <- a[order(-a)]
    a <- a[a > 0]
    a <- cbind(as.integer(names(a)), a)
    rownames(a) <- NULL
    colnames(a) <- c("symbol", "count")
    symbol_count_at_grid <-  a
    for(k in 0:K) { ## haps (1-based)
        ## now - where it goes - 0 based
        if (k > 0) {
            usg[k + 1,] <- usg[k, ]
        }
        usg[k + 1, sx[k]] <- usg[k + 1, sx[k]] + 1
    }
    ## check usg properly made
    expect_equivalent(usg[K + 1, ], a[, 2])
    
    egs <- 100
    n_min_symbols <- 100
    
    usge <- encode_usg(
        usg = usg,
        symbol_count_at_grid = symbol_count_at_grid,
        egs = egs,
        n_min_symbols = n_min_symbols
    )

    v <- 0
    
    for(v in c(0, 1, 2, K + -2:0)) {

        f <- decode_value_of_usge
        
        for(f in c(decode_value_of_usge, Rcpp_decode_value_of_usge)) {
            decoded <- sapply(1:nrow(symbol_count_at_grid), function(s) {

                f(
                    usge = usge,
                    s = s,
                    v = v,
                    egs = egs,
                    n_min_symbols = n_min_symbols
                )
                
            })
            
            expect_equivalent(
                usg[v + 1, ],
                decoded
            )

        }
        
    }

})






test_that("can encode and decode a complete encoding, either list format (usge_all), or super encoding (very large out_mat, out_vec, minimal, etc)", {

    ##
    ## so usge_all is a complete encoding
    ## with a list with one entry per SNP in order
    ## within you have another list, one per symbol
    ## each entry is the encoding for that symbol
    ##
    ## with the super encoding, you have the above turned into 5 objects
    ## the three objects containing the data from the above
    ## plus storage_info, which tells you how to get to them
    ## plus another vector which tells you how to get into storage_info

    usge_all <- list(1:5) ## usge = the below + (e)ncoded
    usg_all <- list(1:5) ## usg = "u plural(s), (g)?"
    all_symbols <- list(1:5)
    egs <- 100
    n_min_symbols <- 100
    
    for(iSNP in 1:5) {
        K <- 10000
        Smax <- sample(7:12, 1)
        prob <- 2 ** (Smax:1)
        prob <- prob / Smax
        sx <- sample(1:Smax, K, prob = prob, replace = TRUE)
        ## ensure no missing values and re-define Sx
        sx <- match(sx, unique(sx))
        Smax <- max(sx)
        if (iSNP == 3) {
            ## here force the max to be just 1 of them
            w <- which(sx == max(sx))
            if (length(w) > 1) {
                sx[w[-1]] <- max(sx) - 1
            }
        }
        ## want to semi-sort s
        sx <- sx[order(sx + rnorm(K) / 5)]
        ## now build
        a <- matrix(0, nrow = Smax, ncol = 2)
        a[, 1] <- 1:Smax
        for(i in 1:Smax) {
            a[i, 2] <- sum(sx == i)
        }
        rownames(a) <- NULL
        colnames(a) <- c("symbol", "count")
        symbol_count_at_grid <-  a
        usg <- array(0L, c(K + 1, Smax)) ## us, for a grid
        usg[] <- 0
        for(k in 0:K) { ## haps (1-based)
            ## now - where it goes - 0 based
            if (k > 0) {
                usg[k + 1,] <- usg[k, ]
            }
            usg[k + 1, sx[k]] <- usg[k + 1, sx[k]] + 1
        }
        ## check usg properly made
        expect_equivalent(usg[K + 1, ], a[, 2])
        ##
        usge <- encode_usg(
            usg = usg,
            symbol_count_at_grid = symbol_count_at_grid,
            egs = egs,
            n_min_symbols = n_min_symbols
        )
        all_symbols[[iSNP]] <- a
        usge_all[[iSNP]] <- usge
        usg_all[[iSNP]] <- usg
    }
        

    ##
    ## build alternate encoding at the end
    ##
    out <- build_super_encoding(
        usge_all = usge_all,
        all_symbols = all_symbols
    )
    out_mat <- out$out_mat
    out_vec <- out$out_vec
    minimal <- out$minimal
    storage_info <- out$storage_info
    cumsum_n_symbols_per_grid <- out$cumsum_n_symbols_per_grid
    
    v <- 0
    iSNP <- 1
    
    for(iSNP in 1:5) {
        for(v in c(0, 1, 2, 10, 20, K + -2:0)) {

            usge <- usge_all[[iSNP]]
            usg <- usg_all[[iSNP]]
            ns <- length(usge)
            
            decoded1 <- sapply(1:ns, function(s) {
                decode_value_of_usge(
                    usge = usge,
                    s = s,
                    v = v,
                    egs = egs,
                    n_min_symbols = n_min_symbols
                )
            })
            
            decoded2 <- sapply(1:ns, function(s) {
                decode_value_from_super_encoding(
                    out_mat = out_mat,
                    out_vec = out_vec,
                    minimal = minimal,
                    storage_info = storage_info,
                    cumsum_n_symbols_per_grid = cumsum_n_symbols_per_grid,
                    iSNP0 = iSNP - 1,
                    s = s,
                    v = v,
                    egs = egs
                )
            })
            
            expect_equivalent(usg[v + 1, ], decoded1)
            expect_equivalent(usg[v + 1, ], decoded2)            

        }
    }


})

