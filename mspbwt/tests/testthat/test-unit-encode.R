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



test_that("can encode and decode columns of u", {

    ## these are vectors with entries 0 through some maximum value
    ## here we focus on the complete ones
    ## the sparse ones are easy to encode separately
    u <- integer(102)
    u[c(2:11, 13:14, 16:17, 50:55)] <- 1L
    uori <- cumsum(u)

    egs <- 20 ## encode grid size
    
    for(i_egs in 1:7) {

        u <- uori

        ## make sure to test exact, minus 1, plus 1, and also 2's on length
        egs <- c(20, 20, 20, 20, 20, 5, 17)[i_egs]
        if (egs == 20) {
            u <- u[1:(97 + i_egs)]
        }

        out <- encode_column_of_u(u, egs, efficient = FALSE)    
        out_mat <- out[["out_mat"]]
        out_vec <- out[["out_vec"]]
        
        ## check all values
        recoded <- sapply(0:(length(u) - 1), function(v) {
            as.integer(decode_value_of_u(out_mat, out_vec, v, egs, do_checks = FALSE))
        })
        
        expect_equal(u, recoded)

    }

})


test_that("can encode and decode a usL", {

    ##
    ## build a usL
    ##
    K <- 10000
    Smax <- 10
    usL <- array(0L, c(K + 1, Smax))    
    usL[] <- 0
    prob <- 2 ** (10:1)
    prob <- prob / 10
    s <- sample(1:10, 10000, prob = prob, replace = TRUE)
    ## want to semi-sort s
    s <- s[order(s + rnorm(10000) / 5)]
    
    a <- table(s, useNA = "always")
    a <- a[order(-a)]
    a <- a[a > 0]
    a <- cbind(as.integer(names(a)), a)
    rownames(a) <- NULL
    colnames(a) <- c("symbol", "count")
    symbol_count_at_grid <-  a
    for(k in 0:(K - 1)) { ## haps (1-based)
        ## now - where it goes - 0 based
        usL[k + 1 + 1,] <- usL[k + 1, ]
        usL[k + 1 + 1, s[k]] <- usL[k + 1 + 1, s[k]] + 1
    }

    
    usL_encoding <- encode_usL(
        usL,
        symbol_count_at_grid,
        egs = 100,
        n_min_symbols = 100
    )

    object.size(usL_encoding)
    object.size(s)
    

})
