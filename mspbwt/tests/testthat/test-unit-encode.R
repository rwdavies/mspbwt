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
    u <- cumsum(u)

    egs <- 20 ## encode grid size
    
    for(egs in c(5, 17, 20)) {
    
        out <- encode_column_of_u(u, egs)    
        out_mat <- out[["out_mat"]]
        out_vec <- out[["out_vec"]]
        
        ## check all values
        recoded <- sapply(0:(length(u) - 1), function(v) {
            as.integer(decode_value_of_u(out_mat, out_vec, v, egs, do_checks = FALSE))
        })
        
        expect_equal(u, recoded)

    }

})
