if (1 == 0) {

    library("crayon")
    library("testthat")
    library("mspbwt")
    dir <- "~/proj/mspbwt/"
    setwd(paste0(dir, "/mspbwt/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)


}

K <- 200 ## number of haps
T <- 50 ## number of SNPs
set.seed(2021)
X <- array(sample(c(0, 1), K * T, replace = TRUE), c(K, T))
Z <- rep(1, T)
## now add things to match to
X[c(10, 18, 25, 15, 20, 40) , ] <- 0
X[c(15, 20), 1:20] <- 1
X[c(10, 18, 25), 11:30] <- 1
X[c(40), 26:50] <- 1



test_that("R Algorithm 5 works with normal implementation", {
    
    indices <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
    top_matches <- MatchZ_Algorithm5(X, indices, Z, verbose = FALSE, do_checks = TRUE)
    check_Algorithm5(X, Z, top_matches, display = FALSE)
    
})


test_that("Rcpp Algorithm 5 works with normal implementation", {
    
    i1 <- BuildIndices_Algorithm5(X, verbose = FALSE, do_checks = TRUE, do_var_check = FALSE)
    i2 <- BuildIndices_Algorithm5_Rcpp(X)
    for(what in c("a", "u", "v", "c", "d")) {
        expect_equivalent(i1[[what]], i2[[what]])
    }

    top_matches <- MatchZ_Algorithm5(X, i1, Z, verbose = FALSE, do_checks = TRUE)
    check_Algorithm5(X, Z, top_matches, display = FALSE)
    
    top_matches_Rcpp <- MatchZ_Algorithm5_Rcpp(X, i2$a, i2$u, i2$v, i2$c, i2$d, Z)

    expect_equal(top_matches_Rcpp, top_matches_Rcpp)    
    
})



