if (1 == 0) {

    library("STITCH")
    library("crayon")
    library("testthat")
    library("mspbwt")
    dir <- "~/proj/STITCH/"
    setwd(paste0(dir, "/STITCH/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    dir <- "~/proj/mspbwt/"
    setwd(paste0(dir, "/mspbwt/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)


}



## bugs found and (should be!) fixed

test_that("bug found", {

    ##
    ## this bug is a simplified version of a larger bug that was found and was reproducible
    ## the issue was that the mspbwt code was using minimal value of X of 0 not 1
    ##

    a1 <- c(0, 0, 0, 0, 0, 256, 0, 0, 65, 134742016, 671252736, 671090727, 268439584, 541148292, 33, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 4194848, 16777216, 0, 0, 0, 0, 5373952, -2147483644, 69892, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 8392864, -1605565856, 67143176, 9541896, 645922864, 114688, 4723202, 1342177280, -2146303946, -2147467110, 4456449, 35651588, 276865088, 559153184, 151553, 1207959584, 160180225, 67112960, 33562626, 136052736, 4344864, 2098176, -1207959550)
    a2 <- c(0, 0, 14942208, 1024, 268960256, 1048576, 8388609, 0, 280641, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 4194848, 16777216, 0, 0, 0, 0, 5373952, -2147483644, 69892, 638586882, 0, 1024, 0, 1073741824, 10240, 33558560, 67109120, 0, 268992521, 1310736, 67117056, 32, 0, 0, 1073741888, 1209008128, 4260096, -2147483623, 2097, 268437536, 16908288, 0, 2097280, 146819328, 147464, 0, 541065216, 256, 0, 536870912, 134217728, 589824, 0, 0, 33554458, 0, 268435456, 33056, 536871456, 524288, 0, 0, 2, 256, 0, 1024, 0, 256)
    a3 <- c(0, 0, 0, 0, 0, 0, 0, 1073743872, 0, 135008256, 537035008, 671090727, 268439552, 536953988, 33, 0, 0, 32, 512, 0, 0, 65536, 67108864, -2147221504, 0, 0, 0, 8, -2147483644, 33556480, 2048, 1024, 0, 268435488, 1480607872, 1048, 1344290880, 537657344, 561153, 0, -1073741760, 134348800, -2145385952, 23069188, 553680896, 67436544, 1745551368, 1107298560, 1476526080, NA, 69892, 638586882, 0, 0, 0, 1073741824, -2147473408, 33558560, 67109120, 0, 268992521, 1310736, 67133440, 0, 64, 0, 0, 1073741824, 0, -2147221488, 2097, 268437536, 16908288, 2, 2097280, 146852096, 1196040, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    a4 <- c(0, 0, 0, 0, 0, 0, 0, 2048, 0, -2012741632, 2097416, -1442838521, 268439552, 197772, 131328, 260, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 4194848, 16777216, 0, 0, 0, 0, 5373952, -2147483644, 69892, 638586882, 0, 1024, 0, 1073741824, 10240, 33558560, 67109120, 0, 268992521, 1310736, 67117056, 32, 0, 0, 1073741888, 1209008128, 4260096, -2147483623, 2097, 268437536, 16908288, 0, 2097280, 146819328, 147464, 0, 541065216, 256, 0, 536870912, 134217728, 589824, 0, 0, 33554458, 0, 268435456, 33056, 536871456, 524288, 0, 0, 2, 256, 0, 1024, 0, 256)

    f <- function(to_keep, recast = FALSE) {

        a <- rbind(a1, a2, a3, a4)
        a <- a[to_keep, ]
        if (recast) {
            a <- apply(a, 2, function(x) match(x, unique(x)) - 1)
        }
        out <- STITCH::make_rhb_t_equality(
            rhb_t = a,
            nSNPs = ncol(a) * 32,
            nMaxDH = 20,
            ref_error = 0.001,
            use_hapMatcherR = FALSE
        )
        hapMatcher <- out[["hapMatcher"]]
        all_symbols <- out[["all_symbols"]]
        ##
        ms_indices <- ms_BuildIndices_Algorithm5(
            X1C = hapMatcher,
            all_symbols = all_symbols,
            check_vs_indices = FALSE,
            indices = NULL,
            egs = 100,
            n_min_symbols = 100,
            do_checks = TRUE
        )
        ms_indices[["all_usg_check"]] <- NULL
        ms_indices_only_Rcpp <- Rcpp_ms_BuildIndices_Algorithm5(
            X1C = hapMatcher,
            all_symbols = all_symbols,
            indices = list(),
            egs = 100,
            n_min_symbols = 100
        )
        expect_equal(ms_indices, ms_indices_only_Rcpp)

        ## can I make the problem smaller
        N <- 15
        ## i1 <- 100
        ## i2 <- 200
        i1 <-1
        i2 <- 2
        Z <- c(hapMatcher[i1, 1:N], hapMatcher[i2, (N + 1):ncol(hapMatcher)])

        ms_top_matches <- ms_MatchZ_Algorithm5(
            X = hapMatcher,
            ms_indices = ms_indices_only_Rcpp,
            Z = Z,
            verbose = FALSE,
            do_checks = TRUE,
            check_vs_indices = FALSE
        )

        ## this isn't always strictly true but is true here
        expect_equivalent(ms_top_matches[1, "indexB0"], 0)
        expect_equivalent(ms_top_matches[1, "start1"], 1)
        expect_equivalent(ms_top_matches[1, "end1"], 51)
        expect_equivalent(ms_top_matches[2, "indexB0"], 1)
        expect_equivalent(ms_top_matches[2, "start1"], 16)
        expect_equivalent(ms_top_matches[2, "end1"], 100)
    }

    f(c(1, 2, 3, 4))
    f(c(1, 2, 4, 3))
    f(c(1, 2, 4), recast = FALSE)
    f(c(1, 2, 4), recast = TRUE)

})
