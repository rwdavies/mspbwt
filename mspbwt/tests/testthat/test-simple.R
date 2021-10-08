
test_that("can properly determine integer is odd using R", {
    expect_equal(R_is_odd(10L), FALSE)
    expect_equal(R_is_odd(11L), TRUE)
})
test_that("can properly determine integer is odd using Rcpp", {
    expect_equal(cpp_is_odd(10L), FALSE)
    expect_equal(cpp_is_odd(11L), TRUE)
})