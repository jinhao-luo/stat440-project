

context("BM")

test_that("BM TMB MakeADFun gives same results as `bm_y_nll()` ",{
  #cycle through test cases
  test_cases <- 1:5
  bm_test_out <- sapply(test_cases, function(ii) {
    t_sigma <- runif(1, ii, ii + 1)
    beta0 <- 10
    beta1 <- 1
    dt <- sample(1:10, 1)
    n_obs <- 99
    ntheta <- 10 # number of parameter sets per test
    # simulate data
    x0 <- rnorm(1, ii)
    Y <- c(NA)
    while (anyNA(Y)) {
      X <- bm_sim(0, t_sigma, dt, n_obs, x0 = x0)
      Y <- y_sim(X, beta0, beta1)
    }
    nll_diff <- replicate(ntheta, expr = {
      sigma <- runif(1, 0.01, ii + 10)

      nll_r <- bm_y_nll(0, sigma, beta0, beta1, X, Y, dt)

      f <- TMB::MakeADFun(data = list(model_type = "bm", dt = dt, Y = Y, beta0 = beta0, beta1 = beta1), parameters = list(X=X, sigma = sigma), DLL="smfret_TMBExports")
      nll_tmb <- f$fn()
      nll_r - nll_tmb
    })
    #   # `nll_diff` should contain a vector of `ntheta` identical values
    #   # the following checks that they are all equal,
    #   # i.e., that the largest difference between any two is very small
    diff_result <-max(abs(diff(nll_diff)))
    expect_equal(diff_result,0,tolerance=1e-3)
  })
  #expect_equal(0,0,tolerance=1e-3)
})

