# Unit tests for the **TMB** model `OUProcess`.

require(TMB)
source("R/smfret-functions.R")

# Compile and load the model.
gr_mod <- "src/TMB/smfret_TMBExports"
dyn.load(dynlib(gr_mod))

context("OU")
# Instead of running unit tests with completely random inputs,
# here's a method to systematically test a few edge cases,
# namely that `p = 1` and `q = 1` work correctly.

# This creates a data frame containing
# each combination of the input arguments.
test_that("OU TMB MakeADFun gives same results as `ou_y_nll()` ",{
  # systematic test cases
  test_cases <- 1:5
# For each test case, we're going to randomly generate a dataset
# and check that the R and TMB negative loglikelihoods are off by
# the identical numerical constant (depending on `y`, `X`, and `Z`)
# for any set of parameters `theta = (beta, gamma)`.

  ntheta <- 20 # number of parameter sets per test case

  # cycle through test cases
  ou_test_out <- sapply(test_cases, function(ii) {
    t_gamma <- runif(1, ii, ii + 1)
    t_mu <- rnorm(1, ii)
    t_sigma <- runif(1, ii, ii + 1)
    beta0 <- 10
    beta1 <- 1
    dt <- sample(1:10, 1)
    n_obs <- 99
    ntheta <- 10 # number of parameter sets per test
    # simulate data
    x0 <- rnorm(1, t_mu, sqrt(t_sigma^2 / 2 / t_gamma))
    Y <- c(NA)
    while (anyNA(Y)) {
      print("in NA")
      X <- ou_sim(t_gamma, t_mu, t_sigma, dt, n_obs, x0)
      Y <- y_sim(X, beta0, beta1)
    }

    nll_diff <- replicate(ntheta, expr = {
      gamma <- runif(1, 0.01, ii + 10)
      mu <- rnorm(1, ii+10)
      sigma <- runif(1, 0.01, ii + 10)

      nll_r <- ou_y_nll(gamma, mu, sigma, beta0, beta1, X, Y, dt)

      f <- MakeADFun(data = list(model_type = "ou", dt = dt, y = Y, beta0 = beta0, beta1 = beta1), parameters = list(X =X, gamma = gamma, mu = mu, sigma = sigma))

      nll_tmb <- f$fn()
      nll_r - nll_tmb
  })
  # `nll_diff` should contain a vector of `ntheta` identical values
  # the following checks that they are all equal,
  # i.e., that the largest difference between any two is very small
  expect_equal(max(abs(diff(nll_diff))),0,tolerance=1e-3)
  })
})

# display the maximum absolute difference next to each test case
#print("OU results")
#cbind(test_cases, max_diff = signif(ou_test_out, 2))

# test_that("BM TMB MakeADFun gives same results as `bm_y_nll()` ",{
#   # cycle through test cases
#   test_cases <- 1:5
#   bm_test_out <- sapply(test_cases, function(ii) {
#     t_sigma <- runif(1, ii, ii + 1)
#     beta0 <- 10
#     beta1 <- 1
#     dt <- sample(1:10, 1)
#     n_obs <- 99
#     ntheta <- 10 # number of parameter sets per test
#     # simulate data
#     x0 <- rnorm(1, ii)
#     Y <- c(NA)
#     while (anyNA(Y)) {
#       print("in NA")
#       X <- bm_sim(0, t_sigma, dt, n_obs, x0 = x0)
#       Y <- y_sim(X, beta0, beta1)
#     }
#
#     nll_diff <- replicate(ntheta, expr = {
#       sigma <- runif(1, 0.01, ii + 10)
#
#       nll_r <- bm_y_nll(0, sigma, beta0, beta1, X, Y, dt)
#
#       f <- MakeADFun(data = list(model_type = "bm", dt = dt, Y = Y, beta0 = beta0, beta1 = beta1), parameters = list(X=X, sigma = sigma))
#       nll_tmb <- f$fn()
#       nll_r - nll_tmb
#   })
#   # `nll_diff` should contain a vector of `ntheta` identical values
#   # the following checks that they are all equal,
#   # i.e., that the largest difference between any two is very small
#   #max(abs(diff(nll_diff)))
#   })
# })

#print("Brownian results")
#cbind(test_cases, max_diff = signif(bm_test_out, 2))
