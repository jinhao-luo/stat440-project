# logliklihood tests for the **TMB** model `OUProcess`.

context("OU")

test_that("OU TMB MakeADFun gives same results as `ou_y_nll()` ",{
  test_cases <- 1:5
  
  ntheta <- 20 # number of parameter sets per test case
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
      X <- ou_sim(t_gamma, t_mu, t_sigma, dt, n_obs, x0)
      Y <- y_sim(X, beta0, beta1)
    }
    
    nll_diff <- replicate(ntheta, expr = {
      gamma <- runif(1, 0.01, ii + 10)
      mu <- rnorm(1, ii+10)
      sigma <- runif(1, 0.01, ii + 10)
      
      nll_r <- ou_y_nll(gamma, mu, sigma, beta0, beta1, X, Y, dt)
      
      f <- TMB::MakeADFun(data = list(model_type = "ou", dt = dt, y = Y, beta0 = beta0, beta1 = beta1), parameters = list(X =X, gamma = gamma, mu = mu, sigma = sigma), DLL="smfret_TMBExports")
      nll_tmb <- f$fn()
      nll_r - nll_tmb
    })
    #   # `nll_diff` should contain a vector of `ntheta` identical values
    #   # the following checks that they are all equal,
    #   # i.e., that the largest difference between any two is very small
    diff_result <-max(abs(diff(nll_diff)))
    expect_equal(diff_result,0,tolerance=1e-3)
  })
  
  ou_omega_test_out <- sapply(test_cases, function(ii) {
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
      omega <- exp(-gamma*dt)
      tau <- sigma/sqrt(2*gamma)
      
      nll_r <- ou_y_nll(gamma, mu, sigma, beta0, beta1, X, Y, dt)
      
      f <- MakeADFun(data = list(model_type = "omega_tau", dt = dt, y = Y, beta0 = beta0, beta1 = beta1), parameters = list(X =X, omega= omega, mu = mu, tau= tau))
      
      nll_tmb <- f$fn()
      print(paste("nll_r=", nll_r, ", nll_tmb=", nll_tmb))
      nll_r - nll_tmb
    })
    diff_result <-max(abs(diff(nll_diff)))
    expect_equal(diff_result,0,tolerance=1e-3)
  })
})