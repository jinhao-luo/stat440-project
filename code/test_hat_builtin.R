# Unit tests for the **TMB** model `OUProcess`.

require(TMB)
require(optimCheck)
require(testthat)
source("multi_start.R")
source("smfret-functions.R")
# graphics.off()

# Compile and load the model.
gr_mod <- "builtin"
dyn.load(dynlib(gr_mod))

# Instead of running unit tests with completely random inputs,
# here's a method to systematically test a few edge cases,
# namely that `p = 1` and `q = 1` work correctly.

# This creates a data frame containing
# each combination of the input arguments.
test_cases <- 1:4
test_cases

# For each test case, we're going to randomly generate a dataset
# and check that the R and TMB negative loglikelihoods are off by
# the identical numerical constant (depending on `y`, `X`, and `Z`)
# for any set of parameters `theta = (beta, gamma)`.

ntheta <- 1 # number of parameter sets per test case

# cycle through test cases
ou_test_out <- sapply(test_cases, function(ii) {
  method <- "Nelder-Mead"
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

    omega <- 0
    multi_start <- TRUE
    num_start <- 5
    start <- 0+0.01
    end <- 1-0.01
    omegas <- seq(start, end, (end-start)/num_start)
    omega <- find_optim_omega(omegas, n_obs, dt, Y, beta0, beta1)

    param <- list(omega= omega, mu = 0, tau= 1, X=rep(0, n_obs))
    data <- list(model_type = "omega_tau", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
    ou_f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE, method=method)

    ou_result <- optim(par = ou_f$par, fn = ou_f$fn, gr = ou_f$gr, control=list(trace=0, maxit=1000, reltol=1e-8), method=method)

    oproj <- optim_proj(fun=ou_f$fn, xsol=ou_result$par, maximize=FALSE, plot=TRUE)

    dproj <- diff(oproj)
    # print(dproj[,c('rel')][["mu"]])
    expect_equal(dproj[,c('rel')][["mu"]], 0, 0.1000001)
    expect_equal(dproj[,c('rel')][["omega"]], 0, 0.1000001)
    expect_equal(dproj[,c('rel')][["tau"]], 0, 0.1000001)

  })
})


bm_test_out <- sapply(test_cases, function(ii) {
  method <- "BFGS"
  # t_sigma <- rexp(1, 1/ii)
  t_sigma <- runif(1, 0, 1/ii)
  beta0 <- 10
  beta1 <- 1
  dt <- 1
  n_obs <- 1000*ii
  ntheta <- 1 # number of parameter sets per test
  # simulate data
  x0 <- rnorm(1, 1+ii)
  Y <- c(NA)
  while (anyNA(Y)) {
    # print("in NA")
    X <- bm_sim(0, t_sigma, dt, n_obs, x0=x0)
    Y <- y_sim(X, beta0, beta1)
  }

  nll_diff <- replicate(ntheta, expr = {
    sigma <- rexp(1, 1/ii)

    nll_r <-  bm_y_nll(0, sigma, beta0, beta1, X, Y, dt)

    param <- list(sigma = ii, X=rep(0, n_obs))
    data <- list(model_type = "bm", dt=dt, Y=Y, beta0 = beta0, beta1 = beta1)
    bm_f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE, method=method)
    bm_result <- optim(par = bm_f$par, fn = bm_f$fn, gr = bm_f$gr, control=list(trace=0, maxit=1000, reltol=1e-8), method=method)
    oproj <- optim_proj(fun=bm_f$fn, xsol=bm_result$par, maximize=FALSE, plot=TRUE)

    dproj <- diff(oproj)
    
    expect_equal(dproj[,c('rel')], 0, 0.1000001)
  })
})
