# Unit tests for the **TMB** model `OUProcess`.

require(TMB)
source("smfret-functions.R")

# Compile and load the model.
gr_mod <- "main"
# compile(paste0(gr_mod, ".cpp"))
dyn.load(dynlib(gr_mod))

# Instead of running unit tests with completely random inputs,
# here's a method to systematically test a few edge cases,
# namely that `p = 1` and `q = 1` work correctly.

# This creates a data frame containing
# each combination of the input arguments.
test_cases <- 1:5
test_cases

# For each test case, we're going to randomly generate a dataset
# and check that the R and TMB negative loglikelihoods are off by
# the identical numerical constant (depending on `y`, `X`, and `Z`)
# for any set of parameters `theta = (beta, gamma)`.

ntests <- nrow(test_cases) # number of test cases
ntheta <- 10 # number of parameter sets per test case

# cycle through test cases
test_out <- sapply(test_cases, function(ii) {
  t_gamma <- rexp(1, 1/ii)
  t_mu <- rexp(1, 1/ii)
  t_sigma <- rexp(1, 1/ii)
  # beta0 <- rnorm(1, ii)
  # beta1 <- rnorm(1, ii)
  beta0 <- 20
  beta1 <- 1
  dt <- sample(1:10, 1)
  n_obs <- sample(50:200, 1)
  ntheta <- 10 # number of parameter sets per test
  # simulate data
  x0 <- rnorm(1, t_mu, sqrt(t_sigma^2/2/t_gamma))
  X <- ou_sim(t_gamma, t_mu, t_sigma, dt, n_obs, x0) 
  Y <- y_sim(X, beta0, beta1)
  
  
  nll_diff <- replicate(ntheta, expr = {
    gamma <- rexp(1, 1/ii)
    mu <- rexp(1, 1/ii)
    sigma <- rexp(1, 1/ii)
    # gamma <- t_gamma
    # mu <- t_mu
    # sigma <- t_sigma

    
    nll_r <- ou_y_nll(gamma, mu, sigma, beta0, beta1, X, Y, dt)
    
    nll_x <- ou_nll(gamma, mu, sigma, X, dt)
    
    f <- MakeADFun(data=list(model_type="ou",x0=x0, dt=dt, y=Y,b0=beta0,b1=beta1, niter=100),parameters=list(gamma=gamma, mu=mu, sigma=sigma))
    nll_tmb <- f$fn()
    # print(paste("nll_tmb is:", nll_tmb))
    # print(paste("nll_r is:", nll_r))
    # print(paste("nll_x is:", nll_x))
    nll_r - nll_tmb
  })
  # `nll_diff` should contain a vector of `ntheta` identical values
  # the following checks that they are all equal,
  # i.e., that the largest difference between any two is very small
  print(paste("nll_diff is:", nll_diff))
  max(abs(diff(nll_diff)))
})

# display the maximum absolute difference next to each test case
cbind(test_cases, max_diff = signif(test_out,2))