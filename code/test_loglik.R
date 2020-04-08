
# Unit tests for the **TMB** model `OUProcess`.

require(TMB)
source("smfret-functions.R")

# Compile and load the model.
gr_mod <- "main"
compile(paste0(gr_mod, ".cpp"))
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
  gamma <- rnorm(ii)
  mu <- rnorm(ii)
  sigma <- rnorm(ii)
  beta0 <- rnorm(ii)
  beta1 <- rnorm(ii)
  dt <- sample(1:10, 1)
  n_obs <- sample(50:200, 1)
  gamma <- rnorm(ii)
  mu <- rnorm(ii)
  sigma <- rnorm(ii)
  beta0 <- rnorm(ii)
  beta1 <- rnorm(ii)
  dt <- sample(1:10, 1)
  n_obs <- sample(50:200, 1)
  # print(X)
  # print(Y)
  
  nll_diff <- replicate(ntheta, expr = {
    gamma <- rnorm(ii)
    mu <- rnorm(ii)
    sigma <- rnorm(ii)
    beta0 <- rnorm(ii)
    beta1 <- rnorm(ii)
    
    data <- list(gamma=gamma, y=Y, b0=beta0, b1=beta1, mu=mu, sigma=sigma, dt=dt)
    data$model_type <- "model1"
    nll_r <- ou_y_nll(gamma, mu, sigma, beta0, beta1, X, Y, dt)
    f1 <-  MakeADFun(data=data,parameters=list(X=X))
    # X_hat <- nlminb(f$par,f$fn,f$gr,lower=c(-10,0.0),upper=c(10.0,10.0))["X"]   
    
    # f2 <-  MakeADFun(data=list(X=X, y=Y, dt=dt, b0=b0, b1=b1),parameters=list(gamma=gamma, mu=mu, sigma=sigma))
    # H <- f2$he()
    
    # f3 <-  MakeADFun(data=list(X=X_hat, y=Y, H=H, dt=dt, b0=b0, b1=b1),parameters=list(gamma=gamma, mu=mu, sigma=sigma),DDL="OUProcess")
    # nll_tmb <- f3$fn()
    nll_tmb <- 0
    
    nll_r - nll_tmb
  })
  # `nll_diff` should contain a vector of `ntheta` identical values
  # the following checks that they are all equal,
  # i.e., that the largest difference between any two is very small
  max(abs(diff(nll_diff)))
})

# display the maximum absolute difference next to each test case
cbind(test_cases, max_diff = signif(test_out,2))