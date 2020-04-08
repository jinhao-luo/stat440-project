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
  t_gamma <- rnorm(ii)
  t_mu <- rnorm(ii)
  t_sigma <- rnorm(ii)
  beta0 <- rnorm(ii)
  beta1 <- rnorm(ii)
  dt <- sample(1:10, 1)
  n_obs <- sample(50:200, 1)
  ntheta <- 10 # number of parameter sets per test
  # simulate data
  x0 <- rnorm(1, mu, sigma^2/2/gamma)
  X <- ou_sim(gamma, mu, sigma, dt, n_obs, x0) 
  Y <- y_sim(X, beta0, beta1)
  
  nll_diff <- replicate(ntheta, expr = {
    # gamma <- rnorm(ii)
    gamma <- rexp(1, 1/ii)
    # mu <- rnorm(ii)
    mu <- rexp(1, 1/ii)
    # sigma <- rnorm(ii)
    sigma <- rexp(1, 1/ii)
    beta0 <- rnorm(ii)
    beta1 <- rnorm(ii)
    
    nll_r <- ou_y_nll(gamma, mu, sigma, beta0, beta1, X, Y, dt)
    # DATA_MATRIX(H); // rate covariate matrix
    
    
    f1=MakeADFun(data=list(model_type="model1",gamma=gamma,y=Y,b0=b0,b1=b1,mu=mu,sigma=sigma,dt=dt),parameters=list(X=X))
    X_hat <- nlminb(f1$par,f1$fn,f1$gr)["X"]
    
    f2 <-  MakeADFun(data=list(model_type="model2",X=X, y=Y, dt=dt, b0=b0, b1=b1),parameters=list(gamma=gamma, mu=mu, sigma=sigma))
    H <- f2$he()

    f3 <-  MakeADFun(data=list(model_type="model3", X=X, y=Y, H=H, dt=dt, b0=b0, b1=b1),parameters=list(gamma=gamma, mu=mu, sigma=sigma))
    nll_tmb <- f3$fn()
    print(paste("nll_tmb is:", nll_tmb))
    print(paste("nll_r is:", nll_r))
    nll_r - nll_tmb
  })
  # `nll_diff` should contain a vector of `ntheta` identical values
  # the following checks that they are all equal,
  # i.e., that the largest difference between any two is very small
  max(abs(diff(nll_diff)))
})

# display the maximum absolute difference next to each test case
cbind(test_cases, max_diff = signif(test_out,2))