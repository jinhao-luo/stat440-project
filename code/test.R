
# Unit tests for the **TMB** model `GammaRegression.cpp`.

require(TMB)
source("smfret-functions.R")

# Compile and load the model.
# Since the model name `GammaRegression` gets used over and over,
# save it to a variable
gr_mod <- "OUProcess"
compile(paste0(gr_mod, ".cpp"))
dyn.load(dynlib(gr_mod))

# Systematic unit testing.


# Instead of running unit tests with completely random inputs,
# here's a method to systematically test a few edge cases,
# namely that `p = 1` and `q = 1` work correctly.

# This creates a data frame containing
# each combination of the input arguments.
test_cases <- expand.grid(p = 1:3, q = 1:3)
test_cases

# For each test case, we're going to randomly generate a dataset
# and check that the R and TMB negative loglikelihoods are off by
# the identical numerical constant (depending on `y`, `X`, and `Z`)
# for any set of parameters `theta = (beta, gamma)`.

ntests <- nrow(test_cases) # number of test cases
ntheta <- 10 # number of parameter sets per test case

# cycle through test cases
test_out <- sapply(1:nrow(test_cases), function(ii) {
  # extract arguments from test_cases
  p <- test_cases$p[ii]
  q <- test_cases$q[ii]
  n <- sample(50:100, 1) # don't need to systematically test n
  ntheta <- 10 # number of parameter sets per test
  # simulate data
  X <- matrix(rnorm(n*p), n, p)
  Z <- matrix(rnorm(n*q), n, q)
  beta0 <- rnorm(p)
  gamma0 <- rnorm(q)
  y <- rgamma(n, shape = exp(X %*% beta0), rate = exp(Z %*% gamma0))
  # Construct the TMB object for this particular dataset
  ##
  ## [code goes here]
  ##
  # check that difference between R and TMB nll is the same constant
  # for any value of theta
  nll_diff <- replicate(ntheta, expr = {
    # randomly generate the parameter values
    beta <- rnorm(p)
    gamma <- rnorm(q)
    # loglikelihood calculation in R
    nll_r <- ou_y_nll(beta, gamma, y, X, Z)
    # loglikelihood calculation in TMB
    ##
    f = MakeADFun(data=list(X=X, dt=dt, b0=b0, ),parameters=list(gamma=gamma, mu=mu, sigma=sigma))
    nll_tmb <- f$fn()
    ##
    nll_r - nll_tmb
  })
  # `nll_diff` should contain a vector of `ntheta` identical values
  # the following checks that they are all equal,
  # i.e., that the largest difference between any two is very small
  max(abs(diff(nll_diff)))
})

# display the maximum absolute difference next to each test case
cbind(test_cases, max_diff = signif(test_out,2))
