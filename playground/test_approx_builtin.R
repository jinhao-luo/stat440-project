# Unit tests for the **TMB** model `OUProcess`.

require(TMB)
source("smfret-functions.R")

gr_mod <- "builtin"
dyn.load(dynlib(gr_mod))

# test theta approx
test_cases <- expand.grid(gamma=seq(1, 10, by=1), mu=10, sigma=1)
beta0 <- 10
beta1 <- 0.5
dt <- 1
n_obs <- 99

ou_test_out <- apply(test_cases, 1, function(tc) {
  t_gamma <- tc[["gamma"]]
  t_mu<- tc[["mu"]]
  t_sigma<- tc[["sigma"]]
  theta <- list(mu=t_mu, sigma=t_sigma, gamma=t_gamma, t = 1 / t_gamma, tau = t_sigma / sqrt(2 * t_gamma))
  # simulate data
  x0 <- rnorm(1, t_mu, sqrt(t_sigma^2 / 2 / t_gamma))
  Y <- c(NA)
  while (anyNA(Y)) {
    print("in NA")
    X <- ou_sim(t_gamma, t_mu, t_sigma, dt, n_obs, x0)
    Y <- y_sim(X, beta0, beta1)
  }

  param <- list(gamma = 1, mu = 10, sigma = 1, X=rep(0, n_obs))
  data <- list(model_type = "ou", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
  f <- MakeADFun(data = data, parameters = param, random = c("X"), inner.control = list(maxit = 100), silent = TRUE)
  result <- optim(par=f$par, fn=f$fn, control = list(maxit=1000, reltol=1e-10))
  theta_hat <- result$par 

  rmse <- sapply(rownames(theta_hat), function (j) {
      sqrt(mean((theta_hat[j,]-theta[[j]])^2, na.rm=FALSE))/theta[[j]]
  })
  test_detail <- list()
  test_detail$true_param <- theta
  test_detail$Y <- Y
  test_detail$rmse <- signif(rmse,2)
  test_detail
})

# display the maximum absolute difference next to each test case
print("OU results")
cbind(test_cases, rmse=t(sapply(ou_test_out, function(tc) {tc$rmse})))
