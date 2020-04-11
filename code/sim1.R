require("TMB")
source("smfret-functions.R")
gr_mod <- "main"
# compile(paste0(gr_mod, ".cpp"))
dyn.load(dynlib(gr_mod))

#' simulation 1
#'
#' @param beta0 the value of beta0
#' @param beta1 the value of beta1
#' @param gamma Scalar mean reversion parameter
#' @param mu Scalar mean parameter
#' @param sigma Scalar diffusion parameter
#' @param dt Interobservation time.
#' @param n_obs Number of observations to generate.
#' @param n_dataset Number of dataset to simulate for each `\beta`
#' @return a vector of rmse: rmse for inference on each simulated dataset
sim_1 <- function(beta0, beta1, gamma = 1, mu = 10, sigma = sqrt(2), dt = 1,
                  n_obs = 99, n_dataset = 100) {
    theta <- list(mu=mu, sigma=sigma, gamma=gamma, t = 1 / gamma, tau = sigma / sqrt(2 * gamma))
    theta_hat <- replicate(n_dataset, expr = {
        Y <- c(NA)
        while (anyNA(Y)) {
            X <- ou_sim(gamma, mu, sigma, dt, n_obs)
            Y <- y_sim(X, beta0, beta1)
        }
        f <- MakeADFun(
            data = list(model_type = "ou", niter = 100, dt = dt, y = Y, beta0 = beta0, beta1 = beta1),
            parameters = list(gamma = gamma, mu = mu, sigma = sigma)
        )
        result <- optim(par = f$par, fn = f$fn, gr = f$gr,
            control = list(maxit = 1000)
        )
        param <- result$par
        param["t"] <- 1/param["gamma"]
        param["tau"] <- param["sigma"] / sqrt(2*param["gamma"])
        param
    })
    # calulate rmse
    sapply(rownames(theta_hat), function (j) {
        sqrt(mean((theta_hat[j,]-theta[[j]])^2))/theta[[j]]
    })
}

test_cases <- data.frame(beta0=c(10,20,30), beta1=c(1,2,3))
result <- apply(test_cases, 1, function(beta) {
    sim_1(beta[["beta0"]], beta[["beta1"]], n_dataset = 2)
})
cbind(test_cases, t(result))