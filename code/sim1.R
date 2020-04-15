require("TMB")
source("smfret-functions.R")
gr_mod <- "builtin"
# gr_mod <- "main"
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
sim_1 <- function(beta0=10, beta1=0.5, gamma = 1, mu = 10, sigma = sqrt(2*gamma), dt = 1,
                  n_obs = 99, n_dataset = 100, method="BFGS") {
    theta <- list(mu=mu, sigma=sigma, gamma=gamma, t = 1 / gamma, tau = sigma / sqrt(2 * gamma))
    test_output <- replicate(n_dataset, expr = {
        test_detail <- list()
        X <- ou_sim(gamma, mu, sigma, dt, n_obs)
        Y <- y_sim(X, beta0, beta1)
        test_detail$Y <- Y
        if (anyNA(Y)) {
            param_names <- c("sigma", "gamma", "mu", "t", "tau")
            param <- rep(NA, length(param_names))
            names(param) <- param_names
            test_detail$param <- param
        } else {
            param <- list(gamma = 1, mu = 0, sigma = 1, X=rep(0, n_obs))
            data <- list(model_type = "ou", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
            f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE)
            # param <- list(gamma = 1, mu = 0, sigma = 1)
            # data <- list(model_type = "ou", niter = 1000, dt = dt, beta0 = beta0, beta1 = beta1, y = Y)
            # f <- MakeADFun(data = data, parameters = param)
            result <- optim(par = f$par, fn = f$fn, gr = f$gr, control=list(maxit=1000,reltol=1e-8))
            param <- result$par
            param["t"] <- 1/param["gamma"]
            param["tau"] <- param["sigma"] / sqrt(2*param["gamma"])
            test_detail$theta_hat <- param
        }
        test_detail
    })

    theta_hat <- apply(test_output, 2, function(tc) {tc$theta_hat})
    # get number of NAs in simulation
    num_na <- sum(is.na(theta_hat["t",]))
    # calulate rmse
    rmse <- sapply(rownames(theta_hat), function (j) {
        sqrt(mean((theta_hat[j,]-theta[[j]])^2, na.rm=TRUE))/theta[[j]]
    })
    sim_output <- list(rmse=signif(rmse,2), num_na=num_na, details=test_output, true_param=theta)
    sim_output
}

sim_1()

test_cases <- expand.grid(beta0=10, beta1=0.5, gamma=c(0.1,1,10), mu=c(1, 10))
result <- apply(test_cases, 1, function(tc) {
    sim_1(tc[["beta0"]], tc[["beta1"]], mu=tc[["mu"]], gamma=tc[["gamma"]], n_dataset = 100, n_obs=99)
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))
