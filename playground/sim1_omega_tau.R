require("TMB")
source("smfret-functions.R")
gr_mod <- "builtin"
# gr_mod <- "main"
dyn.load(dynlib(gr_mod))


#' simulation 1
#'
#' @param beta0 the value of beta0
#' @param beta1 the value of beta1
#' @param omega 
#' @param mu Scalar mean parameter
#' @param tau
#' @param dt Interobservation time.
#' @param n_obs Number of observations to generate.
#' @param n_dataset Number of dataset to simulate for each `\beta`
#' @return a vector of rmse: rmse for inference on each simulated dataset
sim_1 <- function(beta0=10, beta1=0.5, omega= exp(-1), mu = 10, tau= 1, dt = 1,
                  n_obs = 99, n_dataset = 100, method="BFGS") {
    gamma <- -log(omega)/dt
    t <- 1 / gamma
    sigma <- tau*sqrt(2*gamma)
    theta <- list(omega=omega, mu=mu, tau=tau, t=t, gamma= gamma, sigma=sigma)
    test_output <- replicate(n_dataset, expr = {
        test_detail <- list()
        X <- ou_sim(gamma, mu, sigma, dt, n_obs)
        Y <- y_sim(X, beta0, beta1)
        test_detail$Y <- Y
        if (anyNA(Y)) {
            param_names <- c("omega", "mu", "tau", "gamma", "t", "sigma")
            param <- rep(NA, length(param_names))
            names(param) <- param_names
            test_detail$theta_hat <- param
        } else {
            # param <- list(omega= runif(1), mu = -10, tau= 1, X=rep(0, n_obs))
            param <- list(omega= 0.5, mu = 0, tau= 1, X=rep(0, n_obs)) 
            data <- list(model_type = "omega_tau", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
            f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE)
            result <- optim(par = f$par, fn = f$fn, gr = f$gr, control=list(maxit=1000, reltol=1e-8), method=method)
            theta_hat <- result$par
            theta_hat["gamma"] <- -log(theta_hat["omega"])/dt
            theta_hat["t"] <- 1/theta_hat["gamma"]
            theta_hat["sigma"] <- theta_hat["tau"]*sqrt(2*theta_hat["gamma"])
            test_detail$theta_hat <- theta_hat
        }
        test_detail
    })

    theta_hat <- apply(test_output, 2, function(tc) {tc$theta_hat})
    print(theta_hat)
    # get number of NAs in simulation
    num_na <- sum(is.na(theta_hat))/length(theta)
    # calulate rmse
    rmse <- sapply(rownames(theta_hat), function (j) {
        sqrt(mean((theta_hat[j,]-theta[[j]])^2, na.rm=TRUE))/theta[[j]]
    })
    sim_output <- list(rmse=signif(rmse,2), num_na=num_na, details=test_output, true_param=theta)
    sim_output
}


# debug(sim_1)
sim_out0 <-sim_1(10,0.5,n_dataset = 100,n_obs=99)
sim_out1 <-sim_1(10,0.5,n_dataset = 100,n_obs=199)
sim_out2 <-sim_1(10,0.5,n_dataset = 100,n_obs=299)
sim_out3 <-sim_1(10,0.5,n_dataset = 100,n_obs=399)
sim_out4 <-sim_1(10,0.5,n_dataset = 100,n_obs=499)
sim_out0$rmse
sim_out1$rmse
sim_out2$rmse
sim_out3$rmse
sim_out4$rmse
debug(sim_1)


# test_cases <- expand.grid(beta0=10, beta1=0.5, omega=c(exp(-0.1),exp(-1), exp(-10)), mu=c(1, 10), tau=c(10, 1,0.1))
# result <- apply(test_cases, 1, function(tc) {
#     sim_1(tc[["beta0"]], tc[["beta1"]], mu=tc[["mu"]], omega=tc[["omega"]], tau=tc[["tau"]], n_dataset = 10, n_obs=299)
# })
# t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))

test_cases <- expand.grid(beta0=10, beta1=0.5, gamma=c(0.1,1,10), mu=c(1, 10), sigma=c(sqrt(2),sqrt(0.2),sqrt(20)), dt=1)
test_cases <- expand.grid(beta0=10, beta1=0.5, gamma=1, mu=10, sigma=sqrt(2), dt=1)
result <- apply(test_cases, 1, function(tc) {
    omega <- exp(-tc[["gamma"]]*tc[["dt"]])
    tau <- tc[["sigma"]] / sqrt(2*tc[["gamma"]])
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], mu=tc[["mu"]], omega=omega, tau=tau, n_dataset = 1000, n_obs=999)
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))
