require("TMB")
require(NMOF)
source("R/smfret-functions.R")
gr_mod <- "src/TMB/smfret_TMBExports"
dyn.load(dynlib(gr_mod))



#' Simulation 1
#'
#' @param beta0 The value of beta0.
#' @param beta1 The value of beta1.
#' @param omega parameter in Ornstein-Uhlenbeck process.
#' @param mu Scalar mean parameter.
#' @param tau Variance of donor acceptor distance.
#' @param dt Interobservation time.
#' @param n_obs Number of observations to generate.
#' @param n_dataset Number of dataset to simulate for each `beta`.
#' @param num_multistart number of multistart to perform on top of optim
#' @param method passed to `optim()`
#' @return A vector of rmse: rmse for inference on each simulated dataset.
#' @export
#'
sim_1 <- function(beta0=10, beta1=0.5, omega= exp(-1), mu = 10, tau= 1, dt = 1,
                  n_obs = 99, n_dataset = 100, method="BFGS", num_multistart=5) {
    gamma <- -log(omega)/dt
    t <- 1 / gamma
    sigma <- tau*sqrt(2*gamma)
    theta <- list(omega=omega, mu=mu, tau=tau, t=t, gamma= gamma, sigma=sigma)
    test_output <- replicate(n_dataset, expr = {
        test_detail <- list()
        X <- ou_sim(gamma, mu, sigma, dt, n_obs)
        Y <- y_sim(X, beta0, beta1)
        test_detail$Y <- Y
        test_detail$theta_hat <- c(omega=NA, mu=NA, tau=NA, gamma=NA, t=NA, sigma=NA)
        if (!anyNA(Y)) {
            if(num_multistart > 1) {
                omegas <- seq(0+0.01, 1-0.01, length.out=num_multistart)
                test_function <- function(test_omega) {
                    param <- list(omega= test_omega, mu = 0, tau= 1, X=rep(0, n_obs))
                    data <- list(model_type = "omega_tau", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
                    f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE)
                    return(tryCatch({
                        # result <- optim(par = f$par, fn = f$fn, gr = f$gr, control=list(maxit=1000, reltol=1e-8), method=method)
                        result <- constrOptim(f$par, f$fn, f$gr, control=list(maxit=1000, reltol=1e-8), method=method,
                            ui = rbind(c(1,0,0), c(-1,0,0), c(0,0,1)), ci=c(0,-1, 0))
                        theta_hat <- result$par
                        theta_hat["gamma"] <- -log(theta_hat["omega"])/dt
                        theta_hat["t"] <- 1/theta_hat["gamma"]
                        theta_hat["sigma"] <- theta_hat["tau"]*sqrt(2*theta_hat["gamma"])
                        test_detail$theta_hat <- theta_hat
                        f$fn(result$par)
                    }, error= function(cond) {Inf},
                    warning=function(cond) {Inf})) # avoid optim function cannot be evaluated at inital param error TODO: fix this
                    # TODO: maybe remove tryCatch
                }

                sol <- gridSearch(fun=test_function, levels = list(omegas), printDetail=FALSE)
                omega <- sol$minlevel[1] # minlevel could return multiple values if they have same value
            } else {
                omega<- 0.5
            }
            param <- list(omega=omega, mu = 0, tau= 1, X=rep(0, n_obs))
            data <- list(model_type = "omega_tau", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
            f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE)
            # result <- optim(par = f$par, fn = f$fn, gr = f$gr, control=list(maxit=1000, reltol=1e-8), method=method)
            result <- constrOptim(f$par, f$fn, f$gr, control=list(maxit=1000, reltol=1e-8), method=method,
                ui = rbind(c(1,0,0), c(-1,0,0), c(0,0,1)), ci=c(0,-1, 0))
            theta_hat <- result$par
            theta_hat["gamma"] <- log(theta_hat["omega"])/-dt
            theta_hat["t"] <- 1/theta_hat["gamma"]
            theta_hat["sigma"] <- theta_hat["tau"]*sqrt(2*theta_hat["gamma"])
            test_detail$theta_hat <- theta_hat
            }
        test_detail
    })

    theta_hat <- apply(test_output, 2, function(tc) {tc$theta_hat})
    print(theta_hat)
    # calulate rmse
    rmse <- sapply(rownames(theta_hat), function (j) {
        sqrt(mean((theta_hat[j,]-theta[[j]])^2, na.rm=TRUE))/theta[[j]]
    })
    sim_output <- list(rmse=signif(rmse,2), details=test_output, true_param=theta)
    sim_output
}
