require("TMB")
library(NMOF)
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
            rs_num <- 100
            omegas <- seq(0+0.01, 1-0.01, 0.01)

            test_function <- function(test_omega) {
                param <- list(omega= test_omega, mu = 0, tau= 1, X=rep(0, n_obs)) 
                data <- list(model_type = "omega_tau", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
                f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE)
                result <- optim(par = f$par, fn = f$fn, gr = f$gr, control=list(maxit=1000, reltol=1e-8), method=method)
                theta_hat <- result$par
                theta_hat["gamma"] <- -log(theta_hat["omega"])/dt
                theta_hat["t"] <- 1/theta_hat["gamma"]
                theta_hat["sigma"] <- theta_hat["tau"]*sqrt(2*theta_hat["gamma"])
                test_detail$theta_hat <- theta_hat
                return(f$fn(result$par))
            }


            sol <- gridSearch(fun=test_function, levels = list(omegas))
            print(sol$minfun)
            print(sol$minlevel)
            omega <- sol$minlevel
            param <- list(omega= omega, mu = 0, tau= 1, X=rep(0, n_obs)) 
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
    num_na <- NA
    # calulate rmse
    rmse <- sapply(rownames(theta_hat), function (j) {
        sqrt(mean((theta_hat[j,]-theta[[j]])^2, na.rm=TRUE))/theta[[j]]
    })
    sim_output <- list(rmse=signif(rmse,2), num_na=num_na, details=test_output, true_param=theta)
    sim_output
}


# debug(sim_1)
sim_out <-sim_1(10,0.5, omega=exp(-0.1),n_dataset = 5,n_obs=399)
sim_out$rmse
debug(sim_1)


# test_cases <- expand.grid(beta0=10, beta1=0.5, omega=c(exp(-0.1),exp(-1), exp(-10)), mu=c(1, 10), tau=c(10, 1,0.1))
# result <- apply(test_cases, 1, function(tc) {
#     sim_1(tc[["beta0"]], tc[["beta1"]], mu=tc[["mu"]], omega=tc[["omega"]], tau=tc[["tau"]], n_dataset = 10, n_obs=299)
# })
# t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))

# test_cases <- expand.grid(beta0=10, beta1=0.5, gamma=1, mu=c(1, 10), sigma=sqrt(2), dt=1)
# test_cases <- expand.grid(beta0=10, beta1=0.5, gamma=c(0.1,1,10), mu=c(1, 10), sigma=c(sqrt(2),sqrt(0.2),sqrt(20)), dt=1)
result <- apply(test_cases, 1, function(tc) {
    omega <- exp(-tc[["gamma"]]*tc[["dt"]])
    tau <- tc[["sigma"]] / sqrt(2*tc[["gamma"]])
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], mu=tc[["mu"]], omega=omega, tau=tau, n_dataset = 100, n_obs=199)
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))

test_cases <- expand.grid(n_obs=c(99,199,299,399,499,599))
result <- apply(test_cases, 1, function(tc) {
    sim_1(n_dataset=100, n_obs=tc[["n_obs"]])
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))

test_cases <- expand.grid(method=c("BFGS", "Nelder-Mead"))
result <- apply(test_cases, 1, function(tc) {
    sim_1(n_dataset=5, n_obs=99, method=tc[["method"]])
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))


# testing diff
# test_beta0 <- c(10, 20, 1)
test_beta0 <- c(10)
mu <- 10
test_cases <- do.call("rbind", (lapply(test_beta0, function(beta0) {
    data.frame(beta0=beta0, beta1=seq(max(0.1, (beta0-20)/mu), (beta0+1)/mu,length.out=10))
})))
system("osascript -e 'display notification \"beta sim started\" with title \"STAT440\" sound name \"default\"'")
result <- apply(test_cases, 1, function(tc) {
    print(tc)
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=100, n_obs=99)
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))
system("osascript -e 'display notification \"beta sim done\" with title \"STAT440\" sound name \"default\"'")


test_beta0 <- c(10)
mu <- 10
test_cases <- do.call("rbind", (lapply(test_beta0, function(beta0) {
    data.frame(beta0=beta0, beta1=seq(max(0.1, (beta0-20)/mu), (beta0+1)/mu,length.out=10))
})))
system("osascript -e 'display notification \"beta sim started\" with title \"STAT440\" sound name \"default\"'")
result <- apply(test_cases, 1, function(tc) {
    print(tc)
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=100, n_obs=99)
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))
system("osascript -e 'display notification \"beta sim done\" with title \"STAT440\" sound name \"default\"'")