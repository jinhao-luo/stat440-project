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
            omegas <- seq(0+0.01, 1-0.01, length.out=num_multistart) # TODO: might need to change this back?

            test_function <- function(test_omega) {
                param <- list(omega= test_omega, mu = 0, tau= 1, X=rep(0, n_obs)) 
                data <- list(model_type = "omega_tau", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
                f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE)
                return(tryCatch({
                    result <- optim(par = f$par, fn = f$fn, gr = f$gr, control=list(maxit=1000, reltol=1e-8), method=method)
                    theta_hat <- result$par
                    theta_hat["gamma"] <- -log(theta_hat["omega"])/dt
                    theta_hat["t"] <- 1/theta_hat["gamma"]
                    theta_hat["sigma"] <- theta_hat["tau"]*sqrt(2*theta_hat["gamma"])
                    test_detail$theta_hat <- theta_hat
                    f$fn(result$par)
                }, error= function(cond) {Inf}, 
                warning=function(cond) {Inf})) # avoid optim function cannot be evaluated at inital param error TODO: fix this
            }


            sol <- gridSearch(fun=test_function, levels = list(omegas)) # TODO: silent this (printDetail = FALSE)
            print(sol$minlevel)
            omega <- sol$minlevel[1] # minlevel could return multiple values if they have same value
            param <- list(omega=omega, mu = 0, tau= 1, X=rep(0, n_obs)) 
            data <- list(model_type = "omega_tau", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
            f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE)
            if (FALSE) {
            # if (sol$value[1] == Inf) {
                # given method failed to give valid omega \in (0,1)
                test_detail$theta_hat["omega"] <- 0
                warning("Some optimization failed")
            } else {
                print(f$par)
                result <- optim(par = f$par, fn = f$fn, gr = f$gr, control=list(maxit=1000, reltol=1e-8), method=method)
                theta_hat <- result$par
                theta_hat["gamma"] <- log(theta_hat["omega"])/-dt
                theta_hat["t"] <- 1/theta_hat["gamma"]
                theta_hat["sigma"] <- theta_hat["tau"]*sqrt(2*theta_hat["gamma"])
                test_detail$theta_hat <- theta_hat
            }
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


# debug(sim_1)
# sim_out <-sim_1(10,0.5, omega=exp(-0.1))
sim_out0 <-sim_1(10,0.5, num_multistart = 5)
sim_out2 <-sim_1(10,0.5, n_obs=299, num_multistart = 5)
sim_out0$rmse
sim_out2$rmse
system("osascript -e 'display notification \"sim1 done\" with title \"STAT440\" sound name \"default\"'")
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
    sim_1(n_dataset=10, n_obs=tc[["n_obs"]], num_multistart=10, method="BFGS")
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))
system("osascript -e 'display notification \"n_obs done\" with title \"STAT440\" sound name \"default\"'")

test_cases <- expand.grid(method=c("BFGS", "Nelder-Mead"))
result <- apply(test_cases, 1, function(tc) {
    sim_1(n_dataset=5, n_obs=99, method=tc[["method"]])
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))


# testing diff
test_beta0 <- 10
mu <- 10
test_cases <- do.call("rbind", (lapply(test_beta0, function(beta0) {
    data.frame(beta0=beta0, beta1=seq(max(0.1, (beta0-20)/mu), (beta0+1)/mu,length.out=5))
})))
system("osascript -e 'display notification \"beta0=10 sim started\" with title \"STAT440\" sound name \"default\"'")
result <- apply(test_cases, 1, function(tc) {
    print(tc)
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=100, n_obs=99, num_multistart=10)
})
cbind(test_cases, t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))
system("osascript -e 'display notification \"beta0=10 sim done\" with title \"STAT440\" sound name \"default\"'")


test_beta0 <- 5/2^seq(0,2)
mu <- 10
test_cases <- do.call("rbind", (lapply(test_beta0, function(beta0) {
    data.frame(beta0=beta0, beta1=seq(max(0.1, (beta0-20)/mu), (beta0+1)/mu,length.out=5))
})))
system("osascript -e 'display notification \"beta0=5/2^seq(0,2) sim started\" with title \"STAT440\" sound name \"default\"'")
result <- apply(test_cases, 1, function(tc) {
    print(tc)
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=100, n_obs=99, num_multistart=10)
})
cbind(test_cases, t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))
system("osascript -e 'display notification \"beta0=5/2^seq(0,2) sim done\" with title \"STAT440\" sound name \"default\"'")

test_beta0 <- 5/2^seq(3,5)
mu <- 10
test_cases <- do.call("rbind", (lapply(test_beta0, function(beta0) {
    data.frame(beta0=beta0, beta1=seq(max(0.1, (beta0-20)/mu), (beta0+1)/mu,length.out=5))
})))
system("osascript -e 'display notification \"beta0=5/2^seq(3,5) sim started\" with title \"STAT440\" sound name \"default\"'")
result <- apply(test_cases, 1, function(tc) {
    print(tc)
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=100, n_obs=99, num_multistart=10)
})
cbind(test_cases, t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))
system("osascript -e 'display notification \"beta0=5/2^seq(3,5) sim done\" with title \"STAT440\" sound name \"default\"'")

test_beta0 <- 5/2^seq(6,8)
mu <- 10
test_cases <- do.call("rbind", (lapply(test_beta0, function(beta0) {
    data.frame(beta0=beta0, beta1=seq(max(0.1, (beta0-20)/mu), (beta0+1)/mu,length.out=5))
})))
system("osascript -e 'display notification \"beta0=5/2^seq(6,8) sim started\" with title \"STAT440\" sound name \"default\"'")
result <- apply(test_cases, 1, function(tc) {
    print(tc)
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=1, n_obs=99, num_multistart=1)
})
cbind(test_cases, t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))
system("osascript -e 'display notification \"beta0=5/2^seq(6,8) sim done\" with title \"STAT440\" sound name \"default\"'")



test_cases <- expand.grid(n_obs=c(99,199,299,399,499,599), num_multistart=c(1,10), method=c("BFGS", "Nelder-Mead"))
result <- apply(test_cases, 1, function(tc) {
    sim_1(n_dataset=10, n_obs=as.numeric(tc[["n_obs"]]), num_multistart=as.numeric(tc[["num_multistart"]]), method=tc[["method"]])
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))
system("osascript -e 'display notification \"n_obs done\" with title \"STAT440\" sound name \"default\"'")
