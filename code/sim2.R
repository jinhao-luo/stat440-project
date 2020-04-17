# research direction: whether it is possible to distinguish bound donor-acceptor pairs from free diffusion.

require(TMB)
source("smfret-functions.R")
source("multi_start.R")

# Compile and load the model.
gr_mod <- "builtin"
# compile(paste0(gr_mod, ".cpp"))
dyn.load(dynlib(gr_mod))

bm_param_num <- 1
ou_param_num <- 3

#' simulation 2
#'
#' @param from simulate X from which model ('ou' or 'bm')
#' @param beta0 the value of beta0
#' @param beta1 the value of beta1
#' @param gamma Scalar mean reversion parameter
#' @param mu Scalar mean parameter
#' @param sigma Scalar diffusion parameter
#' @param dt Interobservation time.
#' @param n_obs Number of observations to generate.
#' @param n_dataset Number of dataset to simulate for each `\beta`
#' @return a vector of picked model
sim_2 <- function(from, beta0 = 20, beta1 = .5, gamma = 1, mu = 10, sigma, dt = 1,
                  n_obs = 99, n_dataset = 100, multi_start = TRUE, method="Nelder-Mead") {
    theta <- list(mu=mu, sigma=sigma, gamma=gamma, t = 1 / gamma, tau = sigma / sqrt(2 * gamma))
    print("reach0")
    models <- replicate(n_dataset, expr = {
        ou_aic <- NA
        bm_aic <- NA
        while (is.na(ou_aic) || is.na(bm_aic)) {
            Y <- c(NA)
            while (anyNA(Y)) {
                print("NA")
                X <- 0
                if (from == "ou") {
                    X <- ou_sim(gamma, mu, sigma, dt, n_obs)
                } else if (from == "bm") {
                    X <- bm_sim(mu, sigma, dt, n_obs, x0 = mu)
                } else {
                    stop("invalid from")
                }
                Y <- y_sim(X, beta0, beta1)
            }

            omega <- 0
            if (multi_start) {
                num_start <- 5
                start <- 0+0.01
                end <- 1-0.01
                omegas <- seq(start, end, (end-start)/num_start)
                omega <- find_optim_omega(omegas, n_obs, dt, Y, beta0, beta1)
            }
            param <- list(omega= omega, mu = 0, tau= 1, X=rep(0, n_obs)) 
            data <- list(model_type = "omega_tau", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
            print("reach6")
            ou_f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE, method=method)
            print("reach7")
            ou_result <- optim(par = ou_f$par, fn = ou_f$fn, gr = ou_f$gr, control=list(trace=5, maxit=1000, reltol=1e-8), method=method)
            print("reach8")

            param <- list(sigma = 1, X=rep(0, n_obs))
            data <- list(model_type = "bm", dt=dt, Y=Y, beta0 = beta0, beta1 = beta1)
            bm_f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE, method=method)
            bm_result <- optim(par = bm_f$par, fn = bm_f$fn, gr = bm_f$gr, control=list(trace=5, maxit=1000, reltol=1e-8), method=method)
            
            print("reach9")
            # calculate AIC, pick model
            ou_aic <- 2*ou_param_num+2*ou_f$fn(par=ou_result$par)
            bm_aic <- 2*bm_param_num+2*bm_f$fn(par= bm_result$par)
        }

        if (bm_aic < ou_aic) {
            print("reach4")
            "bm"
        } else if (ou_aic < bm_aic) {
            print("reach5")
            "ou"
        } else {
            "bm/ou"
        }
    })
    models
}