# research direction: whether it is possible to distinguish bound donor-acceptor pairs from free diffusion.

require(TMB)
source("smfret-functions.R")

# Compile and load the model.
gr_mod <- "main"
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
                  n_obs = 99, n_dataset = 100) {
    theta <- list(mu=mu, sigma=sigma, gamma=gamma, t = 1 / gamma, tau = sigma / sqrt(2 * gamma))
    print("reach0")
    models <- replicate(n_dataset, expr = {
        bm_aic <- NA
        ou_aic <- NA
        while (is.na(bm_aic) || is.na(ou_aic)) {
            print("in replicate")
            Y <- c(NA)
            while (anyNA(Y)) {
                print("in NA")
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
            print("reach1")
            ou_f <- MakeADFun(
                data = list(model_type = "ou", niter = 100, dt = dt, y = Y, beta0 = beta0, beta1 = beta1),
                parameters = list(gamma = gamma, mu = mu, sigma = sigma)
            )
            print("reach2")
            print(ou_f$par)
            ou_result <- optim(par = ou_f$par, fn = ou_f$fn, gr = ou_f$gr,
                control = list(maxit = 1000)
            )
            bm_f <- MakeADFun(
                data=list(model_type="bm",dt=dt, Y=Y, beta0=beta0, beta1=beta1, niter=100),
                parameters=list(sigma=sigma))
            bm_result <- optim(par = bm_f$par, fn = bm_f$fn, gr = bm_f$gr,
                control = list(maxit = 1000)
            )
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

test_cases <- data.frame(sigma=c(0.001, 0.01, 0.1, 1, 1.1, 1.2, 1.3, 1.4, 1.5))
cur_mu <- 1
cur_gamma <- 0.1
ou_result <- apply(test_cases, 1, function(info) {
    sim_2(from="ou", sigma=info[["sigma"]], n_dataset = 4, beta0=10, beta1=.5, mu=cur_mu, gamma=cur_gamma)
})
print(paste("Simulate from OU with mu=", mu, "gamma=", cur_gamma)
cbind(test_cases, t(ou_result))

test_cases <- data.frame(sigma=c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1))
cur_mu <- 1
cur_gamma <- 0.1
bm_result <- apply(test_cases, 1, function(info) {
    sim_2(from="bm", sigma=info[["sigma"]], n_dataset = 2, n_obs = 200, beta0=10, beta1=0.5, mu=cur_mu, gamma=cur_gamma)
})
print(paste("Simulate from OU with mu=", mu, "gamma=", cur_gamma)
cbind(test_cases, t(bm_result))


