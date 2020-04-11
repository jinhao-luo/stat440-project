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
    models <- replicate(n_dataset, expr = {
        # print("in replicate")
        Y <- c(NA)
        while (anyNA(Y)) {
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
        # print(Y)
        ou_f <- MakeADFun(
            data = list(model_type = "ou", niter = 100, dt = dt, y = Y, beta0 = beta0, beta1 = beta1),
            parameters = list(gamma = gamma, mu = mu, sigma = sigma)
        )
        bm_f <- MakeADFun(
            data=list(model_type="bm",dt=dt, Y=Y, beta0=beta0, beta1=beta1, niter=100),
            parameters=list(sigma=sigma))

        # calculate AIC, pick model
        bm_aic <- 2*bm_param_num+2*bm_f$fn()
        ou_aic <- 2*ou_param_num+2*ou_f$fn()

        print(bm_aic)
        print(ou_aic)
        print(2*bm_f$fn())
        print(2*ou_f$fn())
        if (bm_aic < ou_aic) {
            "bm"
        } else if (ou_aic < bm_aic) {
            "ou"
        } else {
            "bm/ou"
        }
    })
    models
}

test_cases <- data.frame(sigma=c(0.00000001, sqrt(2), sqrt(5)))
ou_result <- apply(test_cases, 1, function(info) {
    sim_2(from="ou", sigma=info[["sigma"]], n_dataset = 20)
})
print("Simulate from OU")
cbind(test_cases, t(ou_result))

test_cases <- data.frame(sigma=c(0.00000001, sqrt(2)))
bm_result <- apply(test_cases, 1, function(info) {
    sim_2(from="bm", sigma=info[["sigma"]], gamma = 0.0000000001, n_dataset = 3, n_obs = 200)
})
print("Simulate from BM")
cbind(test_cases, t(bm_result))

