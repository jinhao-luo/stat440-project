# research direction: whether it is possible to distinguish bound donor-acceptor pairs from free diffusion.

require(TMB)
source("smfret-functions.R")

# Compile and load the model.
gr_mod <- "main"
# compile(paste0(gr_mod, ".cpp"))
dyn.load(dynlib(gr_mod))

# distinguish
correct_num <- 0
bm_param_num <- 1
ou_param_num <- 3
test_num <- 10
NA_num <- 0
for (ii in 1:test_num) {
    print(paste("iteration", ii))
    # simulate from bm
    t_sigma <- rexp(1, 1/ii)
    t_mu <- rexp(1, 1/ii)
    beta0 <- 20
    beta1 <- 1
    dt <- sample(1:10, 1)
    n_obs <- sample(50:200, 1)
    ntheta <- 10 # number of parameter sets per test
    X <- bm_sim(t_mu, t_sigma, dt, n_obs, x0=t_mu)
    Y <- y_sim(X, beta0, beta1)

    # fit bm, OU
    gamma <- rexp(1, 1/ii)
    mu <- rexp(1, 1/ii)
    sigma <- rexp(1, 1/ii)
    bm_f <- MakeADFun(data=list(model_type="bm",dt=dt, Y=Y, beta0=beta0, beta1=beta1, niter=200),parameters=list(sigma=sigma))
    ou_f <- MakeADFun(data=list(model_type="ou",x0=x0, dt=dt, y=Y,beta0=beta0,beta1=beta1, niter=200),parameters=list(gamma=gamma, mu=mu, sigma=sigma))

    # calculate AIC, pick model
    bm_aic <- 2*bm_param_num+2*bm_f$fn()
    ou_aic <- 2*ou_param_num+2*ou_f$fn()
    if (is.na(bm_aic) || is.na(ou_aic)) {
        NA_num <- NA_num + 1
        next
    }
    if (bm_aic < ou_aic) {
        correct_num = correct_num+1
    }
}
print(paste("accuracy is:", correct_num/(test_num-NA_num)))




# simulate from OU
# fit bm, OU
# calculate AIC, pick model


# other model selection criterion
