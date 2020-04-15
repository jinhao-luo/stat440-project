library(NMOF)

find_optim_omega <- function(omegas, n_obs, dt, Y, beta0, beta1, method="Nelder-Mead") {
    test_function <- function(test_omega) {
        param <- list(omega= test_omega, mu = 0, tau= 1, X=rep(0, n_obs)) 
        data <- list(model_type = "omega_tau", dt = dt, y = Y, beta0 = beta0, beta1 = beta1)
        f <- MakeADFun(data = data, parameters = param, random = c("X"), silent = TRUE, method=method)
        # result <- optim(par = f$par, fn = f$fn, gr = f$gr, control=list(trace=5, maxit=1000, reltol=1e-8), method="BFGS")
        result <- optim(par = f$par, fn = f$fn, gr = f$gr, control=list(trace=0, maxit=1000, reltol=1e-8))
        theta_hat <- result$par
        theta_hat["gamma"] <- -log(theta_hat["omega"])/dt
        theta_hat["t"] <- 1/theta_hat["gamma"]
        theta_hat["sigma"] <- theta_hat["tau"]*sqrt(2*theta_hat["gamma"])
        # test_detail$theta_hat <- theta_hat
        return(f$fn(result$par))
    }
    sol <- gridSearch(fun=test_function, levels = list(omegas))
    print(sol$minfun)
    print(sol$minlevel)
    omega <- sol$minlevel
    return(omega)
}