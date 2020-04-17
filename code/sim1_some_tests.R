source("sim1_omega_tau_multiple.R")

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
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=100, n_obs=299, num_multistart=5)
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