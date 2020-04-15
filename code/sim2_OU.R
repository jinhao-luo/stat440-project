source("sim2.R")

cur_beta0 <- 10
cur_beta1 <- 0.5
# cur_omega <- exp(-1)
cur_tau <- 1
cur_mu <- 0

# simulate from OU
# test_cases <- data.frame(gamma=c(0.01, 0.1, 1, 10, 100, 1000, 10000))
test_cases <- data.frame(gamma=c(0.01, 0.05, 1, 2, 3, 4))
# test_cases <- data.frame(gamma=c(100))
# test_cases <- data.frame(gamma=c(1000))

ou_result <- apply(test_cases, 1, function(info) {
    gamma <- info[["gamma"]]
    sigma <- sqrt(2*gamma)*cur_tau
     
    sim_2(from="ou", n_dataset = 6, n_obs = 100, beta0=cur_beta0, beta1=cur_beta1, gamma=gamma, mu=cur_mu, sigma=sigma)
})
print(paste("Simulate from OU with beta0=", cur_beta0, "beta1=", cur_beta1, "tau=", cur_tau, "mu=", cur_mu))
cbind(test_cases, t(ou_result))