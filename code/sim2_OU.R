source("sim2.R")

cur_beta0 <- 10
cur_beta1 <- 0.5
cur_mu <- 10
cur_gamma <- 1
cur_sigma <- sqrt(2)

# simulate from OU
test_cases <- data.frame(gamma=c(0.01, 0.1, 1, 10, 100, 1000, 10000))
# test_cases <- data.frame(gamma=c(100))
# test_cases <- data.frame(gamma=c(1000))

ou_result <- apply(test_cases, 1, function(info) {
    sim_2(from="ou", sigma=cur_sigma, gamma=info[["gamma"]], n_dataset = 6, n_obs = 100, beta0=cur_beta0, beta1=cur_beta1, mu=cur_mu)
})
print(paste("Simulate from OU with beta0=", cur_beta0, "beta1=", cur_beta1, "mu=", cur_mu, "sigma=", cur_sigma))
cbind(test_cases, t(ou_result))