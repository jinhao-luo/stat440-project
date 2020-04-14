source("sim2.R")

cur_beta0 <- 5
cur_beta1 <- 1
cur_mu <- -10
cur_gamma <- 100
cur_sigma <- 100

# simulate from BM
test_cases <- data.frame(sigma=c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1))
# test_cases <- data.frame(sigma=c(0.5))
bm_result <- apply(test_cases, 1, function(info) {
    sim_2(from="bm", sigma=info[["sigma"]], n_dataset = 6, n_obs = 100, beta0=cur_beta0, beta1=cur_beta1, mu=cur_mu, gamma=cur_gamma)
})
print(paste("Simulate from BM with beta0=", cur_beta0, "beta1=", cur_beta1, "mu=", cur_mu, "gamma=", cur_gamma))
cbind(test_cases, t(bm_result))