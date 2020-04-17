source("sim2.R")

all_cases <- read.csv("sim2_OU_cases.csv", row.names=1)
test_cases <- data.frame(gamma=c(0.01, 0.05, 1, 2, 3, 4))

output_file <- "sim2_ou_results.RData"

for(i in 1:nrow(all_cases)) {
    row <- all_cases[i,]
    cur_beta0 <- row[['beta0']]
    cur_beta1 <- row[['beta1']]
    cur_tau <- row[['tau']]
    cur_mu <- row[['mu']]

    ou_result <- apply(test_cases, 1, function(info) {
        gamma <- info[["gamma"]]
        sigma <- sqrt(2*gamma)*cur_tau
        
        sim_2(from="ou", n_dataset = 20, n_obs = 400, beta0=cur_beta0, beta1=cur_beta1, gamma=gamma, mu=cur_mu, sigma=sigma)
    })
    capture.output(print(paste("Simulate from OU with beta0=", cur_beta0, "beta1=", cur_beta1, "tau=", cur_tau, "mu=", cur_mu)), file=output_file, append=TRUE)
    capture.output(print(cbind(test_cases, t(ou_result))), file=output_file, append=TRUE)
}