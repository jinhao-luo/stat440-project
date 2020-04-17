source("sim2.R")

all_cases <- read.csv("sim2_BM_cases.csv", row.names=1)
test_cases <- data.frame(sigma=c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1))

output_file <- "sim2_bm_results.RData"

for(i in 1:nrow(all_cases)) {
    row <- all_cases[i,]
    cur_beta0 <- row[['beta0']]
    cur_beta1 <- row[['beta1']]
    cur_mu <- -row[['mu']]
    cur_gamma <- row[['gamma']]
    cur_method <- toString(row[['method']])
    print(paste("cur_method=", cur_method))

    bm_result <- apply(test_cases, 1, function(info) {
        sim_2(from="bm", sigma=info[["sigma"]], n_dataset = 20, n_obs = 400, beta0=cur_beta0, beta1=cur_beta1, mu=cur_mu, gamma=cur_gamma, method=cur_method)
    })
    capture.output(print(paste("Simulate from BM with beta0=", cur_beta0, "beta1=", cur_beta1, "mu=", cur_mu, "gamma=", cur_gamma, "with method=", cur_method)), file=output_file, append=TRUE)
    capture.output(cbind(test_cases, t(bm_result)), file=output_file, append=TRUE)
}
