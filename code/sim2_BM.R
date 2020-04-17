source("sim2.R")

all_cases_default <- read.csv("sim2_BM_cases.csv", row.names=1)
test_cases_default <- data.frame(sigma=c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1))

output_file_default <- "sim2_bm_results.RData"

#' Generate observations from a Morse SDE.
#'
#' @param all_cases dataframe of cases to test. Each row is a test test case.(see *Details*)
#' @param test_cases dataframe of sigmas to test
#' @param output_file the name of the output file
#' @return 
#' @details 
#' 1. parameter `all_cases` must has column 
#' 	`c("beta0","beta1","gamma","n_dataset","n_obs","num_multistart")`.
#' 2. parameter `test_cases` must have column "sigma"
sim_2_BM_script <- function(all_cases=all_cases_default, test_cases=test_cases_default, output_file=output_file_default) {
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
}

sim_2_BM_script()