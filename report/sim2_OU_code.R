require("smfret")

all_cases_default <- read.csv("sim2_OU_cases.csv", row.names=1)
test_cases_default <- data.frame(gamma=c(0.01, 0.05, 1, 2, 3, 4))
output_file_default <- "sim2_ou_results.RData"


#' Generate observations from a Morse SDE.
#'
#' @param all_cases dataframe of cases to test. Each row is a test test case.(see *Details*)
#' @param test_cases dataframe of sigmas to test
#' @param output_file the name of the output file
#' @return 
#' @details 
#' 1. parameter `all_cases` must has column 
#' 	`c(beta0","beta1", "tau", "mu")`.
#' 2. parameter `test_cases` must have column "gamma"
sim_2_OU_script <- function(all_cases=all_cases_default, test_cases=test_cases_default, 
        output_file=output_file_default) {
    for(i in 1:nrow(all_cases)) {
        row <- all_cases[i,]
        cur_beta0 <- row[['beta0']]
        cur_beta1 <- row[['beta1']]
        cur_tau <- row[['tau']]
        cur_mu <- row[['mu']]
        ou_result <- apply(test_cases, 1, function(info) {
            gamma <- info[["gamma"]]
            sigma <- sqrt(2*gamma)*cur_tau
            
            sim_2(from="ou", n_dataset = 20, n_obs = 400, 
                    beta0=cur_beta0, beta1=cur_beta1, gamma=gamma, mu=cur_mu, sigma=sigma)
        })
        capture.output(print(
            paste("Simulate from OU with beta0=", cur_beta0, "beta1=", cur_beta1, "tau=", 
            cur_tau, "mu=", cur_mu)), file=output_file, append=TRUE)
        capture.output(print(cbind(test_cases, t(ou_result))), file=output_file, append=TRUE)
    }
}

sim2_OU_script()
