# setwd("code")
source("sim1_omega_tau_multiple.R")

all_cases <- read.csv("sim1_cases_new.csv")
# change the line below to your test cases
# my_ci <- 19:33
my_ci <- 31:33

sim_1(10,1,n_dataset = 5)
debug(sim_1)
test_cases <- data.frame(beta0=c(10, 20, 30), beta1=c(1,2,3))
result <- apply(test_cases, 1, function(beta) {
    sim_1(beta[["beta0"]], beta[["beta1"]], n_dataset = 10)
})
cbind(test_cases, t(result))