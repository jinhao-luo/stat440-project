setwd("code")
source("sim1_omega_tau_multiple.R")

all_cases <- read.csv("sim1_cases_new.csv")
# change the line below to your test cases
my_ci <- 16:18

my_cases <- all_cases[my_ci,]
my_cases
test_out_file <- paste0("sim1_case_new_", my_ci[1], "-", my_ci[length(my_ci)], ".RData")

system(paste("osascript -e 'display notification \"", test_out_file,"STARTED\" with title \"STAT440\" sound name \"default\"'"))

result <- apply(my_cases, 1, function(tc) {
    print(tc)
    sim_1(
        beta0=as.numeric(tc[["beta0"]]), beta1=as.numeric(tc[["beta1"]]), 
        omega=exp(-as.numeric(tc[["gamma"]])),
        n_dataset=as.numeric(tc[["n_dataset"]]), n_obs=as.numeric(tc[["n_obs"]]), 
        num_multistart=as.numeric(tc[["num_multistart"]]))
})
# if (file.exists(test_out_file)) {
#     tmp_out_file = tempfile()
#     warning(paste("ERROR IN SIMULATION 1:", test_out_file, "exists... Saving to", tmp_out_file, "Instead. Please manually retrive file if desired."))
#     save(my_cases, result, file=tmp_out_file)
# } else{
    save(my_cases, result, file=test_out_file)
# }

system(paste("osascript -e 'display notification \"", test_out_file,"DONE\" with title \"STAT440\" sound name \"default\"'"))
