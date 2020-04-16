vec_ci <- c("1-3") # TODO: add the remaining

# get data from RData files
all_test_cases <- data.frame()
vec_result <- c()
all_rmse <- data.frame()
for (ci in vec_ci) {
	load(paste0("sim1_case_new_",ci,".RData")) # we get `result` and `my_cases`
	vec_result <- c(vec_result, result)
	all_test_cases <- rbind(all_test_cases, my_cases)
	rmse <- cbind(my_cases, t(sapply(1:nrow(my_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))
	all_rmse <- rbind(all_rmse, rmse)
}
names(vec_result) <- vec_ci
all_rmse