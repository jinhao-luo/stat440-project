# vec_ci <- c("1-5","6-10", "11-15", "16-20","21-25", "26-30") # TODO: add the remaining
# vec_ci <- c("1-3","4-6","7-9", "10-12", "13-15", "16-18")
# vec_ci <- c("19-21", "22-24", "25-27", "28-30", "31-33")
vec_ci <- c("30-31")

# get data from RData files
all_test_cases <- data.frame()
vec_result <- c()
all_rmse <- data.frame()
for (ci in vec_ci) {
	load(paste0("sim1_case_new_",ci,".RData")) # we get `result` and `my_cases`
	vec_result <- c(vec_result, result)
	all_test_cases <- rbind(all_test_cases, my_cases)
	rmse <- cbind(my_cases, t(sapply(1:nrow(my_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))
	print(rmse)
	all_rmse <- rbind(all_rmse, rmse)
}
names(vec_result) <- vec_ci
# all_rmse <- subset(all_rmse, select=-c(X))
all_rmse
# write.csv(apply(all_rmse, 2,  as.numeric), "ohno.csv", row.names=FALSE)
