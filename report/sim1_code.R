require("smfret")

#' Run simulation 1 on all given cases and save RMSE in a csv
#'
#' @param all_cases dataframe of cases to teston. Each row is a test test case.(see **Details**)
#' @param my_ci the indices of cases to test on
#' @param out_file the name of the output file
#' @return 
#' @details 
#' 1. parameter `all_cases` must has column 
#' 	`c("beta0","beta1","gamma","n_dataset","n_obs","num_multistart")`.
#' 2. This function will not overwrite existing files. Instead, a warning will show 
#' 	and output will be saved in a temporary file.
#' 3. To save time, it is recommended that you run cases in different processes and
#' 	merge result afterwards.
#' @examples
#' sim_1_script(out_file="testing.csv", my_ci=1:2) 
sim_1_script <- function(all_cases=read.csv("sim1_cases_new.csv"), 
		my_ci=1:nrow(all_cases), out_file=paste0("sim1_result_",my_ci, ".csv")) {
	my_cases <- all_cases[my_ci,]
	result <- apply(my_cases, 1, function(tc) {
		print(tc)
		sim_1(
			beta0=as.numeric(tc[["beta0"]]), beta1=as.numeric(tc[["beta1"]]), 
			omega=exp(-as.numeric(tc[["gamma"]])),
			n_dataset=as.numeric(tc[["n_dataset"]]), n_obs=as.numeric(tc[["n_obs"]]), 
			num_multistart=as.numeric(tc[["num_multistart"]]))
	})
	rmse <- cbind(my_cases, t(sapply(1:nrow(my_cases), function(i) {
		c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)
		})))
	rmse <- apply(rmse, 2,  as.numeric)
	if (file.exists(out_file)) {
	    tmp_out_file = tempfile()
	    warning(paste("ERROR IN SIMULATION 1:", out_file, "exists... Saving to", 
			tmp_out_file, "Instead. Please manually retrive file if desired."))
	    write.csv(rmse, file=tmp_out_file, row.names=FALSE)
	} else{
	    write.csv(rmse, file=out_file, row.names=FALSE)
	}
}

hi <- sim_1_script(out_file="testing.csv", my_ci=1:2)
