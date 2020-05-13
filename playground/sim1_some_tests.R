source("sim1_omega_tau_multiple.R")

# debug(sim_1)
# sim_out <-sim_1(10,0.5, omega=exp(-0.1))
sim_out0 <-sim_1(10,0.5, num_multistart = 5)
sim_out2 <-sim_1(10,0.5, n_obs=299, num_multistart = 5)
sim_out0$rmse
sim_out2$rmse
system("osascript -e 'display notification \"sim1 done\" with title \"STAT440\" sound name \"default\"'")
debug(sim_1)


# test_cases <- expand.grid(beta0=10, beta1=0.5, omega=c(exp(-0.1),exp(-1), exp(-10)), mu=c(1, 10), tau=c(10, 1,0.1))
# result <- apply(test_cases, 1, function(tc) {
#     sim_1(tc[["beta0"]], tc[["beta1"]], mu=tc[["mu"]], omega=tc[["omega"]], tau=tc[["tau"]], n_dataset = 10, n_obs=299)
# })
# t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))

# test_cases <- expand.grid(beta0=10, beta1=0.5, gamma=1, mu=c(1, 10), sigma=sqrt(2), dt=1)
# test_cases <- expand.grid(beta0=10, beta1=0.5, gamma=c(0.1,1,10), mu=c(1, 10), sigma=c(sqrt(2),sqrt(0.2),sqrt(20)), dt=1)
result <- apply(test_cases, 1, function(tc) {
    omega <- exp(-tc[["gamma"]]*tc[["dt"]])
    tau <- tc[["sigma"]] / sqrt(2*tc[["gamma"]])
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], mu=tc[["mu"]], omega=omega, tau=tau, n_dataset = 100, n_obs=199)
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))

test_cases <- expand.grid(n_obs=c(99,199,299,399,499,599))
result <- apply(test_cases, 1, function(tc) {
    sim_1(n_dataset=10, n_obs=tc[["n_obs"]], num_multistart=10, method="BFGS")
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))
system("osascript -e 'display notification \"n_obs done\" with title \"STAT440\" sound name \"default\"'")

test_cases <- expand.grid(method=c("BFGS", "Nelder-Mead"))
result <- apply(test_cases, 1, function(tc) {
    sim_1(n_dataset=5, n_obs=99, method=tc[["method"]])
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))

new_cases <- rbind(
    c(beta0=15, beta1=1, gamma=1, n_dataset=100, n_obs=299, num_multistart=1),
    c(beta0=15, beta1=1, gamma=1, n_dataset=100, n_obs=299, num_multistart=10),
    c(beta0=11, beta1=1, gamma=1, n_dataset=100, n_obs=299, num_multistart=10),
    c(beta0=9, beta1=1, gamma=1, n_dataset=100, n_obs=299, num_multistart=10),
    c(beta0=15, beta1=0.8, gamma=1, n_dataset=100, n_obs=299, num_multistart=10),
    c(beta0=15, beta1=0.6, gamma=1, n_dataset=100, n_obs=299, num_multistart=10),
    c(beta0=15, beta1=1, gamma=0.1, n_dataset=100, n_obs=299, num_multistart=10),
    c(beta0=15, beta1=1, gamma=0.1, n_dataset=100, n_obs=999, num_multistart=10),
    c(beta0=11, beta1=1, gamma=0.1, n_dataset=100, n_obs=999, num_multistart=10),
    c(beta0=9, beta1=1, gamma=0.1, n_dataset=100, n_obs=999, num_multistart=10),
    c(beta0=15, beta1=0.8, gamma=0.1, n_dataset=100, n_obs=999, num_multistart=10),
    c(beta0=15, beta1=0.6, gamma=0.1, n_dataset=100, n_obs=999, num_multistart=10),
    c(beta0=15, beta1=1, gamma=10, n_dataset=100, n_obs=299, num_multistart=10),
    c(beta0=15, beta1=1, gamma=10, n_dataset=100, n_obs=999, num_multistart=10),
    c(beta0=11, beta1=1, gamma=10, n_dataset=100, n_obs=999, num_multistart=10),
    c(beta0=9, beta1=1, gamma=10, n_dataset=100, n_obs=999, num_multistart=10),
    c(beta0=15, beta1=0.8, gamma=10, n_dataset=100, n_obs=999, num_multistart=10),
    c(beta0=15, beta1=0.6, gamma=10, n_dataset=100, n_obs=999, num_multistart=10)
)
# write.csv(new_cases, "sim1_cases_new.csv", row.names=FALSE)

setwd("code")
source("smfret-functions.R")
tests <- rbind(
    # data.frame(gamma=1, beta0=c(10,15,20), beta1=1),
    # data.frame(gamma=0.1, beta0=5+seq(10,100,by=10), beta1=1:10)
    # data.frame(gamma=0.1, beta0=c(15,16,17,18,19), beta1=1),
    # data.frame(gamma=0.1, beta0=c(23,24,25,26,27), beta1=2)
    # data.frame(gamma=10, beta0=15, beta1=c(0.1,0.2))
    # data.frame(gamma=0.1, beta0=c(10), beta1=0.1)
    # data.frame(gamma=c(0.1,10), beta0=15, beta1=1, n_obs=299),
    # data.frame(gamma=c(0.1,10), beta0=15, beta1=1, n_obs=999)
    data.frame(gamma=1, beta0=15, beta1=c(1.1, 1,0.9,0.8,0.7), n_obs=299)
    )
t(apply(tests, 1, function(tc) {
    # gamma<-1
    # beta0 <- 14
    # beta1 <- 1
    gamma<-tc[["gamma"]]
    beta0 <- tc[["beta0"]]
    beta1 <- tc[["beta1"]]
    n_obs<- tc[["n_obs"]]
    dt<-1
    mu<-10
    tau<-1
    sigma<- tau * sqrt(2*gamma)
    omega<-exp
    result <- replicate(1,  expr={
        X <- ou_sim(gamma, mu, sigma, dt, n_obs)
        Y <- y_sim(X, beta0, beta1)
        # plot(X, log(Y), main=beta1,cex=0.1)
        # plot(beta0-beta1*X, log(Y), main=beta1,cex=0.1, xlim=c(1,10), ylim=c(1,10))
        # plot(exp(beta0-beta1*X), Y, main=beta1,cex=0.1, xlim=c(1,15000), ylim=c(0,15000))
        # c(range=diff(range(Y)), max=range(Y)[2], n_uniq=length(unique(Y)), num_0=sum(Y == 0) + sum(is.na(Y)))
        sapply(1:n_obs, function(i) {
            dpois(Y[i], exp(beta0-beta1*X[i]), log=TRUE)
        })
    })
    hist(result, main=paste(beta1, sum(result)), xlim=c(-11,0), ylim=c(0,n_obs))
    # c(beta0=beta0, beta1=beta1, gamma=gamma, signif(apply(result, 1, mean), 3))
    # print(result)
    # c(gamma=gamma, n_obs=n_obs, mean=apply(result, 1, mean), sd=apply(result, 1 ,sd))
}))
# testing diff
test_beta0 <- 10
mu <- 10
vec_beta1 <- c(1)
vec_beta0 <- c()
for (beta1 in vec_beta1) {
    beta0 <- max(0,mean(beta1*(mu+3*sigma)-1, beta1*(mu-3*sigma)+20))
    print(paste(beta1*(mu+3*sigma)-1, beta1*(mu-3*sigma)+20))
    vec_beta0 <- c(vec_beta0, beta0)
    print(paste("beta0=", beta0, "beta1=",beta1))
    X <- ou_sim(gamma, mu, sigma, dt, n_obs)
    Y <- y_sim(X, beta0, beta1)
    print(paste("mean Y=", mean(Y), log(mean(Y))))
    print(paste(range(Y), sum(Y == 0)))
}

beta0 <- vec_beta0[1]
beta1 <- vec_beta1[1]
test_diff_beta0 <- data.frame(beta0=beta0/2^c(0,1,2,3), beta1=beta1)
test_diff_beta1 <- data.frame(beta0=beta0, beta1=beta1/2^c(0,1,2,3))

all_tests <- data.frame()
for (gamma in c(1, 0.1, 10)) {
    all_tests <- rbind(all_tests, 
        data.frame(beta0=beta0- seq(-2,6,by=2), beta1=beta1, gamma=gamma,
        n_dataset=100, n_obs = 299, num_multistart=10),
        data.frame(beta0=beta0, beta1=beta1-seq(0.1,0.9,by=0.2), gamma=gamma,
        n_dataset=100, n_obs = 299, num_multistart=10)
    )
}
all_tests
test_diff_beta1 <- data.frame(beta0=beta0, beta1=beta1/2^c(0,1,2,3))

}
test_cases <- do.call("rbind", (lapply(test_beta0, function(beta0) {
    data.frame(beta0=beta1, beta1=seq(max(0.1, (beta0-20)/mu), (beta0+1)/mu,length.out=3), num_multistart=2)
})))
test_cases <- data.frame(beta0=vec_beta0, beta1=vec_beta1, gamma=c(0.1,10))
test_cases
system("osascript -e 'display notification \"beta0=10 sim started\" with title \"STAT440\" sound name \"default\"'")
result <- apply(test_cases, 1, function(tc) {
    print(tc)
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=100, n_obs=299, num_multistart=10)
})
cbind(test_cases, t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))
system("osascript -e 'display notification \"beta0=10 sim done\" with title \"STAT440\" sound name \"default\"'")


test_beta0 <- 5/2^seq(0,2)
mu <- 10
test_cases <- do.call("rbind", (lapply(test_beta0, function(beta0) {
    data.frame(beta0=beta0, beta1=seq(max(0.1, (beta0-20)/mu), (beta0+1)/mu,length.out=5))
})))
system("osascript -e 'display notification \"beta0=5/2^seq(0,2) sim started\" with title \"STAT440\" sound name \"default\"'")
result <- apply(test_cases, 1, function(tc) {
    print(tc)
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=100, n_obs=99, num_multistart=10)
})
cbind(test_cases, t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))
system("osascript -e 'display notification \"beta0=5/2^seq(0,2) sim done\" with title \"STAT440\" sound name \"default\"'")

test_beta0 <- 5/2^seq(3,5)
mu <- 10
test_cases <- do.call("rbind", (lapply(test_beta0, function(beta0) {
    data.frame(beta0=beta0, beta1=seq(max(0.1, (beta0-20)/mu), (beta0+1)/mu,length.out=5))
})))
system("osascript -e 'display notification \"beta0=5/2^seq(3,5) sim started\" with title \"STAT440\" sound name \"default\"'")
result <- apply(test_cases, 1, function(tc) {
    print(tc)
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=100, n_obs=99, num_multistart=10)
})
cbind(test_cases, t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))
system("osascript -e 'display notification \"beta0=5/2^seq(3,5) sim done\" with title \"STAT440\" sound name \"default\"'")

test_beta0 <- 5/2^seq(6,8)
mu <- 10
test_cases <- do.call("rbind", (lapply(test_beta0, function(beta0) {
    data.frame(beta0=beta0, beta1=seq(max(0.1, (beta0-20)/mu), (beta0+1)/mu,length.out=5))
})))
system("osascript -e 'display notification \"beta0=5/2^seq(6,8) sim started\" with title \"STAT440\" sound name \"default\"'")
result <- apply(test_cases, 1, function(tc) {
    print(tc)
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], n_dataset=1, n_obs=99, num_multistart=1)
})
cbind(test_cases, t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))
system("osascript -e 'display notification \"beta0=5/2^seq(6,8) sim done\" with title \"STAT440\" sound name \"default\"'")



test_cases <- expand.grid(n_obs=c(99,199,299,399,499,599), num_multistart=c(1,10), method=c("BFGS", "Nelder-Mead"))
result <- apply(test_cases, 1, function(tc) {
    sim_1(n_dataset=10, n_obs=as.numeric(tc[["n_obs"]]), num_multistart=as.numeric(tc[["num_multistart"]]), method=tc[["method"]])
})
t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)}))
system("osascript -e 'display notification \"n_obs done\" with title \"STAT440\" sound name \"default\"'")


# test omega_mu_tau model without multistart
test_cases <- expand.grid(beta0=10, beta1=0.5, tau=c(0.5,1,2), omega=exp(c(-1,-0.1, -10)), mu=c(1, 10), n_obs=c(99,199,299))
result <- apply(test_cases, 1, function(tc) {
    sim_1(beta0=tc[["beta0"]], beta1=tc[["beta1"]], omega=tc[["omega"]], mu=tc[["mu"]], tau=tc[["tau"]], n_dataset = 100, n_obs=tc[["n_obs"]], num_multistart = 1)
})
cbind(test_cases, t(sapply(1:nrow(test_cases), function(i) {c(theta=result[[i]]$true_param, rmse=result[[i]]$rmse)})))