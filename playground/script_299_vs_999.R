source("smfret-functions.R")
tests <- rbind(
    data.frame(gamma=c(0.1,10), beta0=15, beta1=1, n_obs=299),
    data.frame(gamma=c(0.1,10), beta0=15, beta1=1, n_obs=999)
)
result <- t(apply(tests, 1, function(tc) {
    gamma<-tc[["gamma"]]
    beta0 <- tc[["beta0"]]
    beta1 <- tc[["beta1"]]
    n_obs<- tc[["n_obs"]]
    dt<-1
    mu<-10
    tau<-1
    sigma<- tau * sqrt(2*gamma)
    result <- replicate(1000,  expr={
        X <- ou_sim(gamma, mu, sigma, dt, n_obs)
        Y <- y_sim(X, beta0, beta1)
        # c(range=diff(range(Y)), max=range(Y)[2], n_uniq=length(unique(Y)), num_0=sum(Y == 0) + sum(is.na(Y)))
        c(mean.X=mean(X), sd.X=sd(X), mean.Y=mean(Y), sd.Y=sd(Y))
    })
    write.csv(t(result), paste0("effect_n_obs_gamma_",gamma, "_n_obs_", n_obs, ".csv"))
    c(gamma=gamma, n_obs=n_obs, mean=apply(result, 1, mean), sd=apply(result, 1 ,sd))
}))
write.csv(result, "effect_n_obs_summary.csv")