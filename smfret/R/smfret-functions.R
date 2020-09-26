#' Generate observations from an Ornstein-Uhlenbeck process at regular time intervals.
#'
#' @author Martin Lysy \email{mlysy@@uwaterloo.ca}
#' @param gamma Scalar mean reversion parameter (see **Details**).
#' @param mu Scalar mean parameter (see **Details**).
#' @param sigma Scalar diffusion parameter (see **Details**).
#' @param dt Interobservation time.
#' @param n_obs Number of observations to generate.
#' @param x0 Initial value of the process at time `t = 0`.  If missing sampled from the OU stationary distribution.
#' @return A vector of `n_obs` observations of the process at times `t = dt, 2*dt, ..., n_obs * dt`.
#' @details The Ornstein-Uhlenbeck (OU) process satisfies a stochastic differential equation of the form
#' ```
#' dX_t = -gamma * (X_t - mu) dt + sigma dB_t.
#' ```
#' It is a stationary Gaussian Markov process with transition density
#' ```
#' X_s+t | X_s ~ N( rho_t * (X_s - mu) + mu, tau^2 * (1-rho_t^2) ),
#' ```
#' where `rho_t = exp(-gamma * t)` and `tau^2 = sigma^2/(2*gamma)`.  Its stationary distribution is `X_t ~ N(mu, tau^2)`.
#' @export
ou_sim <- function(gamma, mu, sigma, dt, n_obs, x0=mu) {
  tau <- sigma/sqrt(2*gamma) # stationary standard deviation
  if(missing(x0)) x0 <- rnorm(1, mean = mu, sd = tau)
  # generate efficiently using a one-step linear filter
  lrho <- -gamma * dt
  ou_sd <- tau * sqrt((1-exp(2 * lrho))) # conditional sd
  ou_filt <- exp(lrho) # filter coefficients
  z <- rnorm(n_obs, sd = ou_sd) # pre-generate normal draws
  Xt <- filter(x = z, filter = ou_filt,
               method = "recursive", init = x0 - mu)
  return(as.numeric(Xt + mu))
}

#' Negative loglikelihood for the Ornstein-Uhlenbeck model Xn.
#'
#' @author Martin Lysy \email{mlysy@@uwaterloo.ca}
#' @param gamma Scalar mean reversion parameter (see [ou_sim()]).
#' @param mu Scalar mean parameter (see [ou_sim()]).
#' @param sigma Scalar diffusion parameter (see [ou_sim()]).
#' @param Xt Vector of `n_obs` observations from the OU process.
#' @param dt Interobservation time.
#' @return The scalar value of the loglikelihood `loglik(gamma, mu, sigma | Xt)`.
#'
#' @details For simplicity, the loglikelihood contribution of the first observation is ignored.
#' @export
ou_nll <- function(gamma, mu, sigma, Xt, dt) {
  tau <- sigma/sqrt(2*gamma) # stationary standard deviation
  lrho <- -gamma * dt
  ou_sd <- tau * sqrt((1-exp(2 * lrho))) # conditional sd
  ou_filt <- exp(lrho) # filter coefficients
  n <- length(Xt)
  -sum(dnorm(Xt[2:n], log = TRUE,
             mean = ou_filt * (Xt[1:(n-1)] - mu) + mu, sd = ou_sd))
}

#' Negative loglikelihood for the Ornstein-Uhlenbeck model Yn.
#'
#' @param gamma Scalar mean reversion parameter (see [ou_sim()]).
#' @param mu Scalar mean parameter (see [ou_sim()]).
#' @param sigma Scalar diffusion parameter (see [ou_sim()]).
#' @param b0 Scaler intercept parameter for Yt Poisson distribution
#' @param b1 Scaler slope parameter for Yt Poisson distribution
#' @param Xt Vector of `n_obs` observations from the OU process.
#' @param Yt Vector of `n_obs` observations from the OU process.
#' @param dt Interobservation time.
#' @return The scalar value of the loglikelihood `loglik(gamma, mu, sigma | Xt, Yt)`.
#'
#' @export
ou_y_nll <- function(gamma, mu, sigma, b0, b1, Xt, Yt, dt) {
  ou_nll(gamma, mu, sigma, Xt, dt)-sum(dpois(Yt, exp(b0-b1*Xt), log=TRUE))
  # ou_nll(gamma, mu, sigma, Xt, dt)-sum(Yt*(b0-b1*Xt)-exp(b0-b1*Xt))
}


#' Generate observations from a Brownian motion with drift at regular time intervals.
#'
#' @author Martin Lysy \email{mlysy@@uwaterloo.ca}
#' @param mu Scalar drift parameter (see **Details**).
#' @param sigma Scalar diffusion parameter (see **Details**).
#' @param dt Interobservation time.
#' @param n_obs Number of observations to generate.
#' @param x0 Initial value of the process at time `t = 0`.  Default value is `x0 = 0`.
#' @return A vector of `n_obs` observations of the process at times `t = dt, 2*dt, ..., n_obs * dt`.
#' @details Brownian motion with drift is a Gaussian Markov process with transition density
#' ```
#' X_s+t | X_s ~ N(X_s + mu * t, sigma^2 * t).
#' ```
#' @export
bm_sim <- function(mu, sigma, dt, n_obs, x0 = 0) {
  # Brownian increments
  dX <- rnorm(n_obs, mean = mu * dt, sd = sigma * sqrt(dt))
  x0 + cumsum(dX)
}

#' Negative loglikelihood for the Brownian model Xn.
#'
#' @param mu (NOT USED) Scalar mean parameter (see [ou_sim()]). 
#' @param sigma Scalar diffusion parameter (see [ou_sim()]).
#' @param Xt Vector of `n_obs` observations from the OU process.
#' @param dt Interobservation time.
#' @return The scalar value of the loglikelihood `loglik(gamma, mu, sigma | Xt)`.
#'
#' @details For simplicity, the loglikelihood contribution of the first observation is ignored.
#' @export
bm_nll <- function(mu, sigma, Xt, dt) {
  n <- length(Xt)
  -sum(dnorm(Xt[2:n], log = TRUE,mean = Xt[1:(n-1)], sd = sigma*sqrt(dt)))
}

#' Negative loglikelihood for the Brownian model Yn.
#'
#' @param mu Scalar mean parameter (see [ou_sim()]).
#' @param sigma Scalar diffusion parameter (see [ou_sim()]).
#' @param b0 Scaler intercept parameter for Yt Poisson distribution
#' @param b1 Scaler slope parameter for Yt Poisson distribution
#' @param Xt Vector of `n_obs` observations from the OU process.
#' @param Yt Vector of `n_obs` observations from the OU process.
#' @param dt Interobservation time.
#' @return The scalar value of the loglikelihood `loglik(gamma, mu, sigma | Xt, Yt)`.
#'
#' @export
bm_y_nll <- function(mu, sigma, b0, b1, Xt, Yt, dt) {
  bm_nll(mu, sigma, Xt, dt)-sum(dpois(Yt, exp(b0-b1*Xt), log=TRUE))
  # ou_nll(gamma, mu, sigma, Xt, dt)-sum(Yt*(b0-b1*Xt)-exp(b0-b1*Xt))
}

#' Simulate number of photons exchanged at different time
#'
#' @param beta0 Scaler intercept parameter for Yt Poisson distribution
#' @param beta1 Scaler slope parameter for Yt Poisson distribution
#' @param X Vector of `n_obs` observations of donor-acceptor distance.
#' @return Vector of `n_obs` observation of photons exchanged at different time.
#' @export
#'
y_sim <- function(X, beta0, beta1) {
  sapply(X, function(x) rpois(1, exp(beta0-beta1*x)))
}



