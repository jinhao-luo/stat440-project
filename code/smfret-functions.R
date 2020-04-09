#' Generate observations from an Ornstein-Uhlenbeck process at regular time intervals.
#'
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
ou_sim <- function(gamma, mu, sigma, dt, n_obs, x0) {
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
#' @param gamma Scalar mean reversion parameter (see [ou_sim()]).
#' @param mu Scalar mean parameter (see [ou_sim()]).
#' @param sigma Scalar diffusion parameter (see [ou_sim()]).
#' @param Xt Vector of `n_obs` observations from the OU process.
#' @param dt Interobservation time.
#' @return The scalar value of the loglikelihood `loglik(gamma, mu, sigma | Xt)`.
#'
#' @details For simplicity, the loglikelihood contribution of the first observation is ignored.
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
#' @param beta0 Scaler intercept parameter for Yt Poisson distribution
#' @param beta1 Scaler slope parameter for Yt Poisson distribution
#' @param Xt Vector of `n_obs` observations from the OU process.
#' @param Yt Vector of `n_obs` observations from the OU process.
#' @param dt Interobservation time.
#' @return The scalar value of the loglikelihood `loglik(gamma, mu, sigma | Xt, Yt)`.
#'
ou_y_nll <- function(gamma, mu, sigma, b0, b1, Xt, Yt, dt) {
  ou_nll(gamma, mu, sigma, Xt, dt)-sum(dpois(Yt, exp(b0-b1*Xt), log=TRUE))
  # ou_nll(gamma, mu, sigma, Xt, dt)-sum(Yt*(b0-b1*Xt)-exp(b0-b1*Xt))
}


#' Generate observations from a Brownian motion with drift at regular time intervals.
#'
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
bm_sim <- function(mu, sigma, dt, n_obs, x0 = 0) {
  # Brownian increments
  dX <- rnorm(n_obs, mean = mu * dt, sd = sigma * sqrt(dt))
  x0 + cumsum(dX)
}

#' Negative loglikelihood for the Brownian model Xn.
#'
#' @param gamma Scalar mean reversion parameter (see [ou_sim()]).
#' @param mu Scalar mean parameter (see [ou_sim()]).
#' @param sigma Scalar diffusion parameter (see [ou_sim()]).
#' @param Xt Vector of `n_obs` observations from the OU process.
#' @param dt Interobservation time.
#' @return The scalar value of the loglikelihood `loglik(gamma, mu, sigma | Xt)`.
#'
#' @details For simplicity, the loglikelihood contribution of the first observation is ignored.
bm_nll <- function(mu, sigma, Xt, dt) {
  n <- length(Xt)
  -sum(dnorm(Xt[2:n], log = TRUE,mean = Xt[1:(n-1)], sd = sigma*sqrt(dt)))
}

#' Negative loglikelihood for the Brownian model Yn.
#'
#' @param gamma Scalar mean reversion parameter (see [ou_sim()]).
#' @param mu Scalar mean parameter (see [ou_sim()]).
#' @param sigma Scalar diffusion parameter (see [ou_sim()]).
#' @param beta0 Scaler intercept parameter for Yt Poisson distribution
#' @param beta1 Scaler slope parameter for Yt Poisson distribution
#' @param Xt Vector of `n_obs` observations from the OU process.
#' @param Yt Vector of `n_obs` observations from the OU process.
#' @param dt Interobservation time.
#' @return The scalar value of the loglikelihood `loglik(gamma, mu, sigma | Xt, Yt)`.
#'
bm_y_nll <- function(mu, sigma, b0, b1, Xt, Yt, dt) {
  bm_nll(mu, sigma, Xt, dt)-sum(dpois(Yt, exp(b0-b1*Xt), log=TRUE))
  # ou_nll(gamma, mu, sigma, Xt, dt)-sum(Yt*(b0-b1*Xt)-exp(b0-b1*Xt))
}

#' Morse potential function.
#'
#' @param x Vector of nonnegative intermolecular distances.
#' @param gamma Potential well depth parameter (scalar or vector; see [morse_sim()]).
#' @param alpha Potential well width parameter (scalar or vector; see [morse_sim()]).
#' @param mu Equilibrium bond length parameter (scalar or vector; see [morse_sim()]).
#' @return The Morse potential function evaluated at `x`.
#' @details The Morse potential function is
#' ```
#' U(x | gamma, alpha, mu) = gamma * ( 1 - exp( -alpha * (x-mu) ) )^2.
#' ```
morse_pot <- function(x, gamma, alpha, mu) {
  gamma * (1 - exp(-alpha * (x-mu)))^2
}


#' Derivative of the Morse potential function.
#'
#' @param x Vector of nonnegative intermolecular distances.
#' @param gamma Potential well depth parameter (scalar or vector; see [morse_pot()]).
#' @param alpha Potential well width parameter (scalar or vector; see [morse_pot()]).
#' @param mu Equilibrium bond length parameter (scalar or vector; see [morse_pot()]).
#' @return Derivative of Morse potential wrt `x` evaluated at the inputs.
#'
#' @details The Morse potential function is
#' ```
#' U(x | gamma, alpha, mu) = gamma * ( 1 - exp( -alpha * (x-mu) ) )^2.
#' ```
#' Its derivative is
#' ```
#' U'(x | gamma, alpha, mu) = 2 * gamma * ( 1 - exp( -alpha * (x-mu) ) ) * alpha * exp( -alpha * (x-mu) ).
#' ```
morse_dpot <- function(x, gamma, alpha, mu) {
  eax <- exp(-alpha * (x-mu))
  2 * gamma * (1-eax) * alpha * eax
}

#' Generate observations from a Morse SDE.
#'
#' @param gamma Scalar potential well depth parameter (see **Details**).
#' @param alpha Scalar potential well width parameter (see **Details**).
#' @param mu Scalar equilibrium bond length parameter (see **Details**).
#' @param sigma Scalar diffusion parameter (see **Details**).
#' @param dt Interobservation time.
#' @param n_obs Number of observations to generate.
#' @param x0 Initial value of the process at time `t = 0`.
#' @param n_sub Integer greater than zero specifying the sampling interobservation time.  That is, the data is generated with `dt_sim = dt/n_sub`, and only every `n_sub` observation is retained.
#' @return A vector of `n_obs` observations of the process at times `t = dt, 2*dt, ..., n_obs * dt`.
#' @details The Morse SDE is defined as
#' ```
#' dXt = -U'(Xt | gamma, alpha, mu) dt + sigma dBt,
#' ```
#' where `U'(x | gamma, alpha, mu)` is the derivative (wrt `x`) of the Morse potential function (see [morse_pot()] and [morse_dpot()]).  Since the exact transition density `p(X_{t+dt} | X_t)` for this model is not available in closed form, an Euler approximation on the finer timescale `dt_sim = dt/n_sub` is used instead, i.e.,
#' ```
#' X_{t + dt_sim} | X_t ~ N( X_t - U'(X_t | gamma, alpha, mu) * dt_sim, sigma^2 * dt_sim).
#' ```
morse_sim <- function(gamma, alpha, mu, sigma, dt, n_obs, x0, n_sub) {
  n_tot <- n_obs*n_sub # total number of simulated timepoints
  Xt <- rep(NA, n_tot) # preallocate memory
  # pre-generate all normal draws
  dt_sub <- dt/n_sub
  Zt <- rnorm(n_tot, sd = sigma * sqrt(dt_sub))
  x_curr <- x0 # current value
  for(ii in 1:n_tot) {
    # euler approximation
    x_curr <- x_curr - morse_dpot(x_curr, gamma, alpha, mu) * dt_sub + Zt[ii]
    Xt[ii] <- x_curr
  }
  # pruning
  Xt[seq(from = n_sub, to = n_tot, by = n_sub)]
}

y_sim <- function(X, beta0, beta1) {
  sapply(X, function(x) rpois(1, exp(beta0-beta1*x)))
}



