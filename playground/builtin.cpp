/** Implementation using built-in laplace.
 * Usage: MakeADFun(..., random=c("X"))
 */

#include <TMB.hpp>
#include "nll_bm.hpp"
#include "nll_ou.hpp"

template <class Type>
Type objective_function<Type>::operator()()
{
	DATA_STRING(model_type);
	PARAMETER_VECTOR(X);
	if (model_type == "ou")
	{
		// data inputs
		DATA_SCALAR(dt);	//  Time gap
		DATA_SCALAR(beta0); // beta 0
		DATA_SCALAR(beta1); // beta 1
		DATA_VECTOR(y);		// number of photons recorded

		// Parameter inputs
		PARAMETER(gamma); // number of photons recorded
		PARAMETER(mu);	  // parameter estimated
		PARAMETER(sigma); // parameter estimated

		// Construct joint negative log-likelihood
		nll_ou<Type> f_nll(dt, beta0, beta1, y, gamma, mu, sigma);

		Type res = f_nll(X);
		return res;
	}
	else if (model_type == "omega_tau")
	{

		// data inputs
		DATA_SCALAR(dt);	//  Time gap
		DATA_SCALAR(beta0); // beta 0
		DATA_SCALAR(beta1); // beta 1
		DATA_VECTOR(y);		// number of photons recorded

		// Parameter inputs
		PARAMETER(omega); // number of photons recorded
		PARAMETER(mu);	  // parameter estimated
		PARAMETER(tau);	  // parameter estimated

		int N = X.size();
		// computes beta0 - X * beta1
		vector<Type> x_beta = beta0 - beta1 * X;
		// computes (gamma_n*( beta0 - X * beta1 )-exp(beta0 - X * beta1 ))
		Type nll = (exp(x_beta) - y * x_beta).sum();

		// Calculate sigma for OU process
		Type sigma_ou = tau * sqrt(1 - omega * omega);

		for (int i = 1; i < N; i++)
		{
			Type mu_ou = mu + omega * (X[i - 1] - mu);
			nll -= dnorm(X[i], mu_ou, sigma_ou, true); // todo: check later
		}
		return nll;
	}
	else if (model_type == "bm")
	{
		// data inputs
		DATA_SCALAR(dt); // interobservation time
		DATA_SCALAR(beta0);
		DATA_SCALAR(beta1);
		DATA_VECTOR(Y); // number of photons recorded

		// Parameter inputs
		PARAMETER(sigma); // parameter estimated
		int n = Y.size();

		// Construct joint negative log-likelihood
		nll_bm<Type> f_nll(Y, dt, beta0, beta1, sigma);

		Type res = f_nll(X);
		return res;
	}
	else
	{
		error("Unknown model type");
	}
	return 0;
}
