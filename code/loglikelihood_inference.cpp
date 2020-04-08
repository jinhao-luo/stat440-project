#include <TMB.hpp>

template<class Type>

Type objective_function<Type>::operator() () {
  // data inputs
  DATA_VECTOR(X); // donor-acceptor distance,underlying latent random variable
  DATA_INTEGER(dt); //  Time gap
  DATA_INTEGER(b0); // beta 0
  DATA_INTEGER(b1); // beta 1
  DATA_VECTOR(y); // number of photons recorded

  // Parameter inputs
  PARAMETER_VECTOR(gamma); // number of photons recorded
  PARAMETER_VECTOR(mu); // parameter estimated
  PARAMETER_VECTOR(sigma); // parameter estimated


  //This compute the Omega value
  vector<Type> omega = exp(-gamma*dt);
  //This compute the tao value
  vector<Type> tao = sigma/sqrt(2*gamma);

  int N = X.size();
  REPORT(N);
  //this computes X * b1 + b0
  vector<Type> x_beta1 = -X*b1 +b0;
  REPORT(x_beta1);
  // This computes (gamma_n*(X * b1 + b0)-exp(X * b1 + b0))
  Type sum_gamma = (gamma* x_beta1 -exp(x_beta1)).sum();
  REPORT(sum_gamma);

  Type f = sum_gamma;// Define variable that holds the return value

  for(int i=1;i<=N;i++) {
    // Calculate mu and sigma for OU process
    vector<Type> mu_ou = mu+omega*(X[i-1]-mu);
    vector<Type> sigma_ou = tao*sqrt((1-omega*omega));
    f =  f  + dnorm(X[i],mu_ou,sigma_ou,true)[0]; // todo: check later
  }
  return -f; //negative loglikelihood
}
