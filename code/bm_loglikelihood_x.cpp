#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {
  // data inputs
  DATA_VECTOR(gamma); //gamma parameter
  DATA_VECTOR(y); // number of photons recorded
  DATA_INTEGER(b0); // beta 0
  DATA_INTEGER(b1); // beta 1
  DATA_INTEGER(dt); //  Time gap
  DATA_VECTOR(mu); // parameter estimated
  DATA_VECTOR(sigma); // parameter estimated

  //Parameter inputs
  PARAMETER_VECTOR(X); // donor-acceptor distance,underlying latent random variable

  int N = X.size();
  REPORT(N);
  //this computes b0 - X * b1 
  vector<Type> x_beta1 = b0 - X*b1;
  REPORT(x_beta1);
  // This computes (gamma*(b0 - X * b1)-exp(b0 - X * b1))
  Type sum_gamma = (y* x_beta1 -exp(x_beta1)).sum();
  REPORT(sum_gamma);
  Type f = sum_gamma;// Define variable that holds the return value
  
  for(int i=1;i<N;i++) {
    // Calculate mu and sigma for bm process
    vector<Type> mu_bm = X[i-1];
    vector<Type> sigma_bm = sigma*dt;
    f = f+ sum(dnorm(X[i-1],mu_bm,sigma_bm,true));
  }
  return f;
  
}