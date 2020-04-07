
#include <TMB.hpp>

template<class Type>

Type objective_function<Type>::operator() () {
  // data inputs
  DATA_MATRIX(H); // rate covariate matrix
  DATA_VECTOR(X); // donor-acceptor distance,underlying latent random variable
  DATA_INTEGER(dt); //  Time gap
  DATA_INTEGER(b0); // beta 0
  DATA_INTEGER(b1); // beta 1
  DATA_VECTOR(y); // number of photons recorded
  
  // Parameter inputs
  PARAMETER_VECTOR(gamma); // number of photons recorded
  PARAMETER_VECTOR(mu); // parameter estimated
  PARAMETER_VECTOR(sigma); // parameter estimated
  
  
  int N = X.size();
  REPORT(N);
  //this computes b0 - X * b1 
  vector<Type> x_beta1 = b0-b1*X;
  REPORT(x_beta1);
  
  // This computes (gamma *(b0 - X * b1)-exp(b0 - X * b1))
  Type sum_gamma = (gamma* x_beta1 - exp(x_beta1)).sum();
  REPORT(sum_gamma);
  
  Type f = sum_gamma;// Define variable that holds the return value
  
  for(int i=1;i<N;i++) {
    // Calculate mu and sigma for bm process
    vector<Type> mu_bm = X[i-1];
    vector<Type> sigma_bm = sigma*dt;
    f = f+ sum(dnorm(X[i-1],mu_bm,sigma_bm,true));
  }
  //Calculate the norm of the Hessian Matrix
  Type Norm = (H.array()*H.array()).sum();
  REPORT(Norm);
  f =f - 0.5*Norm;
  return f;
  
}




