// # include <TMB.hpp>
/* The negative loglikelihood when x follows brownian motion
    Note: this implementation is inspired by https://github.com/kaskr/adcomp/blob/master/tmb_examples/laplace.cpp
  */

template<class Type>
struct nll_ou {
  
  // data inputs
  Type x0; 
  Type dt; //  Time gap
  Type b0; // beta 0
  Type b1; // beta 1
  vector<Type> y; // number of photons recorded
  
  // Parameter inputs
  Type gamma; // number of photons recorded
  Type mu; // parameter estimated
  Type sigma; // parameter estimated
  
  /* Constructor */
  nll_ou(Type x0_,
         Type dt_,
         Type b0_,
         Type b1_,
         vector<Type> y_,
         Type gamma_,
         Type mu_,
         Type sigma_):
  x0(x0_), dt(dt_), b0(b0_), b1(b1_),
  y(y_), gamma(gamma_), mu(mu_), sigma(sigma_) {}

  template <typename T>
  T operator()(vector<T> X) {

   // compute the Omega value
   T omega = (T) exp(-gamma*dt);
   // compute the tao value
   T tao = (T) (sigma/sqrt(2*gamma));
   
   int N = X.size();
   // computes b0 - X * b1
   vector<T> x_beta = (T) b0 - (T) b1 * X;
   // computes (gamma_n*( b0 - X * b1 )-exp(b0 - X * b1 ))
   T nll = (T) (exp(x_beta) - y.template cast<T>() * x_beta).sum();
   
   for(int i=1;i<N;i++) {
     // Calculate mu and sigma for OU process
     T mu_ou = (T) mu+ omega*(X[i-1] - (T) mu);
     T sigma_ou = tao*sqrt((1-omega*omega));
     nll -= dnorm(X[i], mu_ou, sigma_ou,true); // todo: check later
   }
   return nll;
  }
};
