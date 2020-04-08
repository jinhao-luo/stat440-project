
#include <laplace.hpp>

/* The following is (almost) copy-pasted from the 'spatial' example */
template<class Type>
struct joint_nll {
  /* Data and parameter objects for spatial example: */
  vector<Type> y;
  matrix<Type> X;
  matrix<Type> dd;
  vector<Type> b;
  Type a;
  Type log_sigma;
  
  /* Constructor */
  joint_nll(vector<Type> y_,
            matrix<Type> X_,
            matrix<Type> dd_,
            vector<Type> b_,
            Type a_,
            Type log_sigma_) :
    y(y_), X(X_), dd(dd_), b(b_),
    a(a_), log_sigma(log_sigma_) {}
  
  /* Evaluate the negative joint log-likelihood as function of the
  random effects */
  template <typename T>
  T operator()(vector<T> u) {
    int n = u.size();
    T res=0;
    vector<T> eta = T(exp(log_sigma)) * u;
    vector<Type> mu = X * b;
    eta = eta + mu.template cast<T>();
    matrix<T> cov(n,n); 
    for (int i=0; i<n; i++)
    {
      cov(i,i) = 1.0;
      for (int j=0; j<i; j++)
      {
        // Exponentially decaying correlation
        cov(i,j) = exp(-a * dd(i,j));
        cov(j,i) = cov(i,j);
      }
    }
    density::MVNORM_t<T> neg_log_density(cov);
    res += neg_log_density(u);
    // logdpois = N log lam - lam
    for(int i=0; i<y.size(); i++)
      res -= T(y[i]) * eta[i] - exp(eta[i]);
    return res;
  }
};


template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA_VECTOR(y);
  // DATA_MATRIX(X);
  // DATA_MATRIX(dd);
  // PARAMETER_VECTOR(b);
  // PARAMETER(a);
  // PARAMETER(log_sigma);
  // int n = dd.rows();
  
  // data inputs
  DATA_MATRIX(H); // rate covariate matrix
  DATA_VECTOR(X); // donor-acceptor distance,underlying latent random variable
  DATA_SCALAR(dt); //  Time gap
  DATA_SCALAR(b0); // beta 0
  DATA_SCALAR(b1); // beta 1
  DATA_VECTOR(y); // number of photons recorded
  
  // Parameter inputs
  PARAMETER(gamma); // number of photons recorded
  PARAMETER(mu); // parameter estimated
  PARAMETER(sigma); // parameter estimated
  
  // Construct joint negative log-likelihood
  joint_nll<Type> jnll(y, X, dd, b, a, log_sigma);
  
  // Random effect initial guess
  vector<Type> u(n);
  X.setZero();
  
  // Calculate Laplace approx (updates u)
  DATA_INTEGER(niter);
  Type res = laplace(jnll, u, niter);
  ADREPORT(u)
    
    return res;
}