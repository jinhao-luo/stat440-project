// data inputs
DATA_INTEGER(niter);
DATA_SCALAR(dt); //  Time gap
DATA_SCALAR(beta0); // beta 0
DATA_SCALAR(beta1); // beta 1
DATA_VECTOR(y); // number of photons recorded

// Parameter inputs
PARAMETER(gamma); // number of photons recorded
PARAMETER(mu); // parameter estimated
PARAMETER(sigma); // parameter estimated

// Construct joint negative log-likelihood
nll_ou<Type> f_nll(dt, beta0, beta1, y, gamma, mu, sigma);

int n = y.size();
vector<Type> X(n);
X.setZero();

Type res = laplace(f_nll, X, niter);

  
return res;
