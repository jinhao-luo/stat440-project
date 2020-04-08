// data inputs
DATA_SCALAR(x0); 
DATA_SCALAR(dt); //  Time gap
DATA_SCALAR(b0); // beta 0
DATA_SCALAR(b1); // beta 1
DATA_VECTOR(y); // number of photons recorded

// Parameter inputs
PARAMETER(gamma); // number of photons recorded
PARAMETER(mu); // parameter estimated
PARAMETER(sigma); // parameter estimated

// Construct joint negative log-likelihood
nll_ou<Type> f_nll(x0, dt, b0, b1, y, gamma, mu, sigma);

int n = y.size();
vector<Type> X(n);
X.setZero();

DATA_INTEGER(niter);
Type res = laplace(f_nll, X, niter);

  
return res;
