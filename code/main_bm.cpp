// data inputs
DATA_SCALAR(dt); // interobservation time
DATA_SCALAR(beta0);
DATA_SCALAR(beta1);
DATA_VECTOR(Y); // number of photons recorded

DATA_INTEGER(niter); // number of iterations for laplace

// Parameter inputs
PARAMETER(sigma); // parameter estimated
int n = Y.size();


// Construct joint negative log-likelihood
nll_bm<Type> f_nll(Y, dt, beta0, beta1, sigma);

// Random effect initial guess
vector<Type> u(n);
u.setZero();

// Calculate Laplace approx (updates u)
Type res = laplace(f_nll, u, niter);

return res;
