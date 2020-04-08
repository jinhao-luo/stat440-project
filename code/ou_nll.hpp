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
  
  joint_nll<Type> jnll(x0, dt, b0, b1, y, gamma, mu, sigma);
  
  /* Constructor */
  joint_nll(Type x0_,
            Type dt_,
            Type b0_,
            Type b1_,
            vector<Type> y_,
            Type gamma_,
            Type mu_,
            Type sigma_):
  x0(x0_), dt(dt_), b0(b0_), b1(b1_),
  y(y_), gamma(gamma_), mu(mu_), sigma(sigma_) {}

  
  /* Evaluate the negative joint log-likelihood as function of the
  random effects */
 Type operator()(vector<Type> X) {

   //This compute the Omega value
   Type omega = exp(-gamma*dt);
   //This compute the tao value
   Type tao = sigma/sqrt(2*gamma);
   
   int N = X.size();
   REPORT(N);
   //this computes X * b1 + b0
   vector<Type> x_beta1 = -X*b1 +b0;
   REPORT(x_beta1);
   // This computes (gamma_n*(X * b1 + b0)-exp(X * b1 + b0))
   Type sum_gamma = (y* x_beta1 -exp(x_beta1)).sum();
   REPORT(sum_gamma);
   
   Type f; //= sum_gamma;// Define variable that holds the return value
   
   for(int i=1;i<N;i++) {
     // Calculate mu and sigma for OU process
     Type mu_ou = mu+omega*(X[i-1]-mu);
     Type sigma_ou = tao*sqrt((1-omega*omega));
     f =  f  + dnorm(X[i],mu_ou,sigma_ou,true); // todo: check later
   }
   return -f; //negative loglikelihood
  }
};
