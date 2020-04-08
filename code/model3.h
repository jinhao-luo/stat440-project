
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


  //This compute the Omega value
  Type omega = exp(-gamma*dt);
  //This compute the tao value
  Type tao = sigma/sqrt(2*gamma);

  int N = X.size();
  REPORT(N);
  //this computes X * b1 + b0
  vector<Type> x_beta1 = -X*b1 +b0;
  REPORT(x_beta1);
  // This computes (gamma_n*(-X * b1 + b0)-exp(-X * b1 + b0))
  Type sum_gamma = (gamma* x_beta1 -exp(x_beta1)).sum();
  REPORT(sum_gamma);

  Type f = sum_gamma;// Define variable that holds the return value

  for(int i=1;i<N;i++) {
    // Calculate mu and sigma for OU process
    Type mu_ou = mu+omega*(X[i-1]-mu);
    Type sigma_ou = tao*sqrt((1-omega*omega));
    f =  f  + dnorm(X[i],mu_ou,sigma_ou,true); // todo: check later
  }
  //Calculate the norm of the Hessian Matrix
  Type Norm = (H.array()*H.array()).sum();
  REPORT(Norm);
  f =f - 0.5*Norm;
  return -f; //negative loglikelihood
