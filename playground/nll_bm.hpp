// # include <TMB.hpp>
/* The negative loglikelihood when x follows bm motion
    Note: this implementation is inspired by https://github.com/kaskr/adcomp/blob/master/tmb_examples/laplace.cpp
  */
template<class Type>
struct nll_bm {
    // Data and parameter objects
    vector<Type> Y;
    Type dt;
    Type beta0;
    Type beta1;
    Type sigma;

    /* Constructor */
    nll_bm(vector<Type> Y_,
	    Type dt_,
	    Type beta0_,
	    Type beta1_,
	    Type sigma_ ) :
    Y(Y_), dt(dt_), beta0(beta0_), beta1(beta1_), sigma(sigma_) {}

    template <typename T>
    T operator()(vector<T> X) {
        vector<T> x_beta = (T) beta0 - (T) beta1 * X ;
        T nll = (exp(x_beta) - Y.template cast<T>() * x_beta ).sum();
        int N = X.size();

        T sigma_bm = (T) (sigma*sqrt(dt));

        for(int i=1;i<N;i++) {
            nll -= dnorm(X[i], X[i-1], sigma_bm, true);
        }
        return nll;
    }
};