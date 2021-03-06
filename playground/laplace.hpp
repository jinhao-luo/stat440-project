#include <algorithm>

/* Laplace approximation implementation.
    Note: this implementation is inspired by https://github.com/kaskr/adcomp/blob/master/tmb_examples/laplace.cpp
  */
template<class Type, class Functor>
struct laplace_t {
  Functor f;       // User's implementation of joint nll
  vector<Type>& u; // Random effect vector
  int niter;       // Number of Newton iterations
  laplace_t(Functor f_, vector<Type> &u_, int niter_) :
    f(f_), u(u_), niter(niter_) {}
  Type operator()(){
    // Solve inner problem - Newton iterations

    for (int i=0; i<niter; i++){
      vector<Type> g = autodiff::gradient(f, u);
      
      matrix<Type> H = autodiff::hessian(f, u);
      u = u - atomic::matinv(H) * g;
    }
    // Laplace approximation
    matrix<Type> H = autodiff::hessian(f, u);
    Type ans = .5 * atomic::logdet(H) + f(u);
    return ans;
  }
};
template<class Type, class Functor>
Type laplace(Functor f, vector<Type> &u, int niter){
  laplace_t<Type, Functor> L(f, u, niter);
  return L();
}