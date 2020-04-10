#include <TMB.hpp>
#include "laplace.hpp"
#include "nll_bm.hpp"
#include "nll_ou.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model_type);
  if (model_type == "ou") {
#include "main_ou.cpp"
  } else if (model_type == "bm") {
#include "main_bm.cpp"
  } else {
    error ("Unknown model type");
  }
  return 0;
}
