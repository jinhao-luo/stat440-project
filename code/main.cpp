#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model_type);
  if (model_type == "model1") {
#include "model1.h"
  } else if (model_type == "model2") {
#include "model2.h"
  } else if (model_type == "model3") {
#include "model3.h"
  } else {
    error ("Unknown model type");
  }
  return 0;
}
