# Import necessary library
require(TMB)
compile("tutorial.cpp")
dyn.load(dynlib("tutorial"))

f = MakeADFun(data=list(x=x),parameters=list(mu=0,sigma=1))
