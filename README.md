# mtDBN
Experimental implementation of model tree dynamic Bayesian networks, mtDBN in short

This repository will serve as an initial implementation of the mtDBN model prior to being introduced into the _dbnR_ package. Some of the objectives regarding this model are:

* Test the improvements of the mtDBN model versus a regular DBN model
* Test different kinds of splitting criteria for the tree model
* Test multivariate trees versus univariate
* Test the possibility of a forest model tree
* Implement the mtDBN model as an R6 independent module in order to be able to easily introduce it into _dbnR_
* Test the convergence of the PSOHO algorithm. It seems to range from pretty amazing networks to ones with abysmal forecasting performance  
