# Flexible Modeling of Variable Asymmetries in Cross-covariance Functions for Multivariate Random Fields
This repository provides reproducible code for the manuscript titled *Flexible Modeling of Variable Asymmetries in Cross-covariance Functions for Multivariate Random Fields* by Ghulam A. Qadir, Carolina Éuan, and Ying Sun. The manuscript describes a new approach to introduce flexible asymmetries in the stationary cross-covariance functions through modification of phase spectrum.

## Abstract

The geostatistical analysis of multivariate spatial data for inference as well as joint predictions (co-kriging) ordinarily relies on modeling of the marginal and cross-covariance functions. While the former quantifies the spatial dependence within variables, the latter quantifies the spatial dependence across distinct variables. 
The marginal covariance functions are always symmetric; however, the cross-covariance functions often exhibit asymmetries in the real data. Asymmetric cross-covariance implies change in the value of cross-covariance for interchanged locations on fixed order of variables. Such change of cross-covariance values are often caused due to the spatial delay in effect of the response of one variable on another variable. These spatial delays are common in environmental processes, especially when dynamic phenomena such as prevailing wind, ocean currents, etc., are involved.  Here, we propose a novel approach to introduce flexible asymmetries in the cross-covariances of stationary multivariate covariance functions. The proposed approach involves modeling the phase component of the constrained cross-spectral features to allow for asymmetric cross-covariances. We show the capability of our proposed model to recover the cross-dependence structure and improve spatial predictions against traditionally used models through multiple simulation studies. Additionally, we illustrate our approach on a real trivariate dataset of particulate matter concentration (PM<sub>2.5<sub>), wind speed and relative humidity. The real data example shows that our approach outperforms the traditionally used models, in terms of model fit and spatial predictions.
## Requirements

The codes are written in R, and reproducing would require installing and loading the following R-packages: `fields`,`sp`,`maps`,`maptools`,`geosphere`,`MASS`,`scoringRules`,`doParallel`,`rgdal`,`ggplot2`,`gridExtra`, `RColorBrewer`, and `viridis`. 

##

