
########################################################
#### Simulating symmetric trivariate random field #####
########################################################
setwd("D:/Projects 3/Revision codes - Simulation/Case1")
dir()

##### Loading libraries ######
library(fields)
library(MASS)
library(doParallel)
library(viridis)

registerDoParallel(cores = 50)
getDoParWorkers()

##### Functions #########

my.matern<-function(h,a,sigma,nu)
{
  h[h==0]<-1e-10
  num1<-(sigma^2)*(2^(1-nu))/gamma(nu)
  num2<-(h*a)^nu
  num3<-besselK(x=(h*a), nu=nu)
  return(num1*num2*num3)
}


plotmatrix<-function(C)
{
  n1<-dim(C)[1]
  matx<-rep(1:n1,times=n1)
  maty<-rep(n1:1,each=n1)
  quilt.plot(matx,maty,c(C),nx=n1,ny=n1,axes=F,ylim=c(0,n1+50),xlim=c(-50,n1))
  text(x=(n1/3)*0.5,y=n1+20,"X1")
  text(x=-20,y=(n1/3)*2.5,"X1")
  text(x=(n1/3)*1.5,y=n1+20,"X2")
  text(x=-20,y=(n1/3)*1.5,"X2")
  text(x=(n1/3)*2.5,y=n1+20,"X3")
  text(x=-20,y=(n1/3)*0.5,"X3")
}


########################################################################################
##### parameterizing for the positive definite colocated correlation coefficient #######
########################################################################################




rho_par<-function(a)
{ 
  
  L1<-matrix(c(1,a[1],a[2],0,1,a[3],0,0,1),byrow = T,nrow=3)
  un_norm<-L1%*%t(L1)
  dterm<-1/sqrt(diag(un_norm))
  dmat<-matrix(0,nrow=3,ncol=3)
  diag(dmat)<-dterm
  return(dmat%*%un_norm%*%dmat)
}
##########################################
############# Simulation grid ############
##########################################

n=20   ## number of locations = n^2
x<-y<-seq(0,1,length=n)
grid<-expand.grid(x,y)
par(mfrow=c(1,1))
plot(grid,xlab="x",ylab="y",main="Simulation grid")

######################################################
######### Setting marginal smoothness parameters #####
######################################################

nu1<-1
nu2<-2
nu3<-3

######################################################################################################################
########## Setting the beta matrix in the co-located correlation coefficient of the Parsimonious Matern model ########
######################################################################################################################

beta<-rho_par(c(0.7,0.52,0.50)) #try (0.8,0.6,0.5) condition is invalid but the model is still valid! why?

beta


### For d=2 ### 
rho12<-beta[1,2]*(sqrt(gamma(nu1+1))*sqrt(gamma(nu2+1))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*gamma((nu1+nu2)/2+1))
rho13<-beta[1,3]*(sqrt(gamma(nu1+1))*sqrt(gamma(nu3+1))*gamma((nu1+nu3)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu3))*gamma((nu1+nu3)/2+1))
rho23<-beta[2,3]*(sqrt(gamma(nu2+1))*sqrt(gamma(nu3+1))*gamma((nu2+nu3)/2))/(sqrt(gamma(nu2))*sqrt(gamma(nu3))*gamma((nu2+nu3)/2+1))

rhomat<-matrix(c(1,rho12,rho13,rho12,1,rho23,rho13,rho23,1),nrow = 3)

##############################################################
##### Matrix for the colocated correlation coefficient #######
##############################################################

rhomat

chol(rhomat)

#rhomat[1,2]^2+rhomat[1,3]^2+rhomat[2,3]^2-1<(-2*rhomat[1,2]*rhomat[1,3]*rhomat[2,3])

##########################################################
######### Checking validity conditions ###################
##########################################################

tv<-(beta[1,2]^2+beta[1,3]^2+beta[2,3]^2-1)<=(-2*abs(beta[1,2])*abs(beta[1,3])*abs(beta[2,3]))
tv


##########################################################
##### Setting common spatial scale parameters #############
##########################################################
a1<-a2<-a3<-a12<-a13<-a23<-7.5

##########################################################
##### Setting cross-smoothness parameters ################
##########################################################


nu12<-(nu1+nu2)/2
nu13<-(nu1+nu3)/2
nu23<-(nu2+nu3)/2

#######################################################################
##### Computing marginal and cross-covariance matrices ################
#######################################################################


x.s<-c(grid$Var1)
y.s<-c(grid$Var2)
dist.mat<-rdist(cbind(x.s,y.s))
C11<-my.matern(h=dist.mat,a=a1,sigma = 1,nu=nu1)
C22<-my.matern(h=dist.mat,a=a2,sigma = 1,nu=nu2)
C33<-my.matern(h=dist.mat,a=a3,sigma = 1,nu=nu3)
C13<-rhomat[1,3]*my.matern(h=dist.mat,a=a13,sigma = 1,nu=nu13)
C23<-rhomat[2,3]*my.matern(h=dist.mat,a=a23,sigma = 1,nu=nu23)


#####################################
#### Now introducing lag in C12 #####
#####################################
# C2 is the full covariance matrix
###### Asymmetry parameters set to be 0, therefore symmetry #####
a121=0
a122=0
tempgrid1<-cbind(rep(a121,times=n^2),rep(a122,times=n^2))

shift12grid<-grid-tempgrid1

shift.dist12<-rdist(grid,shift12grid)

C12<-rhomat[1,2]*my.matern(h=shift.dist12,a=a12,sigma = 1,nu=nu12)
C21<-t(C12)

C2<-rbind(cbind(C11,C12,C13),cbind(C21,C22,C23),cbind(t(C13),t(C23),C33))

### Checking nonnegative definiteness
chol(C2)

######################################################################
####### Simulating zero mean trivariate Gaussian random field ########
######################################################################
set.seed(123)
nsim=100
simt<-mvrnorm(n=nsim,mu=rep(0,times=3*(n^2)),Sigma=C2)

########################################
###### Plotting 1st realization ########
########################################
var1<-simt[1,1:(n^2)]
var2<-simt[1,(n^2+1):(2*(n^2))]
var3<-simt[1,(2*(n^2)+1):(3*(n^2))]
par(mfrow=c(1,3))
quilt.plot(x=x.s,y=y.s,var1,nx=n,ny=n,main="Variable 1")
quilt.plot(x=x.s,y=y.s,var2,nx=n,ny=n,main="Variable 2")
quilt.plot(x=x.s,y=y.s,var3,nx=n,ny=n,main="Variable 3")

#### Plotting covariance matrix #######
par(mfrow=c(1,1))

plotmatrix(C2)

######################################################
# mle_marginals_all is a function that takes arguments
# locations (x,y), data (z), n, and parameters (p)
# returns negative loglikelihood value when cross-correlations are ignored
# returns marginal parameters supplied and the corresponding marginal covariance matrices
#######################################################

mle_marginals_all<-function(p,z,x,y,n)
{
  a1<-a2<-a3<-p[1]
  nu1<-p[2]
  sigma1<-p[3]
  nu2<-p[4]
  sigma2<-p[5]
  nu3<-p[6]
  sigma3<-p[7]
  
  if(sum(p<0)!=0)
  {
    nloglikelihood<-10000000
    return(list(mlv=nloglikelihood,params=NULL))
  }
  else
  {
    x.s=x
    y.s=y
    grid<-cbind(x.s,y.s)
    dmat<-rdist(grid)
    C11<-my.matern(h=dmat,a=a1,sigma = sigma1,nu=nu1)
    C22<-my.matern(h=dmat,a=a2,sigma = sigma2,nu=nu2)
    C33<-my.matern(h=dmat,a=a3,sigma = sigma3,nu=nu3)
    zmat<-matrix(0,nrow=nrow(dmat),ncol=ncol(dmat))
    ############## Inverting C11 ##########
    
    C1chol<-chol(C11)
    C11inv<-chol2inv(C1chol)
    C11logdet<-determinant(C11)$modulus
    
    C2chol<-chol(C22)
    C22inv<-chol2inv(C2chol)
    C22logdet<-determinant(C22)$modulus
    
    C3chol<-chol(C33)
    C33inv<-chol2inv(C3chol)
    C33logdet<-determinant(C33)$modulus
    
    logD<-C11logdet+C22logdet+C33logdet
    Cinv<-rbind(cbind(C11inv,zmat,zmat),cbind(zmat,C22inv,zmat),cbind(zmat,zmat,C33inv))
    
    nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi)) 
    if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
    return(list(mlv=nloglikelihood,a1=a1,a2=a2,a3=a3,nu1=nu1,nu2=nu2,nu3=nu3,sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,C11=C11,C22=C22,C33=C33))
  }
  
}

###############################################################################
# mle_marginal_mlv is function that takes pars (parameters) as the arguement
# for the data z=sim, and locations x.s, y.s and n=n
# It returns the negative loglikelihood value
###############################################################################

mle_marginal_mlv<-function(pars)
{
  return(mle_marginals_all(p=pars,z=sim,x=x.s,y=y.s,n=n)$mlv)
}

########################################################################################
# optim_marg_loglik is an optimizer function that minimizes mle_marginal_mlv over pars
#########################################################################################


optim_marg_loglik <- function(par){
  optim(par=par,
        fn = mle_marginal_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}


#### Initializing value #####
marg.init<-c(6.8,1.3,1.1,2.5,1.4,2.8,1.3)

marginal_estimates<-foreach(i=1:nsim) %dopar% {
  library(fields)
  sim<-simt[i,]
  mle_marg<-optim_marg_loglik(par = marg.init)
  mle_marg
}
#testo<-foreach(i=1:50) %dopar% {
#  i
#}
# The marginal estimates are valid for the Symmetric and Asymmetric Parsimonious Matern models #########


###########################################################################################
## Now we will estimate the Marginal estimates for the Separable covariance function
###########################################################################################
# mle marginals_sep_all is a function that takes the argument parameters (p), data (z), location (x,y) and n
# returns negative loglikelihood when ignoring the cross-correlations, parameter values, and marginal covariance matrices
##############################################################################################################################

mle_marginals_sep_all<-function(p,z,x,y,n)
{
  
  a<-p[1]
  nu<-p[2]
  sigma1<-p[3]
  sigma2<-p[4]
  sigma3<-p[5]
  
  if(sum(p<0)!=0)
  {
    nloglikelihood<-10000000
    return(list(mlv=nloglikelihood,params=NULL))
  }
  else
  {
    x.s=x
    y.s=y
    grid<-cbind(x.s,y.s)
    dmat<-rdist(grid)
    C11<-my.matern(h=dmat,a=a,sigma = sigma1,nu=nu)
    C22<-my.matern(h=dmat,a=a,sigma = sigma2,nu=nu)
    C33<-my.matern(h=dmat,a=a,sigma = sigma3,nu=nu)
    zmat<-matrix(0,nrow=nrow(dmat),ncol=ncol(dmat))
    ############## Inverting C11 ##########
    
    C1chol<-chol(C11)
    C11inv<-chol2inv(C1chol)
    C11logdet<-determinant(C11)$modulus
    
    C2chol<-chol(C22)
    C22inv<-chol2inv(C2chol)
    C22logdet<-determinant(C22)$modulus
    
    C3chol<-chol(C33)
    C33inv<-chol2inv(C3chol)
    C33logdet<-determinant(C33)$modulus
    
    logD<-C11logdet+C22logdet+C33logdet
    Cinv<-rbind(cbind(C11inv,zmat,zmat),cbind(zmat,C22inv,zmat),cbind(zmat,zmat,C33inv))
    
    nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi)) 
    if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
    return(list(mlv=nloglikelihood,a=a,nu=nu,sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,C11=C11,C22=C22,C33=C33))
  }
}


##############################################################################################################
# mle_marginals_sep_mlv returns the negative loglikelihood value for parameters (pars) and data z=sim ########
##############################################################################################################
mle_marginals_sep_mlv<-function(pars)
{
  return(mle_marginals_sep_all(p=pars,z=sim,x=x.s,y=y.s,n=n)$mlv)
}
##################################################################################################
# optim_marg_sep_loglik is an optimizer function that minimizes mle_marginal_sep_mlv over pars
##################################################################################################


optim_marg_sep_loglik <- function(par){
  optim(par=par,
        fn = mle_marginals_sep_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}


#### Initializing value #####
marg.sep.init<-c(7,2,1.1,1.1,1.2)
marginal_sep_estimates<-foreach(i=1:nsim) %dopar% {
  library(fields)
  sim<-simt[i,]
  mle_marg<-optim_marg_sep_loglik(par = marg.sep.init)
  mle_marg
}


############################################################################
########### Now we extract the marginal parameters in a matrix #############
############################################################################


margins_par_mat<-matrix(NA,nrow=nsim,ncol=length(marg.init))
margins_sep_mat<-matrix(NA,nrow=nsim,ncol=length(marg.sep.init))

for(i in 1:nsim)
{
  margins_par_mat[i,]<-marginal_estimates[[i]]$par
  margins_sep_mat[i,]<-marginal_sep_estimates[[i]]$par
}

########## saving results till marginal estimates #######
save.image("Marginal_Estimates.RData")


############################################################################
######### Specification of margin matrices #################################
# margins_par_mat: column1=a, column2=nu1 column3=sigma1, column4=nu2, column4=sigma2, column6=nu3, column7=sigma3
# margins_par_mat: column1=a,column2=nu, column3=sigma1, column4=sigma2, column5=sigma3
############################################################################
############################################################################
getDoParWorkers()

############################################################################################################
####### Now we perform cross-parameters and asymmetry parameter estimation for our asymmetric model ########
############################################################################################################

####### parsimonious Matern model ########


############################################################################
## mle_my_asym_pars.mat_all is a function
## that takes argument cross-parameters (p), data (z)
## locations (x,y)
## marginal parameters (margins)
## returns negative loglikelihood value when cross-correlations are considered
## asymmetry and cross-parameters
#############################################################################

mle_my_asym_pars.mat_all<-function(p,z,x,y,n,margins,onlyCOV)
{
  
  beta_p<-p[1:3]                 #### beta paramterization
  a121<-p[4]                     #### asymmetry parameters
  a122<-p[5]
  a131<-p[6]
  a132<-p[7]
  a231<-p[8]
  a232<-p[9]
  a1<-a2<-a3<-a12<-a13<-a23<-margins[1]
  nu1<-margins[2]
  sigma1<-margins[3]
  nu2<-margins[4]
  sigma2<-margins[5]
  nu3<-margins[6]
  sigma3<-margins[7]
  ############ putting hard constraints #####
  nu12<-(nu1+nu2)/2
  nu13<-(nu1+nu3)/2
  nu23<-(nu2+nu3)/2
  
  
  
  beta_mat<-rho_par(beta_p)       ##### beta_(ij)'s #####
  
  
  ####rho_(ij})'s
  rho12<-beta_mat[1,2]*(sqrt(gamma(nu1+1))*sqrt(gamma(nu2+1))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*gamma((nu1+nu2)/2+1))
  rho13<-beta_mat[1,3]*(sqrt(gamma(nu1+1))*sqrt(gamma(nu3+1))*gamma((nu1+nu3)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu3))*gamma((nu1+nu3)/2+1))
  rho23<-beta_mat[2,3]*(sqrt(gamma(nu2+1))*sqrt(gamma(nu3+1))*gamma((nu2+nu3)/2))/(sqrt(gamma(nu2))*sqrt(gamma(nu3))*gamma((nu2+nu3)/2+1))
  test_value<-(beta_mat[1,2]^2+beta_mat[1,3]^2+beta_mat[2,3]^2-1)<=(-2*abs(beta_mat[1,2])*abs(beta_mat[1,3])*abs(beta_mat[2,3]))
  if(!test_value)
  {
    nlikelihood_value<-100000000
    return(list(mlv=nlikelihood_value,params=NULL))
  }
  else{
    x.s<-x
    y.s<-y
    dist.mat<-rdist(cbind(x.s,y.s))
    ###### Marginal covariances
    C11<-my.matern(h=dist.mat,a=a1,sigma = sigma1,nu=nu1)
    C22<-my.matern(h=dist.mat,a=a2,sigma = sigma2,nu=nu2)
    C33<-my.matern(h=dist.mat,a=a3,sigma = sigma3,nu=nu3)
    
    
    ####### Cross covariances
    
    
    grid<-cbind(x,y)
    tempgrid1<-cbind(rep(a121,times=n^2),rep(a122,times=n^2))
    
    shift12grid<-grid-tempgrid1
    shift.dist12<-rdist(grid,shift12grid)
    C12<-rho12*sigma1*sigma2*my.matern(h=shift.dist12,a=a12,sigma = 1,nu=nu12)
    C21<-t(C12)
    
    tempgrid2<-cbind(rep(a131,times=n^2),rep(a132,times=n^2))
    
    shift13grid<-grid-tempgrid2
    shift.dist13<-rdist(grid,shift13grid)
    C13<-rho13*sigma1*sigma3*my.matern(h=shift.dist13,a=a13,sigma = 1,nu=nu13)
    C31<-t(C13)
    tempgrid3<-cbind(rep(a231,times=n^2),rep(a232,times=n^2))
    
    shift23grid<-grid-tempgrid3
    shift.dist23<-rdist(grid,shift23grid)
    C23<-rho23*sigma2*sigma3*my.matern(h=shift.dist23,a=a23,sigma = 1,nu=nu23)
    C32<-t(C23)
    C<-rbind(cbind(C11,C12,C13),cbind(C21,C22,C23),cbind(C31,C32,C33))
    
    if(onlyCOV)
    {
      return(C)
    }
    #Computing log-likelihood
    else
    {
    cholYo <-chol(C)
    nlikelihood_value <-(0.5 * determinant(C)$modulus + 0.5 * t(z) %*% chol2inv(cholYo) %*% z+0.5*length(z)*log(2*pi)) 
    if(abs(nlikelihood_value) == Inf || is.nan(nlikelihood_value)){ nlikelihood_value <- 1e+08}
    return(list(mlv=nlikelihood_value,range=a1,nu1=nu1,nu2=nu2,nu3=nu3,nu12=nu12,nu13=nu13,nu23=nu23,a121=a121,a122=a122,a131=a131,a132=a132,a231=a231,a232=a232,rho12=rho12,rho13=rho13,rho23=rho23))
    }
  } 
}

#################################################################
# mle_my_asym_pars.mat_mlv is a function
# takes argument parameters (pars)
# return lnegative loglikelihood value for 
# data = sim, locations (x.s,y.s), marginal parameters = m.ests
#################################################################
mle_my_asym_pars.mat_mlv<-function(pars)
{
  return(mle_my_asym_pars.mat_all(p=pars,z=sim,x=x.s,y=y.s,n=n,margins=m.ests,onlyCOV = FALSE)$mlv)
}

######## intitalizing value #########
init.my.pars.mat.cross<-c(0.5,0.3,0.3,0,0,0,0,0,0)


###############################################################
# optim_loglik_my_asym_pars.mat is an optimizer for mle_my_asym_pars.mat_mlv over parameters (pars)
###############################################################
optim_loglik_my_asym_pars.mat <- function(par){
  optim(par=par,
        fn = mle_my_asym_pars.mat_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}

cross_my_asym_pars.mat_estimates<-foreach(i=1:nsim) %dopar% {
  library(fields)
  sim<-simt[i,]
  m.ests<-margins_par_mat[i,]
  mle_marg<-optim_loglik_my_asym_pars.mat(par = init.my.pars.mat.cross)
  mle_marg
}
#getDoParWorkers()

########### Separable Model ############


############################################################################
## mle_my_asym_sep.mat_all is a function
## that takes argument cross-parameters (p), data (z)
## locations (x,y)
## marginal parameters (margins)
## returns negative loglikelihood value when cross-correlations are considered
## asymmetry and cross-parameters
#############################################################################



mle_my_asym_sep.mat_all<-function(p,z,x,y,n,margins,onlyCOV)
{
  beta_p<-p[1:3]                 #### beta paramterization
  a121<-p[4]                     #### asymmetry parameters
  a122<-p[5]
  a131<-p[6]
  a132<-p[7]
  a231<-p[8]
  a232<-p[9]
  a<-margins[1]
  nu<-margins[2]
  sigma1<-margins[3]
  sigma2<-margins[4]
  sigma3<-margins[5]
  
  
  beta_mat<-rho_par(beta_p)       ##### beta_(ij)'s #####
  var_mat<-diag(c(sigma1,sigma2,sigma3),nrow=3)
  
  P_mat<-var_mat%*%beta_mat%*%var_mat
  ####rho_(ij})'s (Coherences)######
  rho12<-beta_mat[1,2]
  rho13<-beta_mat[1,3]  
  rho23<-beta_mat[2,3]  
  test_value<-(rho12^2+rho13^2+rho23^2-1)<=(-2*abs(rho12)*abs(rho13)*abs(rho23))
  if(!test_value)
  {
    nlikelihood_value<-100000000
    return(list(mlv=nlikelihood_value,params=NULL))
  }
  else{
    x.s<-x
    y.s<-y
    dist.mat<-rdist(cbind(x.s,y.s))
    ###### Marginal covariances
    C11<-P_mat[1,1]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
    C22<-P_mat[2,2]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
    C33<-P_mat[3,3]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
    
    
    ####### Cross covariances
    
    
    grid<-cbind(x,y)
    tempgrid1<-cbind(rep(a121,times=n^2),rep(a122,times=n^2))
    
    shift12grid<-grid-tempgrid1
    shift.dist12<-rdist(grid,shift12grid)
    C12<-P_mat[1,2]*my.matern(h=shift.dist12,a=a,sigma = 1,nu=nu)
    C21<-t(C12)
    
    tempgrid2<-cbind(rep(a131,times=n^2),rep(a132,times=n^2))
    
    shift13grid<-grid-tempgrid2
    shift.dist13<-rdist(grid,shift13grid)
    C13<-P_mat[1,3]*my.matern(h=shift.dist13,a=a,sigma = 1,nu=nu)
    C31<-t(C13)
    tempgrid3<-cbind(rep(a231,times=n^2),rep(a232,times=n^2))
    
    shift23grid<-grid-tempgrid3
    shift.dist23<-rdist(grid,shift23grid)
    C23<-P_mat[2,3]*my.matern(h=shift.dist23,a=a,sigma = 1,nu=nu)
    C32<-t(C23)
    C<-rbind(cbind(C11,C12,C13),cbind(C21,C22,C23),cbind(C31,C32,C33))
    
    if(onlyCOV)
    {
      return(C)
    }
    
    else{
      
    
    #Computing log-likelihood
    
    cholYo <-chol(C)
        
      nlikelihood_value <-(0.5 * determinant(C)$modulus + 0.5 * t(z) %*% chol2inv(cholYo) %*% z+0.5*length(z)*log(2*pi)) 
      if(abs(nlikelihood_value) == Inf || is.nan(nlikelihood_value)){ nlikelihood_value <- 1e+08}
      return(list(mlv=nlikelihood_value,a=a,nu=nu,a121=a121,a122=a122,a131=a131,a132=a132,a231=a231,a232=a232))
    }
  } 
  
}



#################################################################
# mle_my_asym_sep.mat_mlv is a function
# takes argument parameters (pars)
# return negative loglikelihood value for 
# data = sim, locations (x.s,y.s), marginal parameters = m.ests
#################################################################
mle_my_asym_sep.mat_mlv<-function(pars)
{
  return(mle_my_asym_sep.mat_all(p=pars,z=sim,x=x.s,y=y.s,n=n,margins=m.ests,onlyCOV = FALSE)$mlv)
}


######## intitalizing value #########
init.my.sep.mat.cross<-c(0.5,0.3,0.3,0,0,0,0,0,0)


###############################################################
# optim_loglik_my_asym_sep.mat is an optimizer for mle_my_asym_sep.mat_mlv over parameters (pars)
###############################################################
optim_loglik_my_asym_sep.mat <- function(par){
  optim(par=par,
        fn = mle_my_asym_sep.mat_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}

cross_my_asym_sep.mat_estimates<-foreach(i=1:nsim) %dopar% {
  library(fields)
  sim<-simt[i,]
  m.ests<-margins_sep_mat[i,]
  mle_marg<-optim_loglik_my_asym_sep.mat(par = init.my.sep.mat.cross)
  mle_marg
}


########### Saving till estimates of my asymmetric model ########


save.image("Myasymmetric_cross_Estimates.RData")
###########################################################
######### Extracting the cross parameters #################
###########################################################

####### Specification ######
# my_cross_par.mat_mat is matrix with 9 columns
# C1=rho12, c2=rho13, c3=rho23, c4=a121,c5=a122, c6=a131, c7=a132, c8= a231, c9=a232
# my_cross_sep.mat_mat is a matrix with 12 columns
# C1=P11, C2= P22, C3=P33, C4=P12, C5=P13, C6=P23, C7=a121, C8=122, C9=a131, C10=a132, C11=a231, C12=a232
# my asym_nloglikelihood_pars.mat contains minimized negative log-likelihood values for our asymmetric parsimonious matern model 
# my asym_nloglikelihood_sep.mat contains minimized negative log-likelihood values for our asymmetric separable matern model 
# my asym_aic_pars.mat contains AIC values for our asymmetric parsimonious matern model 
# my asym_aic_sep.mat contains AIC values for our asymmetric separable matern model 

my_cross_par.mat_mat<-matrix(NA,nrow = nsim,ncol=length(init.my.pars.mat.cross))
my_cross_sep.mat_mat<-matrix(NA,nrow = nsim,ncol=(length(init.my.sep.mat.cross)+3))
my_asym_nloglikelihood_pars.mat<-numeric(length = nsim)
my_asym_nloglikelihood_sep.mat<-numeric(length = nsim)
my_asym_aic_pars.mat<-numeric(length=nsim)
my_asym_aic_sep.mat<-numeric(length=nsim)

for(i in 1:nsim)
{
  parval1<-cross_my_asym_pars.mat_estimates[[i]]$par
  parval2<-cross_my_asym_sep.mat_estimates[[i]]$par
  my_asym_nloglikelihood_pars.mat[i]<-cross_my_asym_pars.mat_estimates[[i]]$value
  my_asym_nloglikelihood_sep.mat[i]<-cross_my_asym_sep.mat_estimates[[i]]$value
  my_asym_aic_pars.mat[i]<-2*(length(init.my.pars.mat.cross)+length(marg.init))+2*my_asym_nloglikelihood_pars.mat[i]
  my_asym_aic_sep.mat[i]<-2*(length(init.my.sep.mat.cross)+length(marg.sep.init))+2*my_asym_nloglikelihood_sep.mat[i]
  my_cross_par.mat_mat[i,4:9]<-parval1[4:9]
  my_cross_sep.mat_mat[i,7:12]<-parval2[4:9]
  temp1<-rho_par(parval1[1:3])
  myrho12<-temp1[1,2]
  myrho13<-temp1[1,3]
  myrho23<-temp1[2,3]
  my_cross_par.mat_mat[i,1:3]<-c(myrho12,myrho13,myrho23)
  ###### Now we create a P matrix #####
  temp2<-rho_par(parval2[1:3])
  margin_sigma<-margins_sep_mat[i,3:5]
  dtemp<-diag(margin_sigma,nrow = 3)
  temp3<-dtemp%*%temp2%*%dtemp
  P11<-temp3[1,1]
  P22<-temp3[2,2]
  P33<-temp3[3,3]
  P12<-temp3[1,2]
  P13<-temp3[1,3]
  P23<-temp3[2,3]
  my_cross_sep.mat_mat[i,1:6]<-c(P11,P22,P33,P12,P13,P23)
  
}
par(mfrow=c(1,2))
#########################################################################################
########## Now we do the estimation for the Li-zhangs model #############################
#########################################################################################

############# Parsimonious matern model ############

############################################################################
## mle_liz_asym_pars.mat_all is a function
## that takes argument cross-parameters (p), data (z)
## locations (x,y)
## marginal parameters (margins)
## returns negative loglikelihood value when cross-correlations are considered
## asymmetry and cross-parameters
#############################################################################


mle_liz_asym_pars.mat_all<-function(p,z,x,y,n,margins,onlyCOV)
{
  beta_p<-p[1:3]                 #### beta paramterization
  
  
  as11<-0                     #### asymmetry parameters
  as12<-0
  as21<-p[4]
  as22<-p[5]
  as31<-p[6]
  as32<-p[7]
  
  
  #setting marginal parameters#
  
  a1<-a2<-a3<-a12<-a13<-a23<-margins[1]
  nu1<-margins[2]
  sigma1<-margins[3]
  nu2<-margins[4]
  sigma2<-margins[5]
  nu3<-margins[6]
  sigma3<-margins[7]
  ############ putting hard constraints #####
  nu12<-(nu1+nu2)/2
  nu13<-(nu1+nu3)/2
  nu23<-(nu2+nu3)/2
  
  
  
  beta_mat<-rho_par(beta_p)       ##### beta_(ij)'s
  
  
  ####rho_(ij})'s
  rho12<-beta_mat[1,2]*(sqrt(gamma(nu1+1))*sqrt(gamma(nu2+1))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*gamma((nu1+nu2)/2+1))
  rho13<-beta_mat[1,3]*(sqrt(gamma(nu1+1))*sqrt(gamma(nu3+1))*gamma((nu1+nu3)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu3))*gamma((nu1+nu3)/2+1))
  rho23<-beta_mat[2,3]*(sqrt(gamma(nu2+1))*sqrt(gamma(nu3+1))*gamma((nu2+nu3)/2))/(sqrt(gamma(nu2))*sqrt(gamma(nu3))*gamma((nu2+nu3)/2+1))
  
  
  x.s<-x
  y.s<-y
  dist.mat<-rdist(cbind(x.s,y.s))
  ###### Marginal covariances
  C11<-my.matern(h=dist.mat,a=a1,sigma = sigma1,nu=nu1)
  C22<-my.matern(h=dist.mat,a=a2,sigma = sigma2,nu=nu2)
  C33<-my.matern(h=dist.mat,a=a3,sigma = sigma3,nu=nu3)
  
  
  
  ####### Cross covariances
  
  
  grid<-cbind(x,y)
  tempgrid1<-cbind(rep(as11,times=n^2),rep(as12,times=n^2))
  tempgrid2<-cbind(rep(as21,times=n^2),rep(as22,times=n^2))
  tempgrid3<-cbind(rep(as31,times=n^2),rep(as32,times=n^2))
  
  shift1<-grid-tempgrid1
  shift2<-grid-tempgrid2
  shift3<-grid-tempgrid3
  shift.dist12<-rdist(shift1,shift2)
  shift.dist13<-rdist(shift1,shift3)
  shift.dist23<-rdist(shift2,shift3)
  C12<-rho12*sigma1*sigma2*my.matern(h=shift.dist12,a=a12,sigma = 1,nu=nu12)
  C13<-rho13*sigma1*sigma3*my.matern(h=shift.dist13,a=a13,sigma = 1,nu=nu13)
  C23<-rho23*sigma2*sigma3*my.matern(h=shift.dist23,a=a23,sigma = 1,nu=nu23)
  C<-rbind(cbind(C11,C12,C13),cbind(t(C12),C22,C23),cbind(t(C13),t(C23),C33))
  
  if(onlyCOV)
  {
    return(C)
  }
  #Computing log-likelihood
  else{
    
  cholYo <-chol(C)
  
  
    nlikelihood_value <-(0.5 * determinant(C)$modulus + 0.5 * t(z) %*% chol2inv(cholYo) %*% z+0.5*length(z)*log(2*pi)) 
    if(abs(nlikelihood_value) == Inf || is.nan(nlikelihood_value)){ nlikelihood_value <- 1e+08}
    return(list(mlv=nlikelihood_value,a=a1,nu1=nu1,nu2=nu2,nu3=nu3,nu12=nu12,nu13=nu13,nu23=nu23,as21=as21,as22=as22,as31=as31,as32=as32,rho12=rho12,rho13=rho13,rho23=rho23))
  }
} 




#################################################################
# mle_liz_asym_pars.mat_mlv is a function
# takes argument parameters (pars)
# return negative loglikelihood value for 
# data = sim, locations (x.s,y.s), marginal parameters = m.ests
#################################################################



mle_liz_asym_pars.mat_mlv<-function(pars)
{
  return(mle_liz_asym_pars.mat_all(p=pars,z=sim,x=x.s,y=y.s,n=n,margins=m.ests,onlyCOV = FALSE)$mlv)
}


#########################################################################################################
# optim_loglik_liz_asym_pars.mat is an optimizer for mle_liz_asym_pars.mat_mlv over parameters (pars)####
#########################################################################################################


optim_loglik_liz_asym_pars.mat <- function(par){
  optim(par=par,
        fn = mle_liz_asym_pars.mat_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}

init.liz.pars.mat.cross<-c(0.5,0.3,0.3,0,0,0,0)


cross_liz_asym_pars.mat_estimates<-foreach(i=1:nsim)%dopar%{
  library(fields)
  sim<-simt[i,]
  m.ests<-margins_par_mat[i,]
  mle_marg<-optim_loglik_liz_asym_pars.mat(par=init.liz.pars.mat.cross)
  mle_marg
}

####################### Separable model #################


############################################################################
## mle_liz_asym_sep.mat_all is a function
## that takes argument cross-parameters (p), data (z)
## locations (x,y)
## marginal parameters (margins)
## returns negative loglikelihood value when cross-correlations are considered
## asymmetry and cross-parameters
#############################################################################



mle_liz_asym_sep.mat_all<-function(p,z,x,y,n,margins,onlyCOV)
{
  beta_p<-p[1:3]                 #### beta paramterization
  as11<-0
  as12<-0
  as21<-p[4]                     #### asymmetry parameters
  as22<-p[5]
  as31<-p[6]
  as32<-p[7]
  a<-margins[1]
  nu<-margins[2]
  sigma1<-margins[3]
  sigma2<-margins[4]
  sigma3<-margins[5]
  
  
  beta_mat<-rho_par(beta_p)       ##### beta_(ij)'s #####
  var_mat<-diag(c(sigma1,sigma2,sigma3),nrow=3)
  
  P_mat<-var_mat%*%beta_mat%*%var_mat
  
    x.s<-x
    y.s<-y
    dist.mat<-rdist(cbind(x.s,y.s))
    ###### Marginal covariances
    C11<-P_mat[1,1]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
    C22<-P_mat[2,2]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
    C33<-P_mat[3,3]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
    
    
    ####### Cross covariances
    
    grid<-cbind(x,y)
    tempgrid1<-cbind(rep(as11,times=n^2),rep(as12,times=n^2))
    tempgrid2<-cbind(rep(as21,times=n^2),rep(as22,times=n^2))
    tempgrid3<-cbind(rep(as31,times=n^2),rep(as32,times=n^2))
    
    shift1<-grid-tempgrid1
    shift2<-grid-tempgrid2
    shift3<-grid-tempgrid3
    shift.dist12<-rdist(shift1,shift2)
    shift.dist13<-rdist(shift1,shift3)
    shift.dist23<-rdist(shift2,shift3)
    C12<-P_mat[1,2]*my.matern(h=shift.dist12,a=a,sigma = 1,nu=nu)
    C13<-P_mat[1,3]*my.matern(h=shift.dist13,a=a,sigma = 1,nu=nu)
    C23<-P_mat[2,3]*my.matern(h=shift.dist23,a=a,sigma = 1,nu=nu)
    C<-rbind(cbind(C11,C12,C13),cbind(t(C12),C22,C23),cbind(t(C13),t(C23),C33))
    
    if(onlyCOV)
    {
      return(C)
    }
    #Computing log-likelihood
    else{
    cholYo <-chol(C)
        
      nlikelihood_value <-(0.5 * determinant(C)$modulus + 0.5 * t(z) %*% chol2inv(cholYo) %*% z+0.5*length(z)*log(2*pi)) 
      if(abs(nlikelihood_value) == Inf || is.nan(nlikelihood_value)){ nlikelihood_value <- 1e+08}
      return(list(mlv=nlikelihood_value,range=a,nu=nu,as21=as21,as22=as22,as31=as31,as32=as32,P11=P_mat[1,1],P22=P_mat[2,2],P33=P_mat[3,3],P13=P_mat[1,3],P23=P_mat[2,3],P12=P_mat[1,2]))
    }
}



#################################################################
# mle_liz_asym_sep.mat_mlv is a function
# takes argument parameters (pars)
# return negative loglikelihood value for 
# data = sim, locations (x.s,y.s), marginal parameters = m.ests
#################################################################
mle_liz_asym_sep.mat_mlv<-function(pars)
{
  return(mle_liz_asym_sep.mat_all(p=pars,z=sim,x=x.s,y=y.s,n=n,margins=m.ests,onlyCOV = FALSE)$mlv)
}

######## intitalizing value #########
init.liz.sep.mat.cross<-c(0.5,0.3,0.3,0,0,0,0) 


###############################################################
# optim_loglik_my_asym_sep.mat is an optimizer for mle_my_asym_sep.mat_mlv over parameters (pars)
###############################################################
optim_loglik_liz_asym_sep.mat <- function(par){
  optim(par=par,
        fn = mle_liz_asym_sep.mat_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}

cross_liz_asym_sep.mat_estimates<-foreach(i=1:nsim) %dopar% {
  library(fields)
  sim<-simt[i,]
  m.ests<-margins_sep_mat[i,]
  mle_marg<-optim_loglik_liz_asym_sep.mat(par = init.liz.sep.mat.cross)
  mle_marg
}


###########################################################
######### Extracting the cross parameters #################
###########################################################

####### Specification ######
# liz_cross_par.mat_mat is matrix with 7 columns
# C1=rho12, c2=rho13, c3=rho23, c4=as21,c5=as22, c6=as31, c7=as32, 
# liz_cross_sep.mat_mat is a matrix with 10 columns
# C1=P11, C2= P22, C3=P33, C4=P12, C5=P13, C6=P23, C7=as21, C8=as22, C9=as31, C10=as32
# liz_asym_nloglikelihood_pars.mat contains minimized negative log-likelihood values for our asymmetric parsimonious matern model 
# liz_asym_nloglikelihood_sep.mat contains minimized negative log-likelihood values for our asymmetric separable matern model 
# liz_asym_aic_pars.mat contains AIC values for our asymmetric parsimonious matern model 
# liz_asym_aic_sep.mat contains AIC values for our asymmetric separable matern model 

liz_cross_par.mat_mat<-matrix(NA,nrow = nsim,ncol=length(init.liz.pars.mat.cross))
liz_cross_sep.mat_mat<-matrix(NA,nrow = nsim,ncol=(length(init.liz.sep.mat.cross)+3))
liz_asym_nloglikelihood_pars.mat<-numeric(length = nsim)
liz_asym_nloglikelihood_sep.mat<-numeric(length = nsim)
liz_asym_aic_pars.mat<-numeric(length=nsim)
liz_asym_aic_sep.mat<-numeric(length=nsim)

for(i in 1:nsim)
{
  parval1<-cross_liz_asym_pars.mat_estimates[[i]]$par
  parval2<-cross_liz_asym_sep.mat_estimates[[i]]$par
  liz_asym_nloglikelihood_pars.mat[i]<-cross_liz_asym_pars.mat_estimates[[i]]$value
  liz_asym_nloglikelihood_sep.mat[i]<-cross_liz_asym_sep.mat_estimates[[i]]$value
  liz_asym_aic_pars.mat[i]<-2*(length(init.liz.pars.mat.cross)+length(marg.init))+2*liz_asym_nloglikelihood_pars.mat[i]
  liz_asym_aic_sep.mat[i]<-2*(length(init.liz.sep.mat.cross)+length(marg.sep.init))+2*liz_asym_nloglikelihood_sep.mat[i]
  liz_cross_par.mat_mat[i,4:7]<-parval1[4:7]
  liz_cross_sep.mat_mat[i,7:10]<-parval2[4:7]
  temp1<-rho_par(parval1[1:3])
  lizrho12<-temp1[1,2]
  lizrho13<-temp1[1,3]
  lizrho23<-temp1[2,3]
  liz_cross_par.mat_mat[i,1:3]<-c(lizrho12,lizrho13,lizrho23)
  ###### Now we create a P matrix #####
  temp2<-rho_par(parval2[1:3])
  margin_sigma<-margins_sep_mat[i,3:5]
  dtemp<-diag(margin_sigma,nrow = 3)
  temp3<-dtemp%*%temp2%*%dtemp
  P11<-temp3[1,1]
  P22<-temp3[2,2]
  P33<-temp3[3,3]
  P12<-temp3[1,2]
  P13<-temp3[1,3]
  P23<-temp3[2,3]
  liz_cross_sep.mat_mat[i,1:6]<-c(P11,P22,P33,P12,P13,P23)
  
}

save.image("Lizasymmetry_cross_Estimates.RData")



########################################################
###### Now we do symmetric models ######################
########################################################



############# Parsimonious matern model ####################

############################################################################
## mle_sym_pars.mat_all is a function
## that takes argument cross-parameters (p), data (z)
## locations (x,y)
## marginal parameters (margins)
## returns negative loglikelihood value when cross-correlations are considered
## cross-parameters
#############################################################################


mle_sym_pars.mat_all<-function(p,z,x,y,n,margins,onlyCOV)
{
  
  
  
  beta_p<-p[1:3]                 #### beta paramterization
  #setting marginal parameters#
  
  #setting marginal parameters#
  
  a1<-a2<-a3<-a12<-a13<-a23<-margins[1]
  nu1<-margins[2]
  sigma1<-margins[3]
  nu2<-margins[4]
  sigma2<-margins[5]
  nu3<-margins[6]
  sigma3<-margins[7]
  ############ putting hard constraints #####
  nu12<-(nu1+nu2)/2
  nu13<-(nu1+nu3)/2
  nu23<-(nu2+nu3)/2
  
  
  
  beta_mat<-rho_par(beta_p)       ##### beta_(ij)'s
  
  
  ####rho_(ij})'s
  rho12<-beta_mat[1,2]*(sqrt(gamma(nu1+1))*sqrt(gamma(nu2+1))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*gamma((nu1+nu2)/2+1))
  rho13<-beta_mat[1,3]*(sqrt(gamma(nu1+1))*sqrt(gamma(nu3+1))*gamma((nu1+nu3)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu3))*gamma((nu1+nu3)/2+1))
  rho23<-beta_mat[2,3]*(sqrt(gamma(nu2+1))*sqrt(gamma(nu3+1))*gamma((nu2+nu3)/2))/(sqrt(gamma(nu2))*sqrt(gamma(nu3))*gamma((nu2+nu3)/2+1))
  
  
  x.s<-x
  y.s<-y
  dist.mat<-rdist(cbind(x.s,y.s))
  ###### Marginal covariances
  C11<-my.matern(h=dist.mat,a=a1,sigma = sigma1,nu=nu1)
  C22<-my.matern(h=dist.mat,a=a2,sigma = sigma2,nu=nu2)
  C33<-my.matern(h=dist.mat,a=a3,sigma = sigma3,nu=nu3)
  
  
  
  ####### Cross covariances
  
  C12<-rho12*sigma1*sigma2*my.matern(h=dist.mat,a=a12,sigma = 1,nu=nu12)
  C13<-rho13*sigma1*sigma3*my.matern(h=dist.mat,a=a13,sigma = 1,nu=nu13)
  C23<-rho23*sigma2*sigma3*my.matern(h=dist.mat,a=a23,sigma = 1,nu=nu23)
  C<-rbind(cbind(C11,C12,C13),cbind(t(C12),C22,C23),cbind(t(C13),t(C23),C33))
  
  if(onlyCOV)
  {
    return(C)
  }
  #Computing log-likelihood
  else{
  cholYo <-chol(C)
  
  
    nlikelihood_value <-(0.5 * determinant(C)$modulus + 0.5 * t(z) %*% chol2inv(cholYo) %*% z+0.5*length(z)*log(2*pi)) 
    if(abs(nlikelihood_value) == Inf || is.nan(nlikelihood_value)){ nlikelihood_value <- 1e+08}
    return(list(mlv=nlikelihood_value,a=a1,nu1=nu1,nu2=nu2,nu3=nu3,nu12=nu12,nu13=nu13,nu23=nu23,rho12=rho12,rho13=rho13,rho23=rho23))
  } 
} 
#################################################################
# mle_sym_pars.mat_mlv is a function
# takes argument parameters (pars)
# return negative loglikelihood value for 
# data = sim, locations (x.s,y.s), marginal parameters = m.ests
#################################################################


mle_sym_pars.mat_mlv<-function(pars)
{
  return(mle_sym_pars.mat_all(p=pars,z=sim,x=x.s,y=y.s,n=n,margins=m.ests,onlyCOV = FALSE)$mlv)
}

###############################################################
# optim_loglik_sym_pars.mat is an optimizer for mle_my_asym_sep.mat_mlv over parameters (pars)
###############################################################

optim_loglik_sym_pars.mat <- function(par){
  optim(par=par,
        fn = mle_sym_pars.mat_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}


init.sym.pars.mat.cross<-c(0.5,0.3,0.3)


cross_sym_pars.mat_estimates<-foreach(i=1:nsim) %dopar%{
  library(fields)
  
  sim<-simt[i,]
  m.ests<-margins_par_mat[i,]
  mle_marg<-optim_loglik_sym_pars.mat(init.sym.pars.mat.cross)
  mle_marg
  
}


##################### Separable model #######################




############################################################################
## mle_sym_sep.mat_all is a function
## that takes argument cross-parameters (p), data (z)
## locations (x,y)
## marginal parameters (margins)
## returns negative loglikelihood value when cross-correlations are considered
## and cross-parameters
#############################################################################



mle_sym_sep.mat_all<-function(p,z,x,y,n,margins,onlyCOV)
{
  beta_p<-p[1:3]                 #### beta paramterization
  a<-margins[1]
  nu<-margins[2]
  sigma1<-margins[3]
  sigma2<-margins[4]
  sigma3<-margins[5]
  
  
  beta_mat<-rho_par(beta_p)       ##### beta_(ij)'s #####
  var_mat<-diag(c(sigma1,sigma2,sigma3),nrow=3)
  
  P_mat<-var_mat%*%beta_mat%*%var_mat
  
  x.s<-x
  y.s<-y
  dist.mat<-rdist(cbind(x.s,y.s))
  ###### Marginal covariances
  C11<-P_mat[1,1]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
  C22<-P_mat[2,2]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
  C33<-P_mat[3,3]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
  
  
  ####### Cross covariances
  C12<-P_mat[1,2]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
  C13<-P_mat[1,3]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
  C23<-P_mat[2,3]*my.matern(h=dist.mat,a=a,sigma = 1,nu=nu)
  C<-rbind(cbind(C11,C12,C13),cbind(t(C12),C22,C23),cbind(t(C13),t(C23),C33))
  
  if(onlyCOV)
  {
    return(C)
  }
  #Computing log-likelihood
  else
  {
  cholYo <-chol(C)
  
  
    nlikelihood_value <-(0.5 * determinant(C)$modulus + 0.5 * t(z) %*% chol2inv(cholYo) %*% z+0.5*length(z)*log(2*pi)) 
    if(abs(nlikelihood_value) == Inf || is.nan(nlikelihood_value)){ nlikelihood_value <- 1e+08}
    return(list(mlv=nlikelihood_value,range=a,nu=nu,P11=P_mat[1,1],P22=P_mat[2,2],P33=P_mat[3,3],P13=P_mat[1,3],P23=P_mat[2,3],P12=P_mat[1,2]))
  }
}



#################################################################
# mle_sym_sep.mat_mlv is a function
# takes argument parameters (pars)
# return negative loglikelihood value for 
# data = sim, locations (x.s,y.s), marginal parameters = m.ests
#################################################################
mle_sym_sep.mat_mlv<-function(pars)
{
  return(mle_sym_sep.mat_all(p=pars,z=sim,x=x.s,y=y.s,n=n,margins=m.ests,onlyCOV = FALSE)$mlv)
}

######## intitalizing value #########
init.sym.sep.mat.cross<-c(0.5,0.3,0.3) 


###############################################################
# optim_loglik_my_asym_sep.mat is an optimizer for mle_my_asym_sep.mat_mlv over parameters (pars)
###############################################################
optim_loglik_sym_sep.mat <- function(par){
  optim(par=par,
        fn = mle_sym_sep.mat_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}

cross_sym_sep.mat_estimates<-foreach(i=1:nsim) %dopar% {
  library(fields)
  sim<-simt[i,]
  m.ests<-margins_sep_mat[i,]
  mle_marg<-optim_loglik_sym_sep.mat(par = init.sym.sep.mat.cross)
  mle_marg
}


###########################################################
######### Extracting the cross parameters #################
###########################################################

####### Specification ######
# sym_cross_par.mat_mat is matrix with 3 columns
# C1=rho12, c2=rho13, c3=rho23, 
# sym_cross_sep.mat_mat is a matrix with 6 columns
# C1=P11, C2= P22, C3=P33, C4=P12, C5=P13, C6=P23,
# sym_nloglikelihood_pars.mat contains minimized negative log-likelihood values for our symmetric parsimonious matern model 
# sym_nloglikelihood_sep.mat contains minimized negative log-likelihood values for our symmetric separable matern model 
# sym_aic_pars.mat contains AIC values for our symmetric parsimonious matern model 
# sym_aic_sep.mat contains AIC values for our symmetric separable matern model 

sym_cross_par.mat_mat<-matrix(NA,nrow = nsim,ncol=length(init.sym.pars.mat.cross))
sym_cross_sep.mat_mat<-matrix(NA,nrow = nsim,ncol=(length(init.sym.sep.mat.cross)+3))
sym_nloglikelihood_pars.mat<-numeric(length = nsim)
sym_nloglikelihood_sep.mat<-numeric(length = nsim)
sym_aic_pars.mat<-numeric(length=nsim)
sym_aic_sep.mat<-numeric(length=nsim)

for(i in 1:nsim)
{
  parval1<-cross_sym_pars.mat_estimates[[i]]$par
  parval2<-cross_sym_sep.mat_estimates[[i]]$par
  sym_nloglikelihood_pars.mat[i]<-cross_sym_pars.mat_estimates[[i]]$value
  sym_nloglikelihood_sep.mat[i]<-cross_sym_sep.mat_estimates[[i]]$value
  sym_aic_pars.mat[i]<-2*(length(init.sym.pars.mat.cross)+length(marg.init))+2*sym_nloglikelihood_pars.mat[i]
  sym_aic_sep.mat[i]<-2*(length(init.sym.sep.mat.cross)+length(marg.sep.init))+2*sym_nloglikelihood_sep.mat[i]
  temp1<-rho_par(parval1[1:3])
  symrho12<-temp1[1,2]
  symrho13<-temp1[1,3]
  symrho23<-temp1[2,3]
  sym_cross_par.mat_mat[i,1:3]<-c(symrho12,symrho13,symrho23)
  ###### Now we create a P matrix #####
  temp2<-rho_par(parval2[1:3])
  margin_sigma<-margins_sep_mat[i,3:5]
  dtemp<-diag(margin_sigma,nrow = 3)
  temp3<-dtemp%*%temp2%*%dtemp
  P11<-temp3[1,1]
  P22<-temp3[2,2]
  P33<-temp3[3,3]
  P12<-temp3[1,2]
  P13<-temp3[1,3]
  P23<-temp3[2,3]
  sym_cross_sep.mat_mat[i,1:6]<-c(P11,P22,P33,P12,P13,P23)
  
}


par(mfrow=c(1,2))
boxplot(my_asym_nloglikelihood_pars.mat,my_asym_nloglikelihood_sep.mat,liz_asym_nloglikelihood_pars.mat,liz_asym_nloglikelihood_sep.mat,sym_nloglikelihood_pars.mat,sym_nloglikelihood_sep.mat)
boxplot(my_asym_aic_pars.mat,my_asym_aic_sep.mat,liz_asym_aic_pars.mat,liz_asym_aic_sep.mat,sym_aic_pars.mat,sym_aic_sep.mat)

save.image("Symmetric_cross_Estimates.RData")

library(ggplot2)
nloglikdata<-data.frame(nloglikelihood=c(my_asym_nloglikelihood_pars.mat,my_asym_nloglikelihood_sep.mat,liz_asym_nloglikelihood_pars.mat,liz_asym_nloglikelihood_sep.mat,sym_nloglikelihood_pars.mat,sym_nloglikelihood_sep.mat),
                        Model=c(rep("Our asymmetric pars. Matern",times=nsim),rep("Our asymmetric intrinsic Matern",times=nsim),rep("Li-zhangs asymmetric pars. Matern",times=nsim),rep("Li-zhangs asymmetric intrinsic Matern",times=nsim),rep("Symmetric pars. Matern",times=nsim),rep("Symmetric intrinsic Matern",times=nsim)))

bpdata<-data.frame(AIC=c(my_asym_aic_pars.mat,my_asym_aic_sep.mat,liz_asym_aic_pars.mat,liz_asym_aic_sep.mat,sym_aic_pars.mat,sym_aic_sep.mat),
                   Model=c(rep("Our asymmetric pars. Matern",times=nsim),rep("Our asymmetric intrinsic Matern",times=nsim),rep("Li-zhangs asymmetric pars. Matern",times=nsim),rep("Li-zhangs asymmetric intrinsic Matern",times=nsim),rep("Symmetric pars. Matern",times=nsim),rep("Symmetric intrinsic Matern",times=nsim)))

par(mfrow=c(1,2))
bp<-ggplot(data = bpdata, aes(x=Model, y=AIC)) + geom_boxplot(aes(fill=Model))
plot(bp)
lk<-ggplot(data = nloglikdata, aes(x=Model, y=nloglikelihood)) + geom_boxplot(aes(fill=Model))
plot(lk)


#####################################################################################
############ Now we do the comparison of the prediction performance #################
#####################################################################################
# pred_performance is a function that take argument
# Estimated covariance matrix and corresponding simulated sample
# returns
# MSPE, mLogS
# Predicted values and Prediction variance
#####################################################################################



pred_performance<-function(C,sample)
{
  #### Using Sherman Morrison updating formula from the article 
  ### Relationship between the Inverses of a Matrix and a Submatrix by
  ####E. Ju?rez-Ruiz1
  ####, R. Cort?s-Maldonado2
  ####, F. P?rez-Rodr?guez2
  
  A<-C
  Ainv<-solve(C)
  s.size<-length(sample)
  pred<-pred.var<-numeric(length=s.size)
  for(i in 1:s.size)
  {
    B<-A[-i,-i]
    a<-A[,i]
    e<-rep(0,times=nrow(A))
    e[i]<-1
    u<-a-e
    v<-e
    U<-matrix(u,ncol = 1)
    V<-matrix(v,ncol=1)
    subt<-U%*%t(V)
    num<-(Ainv%*%U)%*%(t(V)%*%Ainv)
    den<-(1-t(V)%*%Ainv%*%U)
    SM_Mat<-Ainv+(c(1/den))*(num)
    
    sigma22inv<-SM_Mat[-i,-i]
    #sigma22<-C[-i,-i]
    z<-matrix(sample[-i],nrow=s.size-1,ncol=1)
    sigma12<-matrix(c(C[i,-i]),nrow=1)
    #sigma22inv<-solve(sigma22)
    weights<-sigma12%*%sigma22inv
    pred[i]<-weights%*%z
    pred.var[i]<-C[i,i]-weights%*%t(sigma12)
    
  }
  
  MSPE<-mean((sample-pred)^2)
  Log_Score<-mean(log(2*pi*pred.var)+(((sample-pred)^2)/pred.var))
  
  return(list(kriged.values=pred,kriged.variance=pred.var,Observed=sample,MSPE=MSPE,Log_Score=Log_Score))  
}



#######################################################################################################
###### Now we compute the #nsim estimated covariance matrix for each of the candidate models ##########
#######################################################################################################
# est.cov is the array of dimension nsim,3*n^2,3*n^2
est.cov1<-array(dim=c(nsim,3*(n^2),3*(n^2)))

for(i in 1:nsim)
{
  tempsim<-simt[i,]
  est.cov1[i,,]<-mle_my_asym_pars.mat_all(p=cross_my_asym_pars.mat_estimates[[i]]$par,z=tempsim,x=x.s,y=y.s,n=n,margins = margins_par_mat[i,],onlyCOV = TRUE)
  rm(tempsim)
  
}
getDoParWorkers()
my_asym_pars.mat_pred_perf<-foreach(i=1:nsim) %dopar%{
  tsample<-simt[i,]
  pred_performance(C=est.cov1[i,,],sample = tsample)
}
rm(est.cov1)
save.image("pred1.RData")



est.cov2<-array(dim=c(nsim,3*(n^2),3*(n^2)))

for(i in 1:nsim)
{
  tempsim<-simt[i,]
  est.cov2[i,,]<-mle_my_asym_sep.mat_all(p=cross_my_asym_sep.mat_estimates[[i]]$par,z=tempsim,x=x.s,y=y.s,n=n,margins = margins_sep_mat[i,],onlyCOV = TRUE)
    rm(tempsim)
}


my_asym_sep.mat_pred_perf<-foreach(i=1:nsim) %dopar%{
  tsample<-simt[i,]
  pred_performance(C=est.cov2[i,,],sample = tsample)
}
rm(est.cov2)
save.image("pred2.RData")


est.cov3<-array(dim=c(nsim,3*(n^2),3*(n^2)))

for(i in 1:nsim)
{
  tempsim<-simt[i,]
  est.cov3[i,,]<-mle_liz_asym_pars.mat_all(p=cross_liz_asym_pars.mat_estimates[[i]]$par,z=tempsim,x=x.s,y=y.s,n=n,margins = margins_par_mat[i,],onlyCOV = TRUE)
  rm(tempsim)
  
}

liz_asym_pars.mat_pred_perf<-foreach(i=1:nsim) %dopar%{
  tsample<-simt[i,]
  pred_performance(C=est.cov3[i,,],sample = tsample)
}
rm(est.cov3)
save.image("pred3.RData")


est.cov4<-array(dim=c(nsim,3*(n^2),3*(n^2)))

for(i in 1:nsim)
{
  tempsim<-simt[i,]
  est.cov4[i,,]<-mle_liz_asym_sep.mat_all(p=cross_liz_asym_sep.mat_estimates[[i]]$par,z=tempsim,x=x.s,y=y.s,n=n,margins = margins_sep_mat[i,],onlyCOV = TRUE)
  rm(tempsim)
  
}

liz_asym_sep.mat_pred_perf<-foreach(i=1:nsim) %dopar%{
  tsample<-simt[i,]
  pred_performance(C=est.cov4[i,,],sample = tsample)
}
rm(est.cov4)
save.image("pred4.RData")

est.cov5<-array(dim=c(nsim,3*(n^2),3*(n^2)))

for(i in 1:nsim)
{
  tempsim<-simt[i,]
  est.cov5[i,,]<-mle_sym_pars.mat_all(p=cross_sym_pars.mat_estimates[[i]]$par,z=tempsim,x=x.s,y=y.s,n=n,margins = margins_par_mat[i,],onlyCOV = TRUE)
  rm(tempsim)
  
}


sym_pars.mat_pred_perf<-foreach(i=1:nsim) %dopar%{
  tsample<-simt[i,]
  pred_performance(C=est.cov5[i,,],sample = tsample)
}
rm(est.cov5)
save.image("pred5.RData")

est.cov6<-array(dim=c(nsim,3*(n^2),3*(n^2)))

for(i in 1:nsim)
{
  tempsim<-simt[i,]
  est.cov6[i,,]<-mle_sym_sep.mat_all(p=cross_sym_sep.mat_estimates[[i]]$par,z=tempsim,x=x.s,y=y.s,n=n,margins = margins_sep_mat[i,],onlyCOV = TRUE)
  rm(tempsim)
  
}


sym_sep.mat_pred_perf<-foreach(i=1:nsim) %dopar%{
  tsample<-simt[i,]
  pred_performance(C=est.cov6[i,,],sample = tsample)
}
rm(est.cov6)
save.image("pred6.RData")




sym_pars.mat_mspe<-sym_sep.mat_mspe<-sym_pars.mat_logs<-sym_sep.mat_logs<-my_asym_pars.mat_mspe<-my_asym_sep.mat_mspe<-my_asym_pars.mat_logs<-my_asym_sep.mat_logs<-liz_asym_pars.mat_mspe<-liz_asym_sep.mat_mspe<-liz_asym_pars.mat_logs<-liz_asym_sep.mat_logs<-numeric(length=nsim)
for(i in 1:nsim)
{
  my_asym_pars.mat_mspe[i]<-my_asym_pars.mat_pred_perf[[i]]$MSPE
  my_asym_pars.mat_logs[i]<-my_asym_pars.mat_pred_perf[[i]]$Log_Score
  my_asym_sep.mat_mspe[i]<-my_asym_sep.mat_pred_perf[[i]]$MSPE
  my_asym_sep.mat_logs[i]<-my_asym_sep.mat_pred_perf[[i]]$Log_Score
  
  liz_asym_pars.mat_mspe[i]<-liz_asym_pars.mat_pred_perf[[i]]$MSPE
  liz_asym_pars.mat_logs[i]<-liz_asym_pars.mat_pred_perf[[i]]$Log_Score
  liz_asym_sep.mat_mspe[i]<-liz_asym_sep.mat_pred_perf[[i]]$MSPE
  liz_asym_sep.mat_logs[i]<-liz_asym_sep.mat_pred_perf[[i]]$Log_Score
  
  sym_pars.mat_mspe[i]<-sym_pars.mat_pred_perf[[i]]$MSPE
  sym_pars.mat_logs[i]<-sym_pars.mat_pred_perf[[i]]$Log_Score
  sym_sep.mat_mspe[i]<-sym_sep.mat_pred_perf[[i]]$MSPE
  sym_sep.mat_logs[i]<-sym_sep.mat_pred_perf[[i]]$Log_Score
}


mspedata<-data.frame(MSPE=c(my_asym_pars.mat_mspe,my_asym_sep.mat_mspe,liz_asym_pars.mat_mspe,liz_asym_sep.mat_mspe,sym_pars.mat_mspe,sym_sep.mat_mspe),
                     Model=c(rep("Our asymmetric pars. Matern",times=nsim),rep("Our asymmetric intrinsic Matern",times=nsim),rep("Li-zhangs asymmetric pars. Matern",times=nsim),rep("Li-zhangs asymmetric intrinsic Matern",times=nsim),rep("Symmetric pars. Matern",times=nsim),rep("Symmetric intrinsic Matern",times=nsim)))


logsdata<-data.frame(Log_Score=c(my_asym_pars.mat_logs,my_asym_sep.mat_logs,liz_asym_pars.mat_logs,liz_asym_sep.mat_logs,sym_pars.mat_logs,sym_sep.mat_logs),
                     Model=c(rep("Our asymmetric pars. Matern",times=nsim),rep("Our asymmetric intrinsic Matern",times=nsim),rep("Li-zhangs asymmetric pars. Matern",times=nsim),rep("Li-zhangs asymmetric intrinsic Matern",times=nsim),rep("Symmetric pars. Matern",times=nsim),rep("Symmetric intrinsic Matern",times=nsim)))

library(ggplot2)
par(mfrow=c(1,2))
mspep<-ggplot(data = mspedata, aes(x=Model, y=MSPE)) + geom_boxplot(aes(fill=Model))
plot(mspep)
logsp<-ggplot(data = logsdata, aes(x=Model, y=Log_Score)) + geom_boxplot(aes(fill=Model))
plot(logsp)

save.image("Fullfinal_Estimates.RData")
nloglikdata<-data.frame(nloglikelihood=c(my_asym_nloglikelihood_pars.mat,my_asym_nloglikelihood_sep.mat,liz_asym_nloglikelihood_pars.mat,liz_asym_nloglikelihood_sep.mat,sym_nloglikelihood_pars.mat,sym_nloglikelihood_sep.mat),
                        Model=c(rep("Our asymmetric pars. Matern",times=nsim),rep("Our asymmetric intrinsic Matern",times=nsim),rep("Li-zhangs asymmetric pars. Matern",times=nsim),rep("Li-zhangs asymmetric intrinsic Matern",times=nsim),rep("Symmetric pars. Matern",times=nsim),rep("Symmetric intrinsic Matern",times=nsim)))

bpdata<-data.frame(AIC=c(my_asym_aic_pars.mat,my_asym_aic_sep.mat,liz_asym_aic_pars.mat,liz_asym_aic_sep.mat,sym_aic_pars.mat,sym_aic_sep.mat),
                   Model=c(rep("Our asymmetric pars. Matern",times=nsim),rep("Our asymmetric intrinsic Matern",times=nsim),rep("Li-zhangs asymmetric pars. Matern",times=nsim),rep("Li-zhangs asymmetric intrinsic Matern",times=nsim),rep("Symmetric pars. Matern",times=nsim),rep("Symmetric intrinsic Matern",times=nsim)))

par(mfrow=c(1,2))
bp<-ggplot(data = bpdata, aes(x=Model, y=AIC)) + geom_boxplot(aes(fill=Model))
plot(bp)
lk<-ggplot(data = nloglikdata, aes(x=Model, y=nloglikelihood)) + geom_boxplot(aes(fill=Model))
plot(lk)


######## Complete ##########



############# Reformatting the boxplots ########
bpdata<-data.frame(AIC=c(my_asym_aic_pars.mat,my_asym_aic_sep.mat,liz_asym_aic_pars.mat,liz_asym_aic_sep.mat,sym_aic_pars.mat,sym_aic_sep.mat),
                   Model=c(rep("M1",times=nsim),rep("M4",times=nsim),rep("M2",times=nsim),rep("M5",times=nsim),rep("M3",times=nsim),rep("M6",times=nsim)))

bp<-ggplot(data = bpdata, aes(x=Model, y=AIC)) + geom_boxplot(aes(fill=Model)) + ylab("AIC")+ theme(axis.text.x=element_text(size=36),axis.text.y=element_text(size=20),
                                                                                                    axis.title.y=element_text(size=26),axis.title.x = element_blank(),legend.position = "none")+scale_x_discrete(labels=c(expression({" "}[1]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{"l.a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"l.a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"})))
plot(bp)

loglikdata<-data.frame(Loglikelihood=c(-my_asym_nloglikelihood_pars.mat,-my_asym_nloglikelihood_sep.mat,-liz_asym_nloglikelihood_pars.mat,-liz_asym_nloglikelihood_sep.mat,-sym_nloglikelihood_pars.mat,-sym_nloglikelihood_sep.mat),
                       Model=c(rep("M1",times=nsim),rep("M4",times=nsim),rep("M2",times=nsim),rep("M5",times=nsim),rep("M3",times=nsim),rep("M6",times=nsim)))

lk<-ggplot(data = loglikdata, aes(x=Model, y=Loglikelihood)) + geom_boxplot(aes(fill=Model)) + ylab("Log likelihood")+ theme(axis.text.x=element_text(size=36),axis.text.y=element_text(size=20),
                                                                                                                             axis.title.y=element_text(size=26),axis.title.x = element_blank(),legend.position = "none")+scale_x_discrete(labels=c(expression({" "}[1]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{"l.a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"l.a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"})))

plot(lk)

mspedata<-data.frame(MSPE=c(my_asym_pars.mat_mspe,my_asym_sep.mat_mspe,liz_asym_pars.mat_mspe,liz_asym_sep.mat_mspe,sym_pars.mat_mspe,sym_sep.mat_mspe),
                     Model=c(rep("M1",times=nsim),rep("M4",times=nsim),rep("M2",times=nsim),rep("M5",times=nsim),rep("M3",times=nsim),rep("M6",times=nsim)))


logsdata<-data.frame(Log_Score=c(my_asym_pars.mat_logs,my_asym_sep.mat_logs,liz_asym_pars.mat_logs,liz_asym_sep.mat_logs,sym_pars.mat_logs,sym_sep.mat_logs),
                     Model=c(rep("M1",times=nsim),rep("M4",times=nsim),rep("M2",times=nsim),rep("M5",times=nsim),rep("M3",times=nsim),rep("M6",times=nsim)))

mspep<-ggplot(data = mspedata, aes(x=Model, y=MSPE)) + geom_boxplot(aes(fill=Model))+ ylab("MSPE")+ theme(axis.text.x=element_text(size=36),axis.text.y=element_text(size=20),
                                                                                                          axis.title.y=element_text(size=26),axis.title.x = element_blank(),legend.position = "none")+scale_x_discrete(labels=c(expression({" "}[1]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{"l.a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"l.a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"})))
plot(mspep)
logsp<-ggplot(data = logsdata, aes(x=Model, y=Log_Score)) + geom_boxplot(aes(fill=Model))+ ylab("LogS")+ theme(axis.text.x=element_text(size=36),axis.text.y=element_text(size=20),
                                                                                                               axis.title.y=element_text(size=26),axis.title.x = element_blank(),legend.position = "none")+scale_x_discrete(labels=c(expression({" "}[1]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{"l.a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"l.a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"})))
plot(logsp)




###### Now visualizing for the manuscript ########

#### creating spatial lags ###
n1<-6
h1lags<-seq(-0.6,0.6,length=3*(2*n1-1))
h2lags<-seq(-0.6,0.6,length=3*(2*n1-1))
##### Creating function for computing norm #####
norm_vec <- function(x) sqrt(sum(x^2))


########## Now reporting the cross-covariances ##########


##### function to visualize cross-covariance function as a function of spatial lags ####
asym.cross<-function(rhoij,sigmai,sigmaj,nuij,a,as11,as12,hlag1,hlag2)
{
  nx<-length(hlag1)
  ny<-length(hlag2)
  cross.cor<-matrix(NA,nrow=nx,ncol=ny)
  for(i in 1:nx)
  {
    for(j in 1:ny)
    {
      htemp<-c(hlag1[i],hlag2[j])
      h1<-htemp-c(as11,as12)
      
      cross.cor[i,j]<-rhoij*sigmai*sigmaj*my.matern(h=norm_vec(h1),a=a,sigma = 1,nu=nuij)
    }
    
  }
  par(cex.axis=1, cex.lab=1.3, cex.main=1, cex.sub=1,mar=c(4,4.2,0.2,0.5)+.1)
  image.plot(hlag1,hlag2,cross.cor,col = viridis(n=200),xlab=expression('h'[1]),ylab=expression('h'[2]),zlim=c(0,max(cross.cor)),las=1,xaxt="n",yaxt="n",legend.width = 0.5)
  contour(hlag1,hlag2,cross.cor,add = T,lty = 1,drawlabels = F,levels = c(0.05,0.15,0.25,0.35))
  abline(v=0,col="grey",lty=1,lwd=1)
  abline(h=0,col="grey",lty=1,lwd=1)
  axis(1,at=c(-0.6,0,0.6),las=1)
  axis(2,at=c(-0.6,0,0.6),las=1)
  points(as11,as12)
  return(list(xlag=hlag1,ylag=hlag2,corr=cross.cor))
  
}

sigma1=1
sigma2=1
sigma3=1
### C12 ###
par(mfrow=c(1,3))
c12.image<-asym.cross(rho12,sigma1,sigma2,nu12,a=a1,as11=a121,a122,hlag1=h1lags,hlag2=h2lags)
c13.image<-asym.cross(rho13,sigma1,sigma3,nu13,a=a1,as11=0,0,hlag1=h1lags,hlag2=h2lags)
c23.image<-asym.cross(rho23,sigma2,sigma2,nu12,a=a1,as11=0,0,hlag1=h1lags,hlag2=h2lags)


###########################################################################################
########### Plotting asymmetry parameters boxplot for the supplementary materials #########
###########################################################################################

l121<-l122<-l131<-l132<-l231<-l232<-numeric(length = 100)
for(i in 1:100)
{
l121[i]<-cross_my_asym_pars.mat_estimates[[i]]$par[4]
l122[i]<-cross_my_asym_pars.mat_estimates[[i]]$par[5]
l131[i]<-cross_my_asym_pars.mat_estimates[[i]]$par[6]
l132[i]<-cross_my_asym_pars.mat_estimates[[i]]$par[7]
l231[i]<-cross_my_asym_pars.mat_estimates[[i]]$par[8]
l232[i]<-cross_my_asym_pars.mat_estimates[[i]]$par[9]

}  
par(cex.axis=1, cex.lab=1, cex.main=1, cex.sub=1,mar=c(3,4,0.5,1)+.1)
par(mfrow=c(1,1))
boxplot(l121,l122,l131,l132,l231,l232,ylim=c(-0.05,0.05),col = "lightsteelblue",names = c(expression("\u2113"[12]['x']),expression("\u2113"[12]['y']),expression("\u2113"[13]['x']),expression("\u2113"[13]['y']),expression("\u2113"[23]['x']),expression("\u2113"[23]['y'])),las=1)
#abline(h=0,col="limegreen",lwd=2)
abline(v=2.5)
abline(v=4.5)
legend("bottomright",lty=1,lwd=2,col="limegreen",c("True value"),cex=1.0)
lines(x=c(0.5,1.5),y=c(a121,a121),col="limegreen",lwd=2)
lines(x=c(1.5,2.5),y=c(a122,a122),col="limegreen",lwd=2)
lines(x=c(2.5,3.5),y=c(0,0),col="limegreen",lwd=2)
lines(x=c(3.5,4.5),y=c(0,0),col="limegreen",lwd=2)
lines(x=c(4.5,5.5),y=c(0,0),col="limegreen",lwd=2)
lines(x=c(5.5,6.5),y=c(0,0),col="limegreen",lwd=2)
