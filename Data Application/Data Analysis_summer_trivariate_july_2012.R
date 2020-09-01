#########################################################################
##### Preparing a trivariate data for PM2.5, WS and other variables #####
#########################################################################
############ Loading required libraries ############
####################################################
library(sp)
library(maps)
library(maptools)
library(geosphere)
library(fields)
library(MASS)
library(doParallel)
library(rgdal)
library(viridis)
registerDoParallel(cores = 6)
library(ggplot2)
library(RColorBrewer)
getDoParWorkers()
#########################################################
################ Data Creation portion ##################
#########################################################

setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Data Analysis codes")
dir()
##### Functions #########

my.matern<-function(h,a,sigma,nu)
{
  h[h==0]<-1e-10
  num1<-(sigma^2)*(2^(1-nu))/gamma(nu)
  num2<-(2*sqrt(nu)*h*a)^nu
  num3<-besselK(x=(2*sqrt(nu)*h*a), nu=nu)
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


rho_par<-function(a)
{ 
  
  L1<-matrix(c(1,a[1],a[2],0,1,a[3],0,0,1),byrow = T,nrow=3)
  un_norm<-L1%*%t(L1)
  dterm<-1/sqrt(diag(un_norm))
  dmat<-matrix(0,nrow=3,ncol=3)
  diag(dmat)<-dterm
  return(dmat%*%un_norm%*%dmat)
}



dir()
####################################
#### Exploratory data analysis #####
####################################
my_work_data<-read.csv("my_data_my_work_data_july2012v2.csv")

####### Transforming Longitude Latitude to Mercator projections####
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}

###### The study region falls into zone 17,18 and 19, of which the zone 18 is the most central, thats why we chose UTM zone 18 for projections
projxy<-LongLatToUTM(my_work_data$Longitude,my_work_data$Latitude,zone = 18) 
mydata<-my_work_data
WS1<-mydata$WS1 #(Wind speed in East-West direction)
WS2<-mydata$WS2 #(Wind speed in North-South direction)
PM25<-mydata$pm25_monthly_mean #(Monthly mean PM2.5 data)
WS<-mydata$WS  #(Magnitude Wind Speed)
TC<-mydata$RH2m  # (TC= Relative Humidity at 2m above ground)
mydata$Longitude<-projxy$X ### Projected locations x-coordinate
mydata$Latitude<-projxy$Y ### Projected locations y-coordinate
my_work_data$Longitude<-projxy$X ### Projected locations x-coordinate
my_work_data$Latitude<-projxy$Y ### Projected locations y-coordinate

hist(WS) #### Histogram for the Wind Speed ###
hist(TC,breaks = 15) #### Histogram for the Relative Humidity ###
hist(log(PM25))  ####### Histogram for the log(PM2.5 monthly mean) ### 
hist((PM25))    ####### Histogram for the PM2.5 monthly mean #####





######################################################################
######### Plotting directional histogram for the windspeed ###########
######################################################################
plot.windrose <- function(data,
                          spd,
                          dir,
                          spdres = 2,
                          dirres = 30,
                          spdmin = 2,
                          spdmax = 20,
                          spdseq = NULL,
                          palette = "YlGnBu",
                          countmax = NA,
                          debug = 0){
  
  
  # Look to see what data was passed in to the function
  if (is.numeric(spd) & is.numeric(dir)){
    # assume that we've been given vectors of the speed and direction vectors
    data <- data.frame(spd = spd,
                       dir = dir)
    spd = "spd"
    dir = "dir"
  } else if (exists("data")){
    # Assume that we've been given a data frame, and the name of the speed 
    # and direction columns. This is the format we want for later use.    
  }  
  
  # Tidy up input data ----
  n.in <- NROW(data)
  dnu <- (is.na(data[[spd]]) | is.na(data[[dir]]))
  data[[spd]][dnu] <- NA
  data[[dir]][dnu] <- NA
  
  # figure out the wind speed bins ----
  if (missing(spdseq)){
    spdseq <- seq(spdmin,spdmax,spdres)
  } else {
    if (debug >0){
      cat("Using custom speed bins \n")
    }
  }
  # get some information about the number of bins, etc.
  n.spd.seq <- length(spdseq)
  n.colors.in.range <- n.spd.seq - 1
  
  # create the color map
  spd.colors <- colorRampPalette(brewer.pal(min(max(3,
                                                    n.colors.in.range),
                                                min(9,
                                                    n.colors.in.range)),                                               
                                            palette))(n.colors.in.range)
  
  if (max(data[[spd]],na.rm = TRUE) > spdmax){    
    spd.breaks <- c(spdseq,
                    max(data[[spd]],na.rm = TRUE))
    spd.labels <- c(paste(c(spdseq[1:n.spd.seq-1]),
                          '-',
                          c(spdseq[2:n.spd.seq])),
                    paste(spdmax,
                          "-",
                          max(data[[spd]],na.rm = TRUE)))
    spd.colors <- c(spd.colors, "grey50")
  } else{
    spd.breaks <- spdseq
    spd.labels <- paste(c(spdseq[1:n.spd.seq-1]),
                        '-',
                        c(spdseq[2:n.spd.seq]))    
  }
  data$spd.binned <- cut(x = data[[spd]],
                         breaks = spd.breaks,
                         labels = spd.labels,
                         ordered_result = TRUE)
  # clean up the data
  data. <- na.omit(data)
  
  # figure out the wind direction bins
  dir.breaks <- c(-dirres/2,
                  seq(dirres/2, 360-dirres/2, by = dirres),
                  360+dirres/2)  
  dir.labels <- c(paste(360-dirres/2,"-",dirres/2),
                  paste(seq(dirres/2, 360-3*dirres/2, by = dirres),
                        "-",
                        seq(3*dirres/2, 360-dirres/2, by = dirres)),
                  paste(360-dirres/2,"-",dirres/2))
  # assign each wind direction to a bin
  dir.binned <- cut(data[[dir]],
                    breaks = dir.breaks,
                    ordered_result = TRUE)
  levels(dir.binned) <- dir.labels
  data$dir.binned <- dir.binned
  
  # Run debug if required ----
  if (debug>0){    
    cat(dir.breaks,"\n")
    cat(dir.labels,"\n")
    cat(levels(dir.binned),"\n")       
  }  
  
  # deal with change in ordering introduced somewhere around version 2.2
  if(packageVersion("ggplot2") > "2.2"){    
    cat("Hadley broke my code\n")
    data$spd.binned = with(data, factor(spd.binned, levels = rev(levels(spd.binned))))
    spd.colors = rev(spd.colors)
  }
  
  # create the plot ----
  p.windrose <- ggplot(data = data,
                       aes(x = dir.binned,
                           fill = spd.binned)) +
    geom_bar() + 
    scale_x_discrete(drop = FALSE,
                     labels = waiver()) +
    coord_polar(start = -((dirres/2)/360) * 2*pi) +
    scale_fill_manual(name = "Wind Speed (m/s)", 
                      values = spd.colors,
                      drop = FALSE) +
    theme(axis.title.x = element_blank())
  
  # adjust axes if required
  if (!is.na(countmax)){
    p.windrose <- p.windrose +
      ylim(c(0,countmax))
  }
  
  # print the plot
  print(p.windrose)  
  
  # return the handle to the wind rose
  return(p.windrose)
}
###############################################
###### Plotting directional histogram #########
###############################################
w.spd<-sqrt(WS1^2+WS2^2)
w.dir<-atan(WS2/WS1)
w.dir2<-w.dir*180/pi
w.dir2[w.dir2<0]<-360+w.dir2[w.dir2<0]
p0 <- plot.windrose(spd = w.spd,
                    dir = w.dir2,spdmin = 0,spdmax = max(w.spd))

###############################################
######### Plotting variables ##################
###############################################

############################################################################################################
###### Now we do regression of variables on spatial locations ##############################################
### We are considering log pm25 as our variable of interest to nealry satisfy the normality assumption######
############################################################################################################
pm25_reg<-lm(log(PM25)~1+my_work_data$Longitude+my_work_data$Latitude)
ws_reg<-lm(WS~1+my_work_data$Longitude+my_work_data$Latitude)
tc_reg<-lm(log(TC)~1+my_work_data$Longitude+my_work_data$Latitude)
#setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Revision work/Data Analysis with standard errors")
#####################################################################################
####### Plotting original variables and their residual after the regression #########
#####################################################################################
par(mfrow=c(2,3))
quilt.plot(x=my_work_data$Longitude,y=my_work_data$Latitude,z=my_work_data$pm25_monthly_mean,main=expression(log(PM[2.5])),xlab="x",ylab="y")
quilt.plot(x=my_work_data$Longitude,y=my_work_data$Latitude,z=my_work_data$WS,main=expression("Wind Speed"),xlab="x",ylab="y")
quilt.plot(x=my_work_data$Longitude,y=my_work_data$Latitude,z=TC,main=expression("Relative Humidity (2m above the groud) (RH2m)"),xlab="x",ylab="y")
quilt.plot(x=my_work_data$Longitude,y=my_work_data$Latitude,z=pm25_reg$residuals,main=expression(paste(log(PM[2.5])," residuals")),xlab="x",ylab="y")
quilt.plot(x=my_work_data$Longitude,y=my_work_data$Latitude,z=ws_reg$residuals,main=expression("Wind Speed Residuals"),xlab="x",ylab="y")
quilt.plot(x=my_work_data$Longitude,y=my_work_data$Latitude,z=tc_reg$residuals,main=expression("RH2m Residuals"),xlab="x",ylab="y")

##################################################
###### Plotting histograms of the residuals ######
##################################################
hist((pm25_reg$residuals),breaks = 15,main=expression(paste(log(PM[2.5])," residuals")))
hist(ws_reg$residuals,main=expression("Wind Speed Residuals"))
hist(tc_reg$residuals,main=expression("RH2m Residuals"),breaks = 15)

############################################
### Now we standardize the the residuals ###
############################################

pm25_reg$residuals<-(pm25_reg$residuals-mean(pm25_reg$residuals))/sd(pm25_reg$residuals)
ws_reg$residuals<-(ws_reg$residuals-mean(ws_reg$residuals))/sd(ws_reg$residuals)
tc_reg$residuals<-(tc_reg$residuals-mean(tc_reg$residuals))/sd(tc_reg$residuals)

my_transform_data<-data.frame(x=mydata$Longitude,y=mydata$Latitude,var1=pm25_reg$residuals,var2=ws_reg$residuals,var3=tc_reg$residuals)


### variable 1= Residual (log pm2.5), Variable 2= Residual (u-component), Variable 3 =Residual (v-component)
hist(my_transform_data$var1,main=expression("Variable 1"))
hist(my_transform_data$var2,main=expression("Variable 2"))
hist(my_transform_data$var3,breaks = 15,main=expression("Variable 3"))

##### Pearson's correlation coefficients #####
cor(my_transform_data$var1,my_transform_data$var2)
cor(my_transform_data$var1,my_transform_data$var3)
cor(my_transform_data$var2,my_transform_data$var3)


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
  nuggetv1<-p[8]
  nuggetv2<-p[9]
  nuggetv3<-p[10]
  if(sum(p<0)!=0) ### Positivity constraint for parameters
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
    C11<-C11+diag(rep(nuggetv1,n))
    C22<-C22+diag(rep(nuggetv2,n))
    C33<-C33+diag(rep(nuggetv3,n))
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
# for the data z, and locations x,, y, and n=n
# It returns the negative loglikelihood value
###############################################################################

mle_marginal_mlv<-function(pars)
{
  return(mle_marginals_all(p=pars,z=c(my_transform_data$var1,my_transform_data$var2,my_transform_data$var3),x=my_transform_data$x,y=my_transform_data$y,n=length(my_transform_data[,1]))$mlv)
}

#########################################################################################
# optim_marg_loglik is an optimizer function that minimizes mle_marginal_mlv over pars###
#########################################################################################


optim_marg_loglik <- function(par){
  optim(par=par,
        fn = mle_marginal_mlv,
        hessian=T,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=5000))
}

#############################
#### Initializing value #####
#############################
## Guessed values for marginal parameters
## Guessed range = 1/160000
## Guessed \nu_i and \sigma_i = 1
## Guessed nuggets = 0.1
marg.init<-c(1/160000,1,1,1,1,1,1,0.1,0.1,0.1)
marg.start<-proc.time()
marginal_pars.mat_estimates<-optim_marg_loglik(par = marg.init)
marg.end<-proc.time()

fit.time.marginal<-marg.end-marg.start
fit.time.marginal[3]/60
# The marginal estimates are valid for the Symmetric and Asymmetric Parsimonious Matern models #########

####Estimated common spatial range paramter in (meters)######
1/marginal_pars.mat_estimates$par[1]


############################################################################
########### Now we extract the marginal parameters in a matrix #############
############################################################################


margins_par_mat<-matrix(NA,nrow=1,ncol=length(marg.init))

for(i in 1:1)
{
  margins_par_mat[i,]<-marginal_pars.mat_estimates$par
}

############################################################################
######### Specification of margin matrices #################################
# margins_par_mat: column1=a, column2=nu1 column3=sigma1, column4=nu2, column4=sigma2, column6=nu3, column7=sigma3, column8=nugget1, column9=nugget2, column10=nugget3
############################################################################
############################################################################
#### Saving the Marginal Estimates #####
save.image("Data_Marginal_Estimates(RH2m).RData")


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
  nuggetv1<-margins[8]
  nuggetv2<-margins[9]
  nuggetv3<-margins[10]
  
  ############ putting hard constraints #####
  nu12<-(nu1+nu2)/2
  nu13<-(nu1+nu3)/2
  nu23<-(nu2+nu3)/2
  
  
  
  beta_mat<-rho_par(beta_p)       ##### beta_(ij)'s #####
  
  
  ####rho_(ij})'s
  rho12<-beta_mat[1,2]*((nu1^(nu1/2))*(nu2^(nu2/2))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*(((nu1+nu2)/2)^((nu1+nu2)/2)))
  rho13<-beta_mat[1,3]*((nu1^(nu1/2))*(nu3^(nu3/2))*gamma((nu1+nu3)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu3))*(((nu1+nu3)/2)^((nu1+nu3)/2)))
  rho23<-beta_mat[2,3]*((nu2^(nu2/2))*(nu3^(nu3/2))*gamma((nu2+nu3)/2))/(sqrt(gamma(nu2))*sqrt(gamma(nu3))*(((nu2+nu3)/2)^((nu2+nu3)/2)))
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
    C11<-C11+diag(rep(nuggetv1,n))
    C22<-C22+diag(rep(nuggetv2,n))
    C33<-C33+diag(rep(nuggetv3,n))
    
    
    ####### Cross covariances
    
    
    grid<-cbind(x,y)
    tempgrid1<-cbind(rep(a121,times=n),rep(a122,times=n))
    
    shift12grid<-grid-tempgrid1
   
    shift.dist12<-rdist(grid,shift12grid)
    
    C12<-rho12*sigma1*sigma2*my.matern(h=shift.dist12,a=a12,sigma = 1,nu=nu12)
    
    
    
    tempgrid2<-cbind(rep(a131,times=n),rep(a132,times=n))
    
    shift13grid<-grid-tempgrid2
    
    shift.dist13<-rdist(grid,shift13grid)
    
    C13<-rho13*sigma1*sigma3*my.matern(h=shift.dist13,a=a13,sigma = 1,nu=nu13)
    
    tempgrid3<-cbind(rep(a231,times=n),rep(a232,times=n))
    
    shift23grid<-grid-tempgrid3
   
    shift.dist23<-rdist(grid,shift23grid)
    
    C23<-rho23*sigma2*sigma3*my.matern(h=shift.dist23,a=a23,sigma = 1,nu=nu23)
    
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
      return(list(mlv=nlikelihood_value,range=a1,nu1=nu1,nu2=nu2,nu3=nu3,nu12=nu12,nu13=nu13,nu23=nu23,a121=a121,a122=a122,a131=a131,a132=a132,a231=a231,a232=a232,rho12=rho12,rho13=rho13,rho23=rho23))
    }
  } 
}


#################################################################
# mle_my_asym_pars.mat_mlv is a function
# takes argument parameters (pars)
# return lnegative loglikelihood value for 
# data = z, locations (x,y), marginal parameters = m.ests
#################################################################
mle_my_asym_pars.mat_mlv<-function(pars)
{
  return(mle_my_asym_pars.mat_all(p=pars,z=c(my_transform_data$var1,my_transform_data$var2,my_transform_data$var3),x=my_transform_data$x,y=my_transform_data$y,n=length(my_transform_data[,1]),margins=m.ests,onlyCOV = FALSE)$mlv)
}


#########################################################################################################
# optim_loglik_my_asym_pars.mat is an optimizer for mle_my_asym_pars.mat_mlv over parameters (pars)######
#########################################################################################################
optim_loglik_my_asym_pars.mat <- function(par){
  optim(par=par,
        fn = mle_my_asym_pars.mat_mlv,
        hessian=T,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}



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
  nuggetv1<-margins[8]
  nuggetv2<-margins[9]
  nuggetv3<-margins[10]
  ############ putting hard constraints #####
  nu12<-(nu1+nu2)/2
  nu13<-(nu1+nu3)/2
  nu23<-(nu2+nu3)/2
  
  
  
  beta_mat<-rho_par(beta_p)       ##### beta_(ij)'s
  
  
  ####rho_(ij})'s
  rho12<-beta_mat[1,2]*((nu1^(nu1/2))*(nu2^(nu2/2))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*(((nu1+nu2)/2)^((nu1+nu2)/2)))
  rho13<-beta_mat[1,3]*((nu1^(nu1/2))*(nu3^(nu3/2))*gamma((nu1+nu3)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu3))*(((nu1+nu3)/2)^((nu1+nu3)/2)))
  rho23<-beta_mat[2,3]*((nu2^(nu2/2))*(nu3^(nu3/2))*gamma((nu2+nu3)/2))/(sqrt(gamma(nu2))*sqrt(gamma(nu3))*(((nu2+nu3)/2)^((nu2+nu3)/2)))
  
  
  x.s<-x
  y.s<-y
  dist.mat<-rdist(cbind(x.s,y.s))
  ###### Marginal covariances
  C11<-my.matern(h=dist.mat,a=a1,sigma = sigma1,nu=nu1)
  C22<-my.matern(h=dist.mat,a=a2,sigma = sigma2,nu=nu2)
  C33<-my.matern(h=dist.mat,a=a3,sigma = sigma3,nu=nu3)
  C11<-C11+diag(rep(nuggetv1,n))
  C22<-C22+diag(rep(nuggetv2,n))
  C33<-C33+diag(rep(nuggetv3,n))
  
  
  
  ####### Cross covariances
  
  
  grid<-cbind(x,y)
  tempgrid1<-cbind(rep(as11,times=n),rep(as12,times=n))
  tempgrid2<-cbind(rep(as21,times=n),rep(as22,times=n))
  tempgrid3<-cbind(rep(as31,times=n),rep(as32,times=n))
  
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
# data = z, locations (x,y), marginal parameters = m.ests
#################################################################



mle_liz_asym_pars.mat_mlv<-function(pars)
{
  return(mle_liz_asym_pars.mat_all(p=pars,z=c(my_transform_data$var1,my_transform_data$var2,my_transform_data$var3),x=my_transform_data$x,y=my_transform_data$y,n=length(my_transform_data[,1]),margins=m.ests,onlyCOV = FALSE)$mlv)
}


#########################################################################################################
# optim_loglik_liz_asym_pars.mat is an optimizer for mle_liz_asym_pars.mat_mlv over parameters (pars)####
#########################################################################################################


optim_loglik_liz_asym_pars.mat <- function(par){
  optim(par=par,
        fn = mle_liz_asym_pars.mat_mlv,
        hessian=T,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}


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
  nuggetv1<-margins[8]
  nuggetv2<-margins[9]
  nuggetv3<-margins[10]
  
  ############ putting hard constraints #####
  nu12<-(nu1+nu2)/2
  nu13<-(nu1+nu3)/2
  nu23<-(nu2+nu3)/2
  
  
  
  beta_mat<-rho_par(beta_p)       ##### beta_(ij)'s
  
  
  ####rho_(ij})'s
  rho12<-beta_mat[1,2]*((nu1^(nu1/2))*(nu2^(nu2/2))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*(((nu1+nu2)/2)^((nu1+nu2)/2)))
  rho13<-beta_mat[1,3]*((nu1^(nu1/2))*(nu3^(nu3/2))*gamma((nu1+nu3)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu3))*(((nu1+nu3)/2)^((nu1+nu3)/2)))
  rho23<-beta_mat[2,3]*((nu2^(nu2/2))*(nu3^(nu3/2))*gamma((nu2+nu3)/2))/(sqrt(gamma(nu2))*sqrt(gamma(nu3))*(((nu2+nu3)/2)^((nu2+nu3)/2)))
  
  
  x.s<-x
  y.s<-y
  dist.mat<-rdist(cbind(x.s,y.s))
  ###### Marginal covariances
  C11<-my.matern(h=dist.mat,a=a1,sigma = sigma1,nu=nu1)
  C22<-my.matern(h=dist.mat,a=a2,sigma = sigma2,nu=nu2)
  C33<-my.matern(h=dist.mat,a=a3,sigma = sigma3,nu=nu3)
  
  C11<-C11+diag(rep(nuggetv1,n))
  C22<-C22+diag(rep(nuggetv2,n))
  C33<-C33+diag(rep(nuggetv3,n))
  
  
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
# data = z, locations (x,y), marginal parameters = m.ests
#################################################################


mle_sym_pars.mat_mlv<-function(pars)
{
  return(mle_sym_pars.mat_all(p=pars,z=c(my_transform_data$var1,my_transform_data$var2,my_transform_data$var3),x=my_transform_data$x,y=my_transform_data$y,n=length(my_transform_data[,1]),margins=m.ests,onlyCOV = FALSE)$mlv)
}

###############################################################
# optim_loglik_sym_pars.mat is an optimizer for mle_my_asym_sep.mat_mlv over parameters (pars)
###############################################################

optim_loglik_sym_pars.mat <- function(par){
  optim(par=par,
        fn = mle_sym_pars.mat_mlv,
        hessian=T,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}


# init.sym.pars.mat.cross<-c(-0.18, 0.07,  0.2) This value is chosen so that the starting point match closely to the empirical cross-correlation coefficients
# The estimates obtained were same as that from the Natural starting point (0,0,0)
# rho_par(c(-0.18, 0.07,  0.2))
init.sym.pars.mat.cross<-c(0, 0,  0) #### Natural starting point for initial value
m.ests<-margins_par_mat[1,]
sym.cross.start<-proc.time()
cross_est_sym<-optim_loglik_sym_pars.mat(init.sym.pars.mat.cross)
sym.cross.end<-proc.time()

fit.cross.sym.time<-sym.cross.end-sym.cross.start+fit.time.marginal

fit.cross.sym.time[3]/60

info.sym<-mle_sym_pars.mat_all(p=cross_est_sym$par,z=c(my_transform_data$var1,my_transform_data$var2,my_transform_data$var3),x=my_transform_data$x,y=my_transform_data$y,n=length(my_transform_data[,1]),margins=m.ests,onlyCOV = FALSE)
##### Estimated colocated correlation coefficients \rho_{ij}
info.sym$rho12
info.sym$rho13
info.sym$rho23
########################################################################################################################################
####### Now we know that the maximum likelihood estimation is tough for the asymmetry parameters #######################################
###### So to get an idea of the cross-parameters, we compute the maximum likelihood estimates for three pair of bivariate process ######
########################################################################################################################################

##################################
#### mle_bivariate_asym is a function that computes negative log likelihod for bivariate process var1 var2 given their fixed marginal parameters
#### The parameters (over which the optimization is to be done) involved are colocated correlation coefficients \rho_{12}, and asymmetric paramter a12<-c(a121,a122)
###################################
mle_bivariate_asym<-function(p,var1,var2,x,y,a.est,nu1.est,nu2.est,sigma1.est,sigma2.est,nug1.est,nug2.est,onlyCOV)
{
  z<-c(var1,var2)
  n<-length(var1)
  beta_p<-p[1]                 #### beta paramterization
  a121<-p[2]                     #### asymmetry parameters
  a122<-p[3]
  a1<-a2<-a12<-a.est
  nu1<-nu1.est
  sigma1<-sigma1.est
  nu2<-nu2.est
  sigma2<-sigma2.est
  nuggetv1<-nug1.est
  nuggetv2<-nug2.est
  
  ############ putting hard constraints #####
  nu12<-(nu1+nu2)/2
  if(abs(beta_p)>1)
  {
    nlikelihood_value<-100000000
    return(list(mlv=nlikelihood_value,params=NULL))
  }
  
  else{
    
    ####rho_(ij})'s
    rho12<-beta_p*((nu1^(nu1/2))*(nu2^(nu2/2))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*(((nu1+nu2)/2)^((nu1+nu2)/2)))
    x.s<-x
    y.s<-y
    dist.mat<-rdist(cbind(x.s,y.s))
    ###### Marginal covariances
    C11<-my.matern(h=dist.mat,a=a1,sigma = sigma1,nu=nu1)
    C22<-my.matern(h=dist.mat,a=a2,sigma = sigma2,nu=nu2)
    C11<-C11+diag(rep(nuggetv1,n))
    C22<-C22+diag(rep(nuggetv2,n))
    
    
    ####### Cross covariances
    
    
    grid<-cbind(x,y)
    tempgrid1<-cbind(rep(a121,times=n),rep(a122,times=n))
    
    shift12grid<-grid-tempgrid1
    
    shift.dist12<-rdist(grid,shift12grid)
    
    C12<-rho12*sigma1*sigma2*my.matern(h=shift.dist12,a=a12,sigma = 1,nu=nu12)
    
    
    
    C<-rbind(cbind(C11,C12),cbind(t(C12),C22))
    
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
      return(list(mlv=nlikelihood_value,range=a1,nu1=nu1,nu2=nu2,nu12=nu12,a121=a121,a122=a122,rho12=rho12))
    }
  } 
}


#########################################
####### For fixed marginal parameters and var1, var2, mle_bi_asym_mlv returns only negative loglikelihood value for the
####### parameters "pars"
#########################################

mle_bi_asym_mlv<-function(pars)
{
  return(mle_bivariate_asym(p=pars,var1=myvar1,var2=myvar2,x=my_transform_data$x,y=my_transform_data$y,a.est=margins_par_mat[1,1],nu1.est=mynu1,nu2.est=mynu2,sigma1.est=mysigma1,sigma2.est=mysigma2,nug1.est=mynug1,nug2.est=mynug2,onlyCOV = FALSE)$mlv)
}
###############################################################
# optim_loglik_my_bi_asym is an optimizer for mle_bi_asym_mlv over parameters (pars)
###############################################################

optim_loglik_my_bi_asym <- function(par){
  optim(par=par,
        fn = mle_bi_asym_mlv,
        hessian=T,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}
#mle_bi_asym_mlv(as.numeric(trial_values[799,]))
#maximum order of correlation#
cor(my_transform_data$var1,my_transform_data$var2)
cor(my_transform_data$var1,my_transform_data$var3)
cor(my_transform_data$var2,my_transform_data$var3)

#### We estimate here C_{12}^a for a bivariate process #######
##### Start with Variable 1 and Variable 2 #######
ini.find.start<-proc.time()

myvar1<-my_transform_data$var1
myvar2<-my_transform_data$var2
mynu1<-margins_par_mat[1,2]
mysigma1<-margins_par_mat[1,3]
mynu2<-margins_par_mat[1,4]
mysigma2<-margins_par_mat[1,5]
mynug1<-margins_par_mat[1,8]
mynug2<-margins_par_mat[1,9]

##############################################################
###### checking on the grid values to get a good start #######
##############################################################
checklength=31 #(must be odd to ensure the presence of 0,0 lags in the grid values)
gridh1<-seq(-(max(my_transform_data$x)-min(my_transform_data$x)),(max(my_transform_data$x)-min(my_transform_data$x)),length=checklength)
gridh2<-seq(-(max(my_transform_data$y)-min(my_transform_data$y)),(max(my_transform_data$y)-min(my_transform_data$y)),length=checklength)
rho_adhoc<-seq(-1,1,by=0.1)
trial_values<-expand.grid(rho_adhoc,gridh1,gridh2)
### Number of iterations required ####
length(trial_values[,1])/getDoParWorkers()
ini.trial12<-foreach(i=1:length(trial_values[,1]))%dopar%
{
  library(fields)
  pars1<-trial_values[i,]
  mle_bi_asym_mlv(as.numeric(pars1))
}
ini.trial12[which.min(ini.trial12)]
trial_values[which.min(ini.trial12),]

tr12<-numeric(length(trial_values[,1]))
for(i in 1:length(trial_values[,1]))
{
  tr12[i]<-c(ini.trial12[[i]])
}
sort(tr12)


####### Based on grid based guess, we estimate C_{12}^{a} or equivalently C_{12}^{l.a}
init.value.12<-optim_loglik_my_bi_asym(as.numeric(trial_values[which.min(ini.trial12),]))
init.value.12
###### init.value.12 is the bivariate asymmetric cross-covariance fit for variable 1 and 2


#### We estimate here C_{13}^a for a bivariate process #######
######## Checking for variable 1 and Variable 3 ###########

myvar1<-my_transform_data$var1
myvar2<-my_transform_data$var3
mynu1<-margins_par_mat[1,2]
mysigma1<-margins_par_mat[1,3]
mynu2<-margins_par_mat[1,6]
mysigma2<-margins_par_mat[1,7]
mynug1<-margins_par_mat[1,8]
mynug2<-margins_par_mat[1,10]

##############################################################
###### checking on the grid values to get a good start #######
##############################################################
checklength=31 #(must be odd to ensure the presence of 0,0 lags in the grid values)
gridh1<-seq(-(max(my_transform_data$x)-min(my_transform_data$x)),(max(my_transform_data$x)-min(my_transform_data$x)),length=checklength)
gridh2<-seq(-(max(my_transform_data$y)-min(my_transform_data$y)),(max(my_transform_data$y)-min(my_transform_data$y)),length=checklength)
rho_adhoc<-seq(-1,1,by=0.1)
trial_values<-expand.grid(rho_adhoc,gridh1,gridh2)
### Number of iterations required ####
length(trial_values[,1])/getDoParWorkers()
ini.trial13<-foreach(i=1:length(trial_values[,1]))%dopar%
{
  library(fields)
  pars1<-trial_values[i,]
  mle_bi_asym_mlv(as.numeric(pars1))
}
ini.trial13[which.min(ini.trial13)]
trial_values[which.min(ini.trial13),]

tr13<-numeric(length(trial_values[,1]))
for(i in 1:length(trial_values[,1]))
{
  tr13[i]<-c(ini.trial13[[i]])
}
sort(tr13) 
init.value.13<-optim_loglik_my_bi_asym(as.numeric(trial_values[which.min(ini.trial13),]))
init.value.13
###### init.value.13 is the bivariate asymmetric cross-covariance fit for variable 1 and 3


#### We estimate here C_{23}^a for a bivariate process #######
#####  Variable 2 and Variable 3 #######
myvar1<-my_transform_data$var2
myvar2<-my_transform_data$var3
mynu1<-margins_par_mat[1,4]
mysigma1<-margins_par_mat[1,5]
mynu2<-margins_par_mat[1,6]
mysigma2<-margins_par_mat[1,7]
mynug1<-margins_par_mat[1,9]
mynug2<-margins_par_mat[1,10]


##############################################################
###### checking on the grid values to get a good start #######
##############################################################
checklength=31#(must be odd to ensure the presence of 0,0 lags in the grid values)
gridh1<-seq(-(max(my_transform_data$x)-min(my_transform_data$x)),(max(my_transform_data$x)-min(my_transform_data$x)),length=checklength)
gridh2<-seq(-(max(my_transform_data$y)-min(my_transform_data$y)),(max(my_transform_data$y)-min(my_transform_data$y)),length=checklength)
rho_adhoc<-seq(-1,1,by=0.1)
trial_values<-expand.grid(rho_adhoc,gridh1,gridh2)

### Number of iterations required ####
length(trial_values[,1])/getDoParWorkers()
ini.trial23<-foreach(i=1:length(trial_values[,1]))%dopar%
{
  library(fields)
  pars1<-trial_values[i,]
  mle_bi_asym_mlv(as.numeric(pars1))
}


ini.trial23[which.min(ini.trial23)]
trial_values[which.min(ini.trial23),]

tr23<-numeric(length(trial_values[,1]))
for(i in 1:length(trial_values[,1]))
{
  tr23[i]<-c(ini.trial23[[i]])
}
sort(tr23) 
init.value.23<-optim_loglik_my_bi_asym(as.numeric(trial_values[which.min(ini.trial23),]))
init.value.23
###### init.value.23 is the bivariate asymmetric cross-covariance fit for variable 1 and 3




######## Based on the bivariate asymmetric fits, we set the initial values to jointly estimate all three cross-covariances in a trivariate process ########

ini.final<-c(init.value.12$par[1],init.value.13$par[1],init.value.23$par[1],init.value.12$par[c(2,3)],init.value.13$par[c(2,3)],init.value.23$par[c(2,3)])
ini.rho<-matrix(c(1,ini.final[1],ini.final[2],ini.final[1],1,ini.final[3],ini.final[2],ini.final[3],1),nrow = 3)
ini.rho
chol(ini.rho)
##### Checking if the estimated beta_{ij}'s from bivariate fits satify our validity condition
testv<-ini.rho[1,2]^2+ini.rho[1,3]^2+ini.rho[2,3]^2<=1-2*abs(ini.rho[1,2]*ini.rho[1,3]*ini.rho[2,3])
testv

### We therefore set the initial values for betas such that it is very close to the matrix ini.rho
#### We choose the values as c(0.18, 0.38, -1.5) since rho_par(c(0.25, 0.38, -1.7)) is very close to ini.rho
ini.rho2<-rho_par(c(-1.1, 0.58, 0.9))
ini.rho2
testv2<-ini.rho2[1,2]^2+ini.rho2[1,3]^2+ini.rho2[2,3]^2<=1-2*abs(ini.rho2[1,2]*ini.rho2[1,3]*ini.rho2[2,3])
testv2
ini.rho
rho_par(c(-1.1, 0.58, 0.9))

ini.final[c(1,2,3)]<-c(-1.1, 0.58, 0.9)
###cross_est_sym

####### Estimating our asymmetric cross-covariances for a trivariate process with initial values based on bivariate fits
ini.find.end<-proc.time()


find.initial.time<-ini.find.end-ini.find.start

find.initial.time[3]/60

my.cross.start<-proc.time()
cross_est_myasym<-optim_loglik_my_asym_pars.mat(par = ini.final)

my.cross.end<-proc.time()

fit.cros.myasym.time<-my.cross.end-my.cross.start+fit.time.marginal

fit.cros.myasym.time[3]/60


save.image("tillmyasym(RH2mv2).RData")



init.my.pars.mat.cross<-ini.final

#################################################
##### Now we estimate Li-Zhangs model ###########
#################################################

liz.cross.start<-proc.time()
#ini.liz1<-ini.final[-c(6,7)]#### Correctly specifying the asymmetry parameters for pairs (1,2) and (2,3) (incorrect for (1,3))
ini.liz1<-ini.final[-c(8,9)]
ini.liz1[c(4,5)]<-ini.final[c(4,5)]
ini.liz1[c(6,7)]<-c(ini.final[4]+ini.final[8],ini.final[5]+ini.final[9])
mle_liz_asym_pars.mat_mlv(ini.liz1)


ini.liz2<-ini.final[-c(8,9)]#### Correctly specifying the asymmetry parameters for pairs (1,2) and (1,3) (incorrect for (2,3))

mle_liz_asym_pars.mat_mlv(ini.liz2)

liz.cross.start1<-proc.time()

ini.liz3<-ini.final[-c(8,9)]#### Correctly specifying the asymmetry parameters for pairs (2,3) and (1,3) (incorrect for (1,2))
ini.liz3[c(4,5)]<-c(ini.final[6]-ini.final[8],ini.final[7]-ini.final[9])
ini.liz3[c(6,7)]<-c(ini.final[6],ini.final[7])
mle_liz_asym_pars.mat_mlv(ini.liz3)

mle_liz_asym_pars.mat_mlv(ini.liz3)

liz.cross.start0<-proc.time()
##### Now we fit Li-Zhangs Model with all the 3 initial values ######
cross_est_lizasym1<-optim_loglik_liz_asym_pars.mat(par=ini.liz1)
liz.cross.start2<-proc.time()
cross_est_lizasym2<-optim_loglik_liz_asym_pars.mat(par=ini.liz2)
liz.cross.start3<-proc.time()
cross_est_lizasym3<-optim_loglik_liz_asym_pars.mat(par=ini.liz3)
liz.cross.start4<-proc.time()

fit.cros.lizasym.time<-liz.cross.start3-liz.cross.start2+fit.time.marginal

fit.cros.lizasym.time[3]/60


cross_est_lizasym<-cross_est_lizasym2 #### We finalize cross_est_lizasym2 because of its minimum negative loglikelihood among the three fits
liz.cross.end<-proc.time()

###### Saving all the cross estimates ######
save.image("tillcrossest(RH2mv2).RData")

########### Computing AIC values

aic.my<-2*cross_est_myasym$value+2*(length(marg.init)+length(init.my.pars.mat.cross))
aic.liz<-2*cross_est_lizasym$value+2*(length(marg.init)+length(ini.liz1))
aic.sym<-2*cross_est_sym$value+2*(length(marg.init)+length(init.sym.pars.mat.cross))
round(aic.my,1)
round(aic.liz,1)
round(aic.sym,1)


round(-cross_est_myasym$value,1)
round(-cross_est_lizasym$value,1)
round(-cross_est_sym$value,1)

######## Reporting time #######
##### Marginal estimation #######

fit.time.marginal[3]

### Symmetric model full estimation time #####
fit.cross.sym.time[3]

### My asymmetric model full estimation time #####
fit.cros.myasym.time[3]


### Lizhang asymmetric model full estimation time #####
fit.cros.lizasym.time[3]

##### Initial value finding #######
find.initial.time[3]


##### Now rporting the parameter estimates, standard errors and fitting time #####
library(HelpersMG)

###Marginal####

margins_par_mat
param.marg<-margins_par_mat
colnames(param.marg)<-c("a","nu1","sigma1","nu2","sigma2","nu3","sigma3","delta1","delta2","delta3")
param.marg[1]<-param.marg[1]*1000 ### multiplying by 100 since we are reporting everything in Kms and coordinate in the projected data are in meters
marg.se<-SEfromHessian(a=marginal_pars.mat_estimates$hessian) 
marg.se[1]<-marg.se[1]*1000
param.marg
marg.se<-matrix(marg.se,nrow=1)
colnames(marg.se)<-c("a","nu1","sigma1","nu2","sigma2","nu3","sigma3","delta1","delta2","delta3")
marg.se
round(param.marg,2)
round(marg.se,2)

##### Symmetric #####

rho_par(cross_est_sym$par)
par.cross.sym<-cross_est_sym$par
cross.sym.se.raw<-SEfromHessian(a=cross_est_sym$hessian)
rho_par(cross.sym.se.raw)

library(numDeriv)

#### This function compute mle for the symmetric model without parametrizing for beta ####
mle_sym_forhessian<-function(p,z,x,y,n,margins,onlyCOV)
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
  nuggetv1<-margins[8]
  nuggetv2<-margins[9]
  nuggetv3<-margins[10]
  
  ############ putting hard constraints #####
  nu12<-(nu1+nu2)/2
  nu13<-(nu1+nu3)/2
  nu23<-(nu2+nu3)/2
  
  
  
  #beta_mat<-rho_par(beta_p)       ##### beta_(ij)'s Removing this parametrization to directly compute the hessian for
                                   #### \beta12, \beta13, \beta23
  
  
  ####rho_(ij})'s
  rho12<-beta_p[1]*((nu1^(nu1/2))*(nu2^(nu2/2))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*(((nu1+nu2)/2)^((nu1+nu2)/2)))
  rho13<-beta_p[2]*((nu1^(nu1/2))*(nu3^(nu3/2))*gamma((nu1+nu3)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu3))*(((nu1+nu3)/2)^((nu1+nu3)/2)))
  rho23<-beta_p[3]*((nu2^(nu2/2))*(nu3^(nu3/2))*gamma((nu2+nu3)/2))/(sqrt(gamma(nu2))*sqrt(gamma(nu3))*(((nu2+nu3)/2)^((nu2+nu3)/2)))
  
  
  x.s<-x
  y.s<-y
  dist.mat<-rdist(cbind(x.s,y.s))
  ###### Marginal covariances
  C11<-my.matern(h=dist.mat,a=a1,sigma = sigma1,nu=nu1)
  C22<-my.matern(h=dist.mat,a=a2,sigma = sigma2,nu=nu2)
  C33<-my.matern(h=dist.mat,a=a3,sigma = sigma3,nu=nu3)
  
  C11<-C11+diag(rep(nuggetv1,n))
  C22<-C22+diag(rep(nuggetv2,n))
  C33<-C33+diag(rep(nuggetv3,n))
  
  
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


mle_sym_pars.mat_mlv.forhess<-function(pars)
{
  return(mle_sym_forhessian(p=pars,z=c(my_transform_data$var1,my_transform_data$var2,my_transform_data$var3),x=my_transform_data$x,y=my_transform_data$y,n=length(my_transform_data[,1]),margins=m.ests,onlyCOV = FALSE)$mlv)
}



sym_p<-rho_par(cross_est_sym$par)

sym.hess<-hessian(mle_sym_pars.mat_mlv.forhess,x=c(sym_p[1,2],sym_p[1,3],sym_p[2,3]))


sym.se<-SEfromHessian(a=sym.hess)

####Cross parameters####
sym.betas<-c(sym_p[1,2],sym_p[1,3],sym_p[2,3])
sym.betas<-matrix(sym.betas,nrow=1)
colnames(sym.betas)<-c("beta12","beta13","beta23")

round(sym.betas,2)
round(sym.se,2)
round(fit.cross.sym.time[3],1)


#### My asymmetric #######

par.cross.myasym<-cross_est_myasym$par
par.cross.myasym<-matrix(par.cross.myasym,nrow=1)

colnames(par.cross.myasym)<-c("beta12","beta13","beta23","l121","l122","l131","l132","l231","l232")
mytemp<-rho_par(par.cross.myasym[1:3])
par.cross.myasym[1:3]<-c(mytemp[1,2],mytemp[1,3],mytemp[2,3])
par.cross.myasym
par.cross.myasym[4:9]<-par.cross.myasym[4:9]/1000 ### Converting in kms
round(par.cross.myasym,2)


######## Computing hessian #######

mle_my_asym_forhessian<-function(p,z,x,y,n,margins,onlyCOV)
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
  nuggetv1<-margins[8]
  nuggetv2<-margins[9]
  nuggetv3<-margins[10]
  
  ############ putting hard constraints #####
  nu12<-(nu1+nu2)/2
  nu13<-(nu1+nu3)/2
  nu23<-(nu2+nu3)/2
  
  
  
  #beta_mat<-rho_par(beta_p)       ##### removing this parameterization
  
  
  ####rho_(ij})'s
  rho12<-beta_p[1]*((nu1^(nu1/2))*(nu2^(nu2/2))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*(((nu1+nu2)/2)^((nu1+nu2)/2)))
  rho13<-beta_p[2]*((nu1^(nu1/2))*(nu3^(nu3/2))*gamma((nu1+nu3)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu3))*(((nu1+nu3)/2)^((nu1+nu3)/2)))
  rho23<-beta_p[3]*((nu2^(nu2/2))*(nu3^(nu3/2))*gamma((nu2+nu3)/2))/(sqrt(gamma(nu2))*sqrt(gamma(nu3))*(((nu2+nu3)/2)^((nu2+nu3)/2)))
  test_value<-(beta_p[1]^2+beta_p[2]^2+beta_p[3]^2-1)<=(-2*abs(beta_p[1])*abs(beta_p[2])*abs(beta_p[3]))
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
    C11<-C11+diag(rep(nuggetv1,n))
    C22<-C22+diag(rep(nuggetv2,n))
    C33<-C33+diag(rep(nuggetv3,n))
    
    
    ####### Cross covariances
    
    
    grid<-cbind(x,y)
    tempgrid1<-cbind(rep(a121,times=n),rep(a122,times=n))
    
    shift12grid<-grid-tempgrid1
    
    shift.dist12<-rdist(grid,shift12grid)
    
    C12<-rho12*sigma1*sigma2*my.matern(h=shift.dist12,a=a12,sigma = 1,nu=nu12)
    
    
    
    tempgrid2<-cbind(rep(a131,times=n),rep(a132,times=n))
    
    shift13grid<-grid-tempgrid2
    
    shift.dist13<-rdist(grid,shift13grid)
    
    C13<-rho13*sigma1*sigma3*my.matern(h=shift.dist13,a=a13,sigma = 1,nu=nu13)
    
    tempgrid3<-cbind(rep(a231,times=n),rep(a232,times=n))
    
    shift23grid<-grid-tempgrid3
    
    shift.dist23<-rdist(grid,shift23grid)
    
    C23<-rho23*sigma2*sigma3*my.matern(h=shift.dist23,a=a23,sigma = 1,nu=nu23)
    
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
      return(list(mlv=nlikelihood_value,range=a1,nu1=nu1,nu2=nu2,nu3=nu3,nu12=nu12,nu13=nu13,nu23=nu23,a121=a121,a122=a122,a131=a131,a132=a132,a231=a231,a232=a232,rho12=rho12,rho13=rho13,rho23=rho23))
    }
  } 
}


#################################################################
# mle_my_asym_pars.mat_mlv is a function
# takes argument parameters (pars)
# return lnegative loglikelihood value for 
# data = z, locations (x,y), marginal parameters = m.ests
#################################################################
mle_my_asym_forhessian_mlv<-function(pars)
{
  return(mle_my_asym_forhessian(p=pars,z=c(my_transform_data$var1,my_transform_data$var2,my_transform_data$var3),x=my_transform_data$x,y=my_transform_data$y,n=length(my_transform_data[,1]),margins=m.ests,onlyCOV = FALSE)$mlv)
}



myasym.hess<-hessian(mle_my_asym_forhessian_mlv,x=c(par.cross.myasym[1:3],cross_est_myasym$par[4:9]))


myasym.se<-SEfromHessian(a=myasym.hess)
myasym.se[4:9]<-myasym.se[4:9]/1000
round(par.cross.myasym,2)
round(myasym.se,2)
#### My asymmetric Time #####
round(fit.cros.myasym.time[3],1)

######## Li-Zhangs asymmetric model #########


par.cross.lizasym<-cross_est_lizasym$par
mytemp2<-rho_par(par.cross.lizasym[1:3])
par.cross.lizasym[1:3]<-c(mytemp2[1,2],mytemp2[1,3],mytemp2[2,3])

par.cross.lizasym[4:7]<-par.cross.lizasym[4:7]/1000
par.cross.lizasym<-matrix(par.cross.lizasym,nrow=1)
colnames(par.cross.lizasym)<-c("beta12","beta13","beta23","l21","l22","l31","l32")
round(par.cross.lizasym,2)


mle_liz_asym_forhessian<-function(p,z,x,y,n,margins,onlyCOV)
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
  nuggetv1<-margins[8]
  nuggetv2<-margins[9]
  nuggetv3<-margins[10]
  ############ putting hard constraints #####
  nu12<-(nu1+nu2)/2
  nu13<-(nu1+nu3)/2
  nu23<-(nu2+nu3)/2
  
  
  
  #beta_mat<-rho_par(beta_p)       ##### removing this parameterization
  
  
  ####rho_(ij})'s
  rho12<-beta_p[1]*((nu1^(nu1/2))*(nu2^(nu2/2))*gamma((nu1+nu2)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu2))*(((nu1+nu2)/2)^((nu1+nu2)/2)))
  rho13<-beta_p[2]*((nu1^(nu1/2))*(nu3^(nu3/2))*gamma((nu1+nu3)/2))/(sqrt(gamma(nu1))*sqrt(gamma(nu3))*(((nu1+nu3)/2)^((nu1+nu3)/2)))
  rho23<-beta_p[3]*((nu2^(nu2/2))*(nu3^(nu3/2))*gamma((nu2+nu3)/2))/(sqrt(gamma(nu2))*sqrt(gamma(nu3))*(((nu2+nu3)/2)^((nu2+nu3)/2)))
  
  
  x.s<-x
  y.s<-y
  dist.mat<-rdist(cbind(x.s,y.s))
  ###### Marginal covariances
  C11<-my.matern(h=dist.mat,a=a1,sigma = sigma1,nu=nu1)
  C22<-my.matern(h=dist.mat,a=a2,sigma = sigma2,nu=nu2)
  C33<-my.matern(h=dist.mat,a=a3,sigma = sigma3,nu=nu3)
  C11<-C11+diag(rep(nuggetv1,n))
  C22<-C22+diag(rep(nuggetv2,n))
  C33<-C33+diag(rep(nuggetv3,n))
  
  
  
  ####### Cross covariances
  
  
  grid<-cbind(x,y)
  tempgrid1<-cbind(rep(as11,times=n),rep(as12,times=n))
  tempgrid2<-cbind(rep(as21,times=n),rep(as22,times=n))
  tempgrid3<-cbind(rep(as31,times=n),rep(as32,times=n))
  
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




mle_liz_asym_forhessian_mlv<-function(pars)
{
  return(mle_liz_asym_forhessian(p=pars,z=c(my_transform_data$var1,my_transform_data$var2,my_transform_data$var3),x=my_transform_data$x,y=my_transform_data$y,n=length(my_transform_data[,1]),margins=m.ests,onlyCOV = FALSE)$mlv)
}


lizasym.hess<-hessian(mle_liz_asym_forhessian_mlv,x=c(par.cross.lizasym[1:3],cross_est_lizasym$par[4:7]))


lizasym.se<-SEfromHessian(a=lizasym.hess)
lizasym.se[4:7]<-lizasym.se[4:7]/1000
round(par.cross.lizasym,2)
round(lizasym.se,2)
#### My asymmetric Time #####
round(fit.cros.lizasym.time[3],1)
round(find.initial.time[3],1)


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
  MSPE_onlyPM25<-mean((sample[1:(s.size/3)]-pred[1:(s.size/3)])^2)
  Log_Score_onlyPM25<-mean(log(2*pi*pred.var[1:(s.size/3)])+(((sample[1:(s.size/3)]-pred[1:(s.size/3)])^2)/pred.var[1:(s.size/3)]))
  
  return(list(kriged.values=pred,kriged.variance=pred.var,Observed=sample,MSPE=MSPE,Log_Score=Log_Score,MSPE_onlyPM25=MSPE_onlyPM25,Log_Score_onlyPM25=Log_Score_onlyPM25))  
}

n1<-length(my_transform_data[,1])
sample<-c(my_transform_data$var1,my_transform_data$var2,my_transform_data$var3)
x.s<-my_transform_data$x
y.s<-my_transform_data$y

####### Computing full covariance matrix using the 3 models (our asymmetric, Lizhangs asymmetric and Symmetric)

C_my_asym<-mle_my_asym_pars.mat_all(p=cross_est_myasym$par,z=sample,x=x.s,y=y.s,n=n1,margins = margins_par_mat[1,],onlyCOV = TRUE)
C_liz_asym<-mle_liz_asym_pars.mat_all(p=cross_est_lizasym$par,z=sample,x=x.s,y=y.s,n=n1,margins = margins_par_mat[1,],onlyCOV = TRUE)
C_sym<-mle_sym_pars.mat_all(p=cross_est_sym$par,z=sample,x=x.s,y=y.s,n=n1,margins = margins_par_mat[1,],onlyCOV = TRUE)


pp_my<-pred_performance(sample = sample,C=C_my_asym)
pp_liz<-pred_performance(sample = sample,C=C_liz_asym)
pp_sym<-pred_performance(sample = sample,C=C_sym)


round(pp_my$MSPE,3)
round(pp_liz$MSPE,3)
round(pp_sym$MSPE,3)


round(pp_my$Log_Score,3)
round(pp_liz$Log_Score,3)
round(pp_sym$Log_Score,3)

dir()
###### saving data till prediction performances ######
save.image("till_pptri(RH2m).RData")

############## Data reporting ##################



#hist(my_transform_data$var1,main="Centered log(PM2.5)",col = "grey")
#hist(my_transform_data$var2,main="Centered Wind speed component u",col = "grey")
#hist(my_transform_data$var3,main="Centered Wind speed component v",col = "grey")

h1max<-1000000
h2max<-1000000

######## Marginal covariance reporting #######
n<-10
h1lags<-seq(-h1max,h1max,length=3*(2*n-1))
h2lags<-seq(-h2max,h2max,length=3*(2*n-1))
##### Creating function for computing norm #####
norm_vec <- function(x) sqrt(sum(x^2))

margins_par_mat

pars.mar.a<-margins_par_mat[1]
pars.mar.nu1<-margins_par_mat[2]
pars.mar.sigma1<-margins_par_mat[3]
pars.mar.nu2<-margins_par_mat[4]
pars.mar.sigma2<-margins_par_mat[5]
pars.mar.nu3<-margins_par_mat[6]
pars.mar.sigma3<-margins_par_mat[7]



mar_cov<-function(a,nu,sigma,hlag1,hlag2,label)
{
  nx<-length(hlag1)
  ny<-length(hlag2)
  cor<-matrix(NA,nrow=nx,ncol=ny)
  for(i in 1:nx)
  {
    for(j in 1:ny)
    {
      htemp<-c(hlag1[i],hlag2[j])
      h1<-htemp
      
      cor[i,j]<-my.matern(h=norm_vec(h1),a=a,sigma = sigma,nu=nu)
    }
    
  }
  par(cex.axis=3, cex.lab=3, cex.main=1, cex.sub=1,mar=c(4,4,1,3)+.1)
  image.plot(hlag1/1000,hlag2/1000,cor,col = viridis(n=100),xlab="Spatial lags (h1)",ylab="Spatial lags (h2)",main=paste(label))
  contour(hlag1,hlag2,cor,add = T,lty = 3,labcex = 2.5)
  abline(v=0,col="grey",lty=2)
  abline(h=0,col="grey",lty=2)
  return(list(xlag=hlag1,ylag=hlag2,corr=cor))
  
  
}

#### Cross-covariances without the nugget effect ####
par(mfrow=c(1,1))

pars.mat.c11<-mar_cov(a=pars.mar.a,nu=pars.mar.nu1,sigma = pars.mar.sigma1,hlag1 =h1lags,hlag2 = h2lags,label = "C11")
pars.mat.c22<-mar_cov(a=pars.mar.a,nu=pars.mar.nu2,sigma = pars.mar.sigma2,hlag1 =h1lags,hlag2 = h2lags,label = "C22")
pars.mat.c33<-mar_cov(a=pars.mar.a,nu=pars.mar.nu3,sigma = pars.mar.sigma3,hlag1 =h1lags,hlag2 = h2lags,label = "C33")






########## Now reporitng the cross-covariances ##########

asym.cross<-function(rhoij,sigmai,sigmaj,nuij,a,as11,as12,hlag1,hlag2,label)
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
  par(cex.axis=1, cex.lab=1, cex.main=1, cex.sub=1,mar=c(4,4,0.1,1)+.1)
  image.plot(hlag1/1000,hlag2/1000,cross.cor,col = viridis(n=100),xlab=expression(paste('h'[1]," (in Kilometers)")),ylab=expression(paste('h'[2]," (in Kilometers)")),las=1,xaxt="n",yaxt="n",legend.width = 0.5)
  #plot(hlag1,hlag2,col="white",xlab="Spatial lags (h1)",ylab="Spatial lags (h2)",main=(label))
  contour(hlag1/1000,hlag2/1000,cross.cor,add = T,nlevels = 6,drawlabels = FALSE,lty = 1)
  #abline(v=0,col="grey",lty=2,lwd=2)
  #abline(h=0,col="grey",lty=2,lwd=2)
  abline(v=0,col="black",lty=2,lwd=2)
  abline(h=0,col="black",lty=2,lwd=2)
  axis(1,at=c(-1000,-500,0,500,1000),las=1)
  axis(2,at=c(-1000,-500,0,500,1000),las=1)
  #innd<-which(cross.cor==max(cross.cor),arr.ind = )
  #points(as11,as12,pch=19,col="red",cex=1)
  return(list(xlag=hlag1,ylag=hlag2,corr=cross.cor))
  
}
#####case1 
n1=481
#n=50
sample<-c(my_transform_data$var1,my_transform_data$var2,my_transform_data$var3)
c1<-mle_my_asym_pars.mat_all(p=cross_est_myasym$par,z=sample,x=x.s,y=y.s,n=n1,margins=margins_par_mat[1,],onlyCOV=F)
par(mfrow=c(1,1))
case1c12<-asym.cross(rhoij = c1$rho12,sigmai = pars.mar.sigma1,sigmaj = pars.mar.sigma2,nuij = c1$nu12,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = c1$a121,as12 = c1$a122,label = expression("Cov(log(PM2.5),WS)"))


case1c13<-asym.cross(rhoij = c1$rho13,sigmai = pars.mar.sigma1,sigmaj = pars.mar.sigma3,nuij = c1$nu13,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = c1$a131,as12 = c1$a132,label = expression("Cov(log(PM2.5),RH2m)"))

case1c23<-asym.cross(rhoij = c1$rho23,sigmai = pars.mar.sigma2,sigmaj = pars.mar.sigma3,nuij = c1$nu23,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = c1$a231,as12 = c1$a232,label = expression("Cov(WS,RH2m)"))






c2<-mle_liz_asym_pars.mat_all(p=cross_est_lizasym$par,z=sample,x=x.s,y=y.s,n=n1,margins=margins_par_mat[1,],onlyCOV=F)


case2c12<-asym.cross(rhoij = c2$rho12,sigmai = pars.mar.sigma1,sigmaj = pars.mar.sigma2,nuij = c2$nu12,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = -(0-c2$as21),as12 = -(0-c2$as22),label = expression("Cov(log(PM2.5),WS)"))


case2c13<-asym.cross(rhoij = c2$rho13,sigmai = pars.mar.sigma1,sigmaj = pars.mar.sigma3,nuij = c2$nu13,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = -(0-c2$as31),as12 = -(0-c2$as32),label = expression("Cov(log(PM2.5),RH2m)"))

case2c23<-asym.cross(rhoij = c2$rho23,sigmai = pars.mar.sigma2,sigmaj = pars.mar.sigma3,nuij = c2$nu23,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = -(c2$as21-c2$as31),as12 = -(c2$as22-c2$as32),label = expression("Cov(WS,RH2m)"))



##### Cross-correlations ######
par(mfrow=c(3,3))
case1c12<-asym.cross(rhoij = c1$rho12,sigmai = 1,sigmaj = 1,nuij = c1$nu12,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = c1$a121,as12 = c1$a122,label = expression("Cor(log(PM2.5),WS)"))


case1c13<-asym.cross(rhoij = c1$rho13,sigmai = 1,sigmaj = 1,nuij = c1$nu13,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = c1$a131,as12 = c1$a132,label = expression("Cor(log(PM2.5),RH2m)"))

case1c23<-asym.cross(rhoij = c1$rho23,sigmai = 1,sigmaj = 1,nuij = c1$nu23,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = c1$a231,as12 = c1$a232,label = expression("Cor(WS,RH2m)"))


c2<-mle_liz_asym_pars.mat_all(p=cross_est_lizasym$par,z=sample,x=x.s,y=y.s,n=n1,margins=margins_par_mat[1,],onlyCOV=F)
#c2$a
#c1$range
case2c12<-asym.cross(rhoij = c2$rho12,sigmai =1,sigmaj = 1,nuij = c2$nu12,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = -(0-c2$as21),as12 = -(0-c2$as22),label = expression("Cor(log(PM2.5),WS)"))


case2c13<-asym.cross(rhoij = c2$rho13,sigmai = 1,sigmaj = 1,nuij = c2$nu13,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = -(0-c2$as31),as12 = -(0-c2$as32),label = expression("Cor(log(PM2.5),RH2m)"))

case2c23<-asym.cross(rhoij = c2$rho23,sigmai = 1,sigmaj = 1,nuij = c2$nu23,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = -(c2$as21-c2$as31),as12 = -(c2$as22-c2$as32),label = expression("Cor(WS,RH2m)"))



##### Plots for the manuscript #######
par(cex.axis=3, cex.lab=3, cex.main=1, cex.sub=1,mar=c(4,4,1,3)+.1)
par(mfrow=c(3,3))

case1c12<-asym.cross(rhoij = c1$rho12,sigmai = 1,sigmaj = 1,nuij = c1$nu12,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = c1$a121,as12 = c1$a122,label = expression(""))


case1c13<-asym.cross(rhoij = c1$rho13,sigmai = 1,sigmaj = 1,nuij = c1$nu13,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = c1$a131,as12 = c1$a132,label = expression(""))

case1c23<-asym.cross(rhoij = c1$rho23,sigmai = 1,sigmaj = 1,nuij = c1$nu23,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = c1$a231,as12 = c1$a232,label = expression(""))



case2c12<-asym.cross(rhoij = c2$rho12,sigmai =1,sigmaj = 1,nuij = c2$nu12,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = -(0-c2$as21),as12 = -(0-c2$as22),label = expression(""))


case2c13<-asym.cross(rhoij = c2$rho13,sigmai = 1,sigmaj = 1,nuij = c2$nu13,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = -(0-c2$as31),as12 = -(0-c2$as32),label = expression(""))

case2c23<-asym.cross(rhoij = c2$rho23,sigmai = 1,sigmaj = 1,nuij = c2$nu23,a=c1$range,hlag1 =h1lags,hlag2 = h2lags,as11 = -(c2$as21-c2$as31),as12 = -(c2$as22-c2$as32),label = expression(""))


###### Plotting bivariate Cross-covariance when estimated independently pairwise ###


##### Start with Variable 1 and Variable 2 #######


myvar1<-my_transform_data$var1
myvar2<-my_transform_data$var2
mynu1<-margins_par_mat[1,2]
mysigma1<-margins_par_mat[1,3]
mynu2<-margins_par_mat[1,4]
mysigma2<-margins_par_mat[1,5]
mynug1<-margins_par_mat[1,8]
mynug2<-margins_par_mat[1,9]


ind.c12<-mle_bivariate_asym(p=init.value.12$par,var1=myvar1,var2=myvar2,x=my_transform_data$x,y=my_transform_data$y,a.est=margins_par_mat[1,1],nu1.est=mynu1,nu2.est=mynu2,sigma1.est=mysigma1,sigma2.est=mysigma2,nug1.est=mynug1,nug2.est=mynu2,onlyCOV=F)
ind.c12.plot<-asym.cross(rhoij = ind.c12$rho12,sigmai =1,sigmaj = 1,nuij = ind.c12$nu12,a=ind.c12$range,hlag1 =h1lags,hlag2 = h2lags,as11 = ind.c12$a121,as12 = ind.c12$a122,label = expression("Cor(log(PM2.5),WS)"))



myvar1<-my_transform_data$var1
myvar2<-my_transform_data$var3
mynu1<-margins_par_mat[1,2]
mysigma1<-margins_par_mat[1,3]
mynu2<-margins_par_mat[1,6]
mysigma2<-margins_par_mat[1,7]
mynug1<-margins_par_mat[1,8]
mynug2<-margins_par_mat[1,10]


ind.c13<-mle_bivariate_asym(p=init.value.13$par,var1=myvar1,var2=myvar2,x=my_transform_data$x,y=my_transform_data$y,a.est=margins_par_mat[1,1],nu1.est=mynu1,nu2.est=mynu2,sigma1.est=mysigma1,sigma2.est=mysigma2,nug1.est=mynug1,nug2.est=mynu2,onlyCOV=F)
ind.c13.plot<-asym.cross(rhoij = ind.c13$rho12,sigmai =1,sigmaj = 1,nuij = ind.c13$nu12,a=ind.c13$range,hlag1 =h1lags,hlag2 = h2lags,as11 = ind.c13$a121,as12 = ind.c13$a122,label = expression(""))


myvar1<-my_transform_data$var2
myvar2<-my_transform_data$var3
mynu1<-margins_par_mat[1,4]
mysigma1<-margins_par_mat[1,5]
mynu2<-margins_par_mat[1,6]
mysigma2<-margins_par_mat[1,7]
mynug1<-margins_par_mat[1,9]
mynug2<-margins_par_mat[1,10]


ind.c23<-mle_bivariate_asym(p=init.value.23$par,var1=myvar1,var2=myvar2,x=my_transform_data$x,y=my_transform_data$y,a.est=margins_par_mat[1,1],nu1.est=mynu1,nu2.est=mynu2,sigma1.est=mysigma1,sigma2.est=mysigma2,nug1.est=mynug1,nug2.est=mynu2,onlyCOV=F)
ind.c23.plot<-asym.cross(rhoij = ind.c23$rho12,sigmai =1,sigmaj = 1,nuij = ind.c23$nu12,a=ind.c23$range,hlag1 =h1lags,hlag2 = h2lags,as11 = ind.c23$a121,as12 = ind.c23$a122,label = expression("Cor(WS,RH2m)"))



######## Percentage reduction in MSPE #######
base.mspe<-pp_sym$MSPE
base.logs<-pp_sym$Log_Score

my.mspe<-pp_my$MSPE
my.logs<-pp_my$Log_Score

liz.mspe<-pp_liz$MSPE
liz.logs<-pp_liz$Log_Score


round(my.mspe,3)
round(my.logs,3)
-round(cross_est_myasym$value,1)
round(aic.my,1)


-round(cross_est_lizasym$value,1)
round(aic.liz,1)
round(liz.mspe,3)
round(liz.logs,3)



-round(cross_est_sym$value,1)
round(aic.sym,1)
round(base.mspe,3)
round(base.logs,3)

####Liz#####
round(((base.mspe-liz.mspe)/base.mspe)*100,2)
round(((base.logs-liz.logs)/abs(base.logs))*100,2)

####My#####
round(((base.mspe-my.mspe)/base.mspe)*100,2)
round(((base.logs-my.logs)/abs(base.logs))*100,2)

aic.my
aic.liz
aic.sym


####### saving all the data ###########
save.image("Fullanalysisv3.RData")





####### Transformed and centered data plot #####
par(cex.axis=3, cex.lab=4, cex.main=1, cex.sub=1,mar=c(3,2.5,0,1)+.1)
par(mfrow=c(1,1))

par(mfrow=c(1,1),mar=c(2,3,2.5,4))
CR_NE <- c("maine", "new hampshire", "vermont", "new york", "massachusetts", "connecticut", "rhode island",
           "pennsylvania", "new jersey", "delaware", "maryland")
x.s<-my_work_data$Longitude
y.s<-my_work_data$Latitude
quilt.plot(x=x.s,y=y.s,z=my_transform_data$var1,main=expression(log(PM[2.5])),xlab="Longitude",ylab="Latitude")
map('state', region = CR_NE, lwd=2, interior = T, add = TRUE)
quilt.plot(x=x.s,y=y.s,z=my_transform_data$var2,main=expression("Wind speed u-component (WSU)"),xlab="Longitude",ylab="Latitude")
map('state', region = CR_NE, lwd=2, interior = T, add = TRUE)

quilt.plot(x=x.s,y=y.s,z=my_transform_data$var3,main=expression("Wind speed v-component (WSV)"),xlab="Longitude",ylab="Latitude")
map('state', region = CR_NE, lwd=2, interior = T, add = TRUE)

