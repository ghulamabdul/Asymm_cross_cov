##############################################
######## Visualisation of parameter space ####
##############################################
setwd("/Users/qadirga/Documents/Project 3/Manuscript/Visualization examples")
resxy<-100
x<-y<-seq(0,1,length=resxy)


#### Based on the inequality obtained in the Appendix #### 
z_full<-function(x,y)
{
  return(c(max(0,x*y-sqrt(1-x^2)*sqrt(1-y^2)),x*y+sqrt(1-x^2)*sqrt(1-y^2)))
}


#### Based on the inequality obtained in the Appendix #### 
z_sub<-function(x,y)
{
  if((-x*y+sqrt(1-x^2)*sqrt(1-y^2))<0)
    return(c(NA,NA))
  else
    return(c(max(0,-x*y-sqrt(1-x^2)*sqrt(1-y^2)),-x*y+sqrt(1-x^2)*sqrt(1-y^2)))
}

xylim<-expand.grid(x,y)
zlim_full<-zlim_sub<-matrix(NA,nrow=length(xylim$Var1),ncol=2)
for(i in 1:length(xylim$Var1))
{
  zlim_full[i,]<-z_full(xylim$Var1[i],xylim$Var2[i])
  zlim_sub[i,]<-z_sub(xylim$Var1[i],xylim$Var2[i])
}


resz<-40
plotdata_full<-plotdata_sub<-matrix(NA,nrow=resz*length(xylim$Var1),ncol=3)
for(i in 1:length(xylim$Var1))
{
  plotdata_full[((i-1)*resz+1):(i*resz),1]<-rep(xylim$Var1[i],times=resz)
  plotdata_full[((i-1)*resz+1):(i*resz),2]<-rep(xylim$Var2[i],times=resz)
  plotdata_full[((i-1)*resz+1):(i*resz),3]<-seq(zlim_full[i,1],zlim_full[i,2],length=resz)
  
  plotdata_sub[((i-1)*resz+1):(i*resz),1]<-rep(xylim$Var1[i],times=resz)
  plotdata_sub[((i-1)*resz+1):(i*resz),2]<-rep(xylim$Var2[i],times=resz)
  if(is.na(zlim_sub[i,1]))
  {
    plotdata_sub[((i-1)*resz+1):(i*resz),3]<-rep(NA,times=resz)
  }
  else
  {
    plotdata_sub[((i-1)*resz+1):(i*resz),3]<-seq(zlim_sub[i,1],zlim_sub[i,2],length=resz)
  }
  
}



library(rgl)


##########################################################
########### Now using persp ##############################
##########################################################
z_full_s3_u<-function(x,y)
{
  t<-x*y+sqrt(1-x^2)*sqrt(1-y^2)
  t[t<0]<-NA
  return(t)
}

z_full_s3_l<-function(x,y)
{ 
  t1<-x*y-sqrt(1-x^2)*sqrt(1-y^2)
  t1[t1<0]<-0
  return(t1)
}

z_sub_s3_u<-function(x,y)
{
  t<--x*y+sqrt(1-x^2)*sqrt(1-y^2)
  t[t<0]<-NA
  return(t)
}

z_sub_s3_l<-function(x,y)
{ t<--x*y+sqrt(1-x^2)*sqrt(1-y^2)
t1<--x*y-sqrt(1-x^2)*sqrt(1-y^2)
t1[t<0]<-NA
t1[t1<0]<-0
return(t1)
}


z_full_high<-function(x,y)
{
  return(x*y+sqrt(1-x^2)*sqrt(1-y^2))
}
x <- y <- seq(0, 1, length= 100)
z <- outer(x, y,z_sub_s3_u)
z2<-outer(x, y,z_sub_s3_l)
z3<-outer(x, y,z_full_s3_u)
z4<-outer(x, y,z_full_s3_l)

plot3d(x=plotdata_sub[,1],y=plotdata_sub[,2],z=plotdata_sub[,3],xlab=expression(paste("|",gamma[12],"(",omega,")","|")),ylab=expression(paste("|",gamma[13],"(",omega,")","|")),zlab=expression(paste("|",gamma[23],"(",omega,")","|")),col = "springgreen")
persp3d(x,y,z3,col="blue",alpha=0.3,add = T)
persp3d(x,y,z4,col="blue",alpha=0.3,add = T)
rgl.snapshot("1.png")


######### Contour plots ########

#plot(xylim,col="grey",pch=19)
#contour(x=x,y=y,z=matrix(zlim_full[,1],nrow=resxy,ncol=resxy),add = T)
#contour(x=x,y=y,z=matrix(zlim_sub[,1],nrow=resxy,ncol=resxy),add = T,col="red")

#zlim_full[,1][zlim_full[,1]>0]
#plot(xylim,col="grey",pch=19)
#contour(x=x,y=y,z=matrix(zlim_full[,2],nrow=resxy,ncol=resxy),add = T,levels = seq(0,1,0.1))
#contour(x=x,y=y,z=matrix(zlim_sub[,2],nrow=resxy,ncol=resxy),add = T,col="red")
resxy<-170
x<-y<-seq(0,1,length=resxy)


#### Based on the inequality obtained in the Appendix #### 
z_full<-function(x,y)
{
  return(c(max(0,x*y-sqrt(1-x^2)*sqrt(1-y^2)),x*y+sqrt(1-x^2)*sqrt(1-y^2)))
}


#### Based on the inequality obtained in the Appendix #### 
z_sub<-function(x,y)
{
  if((-x*y+sqrt(1-x^2)*sqrt(1-y^2))<0)
    return(c(NA,NA))
  else
    return(c(max(0,-x*y-sqrt(1-x^2)*sqrt(1-y^2)),-x*y+sqrt(1-x^2)*sqrt(1-y^2)))
}

xylim<-expand.grid(x,y)
zlim_full<-zlim_sub<-matrix(NA,nrow=length(xylim$Var1),ncol=2)
for(i in 1:length(xylim$Var1))
{
  zlim_full[i,]<-z_full(xylim$Var1[i],xylim$Var2[i])
  zlim_sub[i,]<-z_sub(xylim$Var1[i],xylim$Var2[i])
}

library(fields)
library(viridis)

par(mar=c(4,4.2,0.5,0.5))
image.plot(x=x,y=y,z=matrix(zlim_full[,2],nrow=resxy,ncol=resxy),col=viridis(n=resxy^2),
           xlab=expression(paste("|",gamma[12],"(",omega,")|")),ylab=expression(paste("|",dot(gamma[13]),"(",omega,")|")))
mtext('..', side=1, line=2.3, at=0.450)
mtext('.', side=2, line=3.75, at=0.435)
hyp.angle<-45
lines(x[1:(resxy/2)],y[1:(resxy/2)])
text(x[(resxy/2)]+0.03,y[(resxy/2)]+0.03,bquote("1.0"),srt=hyp.angle)
lines(x[(resxy/2+6):resxy],y[(resxy/2+6):resxy])
contour(x=x,y=y,z=matrix(zlim_full[,2],nrow=resxy,ncol=resxy),add = T,
        levels = seq(0,0.9,0.1),labcex = 1,method = "flattest")
points(x=1,y=1,pch=19,col="red")

image.plot(x=x,y=y,z=matrix(zlim_sub[,2],nrow=resxy,ncol=resxy,byrow = F),col = viridis(n=resxy^2),
           xlab=expression(paste("|",gamma[12],"(",omega,")|")),ylab=expression(paste("|",dot(gamma[13]),"(",omega,")|")))
mtext('..', side=1, line=2.3, at=0.450)
mtext('.', side=2, line=3.75, at=0.435)
contour(x=x,y=y,z=matrix(zlim_sub[,2],nrow=resxy,ncol=resxy),add = T,levels = seq(0,1,0.1),labcex = 1)
points(x=0.5,y=0.5,pch=19,col="red")


image.plot(x=x,y=y,z=matrix(zlim_full[,2],nrow=resxy,ncol=resxy)-matrix(zlim_sub[,2],nrow=resxy,ncol=resxy),col=viridis(n=resxy^2),
           xlab=expression(paste("|",gamma[12],"(",omega,")|")),ylab=expression(paste("|",dot(gamma[13]),"(",omega,")|")))
mtext('..', side=1, line=2.3, at=0.450)
mtext('.', side=2, line=3.75, at=0.435)

contour(x=x,y=y,z=matrix(zlim_full[,2],nrow=resxy,ncol=resxy)-matrix(zlim_sub[,2],nrow=resxy,ncol=resxy),add = T,labcex=1)



image.plot(x=x,y=y,z=matrix(zlim_full[,1],nrow=resxy,ncol=resxy),col = viridis(n=resxy^2),
           xlab=expression(paste("|",gamma[12],"(",omega,")|")),ylab=expression(paste("|",dot(gamma[13]),"(",omega,")|")))
mtext('..', side=1, line=2.3, at=0.450)
mtext('.', side=2, line=3.75, at=0.435)
contour(x=x,y=y,z=matrix(zlim_full[,1],nrow=resxy,ncol=resxy),add = T,labcex = 1)


image.plot(x=x,y=y,z=matrix(zlim_sub[,1],nrow=resxy,ncol=resxy),zlim=c(0,1),col = viridis(n=resxy^2),
           xlab=expression(paste("|",gamma[12],"(",omega,")|")),ylab=expression(paste("|",dot(gamma[13]),"(",omega,")|")))
mtext('..', side=1, line=2.3, at=0.450)
mtext('.', side=2, line=3.75, at=0.435)
contour(x=x,y=y,z=matrix(zlim_sub[,1],nrow=resxy,ncol=resxy),add = T,labcex = 1)

image.plot(x=x,y=y,z=matrix(zlim_sub[,1],nrow=resxy,ncol=resxy)-matrix(zlim_full[,1],nrow=resxy,ncol=resxy),col=viridis(n=resxy^2),zlim=c(0,1),
           xlab=expression(paste("|",gamma[12],"(",omega,")|")),ylab=expression(paste("|",dot(gamma[13]),"(",omega,")|")))
mtext('..', side=1, line=2.3, at=0.450)
mtext('.', side=2, line=3.75, at=0.435)
contour(x=x,y=y,z=matrix(zlim_sub[,1],nrow=resxy,ncol=resxy)-matrix(zlim_full[,1],nrow=resxy,ncol=resxy),add = T,labcex = 1)











#setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Revision work")
#image.plot(x=x,y=y,z=matrix(zlim_full[,2],nrow=resxy,ncol=resxy)-matrix(zlim_sub[,2],nrow=resxy,ncol=resxy),col=viridis(n=resxy^2))
#contour(x=x,y=y,z=matrix(zlim_full[,2],nrow=resxy,ncol=resxy)-matrix(zlim_sub[,2],nrow=resxy,ncol=resxy),add = T)

#abline(h=0.6,v=0.6)

#image.plot(x=x,y=y,z=matrix(zlim_full[,1],nrow=resxy,ncol=resxy),col = viridis(n=resxy^2))
#contour(x=x,y=y,z=matrix(zlim_full[,1],nrow=resxy,ncol=resxy),add = T)


#image.plot(x=x,y=y,z=matrix(zlim_sub[,1],nrow=resxy,ncol=resxy),zlim=c(0,1),col = viridis(n=resxy^2))
#contour(x=x,y=y,z=matrix(zlim_sub[,1],nrow=resxy,ncol=resxy),add = T)

#image.plot(x=x,y=y,z=matrix(zlim_sub[,1],nrow=resxy,ncol=resxy)-matrix(zlim_full[,1],nrow=resxy,ncol=resxy),col=viridis(n=resxy^2),zlim=c(0,1))
#contour(x=x,y=y,z=matrix(zlim_sub[,1],nrow=resxy,ncol=resxy)-matrix(zlim_full[,1],nrow=resxy,ncol=resxy),add = T)




           
