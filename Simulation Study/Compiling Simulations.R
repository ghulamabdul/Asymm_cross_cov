############### Compiling results from all the Scenarios #########
library(ggplot2)
library(dplyr)
library(plyr)
setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Simulation codes/Case1")
load("Fullfinal_Estimates.RData")
#M1 = My asymmetric parsimonious Matern
#M2 = Li Zhang's Asymmetric Parsimonious Matern
#M3 = Symmetric Matern
#M4 = My asymmetric Separable Matern
#M5 = Li Zhangs Separable Matern
#M6 = Symmetric Separable Matern
datas1aic<-data.frame(AIC=c(my_asym_aic_pars.mat,liz_asym_aic_pars.mat,sym_aic_pars.mat,my_asym_aic_sep.mat,liz_asym_aic_sep.mat,sym_aic_sep.mat),Model=c(rep("M1S1",length(my_asym_aic_pars.mat)),rep("M2S1",length(my_asym_aic_pars.mat)),rep("M3S1",length(my_asym_aic_pars.mat)),rep("M4S1",length(my_asym_aic_pars.mat)),rep("M5S1",length(my_asym_aic_pars.mat)),rep("M6S1",length(my_asym_aic_pars.mat))))
datas1loglik<-data.frame(loglik=-c(my_asym_nloglikelihood_pars.mat,liz_asym_nloglikelihood_pars.mat,sym_nloglikelihood_pars.mat,my_asym_nloglikelihood_sep.mat,liz_asym_nloglikelihood_sep.mat,sym_nloglikelihood_sep.mat),Model=c(rep("M1S1",length(my_asym_aic_pars.mat)),rep("M2S1",length(my_asym_aic_pars.mat)),rep("M3S1",length(my_asym_aic_pars.mat)),rep("M4S1",length(my_asym_aic_pars.mat)),rep("M5S1",length(my_asym_aic_pars.mat)),rep("M6S1",length(my_asym_aic_pars.mat))))
datas1mspe<-data.frame(mspe=c(my_asym_pars.mat_mspe,liz_asym_pars.mat_mspe,sym_pars.mat_mspe,my_asym_sep.mat_mspe,liz_asym_sep.mat_mspe,sym_sep.mat_mspe),Model=c(rep("M1S1",length(my_asym_aic_pars.mat)),rep("M2S1",length(my_asym_aic_pars.mat)),rep("M3S1",length(my_asym_aic_pars.mat)),rep("M4S1",length(my_asym_aic_pars.mat)),rep("M5S1",length(my_asym_aic_pars.mat)),rep("M6S1",length(my_asym_aic_pars.mat))))
datas1logs<-data.frame(logs=c(my_asym_pars.mat_logs,liz_asym_pars.mat_logs,sym_pars.mat_logs,my_asym_sep.mat_logs,liz_asym_sep.mat_logs,sym_sep.mat_logs),Model=c(rep("M1S1",length(my_asym_aic_pars.mat)),rep("M2S1",length(my_asym_aic_pars.mat)),rep("M3S1",length(my_asym_aic_pars.mat)),rep("M4S1",length(my_asym_aic_pars.mat)),rep("M5S1",length(my_asym_aic_pars.mat)),rep("M6S1",length(my_asym_aic_pars.mat))))


library(fields)
library(viridis)
#var1<-simt[1,1:(n^2)]
#var2<-simt[1,(n^2+1):(2*(n^2))]
#var3<-simt[1,(2*(n^2)+1):(3*(n^2))]
#par(mfrow=c(1,3))
#x.s<-grid$Var1
#y.s<-grid$Var2
par(mar=c(4,4,0.5,0.5))
quilt.plot(x=x.s,y=y.s,simt[1,1:(n^2)],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(n^2+1):(2*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(2*(n^2)+1):(3*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))





ls()
rmlist=ls()[-c(41,42,43,44)]
rm(list = rmlist)
setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Simulation codes/Case2")
load("Fullfinal_Estimates.RData")
datas2aic<-data.frame(AIC=c(my_asym_aic_pars.mat,liz_asym_aic_pars.mat,sym_aic_pars.mat,my_asym_aic_sep.mat,liz_asym_aic_sep.mat,sym_aic_sep.mat),Model=c(rep("M1S2",length(my_asym_aic_pars.mat)),rep("M2S2",length(my_asym_aic_pars.mat)),rep("M3S2",length(my_asym_aic_pars.mat)),rep("M4S2",length(my_asym_aic_pars.mat)),rep("M5S2",length(my_asym_aic_pars.mat)),rep("M6S2",length(my_asym_aic_pars.mat))))
datas2loglik<-data.frame(loglik=-c(my_asym_nloglikelihood_pars.mat,liz_asym_nloglikelihood_pars.mat,sym_nloglikelihood_pars.mat,my_asym_nloglikelihood_sep.mat,liz_asym_nloglikelihood_sep.mat,sym_nloglikelihood_sep.mat),Model=c(rep("M1S2",length(my_asym_aic_pars.mat)),rep("M2S2",length(my_asym_aic_pars.mat)),rep("M3S2",length(my_asym_aic_pars.mat)),rep("M4S2",length(my_asym_aic_pars.mat)),rep("M5S2",length(my_asym_aic_pars.mat)),rep("M6S2",length(my_asym_aic_pars.mat))))
datas2mspe<-data.frame(mspe=c(my_asym_pars.mat_mspe,liz_asym_pars.mat_mspe,sym_pars.mat_mspe,my_asym_sep.mat_mspe,liz_asym_sep.mat_mspe,sym_sep.mat_mspe),Model=c(rep("M1S2",length(my_asym_aic_pars.mat)),rep("M2S2",length(my_asym_aic_pars.mat)),rep("M3S2",length(my_asym_aic_pars.mat)),rep("M4S2",length(my_asym_aic_pars.mat)),rep("M5S2",length(my_asym_aic_pars.mat)),rep("M6S2",length(my_asym_aic_pars.mat))))
datas2logs<-data.frame(logs=c(my_asym_pars.mat_logs,liz_asym_pars.mat_logs,sym_pars.mat_logs,my_asym_sep.mat_logs,liz_asym_sep.mat_logs,sym_sep.mat_logs),Model=c(rep("M1S2",length(my_asym_aic_pars.mat)),rep("M2S2",length(my_asym_aic_pars.mat)),rep("M3S2",length(my_asym_aic_pars.mat)),rep("M4S2",length(my_asym_aic_pars.mat)),rep("M5S2",length(my_asym_aic_pars.mat)),rep("M6S2",length(my_asym_aic_pars.mat))))


#library(fields)
#library(viridis)
#var1<-simt[1,1:(n^2)]
#var2<-simt[1,(n^2+1):(2*(n^2))]
#var3<-simt[1,(2*(n^2)+1):(3*(n^2))]
#par(mfrow=c(1,3))
#x.s<-grid$Var1
#y.s<-grid$Var2
par(mar=c(4,4,0.5,0.5))
quilt.plot(x=x.s,y=y.s,simt[1,1:(n^2)],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(n^2+1):(2*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(2*(n^2)+1):(3*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))


rmlist=ls()[-(30:37)]
rm(list = rmlist)




setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Simulation codes/Case3")
load("Fullfinal_Estimates.RData")

datas3aic<-data.frame(AIC=c(my_asym_aic_pars.mat,liz_asym_aic_pars.mat,sym_aic_pars.mat,my_asym_aic_sep.mat,liz_asym_aic_sep.mat,sym_aic_sep.mat),Model=c(rep("M1S3",length(my_asym_aic_pars.mat)),rep("M2S3",length(my_asym_aic_pars.mat)),rep("M3S3",length(my_asym_aic_pars.mat)),rep("M4S3",length(my_asym_aic_pars.mat)),rep("M5S3",length(my_asym_aic_pars.mat)),rep("M6S3",length(my_asym_aic_pars.mat))))
datas3loglik<-data.frame(loglik=-c(my_asym_nloglikelihood_pars.mat,liz_asym_nloglikelihood_pars.mat,sym_nloglikelihood_pars.mat,my_asym_nloglikelihood_sep.mat,liz_asym_nloglikelihood_sep.mat,sym_nloglikelihood_sep.mat),Model=c(rep("M1S3",length(my_asym_aic_pars.mat)),rep("M2S3",length(my_asym_aic_pars.mat)),rep("M3S3",length(my_asym_aic_pars.mat)),rep("M4S3",length(my_asym_aic_pars.mat)),rep("M5S3",length(my_asym_aic_pars.mat)),rep("M6S3",length(my_asym_aic_pars.mat))))
datas3mspe<-data.frame(mspe=c(my_asym_pars.mat_mspe,liz_asym_pars.mat_mspe,sym_pars.mat_mspe,my_asym_sep.mat_mspe,liz_asym_sep.mat_mspe,sym_sep.mat_mspe),Model=c(rep("M1S3",length(my_asym_aic_pars.mat)),rep("M2S3",length(my_asym_aic_pars.mat)),rep("M3S3",length(my_asym_aic_pars.mat)),rep("M4S3",length(my_asym_aic_pars.mat)),rep("M5S3",length(my_asym_aic_pars.mat)),rep("M6S3",length(my_asym_aic_pars.mat))))
datas3logs<-data.frame(logs=c(my_asym_pars.mat_logs,liz_asym_pars.mat_logs,sym_pars.mat_logs,my_asym_sep.mat_logs,liz_asym_sep.mat_logs,sym_sep.mat_logs),Model=c(rep("M1S3",length(my_asym_aic_pars.mat)),rep("M2S3",length(my_asym_aic_pars.mat)),rep("M3S3",length(my_asym_aic_pars.mat)),rep("M4S3",length(my_asym_aic_pars.mat)),rep("M5S3",length(my_asym_aic_pars.mat)),rep("M6S3",length(my_asym_aic_pars.mat))))


par(mar=c(4,4,0.5,0.5))
quilt.plot(x=x.s,y=y.s,simt[1,1:(n^2)],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(n^2+1):(2*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(2*(n^2)+1):(3*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))



ls()
rmlist=ls()[-(29:40)]
rm(list = rmlist)


setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Simulation codes/Case4")
load("Fullfinal_Estimates.RData")

datas4aic<-data.frame(AIC=c(my_asym_aic_pars.mat,liz_asym_aic_pars.mat,sym_aic_pars.mat,my_asym_aic_sep.mat,liz_asym_aic_sep.mat,sym_aic_sep.mat),Model=c(rep("M1S4",length(my_asym_aic_pars.mat)),rep("M2S4",length(my_asym_aic_pars.mat)),rep("M3S4",length(my_asym_aic_pars.mat)),rep("M4S4",length(my_asym_aic_pars.mat)),rep("M5S4",length(my_asym_aic_pars.mat)),rep("M6S4",length(my_asym_aic_pars.mat))))
datas4loglik<-data.frame(loglik=-c(my_asym_nloglikelihood_pars.mat,liz_asym_nloglikelihood_pars.mat,sym_nloglikelihood_pars.mat,my_asym_nloglikelihood_sep.mat,liz_asym_nloglikelihood_sep.mat,sym_nloglikelihood_sep.mat),Model=c(rep("M1S4",length(my_asym_aic_pars.mat)),rep("M2S4",length(my_asym_aic_pars.mat)),rep("M3S4",length(my_asym_aic_pars.mat)),rep("M4S4",length(my_asym_aic_pars.mat)),rep("M5S4",length(my_asym_aic_pars.mat)),rep("M6S4",length(my_asym_aic_pars.mat))))
datas4mspe<-data.frame(mspe=c(my_asym_pars.mat_mspe,liz_asym_pars.mat_mspe,sym_pars.mat_mspe,my_asym_sep.mat_mspe,liz_asym_sep.mat_mspe,sym_sep.mat_mspe),Model=c(rep("M1S4",length(my_asym_aic_pars.mat)),rep("M2S4",length(my_asym_aic_pars.mat)),rep("M3S4",length(my_asym_aic_pars.mat)),rep("M4S4",length(my_asym_aic_pars.mat)),rep("M5S4",length(my_asym_aic_pars.mat)),rep("M6S4",length(my_asym_aic_pars.mat))))
datas4logs<-data.frame(logs=c(my_asym_pars.mat_logs,liz_asym_pars.mat_logs,sym_pars.mat_logs,my_asym_sep.mat_logs,liz_asym_sep.mat_logs,sym_sep.mat_logs),Model=c(rep("M1S4",length(my_asym_aic_pars.mat)),rep("M2S4",length(my_asym_aic_pars.mat)),rep("M3S4",length(my_asym_aic_pars.mat)),rep("M4S4",length(my_asym_aic_pars.mat)),rep("M5S4",length(my_asym_aic_pars.mat)),rep("M6S4",length(my_asym_aic_pars.mat))))


par(mar=c(4,4,0.5,0.5))
quilt.plot(x=x.s,y=y.s,simt[1,1:(n^2)],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(n^2+1):(2*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(2*(n^2)+1):(3*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))




ls()
rmlist=ls()[-(29:44)]
rm(list = rmlist)



setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Simulation codes/Case5")
load("Fullfinal_Estimates.RData")

datas5aic<-data.frame(AIC=c(my_asym_aic_pars.mat,liz_asym_aic_pars.mat,sym_aic_pars.mat,my_asym_aic_sep.mat,liz_asym_aic_sep.mat,sym_aic_sep.mat),Model=c(rep("M1S5",length(my_asym_aic_pars.mat)),rep("M2S5",length(my_asym_aic_pars.mat)),rep("M3S5",length(my_asym_aic_pars.mat)),rep("M4S5",length(my_asym_aic_pars.mat)),rep("M5S5",length(my_asym_aic_pars.mat)),rep("M6S5",length(my_asym_aic_pars.mat))))
datas5loglik<-data.frame(loglik=-c(my_asym_nloglikelihood_pars.mat,liz_asym_nloglikelihood_pars.mat,sym_nloglikelihood_pars.mat,my_asym_nloglikelihood_sep.mat,liz_asym_nloglikelihood_sep.mat,sym_nloglikelihood_sep.mat),Model=c(rep("M1S5",length(my_asym_aic_pars.mat)),rep("M2S5",length(my_asym_aic_pars.mat)),rep("M3S5",length(my_asym_aic_pars.mat)),rep("M4S5",length(my_asym_aic_pars.mat)),rep("M5S5",length(my_asym_aic_pars.mat)),rep("M6S5",length(my_asym_aic_pars.mat))))
datas5mspe<-data.frame(mspe=c(my_asym_pars.mat_mspe,liz_asym_pars.mat_mspe,sym_pars.mat_mspe,my_asym_sep.mat_mspe,liz_asym_sep.mat_mspe,sym_sep.mat_mspe),Model=c(rep("M1S5",length(my_asym_aic_pars.mat)),rep("M2S5",length(my_asym_aic_pars.mat)),rep("M3S5",length(my_asym_aic_pars.mat)),rep("M4S5",length(my_asym_aic_pars.mat)),rep("M5S5",length(my_asym_aic_pars.mat)),rep("M6S5",length(my_asym_aic_pars.mat))))
datas5logs<-data.frame(logs=c(my_asym_pars.mat_logs,liz_asym_pars.mat_logs,sym_pars.mat_logs,my_asym_sep.mat_logs,liz_asym_sep.mat_logs,sym_sep.mat_logs),Model=c(rep("M1S5",length(my_asym_aic_pars.mat)),rep("M2S5",length(my_asym_aic_pars.mat)),rep("M3S5",length(my_asym_aic_pars.mat)),rep("M4S5",length(my_asym_aic_pars.mat)),rep("M5S5",length(my_asym_aic_pars.mat)),rep("M6S5",length(my_asym_aic_pars.mat))))

par(mar=c(4,4,0.5,0.5))
quilt.plot(x=x.s,y=y.s,simt[1,1:(n^2)],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(n^2+1):(2*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(2*(n^2)+1):(3*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))


ls()
rmlist=ls()[-(27:46)]
rm(list = rmlist)


setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Simulation codes/Case6")
load("Fullfinal_Estimates.RData")

datas6aic<-data.frame(AIC=c(my_asym_aic_pars.mat,liz_asym_aic_pars.mat,sym_aic_pars.mat,my_asym_aic_sep.mat,liz_asym_aic_sep.mat,sym_aic_sep.mat),Model=c(rep("M1S6",length(my_asym_aic_pars.mat)),rep("M2S6",length(my_asym_aic_pars.mat)),rep("M3S6",length(my_asym_aic_pars.mat)),rep("M4S6",length(my_asym_aic_pars.mat)),rep("M5S6",length(my_asym_aic_pars.mat)),rep("M6S6",length(my_asym_aic_pars.mat))))
datas6loglik<-data.frame(loglik=-c(my_asym_nloglikelihood_pars.mat,liz_asym_nloglikelihood_pars.mat,sym_nloglikelihood_pars.mat,my_asym_nloglikelihood_sep.mat,liz_asym_nloglikelihood_sep.mat,sym_nloglikelihood_sep.mat),Model=c(rep("M1S6",length(my_asym_aic_pars.mat)),rep("M2S6",length(my_asym_aic_pars.mat)),rep("M3S6",length(my_asym_aic_pars.mat)),rep("M4S6",length(my_asym_aic_pars.mat)),rep("M5S6",length(my_asym_aic_pars.mat)),rep("M6S6",length(my_asym_aic_pars.mat))))
datas6mspe<-data.frame(mspe=c(my_asym_pars.mat_mspe,liz_asym_pars.mat_mspe,sym_pars.mat_mspe,my_asym_sep.mat_mspe,liz_asym_sep.mat_mspe,sym_sep.mat_mspe),Model=c(rep("M1S6",length(my_asym_aic_pars.mat)),rep("M2S6",length(my_asym_aic_pars.mat)),rep("M3S6",length(my_asym_aic_pars.mat)),rep("M4S6",length(my_asym_aic_pars.mat)),rep("M5S6",length(my_asym_aic_pars.mat)),rep("M6S6",length(my_asym_aic_pars.mat))))
datas6logs<-data.frame(logs=c(my_asym_pars.mat_logs,liz_asym_pars.mat_logs,sym_pars.mat_logs,my_asym_sep.mat_logs,liz_asym_sep.mat_logs,sym_sep.mat_logs),Model=c(rep("M1S6",length(my_asym_aic_pars.mat)),rep("M2S6",length(my_asym_aic_pars.mat)),rep("M3S6",length(my_asym_aic_pars.mat)),rep("M4S6",length(my_asym_aic_pars.mat)),rep("M5S6",length(my_asym_aic_pars.mat)),rep("M6S6",length(my_asym_aic_pars.mat))))


par(mar=c(4,4,0.5,0.5))
quilt.plot(x=x.s,y=y.s,simt[1,1:(n^2)],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(n^2+1):(2*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))
quilt.plot(x=x.s,y=y.s,simt[1,(2*(n^2)+1):(3*(n^2))],nx=n,ny=n,xlab=bquote(s[x]),ylab=bquote(s[y]),col=viridis(n=400))




ls()
rmlist=ls()[-(35:58)]
rm(list = rmlist)



setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Simulation codes")
save.image("Fullcompiled")

datamspe<-rbind(datas1mspe,datas2mspe,datas3mspe,datas4mspe,datas5mspe,datas6mspe)
max(datamspe$mspe)
sum(is.na(datamspe$mspe))
##### Plotting boxplot for MSPE #####

data<-cbind(datamspe,(c(rep('Case 1',600),rep('Case 2',600),rep('Case 3',600),rep('Case 4',600),rep('Case 5',600),rep('Case 6',600)))
)
colnames(data)<-c("value","individual","group")
data$colorfill<-as.character(rep(rep(1:6,each=100),times=6))

data=data %>% arrange(group)

data$Model<-substr(data$individual,start = 1,stop = 2)
p <- ggplot(data = data, aes(x=as.factor(individual), y=value)) + 
  geom_boxplot(aes(fill=Model),coef=2.7) + 
  ylab("MSPE") + 
  scale_fill_discrete(name = "Model", labels=c(expression({" "}[1]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{"LZ"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"LZ"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"}))) +
  theme(
    legend.title = element_text(color = "black", size = 15),
    legend.text = element_text(color = "black", size = 15),
    legend.position="bottom",
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) 
p<-p + facet_wrap( ~ group, scales="free")


p



############ Plotting for the boxplots for the logarithmic scores ############


datalogs<-rbind(datas1logs,datas2logs,datas3logs,datas4logs,datas5logs,datas6logs)
##### Plotting boxplot for MSPE #####

data<-cbind(datalogs,(c(rep('Case 1',600),rep('Case 2',600),rep('Case 3',600),rep('Case 4',600),rep('Case 5',600),rep('Case 6',600)))
)
colnames(data)<-c("value","individual","group")
data$colorfill<-as.character(rep(rep(1:6,each=100),times=6))

data=data %>% arrange(group)


data$Model<-substr(data$individual,start = 1,stop = 2)
p <- ggplot(data = data, aes(x=as.factor(individual), y=value)) + 
  geom_boxplot(aes(fill=Model),coef=5) + 
  ylab("mLogS") + 
  scale_fill_discrete(name = "Model", labels=c(expression({" "}[1]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{"LZ"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"LZ"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"}))) +
  theme(
    legend.title = element_text(color = "black", size = 15),
    legend.text = element_text(color = "black", size = 15),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) 
p<-p + facet_wrap( ~ group, scales="free")


p





############ Plotting for the boxplots for the AIC ############


dataaic<-rbind(datas1aic,datas2aic,datas3aic,datas4aic,datas5aic,datas6aic)


##### Plotting boxplot for AIC #####

data<-cbind(dataaic,(c(rep('Case 1',600),rep('Case 2',600),rep('Case 3',600),rep('Case 4',600),rep('Case 5',600),rep('Case 6',600)))
)
colnames(data)<-c("value","individual","group")
data$colorfill<-as.character(rep(rep(1:6,each=100),times=6))

data=data %>% arrange(group)
data$Model<-substr(data$individual,start = 1,stop = 2)
p <- ggplot(data = data, aes(x=as.factor(individual), y=value)) + 
  geom_boxplot(aes(fill=Model),coef=3.5) + 
  ylab("AIC") + 
  scale_fill_discrete(name = "Model", labels=c(expression({" "}[1]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{"LZ"}*{"("}*{bold("h")}*{")"}),expression({" "}[1]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"a"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{"LZ"}*{"("}*{bold("h")}*{")"}),expression({" "}[2]*{bold("C")}^{" "}*{"("}*{bold("h")}*{")"}))) +
  theme(
    legend.title = element_text(color = "black", size = 15),
    legend.text = element_text(color = "black", size = 15),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) 
p<-p + facet_wrap( ~ group, scales="free")


p



##### Change the address for different scenario to get the summary of each scenario ####
setwd("/Users/qadirga/Documents/Project 3/Submission to JABES/Simulation codes/Case5")
rm(list = ls())
load("Fullfinal_Estimates.RData")
round(mean(-my_asym_nloglikelihood_pars.mat),1)
round(mean(-liz_asym_nloglikelihood_pars.mat),1)
round(mean(-sym_nloglikelihood_pars.mat),1)
round(mean(-my_asym_nloglikelihood_sep.mat),1)
round(mean(-liz_asym_nloglikelihood_sep.mat),1)
round(mean(-sym_nloglikelihood_sep.mat),1)

round(sd(-my_asym_nloglikelihood_pars.mat),1)
round(sd(-liz_asym_nloglikelihood_pars.mat),1)
round(sd(-sym_nloglikelihood_pars.mat),1)
round(sd(-my_asym_nloglikelihood_sep.mat),1)
round(sd(-liz_asym_nloglikelihood_sep.mat),1)
round(sd(-sym_nloglikelihood_sep.mat),1)



round(mean(my_asym_aic_pars.mat),1)
round(mean(liz_asym_aic_pars.mat),1)
round(mean(sym_aic_pars.mat),1)
round(mean(my_asym_aic_sep.mat),1)
round(mean(liz_asym_aic_sep.mat),1)
round(mean(sym_aic_sep.mat),1)

round(sd(my_asym_aic_pars.mat),1)
round(sd(liz_asym_aic_pars.mat),1)
round(sd(sym_aic_pars.mat),1)
round(sd(my_asym_aic_sep.mat),1)
round(sd(liz_asym_aic_sep.mat),1)
round(sd(sym_aic_sep.mat),1)




#round(sd(my_asym_aic_pars.mat),1)
#round(sd(liz_asym_aic_pars.mat),1)
#round(sd(sym_aic_pars.mat),1)
#round(sd(my_asym_aic_sep.mat),1)
#round(sd(liz_asym_aic_sep.mat),1)
#round(sd(sym_aic_sep.mat),1)






#round(mean(my_asym_pars.mat_mspe),3)
#round(mean(liz_asym_pars.mat_mspe),3)
#round(mean(sym_pars.mat_mspe),3)
#round(mean(my_asym_sep.mat_mspe),3)
#round(mean(liz_asym_sep.mat_mspe),3)
#round(mean(sym_sep.mat_mspe),3)

#round(mean(my_asym_pars.mat_logs),3)
#round(mean(liz_asym_pars.mat_logs),3)
#round(mean(sym_pars.mat_logs),3)
#round(mean(my_asym_sep.mat_logs),3)
#round(mean(liz_asym_sep.mat_logs),3)
#round(mean(sym_sep.mat_logs),3)



####### Percentage Reductions using symmetric parsimonious matern as the base case ###########
mean.MSPE.M1<-(mean(my_asym_pars.mat_mspe))
mean.MSPE.M2<-(mean(liz_asym_pars.mat_mspe))
mean.MSPE.M3<-(mean(sym_pars.mat_mspe))
mean.MSPE.M4<-(mean(my_asym_sep.mat_mspe))
mean.MSPE.M5<-(mean(liz_asym_sep.mat_mspe))
mean.MSPE.M6<-(mean(sym_sep.mat_mspe))

sd.MSPE.M1<-(sd(my_asym_pars.mat_mspe))
sd.MSPE.M2<-(sd(liz_asym_pars.mat_mspe))
sd.MSPE.M3<-(sd(sym_pars.mat_mspe))
sd.MSPE.M4<-(sd(my_asym_sep.mat_mspe))
sd.MSPE.M5<-(sd(liz_asym_sep.mat_mspe))
sd.MSPE.M6<-(sd(sym_sep.mat_mspe))





mean.LogS.M1<-(mean(my_asym_pars.mat_logs))
mean.LogS.M2<-(mean(liz_asym_pars.mat_logs))
mean.LogS.M3<-(mean(sym_pars.mat_logs))
mean.LogS.M4<-(mean(my_asym_sep.mat_logs))
mean.LogS.M5<-(mean(liz_asym_sep.mat_logs))
mean.LogS.M6<-(mean(sym_sep.mat_logs))

sd.LogS.M1<-(sd(my_asym_pars.mat_logs))
sd.LogS.M2<-(sd(liz_asym_pars.mat_logs))
sd.LogS.M3<-(sd(sym_pars.mat_logs))
sd.LogS.M4<-(sd(my_asym_sep.mat_logs))
sd.LogS.M5<-(sd(liz_asym_sep.mat_logs))
sd.LogS.M6<-(sd(sym_sep.mat_logs))




base.mspe<-mean.MSPE.M3
base.LogS<-mean.LogS.M3

base.mspe.sd<-sd.MSPE.M3
base.LogS.sd<-sd.LogS.M3

##### percentage reductions in mean MSPE #####
round(((base.mspe-mean.MSPE.M1)/(abs(base.mspe)))*100,2)
round(((base.mspe-mean.MSPE.M2)/(abs(base.mspe)))*100,2)
round(((base.mspe-mean.MSPE.M3)/(abs(base.mspe)))*100,2)
round(((base.mspe-mean.MSPE.M4)/(abs(base.mspe)))*100,2)
round(((base.mspe-mean.MSPE.M5)/(abs(base.mspe)))*100,2)
round(((base.mspe-mean.MSPE.M6)/(abs(base.mspe)))*100,2)


round(((base.mspe.sd-sd.MSPE.M1)/(abs(base.mspe.sd)))*100,2)
round(((base.mspe.sd-sd.MSPE.M2)/(abs(base.mspe.sd)))*100,2)
round(((base.mspe.sd-sd.MSPE.M3)/(abs(base.mspe.sd)))*100,2)
round(((base.mspe.sd-sd.MSPE.M4)/(abs(base.mspe.sd)))*100,2)
round(((base.mspe.sd-sd.MSPE.M5)/(abs(base.mspe.sd)))*100,2)
round(((base.mspe.sd-sd.MSPE.M6)/(abs(base.mspe.sd)))*100,2)






##### percentage reductions in mean LogS #####
round(((base.LogS-mean.LogS.M1)/(abs(base.LogS)))*100,2)
round(((base.LogS-mean.LogS.M2)/(abs(base.LogS)))*100,2)
round(((base.LogS-mean.LogS.M3)/(abs(base.LogS)))*100,2)
round(((base.LogS-mean.LogS.M4)/(abs(base.LogS)))*100,2)
round(((base.LogS-mean.LogS.M5)/(abs(base.LogS)))*100,2)
round(((base.LogS-mean.LogS.M6)/(abs(base.LogS)))*100,2)


round(((base.LogS.sd-sd.LogS.M1)/(abs(base.LogS.sd)))*100,2)
round(((base.LogS.sd-sd.LogS.M2)/(abs(base.LogS.sd)))*100,2)
round(((base.LogS.sd-sd.LogS.M3)/(abs(base.LogS.sd)))*100,2)
round(((base.LogS.sd-sd.LogS.M4)/(abs(base.LogS.sd)))*100,2)
round(((base.LogS.sd-sd.LogS.M5)/(abs(base.LogS.sd)))*100,2)
round(((base.LogS.sd-sd.LogS.M6)/(abs(base.LogS.sd)))*100,2)






