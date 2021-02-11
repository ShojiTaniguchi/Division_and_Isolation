library(parallel)
library(abc)

source("summary_functions.R")
DD<-read.csv("dist_pointgroup8.csv")
leng<-nrow(DD)
leneage<-read.csv("dat_in_japan.csv")
place<-as.vector(leneage$clade)
place[place=="c"]<-1; place[place=="f"]<-2; place[place=="g"]<-4
place<-as.numeric(place)
leng<-length(place)

DD <- read.csv("kawamutsu_dist_samplingpoints2.csv")[,2:(1+leng)]
DD <- DD*100
stat.obs <-c(summary_quant(p=place,l=leneage[,5]),summary_autocorrelation(p=place,d=DD,c=c(0.005,0.01,0.02)),
            +             summary_templeton(p=place,d=DD),summary_concordance(p=place,obs=place))

#############
# model S
#############
para_stat<-read.table("modelG_stat.txt")
dim(para_stat)
param.set<-para_stat[,1:4]
colnames(param.set)<-c("m","s","alpha","r")
stati.set<-para_stat[,7:13]
abc.result<-abc(target=stat.obs[c(7,8,9,10,11,12,13)], param=param.set,sumstat=stati.set,tol=.025,method="neuralnet")
abc.summary<-summary(abc.result)

# Timing of the migration (t Ma)
timing <- 1.313 * (40 - abc.summary[,4])/40
timing