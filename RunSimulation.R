# test run

################
# data generation for ABC (Model G: point group)
################

source("simulation_model.R")
source("summary_functions.R")
DD<-read.csv("dist_pointgroup8.csv")

leng<-nrow(DD)
library(parallel)

cl  <- makeCluster(72)
ncl <- length(cl)
nsimu <-72 # number of parallel computation
nbl <- as.list(rep(nsimu %/% ncl, times=ncl)) 
if((rem <- nsimu %% ncl) > 0) nbl[1:rem] <- lapply(nbl[1:rem], "+", 1) # correction of the remainder

result.data <- NULL

for(n1 in 1:5000){
  results_simu=parLapply(cl,nbl,simulate,leng=leng,DD=DD)
  s_star<-matrix(0,nrow=sum(unlist(nbl)),ncol=leng+4);k<-1
  for(i in 1:length(nbl)){
    for(j in 1:nbl[[i]]){
      s_star[k,]<-results_simu[[i]][[j]]
      k<-k+1
    }
  }
  result.data <-rbind(result.data,s_star)
  if(n1 %% 10 == 0){
  	print(n1)
  	write.table(result.data,file="modelG_result.txt")
  }
}

last.col<-ncol(result.data)
result_a <- matrix(parRapply(cl=cl,x=result.data[,5:last.col],FUN=summary_autocorrelation,d=DD,c=c(0.005,0.01,0.02)),ncol=3,byrow=T)
result_t <- matrix(parRapply(cl=cl,x=result.data[,5:last.col],FUN=summary_templeton,d=DD),ncol=6,byrow=T)
para_stat<-cbind(result.data[,1:4],result_a,result_t)
write.table(para_stat,"modelG_stat.txt")
