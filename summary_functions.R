##############
# define summary function
##############

# p; place, l; longitude
summary_quant<-function(p,l){
  e_quanti1<-mean(l[p==1]); e_quanti2<-mean(l[p==2]); e_quanti4<-mean(l[p==4])
  if(is.na(e_quanti1)){e_quanti1<-min(l)}
  if(is.na(e_quanti2)){e_quanti2<-min(l)}
  if(is.na(e_quanti4)){e_quanti4<-min(l)}
  
  between_class<-c(var(l[p==1]),var(l[p==2]),var(l[p==4]))
  within_class<-c(mean(l[p==1]),mean(l[p==2]),mean(l[p==4]))
  between_class[is.na(between_class)]<-0
  within_class[is.na(within_class)]<-0
  between_within_var<-log(var(within_class)/var(between_class))
  quant<-as.vector(c(e_quanti1,e_quanti2,e_quanti4,between_within_var))
  names(quant)<-c("mean1","mean2","mean4","overlap")
  return(quant)
}

# p; place, d; geographic distance, c; variable
summary_autocorrelation<-function(p,d,c){
  genetic_dist<-matrix(1,nrow=length(p),ncol=length(p))
  genetic_dist[p==1,p==1]<-0; genetic_dist[p==2,p==2]<-0; genetic_dist[p==4,p==4]<-0
  auto<-NULL
  for(c1 in c){
    auto1 <- sum( exp(-c1*as.vector(d)) * as.vector(genetic_dist) ) * 1/sum(exp(-c1*as.vector(d)))
    auto<-append(auto,auto1)
  }
  names(auto) <- paste("auto",as.character(c),sep="")
  return(auto)
}

# p; place, d; geographic distance
summary_templeton<-function(p,d){
  n_1<-sum(p==1); n_2<-sum(p==2); n_4<-sum(p==4)
  d1<-d[p==1,p==1];d2<-d[p==2,p==2];d4<-d[p==4,p==4]
  pnumber<-1:length(p)
  
  # center_0
  eachdist<-apply(d,FUN=sum,MARGIN=1)
  mindist<-eachdist==min(eachdist)
  center_0<-pnumber[mindist]
  if(sum(mindist)>1){
    crow<-sample(1:sum(mindist),1)
    center_0<-center_0[crow]
  }
  
  # center_1
  if(n_1>0){
    if(n_1>1){
      eachdist<-apply(d1,FUN=sum,MARGIN=1)
    }else{
      eachdist<-d1
    }
    
    mindist<-eachdist==min(eachdist)
    center_1<-(1:n_1)[mindist]
    if(sum(mindist)>1){
      crow<-sample(1:sum(mindist),1)
      center_1<-center_1[crow]
    }
  }else{
    center_1<-NA
  }
  
  # center_2
  if(n_2>0){
    if(n_2>1){
      eachdist<-apply(d2,FUN=sum,MARGIN=1)
    }else{
      eachdist<-d2
    }
    
    mindist<-eachdist==min(eachdist)
    center_2<-(1:n_2)[mindist]
    if(sum(mindist)>1){
      crow<-sample(1:sum(mindist),1)
      center_2<-center_2[crow]
    }
  }else{
    center_2<-NA
  }
  
  # center_4
  if(n_4>0){
    if(n_4>1){
      eachdist<-apply(d4,FUN=sum,MARGIN=1)
    }else{
      eachdist<-d4
    }
    
    mindist<-eachdist==min(eachdist)
    center_4<-(1:n_4)[mindist]
    if(sum(mindist)>1){
      crow<-sample(1:sum(mindist),1)
      center_4<-center_4[crow]
    }
  }else{
    center_4<-NA
  }
  
  # clade distance D_c1, D_c2, D_c4
  d1<-as.matrix(d1)
  if(n_1>0){
    D_c1<-mean(d1[center_1,])
  }else{
    D_c1<-0
  }
  
  d2<-as.matrix(d2)
  if(n_2>0){
    D_c2<-mean(d2[center_2,])
  }else{
    D_c2<-0
  }
  
  d4<-as.matrix(d4)
  if(n_4>0){
    D_c4<-mean(d4[center_4,])
  }else{
    D_c4<-0
  }
  
  # nested distance D_n1, D_n2, D_n4
  center_1<-(pnumber[p==1])[center_1]; center_2<-(pnumber[p==2])[center_2]; center_4<-(pnumber[p==4])[center_4]
  D_n1<-d[center_1,center_0]; D_n2<-d[center_2,center_0]; D_n4<-d[center_4,center_0]
  if(is.na(D_n1)==TRUE){D_n1<-0}; if(is.na(D_n2)==TRUE){D_n2<-0}; if(is.na(D_n4)==TRUE){D_n4<-0}
  
  templeton<-(c(D_c1,D_n1,D_c2,D_n2,D_c4,D_n4))
  names(templeton)<-c("D_c1","D_n1","D_c2","D_n2","D_c4","D_n4")
  return(templeton)
}

# p; simukation distribution, obs; observed distribution
summary_concordance<-function(p,obs){
  concordance_rate<-sum(p==obs)/length(p)
  names(concordance_rate)<-"conc_rate"
  return(concordance_rate)
}