# simulation with pointgroup8
####################
#  define the function "simulate"
####################

simulate<-function(nb,leng,DD){
  
  m1<- 0; while(m1 == 0){m1 <- runif(n=1, 0, 5)}
  s <- 0; while(s  == 0){s  <- runif(n=1, 0, 50)}
  alpha<-0; while(alpha<=0.5|alpha>=1){alpha<-runif(1,0.5,1)}
  inv<-sample(1:39,1)
  t.inv<-1.313/40*1000
  D1 <- pgamma(q=as.matrix(DD), shape = t.inv * m1 / s, scale = s, lower.tail = F)

  Result<-NULL
  
  #####################
  #   define the function "sample_simu"  "sample_lineage"
  #####################
  
  sample_simu<-function(x, size, replace = FALSE, prob = NULL){
    if(any(prob %in% 1)){
      rep(x[prob==1],size)
    }else{
      sample(x,size,replace=replace,prob=prob)
    }
  }
  
  sample_lineage<-function(x,size=1,replace=F,prob=NULL){
    if(sum(!is.na(x))==0){
      return(NA)
    }else{
      if(sum(!is.na(x))==1){
        return(x[!is.na(x)])
      }
      else{
        x<-x[!is.na(x)]
        sample(x,size=1)
      }
    }
  }
  
  samplem <- function(p){
    y<-c(1,0)
    sample(y,1,prob=c(p,1-p))
  }
  
  #1,2 
  p1<-(1-alpha) ;p2<-alpha
  #1,4
  const<-((1-alpha)^2+alpha^2)
  p3<-(1-alpha)^2/const ;p4<-alpha^2/const
  #2,4
  p5<-(1-alpha);p6<-alpha
  #1,2,4
  const<-((1-alpha)^2 + alpha*(1-alpha) + alpha^2)
  p7<-(1-alpha)^2/const; p8<-alpha*(1-alpha)/const; p9<-alpha^2/const
  
  for(l in 1:nb){
    place<-rep(1,leng); place[c(200,201,202)]<-2

    for(x in 1:inv){
      place1<-numeric(leng) #アルゴリズムのため、ベクトル作成
      place2<-numeric(leng) #アルゴリズムのため、ベクトル作成
      prob_mat1=matrix(0,leng,leng)
      prob_mat1<-apply(D1,c(1,2),FUN=samplem)

      iuse1=place==1 ##1に占拠されている場所
      juse1=apply(as.matrix(prob_mat1[,iuse1],c(leng,sum(iuse1))),MARGIN=1,FUN=sum) >= 1 ##少なくとも一か所から1が流れてきた場所
      place1[juse1]=1
      iuse2=place==2  ##2に占拠されている場所
      juse2=apply(as.matrix(prob_mat1[,iuse2],c(leng,sum(iuse2))),MARGIN=1,FUN=sum) >= 1 ##少なくとも一か所から2が流れてきた場所
      place2[juse2]=2
      
      #置きかわりを考える
      place=place1+place2
      pos_12 = (place1-1)+(place2-2)==0 ##1と2が食い合う場所（論理積？）
      s12=sample(c(1,2),sum(pos_12),prob=c(p1,p2),replace=T)
      place[pos_12]=s12
    }
    
    place[c(200,201,202)]<-4

    for(x in 1:(40-inv)){
      place1<-numeric(leng) #アルゴリズムのため、ベクトル作成
      place2<-numeric(leng) #アルゴリズムのため、ベクトル作成
      place4<-numeric(leng) #アルゴリズムのため、ベクトル作成
      
      prob_mat1=matrix(0,leng,leng)
      prob_mat1<-apply(D1,c(1,2),FUN=samplem)
      
      iuse1=place==1 ##1に占拠されている場所
      juse1=apply(as.matrix(prob_mat1[,iuse1],c(leng,sum(iuse1))),MARGIN=1,FUN=sum) >= 1 ##少なくとも一か所から1が流れてきた場所
      place1[juse1]=1
      
      iuse2=place==2  ##2に占拠されている場所
      juse2=apply(as.matrix(prob_mat1[,iuse2],c(leng,sum(iuse2))),MARGIN=1,FUN=sum) >= 1 ##少なくとも一か所から2が流れてきた場所
      place2[juse2]=2
      
      iuse4=place==4  ##2に占拠されている場所
      juse4=apply(as.matrix(prob_mat1[,iuse4],c(leng,sum(iuse4))),MARGIN=1,FUN=sum) >= 1 ##少なくとも一か所から2が流れてきた場所
      place4[juse4]=4
      
      place=place1+place2+place4
      pos_12=place==3
      pos_14=place==5
      pos_24=place==6
      pos_124=place==7
      
      s12=sample(c(1,2),sum(pos_12),prob=c(p1,p2),replace=T)
      s24=sample(c(2,4),sum(pos_24),prob=c(p5,p6),replace=T)
      s14=sample(c(1,4),sum(pos_14),prob=c(p3,p4),replace=T)
      s124=sample(c(1,2,4),sum(pos_124),prob=c(p7,p8,p9),replace=T)
      
      place[pos_12]=s12
      place[pos_14]=s14
      place[pos_24]=s24
      place[pos_124]=s124
    }
    
    
    ###############
    # funal output
    ###############
    
    Result[[l]]<-c(m1, s, alpha, inv, place)
  }
  Result #return statistics
  
}
