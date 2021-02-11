read.csv("N_temminckii.csv") -> leneage2

head(leneage2)
#################
# Delete data from Korea
#################

a<-(35.0519-33.3361)/(129.418-127.466)	
b<-35.0519-129.418*a					
leneage2_sub <-subset(leneage2,leneage2$lat<=a*leneage2$lon+b)
# leneage1$V4 <- paste(leneage2$number,"_",leneage1$V2,sep="")
# leng<-nrow(leneage1)

library("geosphere")
distm(as.matrix(leneage2_sub[,c(5,4)])) -> dist_point
dist_point <- dist_point/100000

write.csv(dist_point, "kawamutsu_dist_samplingpoints2.csv")