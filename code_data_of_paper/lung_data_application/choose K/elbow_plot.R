load("lung_disease_our_top2000_chooseK.Rdata")

RMSE.all<-matrix(nrow=length(K.all),ncol=length(lambda.all))
R2.all<-matrix(nrow=length(K.all),ncol=length(lambda.all))
rownames(RMSE.all)=rownames(R2.all)=K.all
colnames(RMSE.all)=colnames(R2.all)=lambda.all
l=1
for(i in 1:length(K.all)){
  for(j in 1:length(lambda.all)){
    RMSE.all[i,j]<-res.all[[l]][["RMSE"]]
    R2.all[i,j]<-res.all[[l]][["R2"]]
    l=l+1
  }
}
RMSE.all
R2.all

par(mfrow=c(1,2))
for(i in 8){
  plot(x=rownames(R2.all),y=R2.all[,i],ylim=c(0,1),ylab="R2",xlab="K")
}
for(i in 8){
  plot(x=rownames(RMSE.all),y=RMSE.all[,i],ylim=c(0,0.25),ylab="RMSE",xlab="K")
}

