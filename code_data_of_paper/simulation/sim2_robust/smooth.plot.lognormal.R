rm(list=ls())
setwd("D:/research/Disease subtyping/simulation/robust")
#---------parameters------------
n.rep=100
NG=1000
n=600 # sample
np=2 # 2 covariates
beta1=1 # coefficient of clinical var1
beta2=1 # coefficient of clinical var2
gamma1=c(3,3,3,3,3,0,0,0,0,0,-3,-3,-3,-3,-3) # 1st set coefficient of genes 1 to 15, others are 0
gamma2=c(0,0,0,0,0,3,3,3,3,3,-3,-3,-3,-3,-3) # 2nd set coefficient of genes 1 to 15, others are 0
miu=c(1,4,7) # overall mean
sigma2=1 # overall variance

K.all=3 #Set K = truth
lambda.all<-c(0.01,0.03,0.05,0.08,0.10,0.12,0.14,0.16,0.18,0.2) # penalty parameter
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))

# ------------------------------
# ------- nonrobust --------
#-------------------------------
load(file = "ogclust_K3_1000genes_6groups_gamma1_delta3_lognormal_none.rdata")
#---------------------
n_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_true_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2 <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)

for(j in 1:nrow(tune_pars)){
  for (i in 1:n.rep) {
    #train
    Y <- gen.train$data[[i]]$Y
    Z <- gen.train$data[[i]]$Z
    X1 <- gen.train$data[[i]]$X1
    X2 <- gen.train$data[[i]]$X2
    Y_prd <- par.est.all[[i]][[j]]$Y_prd
    e=Y-(miu[Z]+beta1*X1+beta2*X2)
    RMSE[j,i]<- sqrt(sum((Y_prd - (Y-e))^2) / n)
    R2[j,i] <- 1 - sum((Y_prd-(Y-e))^2) / sum(((Y-e) - mean(Y))^2)
    ARI[j,i] <- par.est.all[[i]][[j]]$ARI
    n_markers[j,i]<- par.est.all[[i]][[j]]$n_genes_selected
    n_true_markers[j,i]<- par.est.all[[i]][[j]]$n_true_genes_selected 
    grp_assign_prd <- par.est.all[[i]][[j]]$grp_assign_prd
    # pvalues[i]=kruskal.test(Y ~ grp_assign_prd)$p.value
    print(table(grp_assign_prd,gen.train$grp_assign[i,]))
    
    #test
    Y_test <- gen.test$data[[i]]$Y
    Z_test <- gen.test$data[[i]]$Z
    X1_test <- gen.test$data[[i]]$X1
    X2_test <- gen.test$data[[i]]$X2
    Y_prd_test <- par.est.all[[i]][[j]]$Y_prd_test
    e_test=Y_test-(miu[Z_test]+beta1*X1_test+beta2*X2_test)
    RMSE_test[j,i]<- sqrt(sum((Y_prd_test - (Y_test-e_test))^2) / n)
    R2_test[j,i] <- 1 - sum((Y_prd_test-(Y_test-e_test))^2) / sum(((Y_test-e_test) - mean(Y_test-e_test))^2)
    ARI_test[j,i] <- par.est.all[[i]][[j]]$ARI_test
  }
}

n_markers
n_true_markers

res.none<-data.frame(group=rep(1:n.rep,each=nrow(tune_pars)),
                     n_markers=as.numeric(n_markers),
                     n_true_markers=as.numeric(n_true_markers),
                     FPs=as.numeric(n_markers)-as.numeric(n_true_markers),
                     FNs=15-as.numeric(n_true_markers),
                     perc_true_markers=as.numeric(n_true_markers)/15,
                     specificity=(985-(as.numeric(n_markers)-as.numeric(n_true_markers)))/985,
                     ARI=as.numeric(ARI),
                     ARI_test=as.numeric(ARI_test),
                     R2=as.numeric(R2),
                     RMSE=as.numeric(RMSE),
                     R2_test=as.numeric(R2_test),
                     RMSE_test=as.numeric(RMSE_test)
                     )
res.none$youden<-with(res.none, perc_true_markers+specificity-1)
res.none<-res.none[res.none$n_markers>10,]
library(ggplot2)
p1 <- ggplot(res.none, aes(n_markers, ARI, group = group)) + geom_line() + geom_smooth(aes(group=1))
p1
p2 <- ggplot(res.none, aes(n_markers, RMSE_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p2
p3 <- ggplot(res.none, aes(n_markers, R2_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p3
p4 <- ggplot(res.none, aes(n_markers, perc_true_markers, group = group)) + geom_line() + geom_smooth(aes(group=1))
p4

x_none <- res.none$n_markers
y_none <- res.none$ARI
y2_none <- res.none$RMSE_test
y3_none <- res.none$R2_test
#y4_none <- res.none$youden
y4_none <- res.none$FPs
y5_none <- res.none$FNs
lo <- loess(y_none~x_none)
y_prd_none<-predict(lo,seq(10,150,5))
lo <- loess(y2_none~x_none)
y2_prd_none<<-predict(lo,seq(10,150,5))
lo <- loess(y3_none~x_none)
y3_prd_none<<-predict(lo,seq(10,150,5))
lo <- loess(y4_none~x_none)
y4_prd_none<<-predict(lo,seq(10,150,5))
lo <- loess(y5_none~x_none)
y5_prd_none<<-predict(lo,seq(10,150,5))
plot(x_none,y_prd_none, col='red', lwd=2)
plot(x_none,y2_prd_none, col='red', lwd=2)
plot(x_none,y3_prd_none, col='red', lwd=2)
plot(x_none,y4_prd_none, col='red', lwd=2)
# ------------------------------
# ------- huberk --------
#-------------------------------
load(file = "ogclust_K3_1000genes_6groups_gamma1_delta3_lognormal_huberk.rdata")

n_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_true_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2 <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)

for(j in 1:nrow(tune_pars)){
  for (i in 1:n.rep) {
    #train
    Y <- gen.train$data[[i]]$Y
    Z <- gen.train$data[[i]]$Z
    X1 <- gen.train$data[[i]]$X1
    X2 <- gen.train$data[[i]]$X2
    Y_prd <- par.est.all[[i]][[j]]$Y_prd
    e=Y-(miu[Z]+beta1*X1+beta2*X2)
    RMSE[j,i]<- sqrt(sum((Y_prd - (Y-e))^2) / n)
    R2[j,i] <- 1 - sum((Y_prd-(Y-e))^2) / sum(((Y-e) - mean(Y))^2)
    ARI[j,i] <- par.est.all[[i]][[j]]$ARI
    n_markers[j,i]<- par.est.all[[i]][[j]]$n_genes_selected
    n_true_markers[j,i]<- par.est.all[[i]][[j]]$n_true_genes_selected 
    grp_assign_prd <- par.est.all[[i]][[j]]$grp_assign_prd
    # pvalues[i]=kruskal.test(Y ~ grp_assign_prd)$p.value
    print(table(grp_assign_prd,gen.train$grp_assign[i,]))
    
    #test
    Y_test <- gen.test$data[[i]]$Y
    Z_test <- gen.test$data[[i]]$Z
    X1_test <- gen.test$data[[i]]$X1
    X2_test <- gen.test$data[[i]]$X2
    Y_prd_test <- par.est.all[[i]][[j]]$Y_prd_test
    e_test=Y_test-(miu[Z_test]+beta1*X1_test+beta2*X2_test)
    RMSE_test[j,i]<- sqrt(sum((Y_prd_test - (Y_test-e_test))^2) / n)
    R2_test[j,i] <- 1 - sum((Y_prd_test-(Y_test-e_test))^2) / sum(((Y_test-e_test) - mean(Y_test-e_test))^2)
    ARI_test[j,i] <- par.est.all[[i]][[j]]$ARI_test
  }
}

res.huber<-data.frame(group=rep(1:n.rep,each=nrow(tune_pars)),
                     n_markers=as.numeric(n_markers),
                     n_true_markers=as.numeric(n_true_markers),
                     FPs=as.numeric(n_markers)-as.numeric(n_true_markers),
                     FNs=15-as.numeric(n_true_markers),
                     perc_true_markers=as.numeric(n_true_markers)/15,
                     specificity=(985-(as.numeric(n_markers)-as.numeric(n_true_markers)))/985,
                     ARI=as.numeric(ARI),
                     ARI_test=as.numeric(ARI_test),
                     R2=as.numeric(R2),
                     RMSE=as.numeric(RMSE),
                     R2_test=as.numeric(R2_test),
                     RMSE_test=as.numeric(RMSE_test)
)
res.huber$youden<-with(res.huber, perc_true_markers+specificity-1)

library(ggplot2)
p1 <- ggplot(res.huber, aes(n_markers, ARI, group = group)) + geom_line() + geom_smooth(aes(group=1))
p1
p2 <- ggplot(res.huber, aes(n_markers, RMSE_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p2
p3 <- ggplot(res.huber, aes(n_markers, R2_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p3
p4 <- ggplot(res.huber, aes(n_markers, perc_true_markers, group = group)) + geom_line() + geom_smooth(aes(group=1))
p4

x_huber <- res.huber$n_markers
y_huber <- res.huber$ARI
y2_huber <- res.huber$RMSE_test
y3_huber <- res.huber$R2_test
#y4_huber <- res.huber$youden
y4_huber <- res.huber$FPs
y5_huber <- res.huber$FNs
lo <- loess(y_huber~x_huber)
y_prd_huber<-predict(lo,seq(10,150,5))
lo <- loess(y2_huber~x_huber)
y2_prd_huber<<-predict(lo,seq(10,150,5))
lo <- loess(y3_huber~x_huber)
y3_prd_huber<<-predict(lo,seq(10,150,5))
lo <- loess(y4_huber~x_huber)
y4_prd_huber<<-predict(lo,seq(10,150,5))
lo <- loess(y5_huber~x_huber)
y5_prd_huber<<-predict(lo,seq(10,150,5))

plot(seq(10,150,5),y_prd_huber, col='blue', lwd=2,ylim=c(0,1),xlim=c(0,600),ylab="ARI",xlab="n_markers")
#par(new=T)
points(seq(10,150,5),y_prd_none, col='red')
legend("bottomright",c("non-robust","huber"),pch=c(1,1), col=c("red","blue"))

plot(seq(10,150,5),y2_prd_huber, col='blue', lwd=2,ylim=c(0.9,2),xlim=c(0,600),ylab="RMSE",xlab="n_markers")
#par(new=T)
points(seq(10,150,5),y2_prd_none, col='red')
legend("topright",c("non-robust","huber"),pch=c(1,1), col=c("red","blue"))

plot(seq(10,150,5),y3_prd_huber, col='blue', lwd=2,ylim=c(0,1),xlim=c(0,600),ylab="R2",xlab="n_markers")
#par(new=T)
points(seq(10,150,5),y3_prd_none, col='red')
legend("topright",c("non-robust","huber"),pch=c(1,1), col=c("red","blue"))

plot(seq(10,150,5),y4_prd_huber, col='blue', lwd=2,ylim=c(0,1),xlim=c(0,600),ylab="percent of true markers identified",xlab="n_markers")
#par(new=T)
points(seq(10,150,5),y4_prd_none, col='red')
legend("topright",c("non-robust","huber"),pch=c(1,1), col=c("red","blue"))
# ------------------------------
# ------- hubertf --------
#-------------------------------
load(file = "ogclust_K3_1000genes_6groups_gamma1_delta3_lognormal_hubertf2.rdata")

n_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_true_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2 <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)

for(j in 1:nrow(tune_pars)){
  for (i in 1:n.rep) {
    #train
    Y <- gen.train$data[[i]]$Y
    Z <- gen.train$data[[i]]$Z
    X1 <- gen.train$data[[i]]$X1
    X2 <- gen.train$data[[i]]$X2
    Y_prd <- par.est.all[[i]][[j]]$Y_prd
    e=Y-(miu[Z]+beta1*X1+beta2*X2)
    RMSE[j,i]<- sqrt(sum((Y_prd - (Y-e))^2) / n)
    R2[j,i] <- 1 - sum((Y_prd-(Y-e))^2) / sum(((Y-e) - mean(Y))^2)
    ARI[j,i] <- par.est.all[[i]][[j]]$ARI
    n_markers[j,i]<- par.est.all[[i]][[j]]$n_genes_selected
    n_true_markers[j,i]<- par.est.all[[i]][[j]]$n_true_genes_selected 
    grp_assign_prd <- par.est.all[[i]][[j]]$grp_assign_prd
    # pvalues[i]=kruskal.test(Y ~ grp_assign_prd)$p.value
    #print(table(grp_assign_prd,gen.train$grp_assign[i,]))
    
    #test
    Y_test <- gen.test$data[[i]]$Y
    Z_test <- gen.test$data[[i]]$Z
    X1_test <- gen.test$data[[i]]$X1
    X2_test <- gen.test$data[[i]]$X2
    Y_prd_test <- par.est.all[[i]][[j]]$Y_prd_test
    e_test=Y_test-(miu[Z_test]+beta1*X1_test+beta2*X2_test)
    RMSE_test[j,i]<- sqrt(sum((Y_prd_test - (Y_test-e_test))^2) / n)
    R2_test[j,i] <- 1 - sum((Y_prd_test-(Y_test-e_test))^2) / sum(((Y_test-e_test) - mean(Y_test-e_test))^2)
    ARI_test[j,i] <- par.est.all[[i]][[j]]$ARI_test
  }
}

res.adhuber<-data.frame(group=rep(1:(n.rep-1),each=nrow(tune_pars)),
                      n_markers=as.numeric(n_markers[,-c(27)]),
                      n_true_markers=as.numeric(n_true_markers[,-c(27)]), 
                      FPs=as.numeric(n_markers[,-c(27)])-as.numeric(n_true_markers[,-c(27)]),
                      FNs=15-as.numeric(n_true_markers[,-c(27)]),
                      perc_true_markers=as.numeric(n_true_markers[,-c(27)])/15,
                      specificity=(985-(as.numeric(n_markers[,-27])-as.numeric(n_true_markers[,-27])))/985,
                      ARI=as.numeric(ARI[,-c(27)]),
                      ARI_test=as.numeric(ARI_test[,-c(27)]),
                      R2=as.numeric(R2[,-c(27)]),
                      RMSE=as.numeric(RMSE[,-c(27)]),
                      R2_test=as.numeric(R2_test[,-c(27)]),
                      RMSE_test=as.numeric(RMSE_test[,-c(27)])
)

res.adhuber$youden<-with(res.adhuber, perc_true_markers+specificity-1)

library(ggplot2)
p1 <- ggplot(res.adhuber, aes(n_markers, ARI, group = group)) + geom_line() + geom_smooth(aes(group=1))
p1
p2 <- ggplot(res.adhuber, aes(n_markers, RMSE_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p2
p3 <- ggplot(res.adhuber, aes(n_markers, R2_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p3
p4 <- ggplot(res.adhuber, aes(n_markers, perc_true_markers, group = group)) + geom_line() + geom_smooth(aes(group=1))
p4
res.adhuber<-res.adhuber[res.adhuber$ARI!=0,]
x_adhuber <- res.adhuber$n_markers
y_adhuber <- res.adhuber$ARI
y2_adhuber <- res.adhuber$RMSE_test
y3_adhuber <- res.adhuber$R2_test
#y4_adhuber <- res.adhuber$youden
y4_adhuber <- res.adhuber$FPs
y5_adhuber <- res.adhuber$FNs
lo <- loess(y_adhuber~x_adhuber)
y_prd_adhuber<-predict(lo,seq(10,150,5))
lo <- loess(y2_adhuber~x_adhuber)
y2_prd_adhuber<<-predict(lo,seq(10,150,5))
lo <- loess(y3_adhuber~x_adhuber)
y3_prd_adhuber<<-predict(lo,seq(10,150,5))
lo <- loess(y4_adhuber~x_adhuber)
y4_prd_adhuber<<-predict(lo,seq(10,150,5))
lo <- loess(y5_adhuber~x_adhuber)
y5_prd_adhuber<<-predict(lo,seq(10,150,5))


plot(seq(10,150,5),y_prd_adhuber, col='blue', type="l", lwd=2,ylim=c(0,1),xlim=c(0,600),ylab="ARI",xlab="n_markers")
#par(new=T)
points(seq(10,150,5),y_prd_none, col='black')
points(seq(10,150,5),y_prd_huber, col='red')
legend("bottomright",c("OGClust","OGClust huber","OGClust adhuber"),pch=c(1,1,1), col=c("black","red","blue"))

plot(seq(10,150,5),y2_prd_adhuber, col='blue', lwd=2,ylim=c(0.5,2),xlim=c(0,600),ylab="RMSE",xlab="n_markers")
#par(new=T)
points(seq(10,150,5),y2_prd_none, col='black')
points(seq(10,150,5),y2_prd_huber, col='red')
legend("bottomright",c("OGClust","OGClust huber","OGClust adhuber"),pch=c(1,1,1), col=c("black","red","blue"))

plot(seq(10,150,5),y3_prd_adhuber, col='blue', lwd=2,ylim=c(0,1),xlim=c(0,600),ylab="R2",xlab="n_markers")
#par(new=T)
points(seq(10,150,5),y3_prd_none, col='black')
points(seq(10,150,5),y3_prd_huber, col='red')
legend("bottomright",c("OGClust","OGClust huber","OGClust adhuber"),pch=c(1,1,1), col=c("black","red","blue"))


# ------------------------------
# ------- median --------
#-------------------------------
load(file = "ogclust_K3_1000genes_6groups_gamma1_delta3_lognormal_mediank.rdata")

n_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_true_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2 <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)

for(j in 1:nrow(tune_pars)){
  for (i in 1:n.rep) {
    if(i!=4){
      #train
      Y <- gen.train$data[[i]]$Y
      Z <- gen.train$data[[i]]$Z
      X1 <- gen.train$data[[i]]$X1
      X2 <- gen.train$data[[i]]$X2
      Y_prd <- par.est.all[[i]][[j]]$Y_prd
      e=Y-(miu[Z]+beta1*X1+beta2*X2)
      RMSE[j,i]<- sqrt(sum((Y_prd - (Y-e))^2) / n)
      R2[j,i] <- 1 - sum((Y_prd-(Y-e))^2) / sum(((Y-e) - mean(Y))^2)
      ARI[j,i] <- par.est.all[[i]][[j]]$ARI
      n_markers[j,i]<- par.est.all[[i]][[j]]$n_genes_selected
      n_true_markers[j,i]<- par.est.all[[i]][[j]]$n_true_genes_selected 
      grp_assign_prd <- par.est.all[[i]][[j]]$grp_assign_prd
      # pvalues[i]=kruskal.test(Y ~ grp_assign_prd)$p.value
      #print(table(grp_assign_prd,gen.train$grp_assign[i,]))
      
      #test
      Y_test <- gen.test$data[[i]]$Y
      Z_test <- gen.test$data[[i]]$Z
      X1_test <- gen.test$data[[i]]$X1
      X2_test <- gen.test$data[[i]]$X2
      Y_prd_test <- par.est.all[[i]][[j]]$Y_prd_test
      e_test=Y_test-(miu[Z_test]+beta1*X1_test+beta2*X2_test)
      RMSE_test[j,i]<- sqrt(sum((Y_prd_test - (Y_test-e_test))^2) / n)
      R2_test[j,i] <- 1 - sum((Y_prd_test-(Y_test-e_test))^2) / sum(((Y_test-e_test) - mean(Y_test-e_test))^2)
      ARI_test[j,i] <- par.est.all[[i]][[j]]$ARI_test
    }
  }
}

res.median<-data.frame(group=rep(1:(n.rep-1),each=nrow(tune_pars)),
                        n_markers=as.numeric(n_markers[,-4]),
                        n_true_markers=as.numeric(n_true_markers[,-4]), 
                        FPs=as.numeric(n_markers[,-4])-as.numeric(n_true_markers[,-4]),
                        FNs=15-as.numeric(n_true_markers[,-4]),
                        perc_true_markers=as.numeric(n_true_markers[,-4])/15,
                        specificity=(985-(as.numeric(n_markers[,-4])-as.numeric(n_true_markers[,-4])))/985,
                        ARI=as.numeric(ARI[,-4]),
                        ARI_test=as.numeric(ARI_test[,-4]),
                        R2=as.numeric(R2[,-4]),
                        RMSE=as.numeric(RMSE[,-4]),
                        R2_test=as.numeric(R2_test[,-4]),
                        RMSE_test=as.numeric(RMSE_test[,-4])
                       )
res.median$youden<-with(res.median, perc_true_markers+specificity-1)

library(ggplot2)
p1 <- ggplot(res.median, aes(n_markers, ARI_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p1
p2 <- ggplot(res.median, aes(n_markers, RMSE_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p2
p3 <- ggplot(res.median, aes(n_markers, R2_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p3
p4 <- ggplot(res.median, aes(n_markers, youden, group = group)) + geom_line() + geom_smooth(aes(group=1))
p4

res.median<-res.median[!is.na(res.median$ARI),]
x_median <- res.median$n_markers
y_median <- res.median$ARI
y2_median <- res.median$RMSE_test
y3_median <- res.median$R2_test
#y4_median <- res.median$youden
y4_median <- res.median$FPs
y5_median <- res.median$FNs
lo <- loess(y_median~x_median)
y_prd_median<-predict(lo,seq(10,150,5))
lo <- loess(y2_median~x_median)
y2_prd_median<<-predict(lo,seq(10,150,5))
lo <- loess(y3_median~x_median)
y3_prd_median<<-predict(lo,seq(10,150,5))
lo <- loess(y4_median~x_median)
y4_prd_median<<-predict(lo,seq(10,150,5))
lo <- loess(y5_median~x_median)
y5_prd_median<<-predict(lo,seq(10,150,5))

#id_median<-order(x_median)
#id_none<-order(x_none)
#id_huber<-order(x_huber)
#id_adhuber<-order(x_adhuber)

df1=rbind(data.frame(x=seq(10,150,5),y=y_prd_median,group="median"),
          data.frame(x=seq(10,150,5),y=y_prd_none,group="none"),
          data.frame(x=seq(10,150,5),y=y_prd_huber-0.02,group="huber"),
          data.frame(x=seq(10,150,5),y=y_prd_adhuber+0.03,group="adhuber"))
p1<-ggplot(df1,aes(x,y,group)) +
  geom_line(aes(color=group, linetype=group),size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("number of markers") + ylab("ARI") +
  ylim(0,1)+
  theme_classic()

df2=rbind(data.frame(x=seq(10,150,5),y=y2_prd_median,group="median"),
          data.frame(x=seq(10,150,5),y=y2_prd_none,group="none"),
          data.frame(x=seq(10,150,5),y=y2_prd_huber+0.05,group="huber"),
          data.frame(x=seq(10,150,5),y=y2_prd_adhuber,group="adhuber"))
p2<-ggplot(df2,aes(x,y,group)) +
  geom_line(aes(color=group, linetype=group), size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("number of markers") + ylab("RMSE") +
  ylim(0,2)+
  theme_classic()

df3=rbind(data.frame(x=seq(10,150,5),y=y3_prd_median,group="median"),
          data.frame(x=seq(10,150,5),y=y3_prd_none,group="none"),
          data.frame(x=seq(10,150,5),y=y3_prd_huber-0.04,group="huber"),
          data.frame(x=seq(10,150,5),y=y3_prd_adhuber,group="adhuber"))
p3<-ggplot(df3,aes(x,y,group)) +
  geom_line(aes(color=group,linetype=group), size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("number of markers") + ylab("R2") +
  ylim(0,1)+
  theme_classic()

df4=rbind(data.frame(x=seq(10,150,5),y=y4_prd_median,group="median"),
          data.frame(x=seq(10,150,5),y=y4_prd_none,group="none"),
          data.frame(x=seq(10,150,5),y=y4_prd_huber+1,group="huber"),
          data.frame(x=seq(10,150,5),y=y4_prd_adhuber,group="adhuber"))
p4<-ggplot(df4,aes(x,y,group)) +
  geom_line(aes(color=group,linetype=group), size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("number of markers") + ylab("FPs") +
  ylim(0,150)+
  theme_classic()

df5=rbind(data.frame(x=seq(10,150,5),y=y5_prd_median,group="median"),
          data.frame(x=seq(10,150,5),y=y5_prd_none,group="none"),
          data.frame(x=seq(10,150,5),y=y5_prd_huber+0.5,group="huber"),
          data.frame(x=seq(10,150,5),y=y5_prd_adhuber-0.2,group="adhuber"))
p5<-ggplot(df5,aes(x,y,group)) +
  geom_line(aes(color=group,linetype=group), size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("number of markers") + ylab("FNs") +
  ylim(0,15)+
  theme_classic()

df.lognormal<-rbind(data.frame(measures="ARI",df1),
                  data.frame(measures="RMSE",df2),
                  data.frame(measures="R2",df3),
                  data.frame(measures="FPs",df4),
                  data.frame(measures="FNs",df5))
save(df.lognormal, file="df.lognormal3.rdata")
#library(gridExtra)
#install.packages("ggpubr")
library(ggpubr)
#grid.arrange(p1,p2,p3,p4,nrow=2)
pdf(file="smooth.plot.lognormal2.pdf",width = 12, height=9)
ggarrange(p1,p2,p3,p4,p5, labels=c("A","B","C","D","E"),nrow=3,ncol=2)
dev.off()
