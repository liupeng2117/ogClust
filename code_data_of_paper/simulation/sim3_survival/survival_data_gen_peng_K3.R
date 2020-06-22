rho_loglogistic <- function(time, z, zCoef, p) {
  lambda <- 1 / parametrize(z, zCoef, FUN = "exponential")
  lambda * p * (lambda * time) ^ (p - 1) / (1 + (lambda * time) ^ p)
}

#rho_loglogistic <- function(time, mu, p) {
#  lambda <- exp(mu)
#  lambda * p * (lambda * time) ^ (p - 1) / (1 + (lambda * time) ^ p)
#}

surv_gen=function(nProcess,mu,beta,X1,X2,sigma){
  #zcoe=c(mu,beta)
  res=simEventData(z = cbind(rep(1,nProcess), X1,X2),
                   zCoef = c(mu,beta), end=rep(100,nProcess),
                   recurrent = FALSE, relativeRisk="none", rho = rho_loglogistic, 
                   arguments = list(rho = list(p=1/sigma)))
  return(res)
}




survival_data_gen=function(n.rep,NG,n,np,theta,cluster_index,exp=1){
  #params
  beta1 = theta$beta1# coefficient of clinical var1
  beta2 = theta$beta2 # coefficient of clinical var2
  gamma1 = theta$gamma1 # 1st set of coefficients of gene 1 to 15
  gamma2 = theta$gamma2 # 2nd set of coefficient of gene 1 to 15
  miu = theta$miu # overall mean
  sigma2 = theta$sigma2 # overall variance
  
  #cluster index
  c1_index=cluster_index$c1_index
  c2_index=cluster_index$c2_index
  c3_index=cluster_index$c3_index
  c4_index=cluster_index$c4_index
  c5_index=cluster_index$c5_index
  c6_index=cluster_index$c6_index
  
  data<-list()
  grp_assign=matrix(,ncol=n,nrow=n.rep)
  delta_assign=matrix(,ncol=n,nrow=n.rep)
  for(rep in 1:n.rep){
    # T=c(rep(0,n/2),rep(1,n/2)) #treatment
    X1=rnorm(n,mean=1,sd=1) # clinical var1
    X2=rnorm(n,mean=0,sd=1) # clinical var2
    G1_5=matrix(rnorm(5*n,0,1),ncol=5) # genes 1 to 5
    G6_10=matrix(rnorm(5*n,0,1),ncol=5) # genes 6 to 10
    G11_15=matrix(rnorm(5*n,0,1),ncol=5) # genes 11 to 15
    G16_20=matrix(rnorm(5*n,0,1),ncol=5) # genes 1 to 5
    G21_25=matrix(rnorm(5*n,0,1),ncol=5) # genes 6 to 10
    G26_30=matrix(rnorm(5*n,0,1),ncol=5) # genes 11 to 15
    G1_5[c1_index,]=rnorm(5*n/3, 1, 1) # cluster 1 gene 1 to 5
    G6_10[c2_index,]=rnorm(5*n/3, 1, 1) # cluster 2 gene 6 to 10
    G11_15[c3_index,]=rnorm(5*n/3, 1, 1) # cluster 3 gene 11 to 15
    G16_20[c4_index,]=matrix(rnorm(5*n/3, exp, 1),ncol=5) # cluster 4 gene 11 to 15
    G21_25[c5_index,]=matrix(rnorm(5*n/3, exp, 1),ncol=5) # cluster 5 gene 16 to 20
    G26_30[c6_index,]=matrix(rnorm(5*n/3, exp, 1),ncol=5) # cluster 6 gene 21 to 25
    G31_100=matrix(rnorm((NG-30)*n,0,1),ncol=(NG-30))
    
    G=cbind(G1_5,G6_10,G11_15,G16_20,G21_25,G26_30,G31_100)
    vProb = cbind(exp(G[,1:15] %*% rep(0,15)), exp(G[,1:15] %*% gamma1), exp(G[,1:15] %*% gamma2))
    mChoices=t(apply(vProb, 1, rmultinom, n=1, size=1))
    Z=sapply(1:n,function(x) which(mChoices[x,]==1))
    grp_assign[rep,]=Z
    #calculate outcome by regression function
    colnames(G)<-paste0("G",1:NG)
    beta=c(beta1,beta2)
    
    idx1=which(Z==1)
    idx2=which(Z==2)
    idx3=which(Z==3)
    n1=length(idx1)
    Y1=surv_gen(n1,miu[1],beta,X1[idx1],X2[idx1],sigma)
    n2=length(idx2)
    Y2=surv_gen(n2,miu[2],beta,X1[idx2],X2[idx2],sigma)
    n3=length(idx3)
    Y3=surv_gen(n3,miu[3],beta,X1[idx3],X2[idx3],sigma)
    Y=rbind(Y1, Y2, Y3)
    Y[idx1,]<-Y1
    Y[idx2,]<-Y2
    Y[idx3,]<-Y3
    data[[rep]]<-data.frame(X1,X2,G,Z,Y)
    delta_assign[rep,]=Y$event
  }
  return(list(data=data,grp_assign=grp_assign,delta_assign=delta_assign))
}
