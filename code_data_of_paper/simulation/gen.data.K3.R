gen.data.K3<-function(n.rep,NG,n,theta,cluster_index, exp=3, error="normal", outlier_typeI=FALSE, outlier_typeII=FALSE, propI=0.1, propII=0.1){
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
  grp_assign=matrix(ncol=n,nrow=n.rep)
  outliers.id<-matrix(ncol=n/10,nrow=n.rep)
  for(rep in 1:n.rep){
    X1=rnorm(n,1,1) # clinical var1
    X2=rnorm(n,1,0.5) # clinical var2
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
    if(error=="lognormal"){
      e=rlnorm(n,meanlog=0,sigma2)
      #e=ifelse(runif(n)>0.5,e,-e)
      e=e-mean(e)
    } else if(error == "pareto") {
      library(EnvStats)
      e=rpareto(n, location=1, shape=1)
      #e=ifelse(runif(n)>0.5,e,-e)
      e=e-mean(e)
    } else if(error == "t2"){
      e=rt(n,df=2)
    }else if(error == "normal"){
      e=rnorm(n,0,sigma2) 
    }
    Y=miu[Z]+beta1*X1+beta2*X2+e
    #Add outliers
    if(outlier_typeI==TRUE){
      # 10% subjects exchange their subgroup labels
      outliers.pos<-sample(1:n,n*propI)
      Z[outliers.pos]<-ifelse(Z[outliers.pos]==1, 2, 1)
      Y[outliers.pos]=miu[Z[outliers.pos]]+beta1*X1[outliers.pos]+beta2*X2[outliers.pos]+e[outliers.pos]
      outliers.id[rep,]<-outliers.pos
    }
    if(outlier_typeII==TRUE){
      #Add 10% outliers on Y only, in this case, they are outliers but not high leverage points
      outliers.pos<-sample(1:n,n*propII)
      Y[outliers.pos]<-runif(n*propII,range(Y)[1]-10,range(Y)[2]+10)
      outliers.id[rep,]<-outliers.pos
    }
    if(outlier_typeII=="outliers and high leverage"){
      #Add 10% outliers on both X and Y, in this case, they are both outliers and high leverage points
      outliers.pos<-sample(1:n,n*propII)
      X1[outliers.pos]<-runif(n*propII, -10,10)
      X2[outliers.pos]<-runif(n*propII,-10,10)
      Y[outliers.pos]<-runif(n*propII,-20,20)
      outliers.id[rep,]<-outliers.pos
    }
    data[[rep]]<-data.frame(X1,X2,G,Z,Y)
  }
  return(list(data=data, grp_assign=grp_assign, outliers_id=outliers.id))
}