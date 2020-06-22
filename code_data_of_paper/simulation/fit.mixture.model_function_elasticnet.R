fit.mixture.model<-function(n, K, np, NG, lambda, alpha,df, G, Y, X, robust="none", tau=1.345){
  #inital values of theta
  beta_int=runif(np,0,3)
  gamma_int=runif((K-1)*(NG+1),0,1)
  beta0_int=runif(K,0,3)
  sigma2_int=runif(1,1,3)
  
  theta_int=c(beta_int, 
              gamma_int, 
              beta0_int, sigma2_int)
  theta_est<-EM(theta_int, lambda=lambda, df, np=np, K=K, NG=NG, alpha=alpha, robust=robust,tau=tau)

  #estimated parameters
  beta_est=theta_est[1:np]
  gamma_est=theta_est[(np+1):((K-1)*(NG+1)+np)]
  beta0_est=theta_est[((K-1)*(NG+1)+np+1):((K-1)*(NG+1)+np+K)]
  sigma2_est=theta_est[((K-1)*(NG+1)+np+K+1)]
  
  gamma_est_matrix=matrix(gamma_est,ncol=K-1,byrow=T)
  gamma_est_matrix=cbind(gamma_est_matrix,0)
  G=cbind(1,G)
  pai_est=sapply(1:K, function(k) 
    exp(G %*% gamma_est_matrix[,k,drop=F])/rowSums(exp(G %*% gamma_est_matrix)))
  f_est<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_est))*exp(-(Y-beta0_est[x]-X %*% beta_est)^2/(2*sigma2_est)))
  (ll=sum(log(diag(pai_est %*% t(f_est)))))
  #calculate the expected value of Y and R2
  Y_prd=apply(sapply(1:K, function(x) pai_est[,x]*(beta0_est[x]+X %*% beta_est)),1,sum)
  R2=1-sum((Y-Y_prd)^2)/sum((Y-mean(Y))^2)
  
  #Calculate AIC BIC
  AIC=2*sum(theta_est!=0)-2*ll
  BIC=log(n)*sum(theta_est!=0)-2*ll

  #EBIC=BIC+(1-1/(2*log(length(theta_est),base=n)))*choose(length(theta_est), sum(theta_est!=0))
  
  #prosterior prob
  w_est=sapply(1:K, function(k) 
    (pai_est[,k]*f_est[,k])/diag(pai_est %*% t(f_est)))
  cl.assign<-apply(w_est,1,which.max)
  
  return(list(res=c(theta_est,ll,R2,AIC,BIC,lambda),
              theta_est=theta_est,
              Y_prd=Y_prd,
              grp_assign=cl.assign))
}
