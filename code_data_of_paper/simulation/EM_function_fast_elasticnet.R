EM<-function(theta, lambda, data, np, K, NG, robust, alpha, tau=1.345){
  X=as.matrix(data[,1:np])
  G=cbind(1,as.matrix(data[,(np+1):(np+NG)]))
  Y=data$Y
  
  l=1
  repeat{
    #print(l)
    beta_old=theta[1:np]
    gamma_old=theta[(1+np):((K-1)*(NG+1)+np)]
    miu_old=theta[((K-1)*(NG+1)+np+1):((K-1)*(NG+1)+np+K)]
    sigma2_old=theta[((K-1)*(NG+1)+np+K+1)]
    
    #==E-STEP==#
    gamma_old_matrix=matrix(gamma_old,ncol=K-1,byrow=T)
    gamma_old_matrix=cbind(gamma_old_matrix,0)
    pai_old=sapply(1:K, function(k) 
      exp(G %*% gamma_old_matrix[,k,drop=F])/rowSums(exp(G %*% gamma_old_matrix)))
    #pai_old=cbind(pai_old,1-rowSums(pai_old))
#    if(robust=="none"){
      f_old<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_old))*exp(-(Y-miu_old[x] - X %*% beta_old)^2/(2*sigma2_old)))
#    } else if(robust %in% c("huberk","hubertf","hubertf2","huber")){
#      e=sapply(1:K, function(x) Y-miu_old[x]-X%*%beta_old)
#      f_old<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_old))*exp(-ifelse(abs(e[,x])<=tau,e[,k]^2,tau*abs(e[,k])-0.5*tau^2)/(2*sigma2_old)))
      
#    } else if(robust == "mediank"){
#      e=sapply(1:K, function(x) Y-miu_old[x]-X%*%beta_old)
#      f_old<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_old))*exp(-ifelse(abs(e[,x])<=median(e[,x]),e[,k]^2,0)/(2*sigma2_old)))
      
#    }
    
    
    #calculate the expected value of Z
    w_old=sapply(1:K, function(k) 
      (pai_old[,k]*f_old[,k])/diag(pai_old %*% t(f_old)))
    #==M-STEP==#
    library(glmnet)
    gamma_new_matrix = tryCatch({
      fit<-glmnet(x=G[,-1], y=w_old, lambda=lambda, family = "multinomial", alpha = alpha, type.multinomial = "grouped")
      gamma_new_matrix=rbind(t(fit$a0),sapply(1:K, function(x) as.numeric(fit$beta[[x]])))
      gamma_new_matrix=sapply(1:K, function(x) gamma_new_matrix[,x]-gamma_new_matrix[,K])
    }, error =function(e) {
      return(gamma_old_matrix)
    })
    

    
    if(robust=="none"){
      #update miu, alpha
      miu_new=sapply(1:K, function(k) {sum(w_old[,k]*(Y-X%*%beta_old),na.rm=T)/sum(w_old[,k],na.rm=T)})
      #update beta
      beta_new<-beta_old
      for(i in 1:np){
        beta_new[i]<-sum(sapply(1:K, function(k) w_old[,k]*X[,i]*(Y-miu_new[k]-X[,-i,drop=F] %*% beta_new[-i])),na.rm=T)/sum(w_old*X[,i]^2,na.rm=T)
      }
      #update sigma2
      sigma2_new=sum(sapply(1:K, function(k) w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2),na.rm=T)/n
    }
    
    
    if(robust=="median1"){
      e=Y - apply(sapply(1:K, function(k) miu_old[k]+X%*%beta_old) * pai_old,1,sum)
      #update miu, alpha
      miu_new=sapply(1:K, function(k) sum((w_old[,k]*(Y-X%*%beta_old))[abs(e)<=median(abs(e))],na.rm=T)/sum(w_old[abs(e)<=median(abs(e)),k],na.rm=T))
      #update beta1 beta2
      beta_new<-beta_old
      for(i in 1:np){
        beta_new[i]=sum(sapply(1:K, function(k) (w_old[,k]*X[,i]*(Y-miu_new[k]-X[,-i, drop=F]%*%beta_new[-i]))[abs(e)<=median(abs(e))]),na.rm=T)/sum(sapply(1:K, function(k) (w_old[,k]*X[,i]^2)[abs(e)<=median(abs(e))]),na.rm=T)
      }
      #update sigma2
      sigma2_new=sum(sapply(1:K, function(k) (w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2)[abs(e)<=median(abs(e))]),na.rm=T)/sum(w_old[abs(e)<=median(abs(e)),],na.rm=T)
    }
    if(robust=="mediank"){
      e=sapply(1:K, function(k) Y-miu_old[k]-X%*%beta_old)
      #update miu, alpha
      miu_new=sapply(1:K, function(k) sum((w_old[,k]*(Y-X%*%beta_old))[abs(e[,k])<=median(abs(e[,k]))],na.rm=T)/sum(w_old[abs(e[,k])<=median(abs(e[,k])),k],na.rm=T))
      #update beta1 beta2
      beta_new<-beta_old
      for(i in 1:np){
        beta_new[i]=sum(sapply(1:K, function(k) (w_old[,k]*X[,i]*(Y-miu_new[k]-X[,-i,drop=F]%*%beta_new[-i]))[abs(e[,k])<=median(abs(e[,k]))]),na.rm=T)/sum(sapply(1:K, function(k) (w_old[,k]*X[,i]^2)[abs(e[,k])<=median(abs(e[,k]))]),na.rm=T)
      }
      #update sigma2
      sigma2_new=sum(sapply(1:K, function(k) (w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2)[abs(e[,k])<=median(abs(e[,k]))]),na.rm=T)/sum(sapply(1:K, function(k) w_old[abs(e[,k])<=median(abs(e[,k])),k]),na.rm=T)
      
    }
    if(robust=="hubertf"){
      beta_new=beta_old
      miu_new=miu_old
      grp.id<-apply(pai_old,1,which.max)
      for(k in 1:K){
        if(sum(grp.id==k)>1){
          Y_k=Y[grp.id==k]
          X_k=X[grp.id==k,]
          library(tfHuber)
          listHuber=huberReg(X_k, Y_k)
          miu_new[k]<-listHuber$theta[1]
          beta_new<-listHuber$theta[-1]
        }
      }
      tau=listHuber$tauCoef
      e=sapply(1:K, function(k) Y-miu_new[k]-X%*%beta_new)
      sigma2_new=sum(sapply(1:K, function(k) c((w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2)[abs(e[,k])<=tau],(w_old[,k]*(2*tau*abs(Y-miu_new[k]-X%*%beta_new)-tau^2))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[,k])<=tau,],na.rm=T)),na.rm=T)
      #e=Y-rowSums(sapply(1:K,function(k) pai_old[,k]*(miu_new[k]+X%*%beta_new)))
      #sigma2_new=sum(ifelse(abs(e)<tau, e^2,2*abs(e)*tau-tau^2))/(n-(np+K+(K-1)))
     
       #sigma2_new=sum(e^2)/(n-(np+K+(K-1)))
      #sigma2_new=(median(abs(e-median(e)))/qnorm(0.75))^2
    }
    if(robust=="hubertf2"){
      beta_new=beta_old
      miu_new=miu_old
      grp.id<-apply(pai_old,1,which.max)-1
      X_tf<-cbind(grp.id,X)
      library(tfHuber)
      listHuber=huberReg(X_tf, Y)
      mm<-diag(K)
      mm[1,]<-1
      miu_new<-listHuber$theta[1:K] %*% mm
      beta_new<-listHuber$theta[-c(1:K)]
      tau=listHuber$tauCoef
      e=sapply(1:K, function(k) Y-miu_new[k]-X%*%beta_new)
      sigma2_new=sum(sapply(1:K, function(k) c((w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2)[abs(e[,k])<=tau],(w_old[,k]*(2*tau*abs(Y-miu_new[k]-X%*%beta_new)-tau^2))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[,k])<=tau,],na.rm=T)),na.rm=T)
      #e=Y-rowSums(sapply(1:K,function(k) pai_old[,k]*(miu_new[k]+X%*% beta_new)))
      #sigma2_new=sum(ifelse(abs(e)<tau, e^2,2*abs(e)*tau-tau^2))/(n-(np+K+(K-1)))
     
      #sigma2_new=sum(e^2)/(n-(np+K+(K-1)))
      #sigma2_new=(median(abs(e-median(e)))/qnorm(0.75))^2
    }
    
    if(robust=="huberad"){ #tau is the same in k clusters
      ll=1
      repeat{
        par_old=c(beta_old,miu_old,sigma2_old,tau)
        #residuals
        e=sapply(1:K, function(k) Y-miu_old[k]-X%*%beta_old)
        rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        #update tau
        tau=sqrt(sum(ifelse(rr^2>tau^2,tau^2,rr^2))/(np+log(n)))
        #update miu, alpha
        miu_new=sapply(1:K, function(k) sum(c((w_old[,k]*(Y-X%*%beta_old))[abs(e[,k])<=tau],(w_old[,k]*tau*sign(e[,k]))[abs(e[,k])>tau]),na.rm=T)/sum(w_old[abs(e[,k])<=tau,k],na.rm=T))
        #update beta
        beta_new=beta_old
        for(i in 1:np){
          beta_new[i]=sum(sapply(1:K, function(k) c((w_old[,k]*X[,i]*(Y-miu_new[k]-X[,-i,drop=F]%*%beta_new[-i]))[abs(e[,k])<=tau],(w_old[,k]*tau*sign(e[,k]))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum((w_old[,k]*X[,i]^2)[abs(e[,k])<=tau],na.rm=T)),na.rm=T)
        }
        #update sigma2
        #e=sapply(1:K, function(k) Y-miu_new[k]-X%*%beta_new)
        #rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        #sigma2_new=(median(abs(rr-median(rr)))/qnorm(0.75))^2
        sigma2_new=sum(sapply(1:K, function(k) c((w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2)[abs(e[,k])<=tau],(w_old[,k]*(2*tau*abs(Y-miu_new[k]-X%*%beta_new)-tau^2))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[,k])<=tau,],na.rm=T)),na.rm=T)
        #e=sapply(1:K, function(k) Y-miu_new[k]-X%*%beta_new)
        #tau=sqrt(sum(sapply(1:K, function(k) w_old[,k]*ifelse(e[,k]^2>tau^2,tau^2,e[,k]^2)))/(np+log(n)))
        #rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        #tau=sqrt(sum(ifelse(rr^2>tau^2,tau^2,rr^2))/(np+log(n)))
        beta_old=beta_new
        miu_old=miu_new
        sigma2_old=sigma2_new
        par_new=c(beta_new,miu_new,sigma2_new,tau)
        dis_par=sqrt(sum((par_new-par_old)^2,na.rm=T))
        ll=ll+1
        if(dis_par<1*10^(-6) |ll>100){
          break
        }
      }
    }
    
    if(robust=="huberad3"){ #tau is the same in k clusters
      ll=1
      repeat{
        par_old=c(beta_old,miu_old,sigma2_old,tau)
        #residuals
        e=sapply(1:K, function(k) Y-miu_old[k]-X%*%beta_old)
        rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        #update tau
        tau=sqrt(sum(ifelse(rr^2>tau^2,tau^2,rr^2))/(np+log(n)))
        #update miu, alpha
        group.id=apply(w_old,1,which.max)
        miu2_old=miu_old[2]-miu_old[1]
        miu_new=miu_old
        #miu_new=sapply(1:K, function(k) sum(c((Y[group.id==k]-X[group.id==k,] %*% beta_old)[abs(rr[group.id==k])<=tau],(tau*sign(rr[group.id==k]))[abs(rr[group.id==k])>tau]),na.rm=T)/sum(abs(rr[group.id==k])<=tau))
        miu_new[1]=sum(c((Y- miu2_old * (group.id-1) -X %*% beta_old)[abs(rr)<=tau],(tau*sign(rr))[abs(rr)>tau]),na.rm=T)/sum(abs(rr<=tau))
        miu2_new=sum(c((Y[group.id==2]-miu_new[1]-X[group.id==2,] %*% beta_old)[abs(rr[group.id==2])<=tau],(tau*sign(rr[group.id==2]))[abs(rr[group.id==2])>tau]),na.rm=T)/sum(abs(rr[group.id==2])<=tau)
        miu_new[2]=miu_new[1]+miu2_new#update beta
        beta_new=beta_old
        for(i in 1:np){
          beta_new[i]=sum(unlist(sapply(1:K, function(k) 
            c((X[group.id==k,i]*(Y[group.id==k]-miu_new[k]-X[group.id==k,-i,drop=F]%*%beta_new[-i]))[abs(rr[group.id==k])<=tau],(tau*sign(rr[group.id==k]))[abs(rr[group.id==k])>tau]))),na.rm=T)/sum(sapply(1:K, function(k) sum((X[group.id==k,i]^2)[abs(rr[group.id==k])<=tau],na.rm=T)),na.rm=T)
        }
        #update sigma2
        #e=sapply(1:K, function(k) Y-miu_new[k]-X%*%beta_new)
        #rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        #sigma2_new=(median(abs(rr-median(rr)))/qnorm(0.75))^2
        #sigma2_new=sum(sapply(1:K, function(k) c((w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2)[abs(e[,k])<=tau],(w_old[,k]*(2*tau*abs(Y-miu_new[k]-X%*%beta_new)-tau^2))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[,k])<=tau,],na.rm=T)),na.rm=T)
        e=sapply(1:K, function(k) Y-miu_new[k]-X %*% beta_new)
        rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        sigma2_new=sum(rr^2)/(n-(np+K+(K-1)))
        #e=sapply(1:K, function(k) Y-miu_new[k]-X%*%beta_new)
        #tau=sqrt(sum(sapply(1:K, function(k) w_old[,k]*ifelse(e[,k]^2>tau^2,tau^2,e[,k]^2)))/(np+log(n)))
        #rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        #tau=sqrt(sum(ifelse(rr^2>tau^2,tau^2,rr^2))/(np+log(n)))
        beta_old=beta_new
        miu_old=miu_new
        sigma2_old=sigma2_new
        par_new=c(beta_new,miu_new,sigma2_new,tau)
        dis_par=sqrt(sum((par_new-par_old)^2,na.rm=T))
        ll=ll+1
        if(dis_par<1*10^(-6) |ll>100){
          break
        }
      }
    }
    
    if(robust=="huberad2"){ #tau is the same in k clusters
      ll1=1
      repeat{
        par_old=c(beta_old,miu_old,sigma2_old,tau)
        #residual
        e=sapply(1:K, function(k) Y-miu_old[k]-X%*%beta_old)
        rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        #update tau
        tau=1.345*median(abs(rr-median(rr)))/qnorm(0.75)
        #update miu, alpha
        miu_new=sapply(1:K, function(k) sum(c((w_old[,k]*(Y-X%*%beta_old))[abs(e[,k])<=tau],(w_old[,k]*tau*sign(e[,k]))[abs(e[,k])>tau]),na.rm=T)/sum(w_old[abs(e[,k])<=tau,k],na.rm=T))
        #update beta
        beta_new=beta_old
        for(i in 1:np){
          beta_new[i]=sum(sapply(1:K, function(k) c((w_old[,k]*X[,i]*(Y-miu_new[k]-X[,-i,drop=F]%*%beta_new[-i]))[abs(e[,k])<=tau],(w_old[,k]*tau*sign(e[,k]))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum((w_old[,k]*X[,i]^2)[abs(e[,k])<=tau],na.rm=T)),na.rm=T)
        }
        #update sigma2
        #e=sapply(1:K, function(k) Y-miu_new[k]-X%*%beta_new)
        #rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        #sigma2_new=(median(abs(rr-median(rr)))/qnorm(0.75))^2
        sigma2_new=sum(sapply(1:K, function(k) c((w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2)[abs(e[,k])<=tau],(w_old[,k]*(2*tau*abs(Y-miu_new[k]-X%*%beta_new)-tau^2))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[,k])<=tau,],na.rm=T)),na.rm=T)
        beta_old=beta_new
        miu_old=miu_new
        sigma2_old=sigma2_new
        par_new=c(beta_new,miu_new,sigma2_new,tau)
        dis_par=sqrt(sum((par_new-par_old)^2,na.rm=T))
        ll1=ll1+1
        if(dis_par<1*10^(-6) |ll1>100){
          break
        }
      }
      
      ll2=1
      tau2=tau
      repeat{
        par_old=c(miu_old,tau2)
        delta=Y-X%*%beta_new
        #residuals
        e=sapply(1:K, function(k) delta-miu_old[k])
        rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        #update tau2
        tau2=sqrt(sum(ifelse(rr^2>tau2^2,tau2^2,rr^2))/(log(n)))
        #update miu
        miu_new=sapply(1:K, function(k) sum(c((w_old[,k]*delta)[abs(e[,k])<=tau2],(w_old[,k]*tau2*sign(e[,k]))[abs(e[,k])>tau2]),na.rm=T)/sum(w_old[abs(e[,k])<=tau2,k],na.rm=T))
        #update beta
        #beta_new=beta_old
        #for(i in 1:np){
        #  beta_new[i]=sum(sapply(1:K, function(k) c((w_old[,k]*X[,i]*(Y-miu_new[k]-X[,-i,drop=F]%*%beta_new[-i]))[abs(e[,k])<=tau],(w_old[,k]*tau*sign(e[,k]))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum((w_old[,k]*X[,i]^2)[abs(e[,k])<=tau],na.rm=T)),na.rm=T)
        #}
        #update sigma2
        #sigma2_new=sum(sapply(1:K, function(k) c((w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2)[abs(e[,k])<=tau],(w_old[,k]*(2*tau*abs(Y-miu_new[k]-X%*%beta_new)-tau^2))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[,k])<=tau,],na.rm=T)),na.rm=T)
        #e=sapply(1:K, function(k) Y-miu_new[k]-X%*%beta_new)
        #tau=sqrt(sum(sapply(1:K, function(k) w_old[,k]*ifelse(e[,k]^2>tau^2,tau^2,e[,k]^2)))/(np+log(n)))
        #rr=apply(sapply(1:K, function(k) pai_old[,k]*e[,k]),1,sum)
        #tau=sqrt(sum(ifelse(rr^2>tau^2,tau^2,rr^2))/(np+log(n)))
        #beta_old=beta_new
        miu_old=miu_new
        #sigma2_old=sigma2_new
        par_new=c(miu_new,tau2)
        dis_par=sqrt(sum((par_new-par_old)^2,na.rm=T))
        ll2=ll2+1
        if(dis_par<1*10^(-6) |ll2>100){
          break
        }
      }
      print(paste0("tau1=",tau,"; tau2=",tau2))
    }
    
    
    if(robust=="huberk"){
      e=sapply(1:K, function(k) Y-miu_old[k]-X%*%beta_old)
      #update miu, alpha
      miu_new=sapply(1:K, function(k) sum(c((w_old[,k]*(Y-X%*%beta_old))[abs(e[,k])<=tau],(w_old[,k]*tau*sign(e[,k]))[abs(e[,k])>tau]),na.rm=T)/sum(w_old[abs(e[,k])<=tau,k],na.rm=T))
      #update beta
      beta_new=beta_old
      for(i in 1:np){
        beta_new[i]=sum(sapply(1:K, function(k) c((w_old[,k]*X[,i]*(Y-miu_new[k]-X[,-i,drop=F]%*%beta_new[-i]))[abs(e[,k])<=tau],(w_old[,k]*tau*sign(e[,k]))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum((w_old[,k]*X[,i]^2)[abs(e[,k])<=tau],na.rm=T)),na.rm=T)
      }
      #update sigma2
      sigma2_new=sum(sapply(1:K, function(k) c((w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2)[abs(e[,k])<=tau],(w_old[,k]*(2*tau*abs(Y-miu_new[k]-X%*%beta_new)-tau^2))[abs(e[,k])>tau])),na.rm=T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[,k])<=tau,],na.rm=T)),na.rm=T)
      #sigma2_new=sum(sapply(1:K, function(k) w_old[,k]*(Y-miu_new[k]-X%*%beta_new)^2),na.rm=T)/n      
    }
    theta_new=c(beta_new, 
                as.numeric(t(gamma_new_matrix[,-K])), 
                miu_new, sigma2_new)
    dis=sqrt(sum((theta_new-theta)^2,na.rm=T))
    theta=theta_new
    #print(theta)
    l=l+1
    if(dis<1*10^(-6) |l>500){
      break
    }
  }
  return(theta)
}

