#measure
beta_distance <- function(beta1,beta2){
  # two levels list, the length of first level is rank, the length of second level is dimension
  rank1=length(beta1)
  rank2=length(beta2)
  mp1=sapply(beta1[[1]],length)
  mp2=sapply(beta2[[1]],length)
  if(rank1 !=rank2) stop("rank1 =",rank1,"is not equal to rank2 =",rank2)
  if(mp1!=mp2) stop("mp1 =",mp1,"is not equal to mp2 =",mp2)
  dis.r=lapply(1:rank1,function(r1){
    sapply(1:rank2,function(r2) 
      sapply(1:length(beta1[[1]]),function(p) 
        min(norm(beta1[[r1]][[p]] - beta2[[r2]][[p]],type="2"), norm(beta1[[r1]][[p]] + beta2[[r2]][[p]],type="2")) 
      ) 
    )
  })
  dis.r0=sapply(1:rank1,function(r1)min(colMeans(dis.r[[r1]])))
  sum(dis.r0)
}
myerror_betaw = function(beta.est, beta.true, w.est, w.true){
  #compute the final estimation error of components and weight w.
  K=length(beta.true)
  K.est=length(beta.est)
  rank=length(beta.est[[1]])
  P=length(beta.est[[1]][[1]])
  mp=sapply(beta.est[[1]][[1]],length)
  
  if(K.est > K){
    # more latent decomposition. truncate it to first K columns
    beta.est=beta.est[1:K]
    w.est=w.est[1:K,]
    warning("K.est = ", K.est, " is larger than true K = ", K, ". Use first K latent estimations!")
  }else if(K.est <K){
    # less latent decomposition. add columns with zeros to fill out K columns
    beta.plus=lapply(1:P,function(p)c(1,rep(0,mp[p]-1)))
    beta.est[(K.est+1):K]=beta.plus
    w.est[(K.est+1):K,1:Rank]=0
    warning("K.est = ", K.est, " is smaller than true K = ", K, ". Fill up zeros to get K latent estimations!")
  }
 
  ## estimation error 
  ERROR = matrix(NA,K,K)
  j.best = rep(0,K)
  
  for(i in 1:K){
    for(j in 1:K){	
      ERROR[i,j] = beta_distance(beta.est[[i]],beta.true[[j]])
    }
    j.best[i] = which.min(ERROR[i,])	
  }
  
  error.k = rep(0, K)
  error.w = rep(0,K)
  for(ii in 1:K){
    error.k[ii] = ERROR[ii,j.best[ii]]
    error.w[ii] = sqrt(sum((abs(w.est[ii]) - abs(w.true[j.best[ii]]))^2))/sqrt(sum(w.true[j.best[ii]]^2))
  }
  
  P = list()
  P$error.beta = mean(error.k)
  P$error.w = mean(error.w)
  P
}


myerror = function(mu.est,mu.true,Sigma.est,Sigma.true){
  #compute the tensor recovery error
  K=length(mu.true)
  K.est=length(mu.est)
  ERROR=sapply(1:K,function(k1){
    sapply(1:K.est,function(k2){
      fnorm(mu.est[[k1]] - mu.true[[k2]])/fnorm(mu.true[[k2]])
    })
  })
  error.k=sapply(1:K,function(k1) min(ERROR[,k1]))
  k.est_id=sapply(1:K,function(k1) which.min(ERROR[,k1]))
  mu_error=mean(error.k)
  Sigma_error0=sapply(1:K,function(k1) sqrt(sum((Sigma.est[[k1]]-Sigma.true[[k.est_id[k1]]])^2))/sqrt(sum(Sigma.true[[k.est_id[k1]]])^2)) 
  list(mu_error=mu_error,Sigma_error=mean(Sigma_error0),k.est_id=k.est_id)
}

tpr <-function(sigma1,sigma2){
  d_trans1=sum((sigma1-t(sigma1))^2)
  d_trans2=sum((sigma2-t(sigma2))^2)
  if(d_trans1>10^(-6)) stop("Sigma1 is not a symmetric matrix") 
  if(d_trans2>10^(-6)) stop("Sigma2 is not a symmetric matrix") 
  diag(sigma1)=diag(sigma2)=0
  upper=sum(sigma1*sigma2!=0)/2
  down=sum(sigma2!=0)/2
  tpr=upper/down
  tpr
}

tpr_P <- function(omega1,omega2){
  P=length(omega1)
  mp=sapply(1:P,function(p) nrow(omega1[[p]]) )
  omega1_offdiag=lapply(1:P,function(p){
    omega1p=omega1[[p]]
    diag(omega1p)=0
    omega1p
  })
  omega2_offdiag=lapply(1:P,function(p){
    omega2p=omega2[[p]]
    diag(omega2p)=0
    omega2p
  })
  upper_p=sapply(1:P,function(p) sum(omega1_offdiag[[p]]*omega2_offdiag[[p]]!=0)/2 )
  down_p=sapply(1:P,function(p) sum(omega2_offdiag[[p]]!=0)/2 )
  if(sum(down_p)!=0) return(sum(upper_p)/sum(down_p))
  if(sum(down_p)==0){
    upper_p=sapply(1:P,function(p) sum(omega1[[p]]*omega2[[p]]!=0) )
    down_p=sapply(1:P,function(p) sum(omega2[[p]]!=0) )
    return(sum(upper_p)/sum(down_p))
  }
}

fpr_P <- function(omega1,omega2){
  P=length(omega1)
  mp=sapply(1:P,function(p) nrow(omega1[[p]]) )
  
  omega1_offdiag=lapply(1:P,function(p){
    omega1p=omega1[[p]]
    diag(omega1p)=0
    omega1p
  })
  omega2_offdiag_change=lapply(1:P,function(p){ #0->1,not 0->0
    omega2p=omega2[[p]]
    id0=which(omega2p==0)
    omega20=matrix(0,mp[p],mp[p])
    omega20[id0]=1
    diag(omega20)=0
    omega20
  })
  upper_p=sapply(1:P,function(p) sum(omega1_offdiag[[p]]*omega2_offdiag_change[[p]]!=0)/2 )
  down_p=sapply(1:P,function(p) sum(omega2_offdiag_change[[p]]!=0)/2 )
  sum(upper_p)/sum(down_p)
}

myrecovery <- function(mu.est,mu.true,Sigma.est,Sigma.true,omega.est,omega.true){
  error.est=myerror(mu.est,mu.true,Sigma.est,Sigma.true)
  k.est_id=error.est$k.est_id
  omega.true0=omega.true[k.est_id]
  K=length(omega.true)
  P=length(omega.true[[1]])
  TPR.k=sapply(1:K,function(k){
    tpr_P(omega1=omega.est[[k]],omega2=omega.true0[[k]])
  })
  FPR.k=sapply(1:K,function(k){
    fpr_P(omega1=omega.est[[k]],omega2=omega.true0[[k]])
  })
  list(TPR=mean(TPR.k),FPR=mean(FPR.k))
}

mycluster_error <- function(l1,l2){
  n=length(l1)
  error=0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      I1=(l1[i]==l1[j])*1
      I2=(l2[i]==l2[j])*1
      error=error+(I1!=I2)*1
    }
  }
  error*2/(n*(n-1))
}
