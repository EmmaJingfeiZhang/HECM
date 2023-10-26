library(rTensor)
library(MASS)
library(glasso)
#library(mvtnorm)
library(clusterGeneration)
k_list <- function(x,P){
  if(P==1) Sigma=x[[1]]
  if(P>1) Sigma=kronecker_list(x)
  return(Sigma)
}
outer_list <- function(x,P){
  if(P==1){
    x0=as.tensor(x[[1]])
    return(x0)
  } 
  if(P>1){
    x0=x[[1]]
    for(p in 2:P) x0=x0%o%x[[p]]
    x0=as.tensor(x0)
    return(x0)
  }
}

# CP decomposition of mu
fac <- function(mu,Omega,R,nk=1,n=1,lambda1=0,nt=20){
  mp=dim(mu)
  Omega1=Omega
  P=length(mp)
  Omega=lapply(1:P,function(p)diag(mp[p]))
  weight.new=weight=rep(1,R)
  beta.new=beta=lapply(1:R,function(r){
    lapply(1:P,function(p) {
      x0=sapply(1:mp[p],function(p0)runif(1))
      x0/sqrt(sum(x0^2))
    } )
  })
  mu.r=lapply(1:R,function(r){
    outer_list(beta[[r]],P)
  })
  
  for(r in 1:R){
    delta=10
    t=1
    mu.rs=-weight.new[r]*mu.r[[r]]
    for(r0 in 1:r) mu.rs=mu.rs+weight.new[r0]*mu.r[[r0]]
    mu.new=mu-mu.rs
    while(delta>10^(-4) & t<nt){
      beta[[r]]=beta.new[[r]]
      for(p in 1:P){
        mu.new.p0=k_unfold(mu.new,p)
        dims=dim(mu.new.p0)
        mu.new.p=array(vec(mu.new.p0),dim=dims)
        Omega.dp=Omega[-p]
        Omega.dpi=rev(Omega.dp)
        Omega.dpi.kp=k_list(Omega.dpi,P-1)
        betar.dp=beta.new[[r]][-p]
        betar.dp.out=outer_list(betar.dp,P-1)
        betar.dp.out.v=as.matrix(vec(betar.dp.out))
        
        C_krp=(t(betar.dp.out.v)%*%Omega.dpi.kp%*%betar.dp.out.v)[1,1]*weight.new[r]^2
        t1=Omega[[p]]%*%((mu.new.p%*%Omega.dpi.kp%*%betar.dp.out.v)*weight.new[r]/C_krp)
        t2=Omega[[p]]%*%beta.new[[r]][[p]]
        t12=(t1-t2)/diag(Omega[[p]])
        g1=t12+beta.new[[r]][[p]]
        betarp=g1
        beta.new[[r]][[p]]=as.vector(betarp/sqrt(sum(betarp^2)))
      }
      delta=sum(sapply(1:P,function(p0)sum((beta.new[[r]][[p0]]-beta[[r]][[p0]])^2)))
      t=t+1
    }
    mu.r[[r]]=outer_list(beta.new[[r]],P)
    Omegai=rev(Omega)
    Omegai.kp=k_list(Omegai,P)
    mur.v=as.matrix(vec(mu.r[[r]]))
    mu.new.v=as.matrix(vec(mu.new))
    upper=t(mu.new.v)%*%Omegai.kp%*%mur.v
    down=t(mur.v)%*%Omegai.kp%*%mur.v
    weight.new[r]=upper/down
  }
  mu.hat=weight.new[1]*mu.r[[1]]
  if(R>1){
    for(r0 in 2:R) mu.hat=mu.hat+weight.new[r0]*mu.r[[r0]]
  }
  Omega=Omega1
  for(r in 1:R){
    delta=10
    t=1
    mu.rs=-weight.new[r]*mu.r[[r]]
    for(r0 in 1:r) mu.rs=mu.rs+weight.new[r0]*mu.r[[r0]]
    mu.new=mu-mu.rs
    while(delta>10^(-4) & t<nt){
      beta[[r]]=beta.new[[r]]
      for(p in 1:P){
        mu.new.p0=k_unfold(mu.new,p)
        dims=dim(mu.new.p0)
        mu.new.p=array(vec(mu.new.p0),dim=dims)
        Omega.dp=Omega[-p]
        Omega.dpi=rev(Omega.dp)
        Omega.dpi.kp=k_list(Omega.dpi,P-1)
        betar.dp=beta.new[[r]][-p]
        betar.dp.out=outer_list(betar.dp,P-1)
        betar.dp.out.v=as.matrix(vec(betar.dp.out))
        
        C_krp=(t(betar.dp.out.v)%*%Omega.dpi.kp%*%betar.dp.out.v)[1,1]*weight.new[r]^2
        t1=Omega[[p]]%*%((mu.new.p%*%Omega.dpi.kp%*%betar.dp.out.v)*weight.new[r]/C_krp)
        t2=Omega[[p]]%*%beta.new[[r]][[p]]
        t12=(t1-t2)/diag(Omega[[p]])
        g1=t12+beta.new[[r]][[p]]
        g2=C_krp*nk*g1*diag(Omega[[p]]) #nk*g1*diag(Omega[[p]])*den*weight.new[r]
        beta_zero=which(abs(g2)<=lambda1)
        betarp=g1-n*lambda1/(nk*C_krp)/diag(Omega[[p]])*sign(beta.new[[r]][[p]])
        betarp[beta_zero]=0
        if(sum(abs(betarp))==0){
          return(list(mu=mu.hat,weight=weight.new,beta=beta.new,t=t))
          stop("penalty is too big")
        }  
        beta.new[[r]][[p]]=as.vector(betarp/sqrt(sum(betarp^2)))
      }
      delta=sum(sapply(1:P,function(p0)sum((beta.new[[r]][[p0]]-beta[[r]][[p0]])^2)))
      t=t+1
    }
    mu.r[[r]]=outer_list(beta.new[[r]],P)
    Omegai=rev(Omega)
    Omegai.kp=k_list(Omegai,P)
    mur.v=as.matrix(vec(mu.r[[r]]))
    mu.new.v=as.matrix(vec(mu.new))
    upper=t(mu.new.v)%*%Omegai.kp%*%mur.v
    down=t(mur.v)%*%Omegai.kp%*%mur.v
    weight.new[r]=upper/down
  }
  mu.hat=weight.new[1]*mu.r[[1]]
  if(R>1){
    for(r0 in 2:R) mu.hat=mu.hat+weight.new[r0]*mu.r[[r0]]
  }
  list(mu=mu.hat,weight=weight.new,beta=beta.new,t=t)
}

# update omega
stgm <- function(x,Omega,lambda,Lk,tmax){ #sparse tensor graphical model 
  mp=dim(x[[1]])
  P=length(mp)
  m=prod(mp)
  n=length(x)
  lambda.p=lapply(1:P,function(p0){
    lambda0=matrix(lambda[p0],mp[p0],mp[p0])
    diag(lambda0)=0
    lambda0
  })
  Omega.new=Omega
  sigma=Omega
  t=1
  delta=10
  while(delta>10^(-3) & t<tmax){
    Omega=Omega.new
    for(p in 1:P){
      tensorp_matrix <- function(x){
        xp=k_unfold(x,p)
        dim0=dim(xp)
        array(vec(xp),dim=dim0)
      }
      x.p=lapply(x,tensorp_matrix)
      Omega.p=Omega.new[-p]
      Omega.pv=rev(Omega.p)
      Omega.pv.k=k_list(Omega.pv,P-1)
      va.m <- function(x) x%*%Omega.pv.k%*%t(x)
      VVn0=lapply(x.p,va.m)
      VVn=Map("*",Lk,VVn0)
      
      Sp0=Reduce("+",VVn)
      Sp=Sp0*mp[p]/m/sum(Lk)
      Omegap.up=glasso(Sp,rho=lambda.p[[p]])#$wi
      wi0=(Omegap.up$wi+t(Omegap.up$wi))/2
      Omega.new[[p]]=wi0/sqrt(sum(wi0^2))#*sqrt(mp[p])
      sigma[[p]]=Omegap.up$w*sqrt(sum(wi0^2))#/sqrt(mp[p])#/sqrt(sum(Omegap.up$w^2))
    }
    delta=sum(sapply(1:P,function(p0) sum((Omega.new[[p0]]-Omega[[p0]])^2) ))
    t=t+1
  }
  x0=t(sapply(x,vec))
  #sigma1=t(x0)%*%diag(Lk)%*%x0/sum(Lk)
  sigma.hat0=k_list(rev(sigma),P)
  #weight=sum(sigma.hat0*sigma1)/sum(sigma.hat0^2)
  omega.hat0=k_list(rev(Omega.new),P)
  weight01=diag(x0%*%omega.hat0%*%t(x0))
  weight0=sum(weight01*Lk)
  weight.omega=sum(Lk)*m/weight0
  weight=1/weight.omega
  sigma.kp=weight*sigma.hat0
  list(Omega=Omega.new,sigma=sigma,sigma.kp=sigma.kp,weight=weight,t=t)
}

# EM algorithm for mixture gaussian model
mve <- function(x,K,lambda,R,lambda1,delta0,tpmax,tmax){
  mp=dim(x[[1]])
  P=length(mp)
  m=prod(mp)
  n=length(x)
  pri=rep(1/K,K)
  x0=t(sapply(x,vec))
  t0=proc.time()
  #initialization
  sigma0=t(x0)%*%x0/n
  scale=sqrt(sum(sigma0^2))
  Omega.new=Omega=lapply(1:K,function(k){
    lapply(1:P,function(p) diag(1/sqrt(scale/mp[p]),mp[p]) )
  })
  t1=proc.time()#24.77
  cl=kmeans(x0,K,iter.max=100,nstart=100)
  centers=cl$centers
  centers=lapply(1:K,function(k)as.tensor(array(centers[k,],dim=mp)))
  mu.new=vector("list",K)
  beta=vector("list",K)
  weight=matrix(0,K,R)
  for(k in 1:K){
    dec = fac(centers[[k]],Omega[[1]],R,lambda1=0)
    beta[[k]]=dec$beta
    weight[k,]=dec$weight
    mu.new[[k]]=dec$mu
  }
  weight.new=weight
  beta.new=beta
  t2=proc.time()#106.98
  weight.sig=rep(1,K)
  
  sigma.new=lapply(1:K,function(k){
    lapply(1:P,function(p) solve(Omega[[k]][[p]]))
  })
  sigma.new.kp=lapply(1:K,function(k)k_list(rev(sigma.new[[k]]),P))
  omega.new.kp=lapply(1:K,function(k)k_list(rev(Omega.new[[k]]),P))
  sigma.new.det.com=t(sapply(1:K,function(k){
    sapply(1:P,function(p)(det(sigma.new[[k]][[p]])/det(sigma.new[[1]][[p]]))^(-m/(2*mp[p])))
  }))
  t3=proc.time() #9.14
  t=1
  delta=10
  while(delta>delta0 & t<tmax){
    #update L and pi(pri)
   
    term=sapply(1:K,function(k){
      mu_k=matrix(rep(vec(mu.new[[k]]),n),nrow=n,ncol=m,byrow=TRUE)
      xk0= x0-mu_k
      term0=diag(xk0%*%omega.new.kp[[k]]%*%t(xk0))
    })
    L0=t(sapply(1:n,function(i){
      term=term[i,]
      k.min=which.min(term)
      sapply(1:K,function(k1){
        log_term=log(prod(sigma.new.det.com[k1,]))-log(prod(sigma.new.det.com[k.min,]))-m/2*log(weight.sig[k1]/weight.sig[k.min])+log(pri[k1]/pri[k.min])-0.5*(term[k1]-term[k.min])
        #prod(sigma.new.det.com[k1,])/prod(sigma.new.det.com[k.min,])*exp(-0.5*(term[k1]-term[k.min]))*pri[k1]/pri[k.min]*(weight.sig[k1]/weight.sig[k.min])^(-m/2)
        exp(log_term)
      })
    }))
   
    L0.rs=sapply(1:K,function(k) rowSums(L0)) #row sums of L0
    L=L0/L0.rs
    
    nk=colSums(L)
    if(sum(nk>10,na.rm=TRUE)!=K){
      print("the sample for one cluster is too small")
      return(list(L=L,weight=weight.new,beta=beta.new,mu=mu.new,Omega=Omega.new,sigma=sigma.new,sigma.kp=sigma.new.kp,omega.kp=omega.new.kp,weight.kp=weight.sig,t=t))
      stop("the sample for one cluster is too small")
    } 
    pri=nk/n
    #update mu and omega for each k
    beta=beta.new
    weight=weight.new
    Omega=Omega.new
  
    for(k in 1:K){
      x.k0=t(L[,k])%*%x0
      x.k=as.tensor(array(x.k0,dim=mp))/nk[k]

      muk.up=fac(x.k,Omega=Omega.new[[k]],R,nk=nk[k],n=n,lambda1,tpmax)
      
      beta.new[[k]]=muk.up$beta
      weight.new[k,]=muk.up$weight
      mu.new[[k]]=muk.up$mu
      xk.dm=lapply(1:n,function(i)x[[i]]-mu.new[[k]])

      Omegak.up=stgm(xk.dm,Omega=Omega.new[[k]],lambda,Lk=L[,k],tmax)

      Omega.new[[k]]=Omegak.up$Omega
      sigma.new.kp[[k]]=Omegak.up$sigma.kp
      sigma.new[[k]]=Omegak.up$sigma
      sigma.new.det.com[k,]=sapply(1:P,function(p)(det(sigma.new[[k]][[p]])/det(sigma.new[[1]][[p]]))^(-m/(2*mp[p])))
      weight.sig[k]=Omegak.up$weight
      omega.new.kp[[k]]=weight.sig[k]^(-1)*k_list(rev(lapply(sigma.new[[k]],solve)),P)#solve(sigma.new.kp[[k]])
    }
    
    delta1=sum(sapply(1:K,function(k){
      sum(sapply(1:R,function(r) sapply(1:P,function(p) sqrt(sum((beta.new[[k]][[r]][[p]]-beta[[k]][[r]][[p]])^2))  ) ))
    }))
    delta2=sum(sapply(1:K,function(k){
      sum(sapply(1:P,function(p0) sqrt(sum((Omega.new[[k]][[p0]]-Omega[[k]][[p0]])^2)) ))
    }))
    delta=delta1+delta2
    t=t+1
  }
  return(list(L=L,weight=weight.new,beta=beta.new,mu=mu.new,Omega=Omega.new,sigma=sigma.new,sigma.kp=sigma.new.kp,omega.kp=omega.new.kp,weight.kp=weight.sig,t=t))
}
# BIC
myBIC = function(L,mu,sigma,Omega,sigma.kp,omega.kp,beta,weight.kp,x,method="BIC"){
  ## L: the matrix of the probability in each cluser
  ## mu,Omega,sigma,beta: the output from mve function
  ## x: the original tensor
  ## output: BIC value
  ## method: bic(default) or ebic
  n=length(x)
  K=length(beta)
  R=length(beta[[1]])
  P=length(beta[[1]][[1]])
  mp=dim(mu[[1]])
  m=prod(mp)
  pri=colSums(L)/n
  s1=sapply(1:K,function(k){
    nonzero=sapply(1:R,function(r)sapply(1:P,function(p) sum(beta[[k]][[r]][[p]]!=0) ))
    sum(nonzero)
  })
  s2=sapply(1:K,function(k){
    nonzero=sapply(1:P,function(p) {
      Omega.minus=Omega[[k]][[p]]
      diag(Omega.minus)=0
      sum(Omega.minus!=0)/2
    } )
    sum(nonzero)
  })
  sigma_eigen=lapply(1:K,function(k){
    all_values=lapply(1:P,function(p)eigen(sigma[[k]][[p]])$values)
    eigen0=outer_list(all_values,P)
    eigen0=vec(eigen0)
    eigen0[order(eigen0,decreasing=TRUE)]*weight.kp[k]
  }
  )
  
  term=t(sapply(1:K,function(k1){
    sigma_eigen0=sigma_eigen[[k1]]
    fk <- function(x1){
      log(pri[k1])-0.5*sum(log(sigma_eigen0*2*pi))-t(vec(x1-mu[[k1]]))%*%omega.kp[[k1]]%*%vec(x1-mu[[k1]])/2
    }  
    sapply(x,fk)
  }))
  
  xpdf_id0=sapply(1:n,function(i){
    term0=term[,i]
    term1=term0[term0>=max(term)-1000]
    scale=mean(term1)
    term2=exp(term1-scale)
    log(sum(term2))+scale
  })
  #scale=mean(term)
  #term0=exp(term-scale)
  #xpdf_id0=log(colSums(term0))+scale
  
  likelihood=-2*sum(xpdf_id0)
  
  
  parameter=sum(sapply(1:K,function(k) (s1[k]+s2[k]) ))
  if(method=="BIC"){
    n0=K*R*sum(mp)
    n1=K*sum(sapply(1:P,function(p0) mp[p0]*(mp[p0]-1) ))
    bic=likelihood+log(n)*sum(s1)+sum(s2)*2*log(n)+0.5*(sum(s1)+2*sum(s2))*log(n0+n1)
  }else if(method=="eBIC"){
    n0=K*R*sum(mp)
    n1=K*sum(sapply(1:P,function(p0) mp[p0]*(mp[p0]-1)/2 ))
    bic=likelihood+log(n)*parameter
  }
  bic
}

# tunning
mytune <- function(x,K,rank_list=seq(1,5,1),lambda1_list=seq(0,.01,.002),lambda_list=seq(0,.1,0.01),tune_method = "together",bic_method="BIC"){  
  # x: tensor data
  # rank_list: list of possible ranks for tuning. default rank_list = seq(1,5,1)
  # sparse_list: list of possible sparsity parametters. default sparse_list = seq(0.1, 0.9, 0.1).
  # lambda_list: list of possible lambda parametters. default lambda_list = seq(0.1, 0.9, 0.1).
  # tune_method: {"sequential", "together"}. "sequential" tunes Rank --> lambda --> Sparse. "together" tunes all three parameters together. default = "sequential"
  #bic_method: bic (default) or ebic
  mp = dim(x[[1]])
  P = length(mp) 
  n = length(x)
  len_rank = length(rank_list)
  len_tune_lambda1 = length(lambda1_list)
  len_tune_lambda = length(lambda_list)
  TIME = array(0, c(len_rank, len_tune_lambda1, len_tune_lambda))
  if(tune_method == "together"){
    ## tune_method = "together"
    bic_value_mat = rep(NA, 4)
    for(irank in 1:len_rank){
      for(ilambda1 in 1:len_tune_lambda1){
        for(ilambda in 1:len_tune_lambda){
          tt2 = proc.time()
          Rank = rank_list[irank]
          lambda1_para = lambda1_list[ilambda1]
          lambda_para = lambda_list[ilambda]
          
          out_est = mve(x,K,lambda=lambda_para,R=Rank,lambda1=lambda1_para,delta0=10^(-3),tpmax=20,tmax=20)
          BIC_OUT = myBIC(out_est$L,out_est$mu,out_est$Omega,sigma=out_est$sigma.kp,omega.kp=out_est$omega.kp,out_est$beta,x,method=bic_method)
          bic_value_mat = rbind(bic_value_mat, c(BIC_OUT, Rank, lambda1_list[ilambda1], lambda_list[ilambda]))
          TIME[irank, ilambda1, ilambda] = ((proc.time() - tt2)[3])
        }
      }
    }
    time.rank = sum(apply(TIME, 1, mean))
    time.lambda1 = sum(apply(TIME, c(1,2), mean)) - time.rank
    time.lambda = sum(apply(TIME, c(1,3), mean)) - time.rank
    bic_value_mat = bic_value_mat[-1,]
    colnames(bic_value_mat) = c("bic", "Rank", "Lambda1", "Lambda")
    min_BIC_index = which.min(bic_value_mat[,1])
    Rank.opt = bic_value_mat[min_BIC_index, 2]
    lambda1.opt = bic_value_mat[min_BIC_index,3]
    lambda.opt = bic_value_mat[min_BIC_index,4]
    
  }else{
    
    ## tune_method = "sequential"
    tt1 = proc.time()
    bic_value_mat_rank = rep(NA, 2)
    for(irank in 1:len_rank){
      Rank = rank_list[irank]
      
      out_est = mve(x,K,lambda=rep(0,P),R=Rank,lambda1=0,delta0=10^(-3),tpmax=20,tmax=20)
      BIC_OUT = myBIC(out_est$L,out_est$mu,sigma=out_est$sigma,Omega=out_est$Omega,sigma.kp=out_est$sigma.kp,omega.kp=out_est$omega.kp,beta=out_est$beta,weight.kp=out_est$weight.kp,x,method=bic_method)
      bic_value_mat_rank = rbind(bic_value_mat_rank, c(BIC_OUT, Rank))
    }
    bic_value_mat_rank = bic_value_mat_rank[-1,]
    min_BIC_index = which.min(bic_value_mat_rank[,1])
    Rank.opt = bic_value_mat_rank[min_BIC_index, 2]
    
    time.rank = (proc.time() - tt1)[3]
    
    bic_value_mat_lambda_all = rep(NA, 2)
    lambda=rep(0,P)
    for(p in 1:P){
      tt1 = proc.time()		
      bic_value_mat_lambda = rep(NA, 2)
      for(ilambda in 1:len_tune_lambda){
        
        lambda[p] = lambda_list[ilambda]
        out_est = mve(x,K,lambda=lambda,R=Rank.opt,lambda1=0,delta0=10^(-3),tpmax=20,tmax=20)
        BIC_OUT = myBIC(out_est$L,out_est$mu,sigma=out_est$sigma,Omega=out_est$Omega,sigma.kp=out_est$sigma.kp,omega.kp=out_est$omega.kp,beta=out_est$beta,weight.kp=out_est$weight.kp,x,method=bic_method)
        bic_value_mat_lambda = rbind(bic_value_mat_lambda, c(BIC_OUT, lambda_list[ilambda]))
        rm(out_est)
      }
      bic_value_mat_lambda = bic_value_mat_lambda[-1,]
      bic_value_mat_lambda_all = rbind(bic_value_mat_lambda_all,bic_value_mat_lambda)
      min_BIC_index = which.min(bic_value_mat_lambda[,1])
      lambda.opt = lambda_list[min_BIC_index]
      lambda[p]=lambda.opt
      time.lambda = (proc.time() - tt1)[3]
    }
    bic_value_mat_lambda_all = bic_value_mat_lambda_all[-1,]
    
    
    
    tt1 = proc.time()
    bic_value_mat_lambda1 = rep(NA, 2)
    for(ilambda1 in 1:len_tune_lambda1){
      lambda1_para = lambda1_list[ilambda1]
      
      out_est = mve(x,K,lambda=lambda,R=Rank.opt,lambda1=lambda1_para,delta0=10^(-3),tpmax=20,tmax=20)
      BIC_OUT = myBIC(out_est$L,out_est$mu,sigma=out_est$sigma,Omega=out_est$Omega,sigma.kp=out_est$sigma.kp,omega.kp=out_est$omega.kp,beta=out_est$beta,weight.kp=out_est$weight.kp,x,method=bic_method)
      bic_value_mat_lambda1 = rbind(bic_value_mat_lambda1, c(BIC_OUT, lambda1_list[ilambda1]))
      rm(out_est)
    }
    bic_value_mat_lambda1 = bic_value_mat_lambda1[-1,]
    min_BIC_index = which.min(bic_value_mat_lambda1[,1])
    lambda1.opt = lambda1_list[min_BIC_index]
    
    time.lambda1 = (proc.time() - tt1)[3]
    
    bic_value_mat = rbind(bic_value_mat_rank, bic_value_mat_lambda_all, bic_value_mat_lambda1)
  }
  
  P <- list()
  P$Rank.opt = Rank.opt
  P$lambda1.opt = lambda1.opt
  P$lambda.opt = lambda
  P$time.rank = time.rank
  P$time.lambda = time.lambda
  P$time.lambda1 = time.lambda1
  P$bic_value_mat = bic_value_mat
  
  return(P)
}