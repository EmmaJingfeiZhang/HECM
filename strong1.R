library("iterators")
library("parallel")
library("foreach")
library("doParallel")
rm(fnorm)
source("/nethome/jxz280/tensor/all31.R")
source("/nethome/jxz280/tensor/data.R")
source("/nethome/jxz280/tensor/measure.R")

mytune <- function(x,K,rank_list=seq(1,10,1),lambda1_list=seq(0,.01,.002),lambda_list=seq(0,.1,0.01),tune_method = "together",bic_method="BIC"){  
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
    #out_est_all=vector("list",33)
    lambda.opt=rep(0,P)
    lambda=rep(0,P)
    for(p in 1:P){
      tt1 = proc.time()		
      bic_value_mat_lambda = rep(NA, 2)
      for(ilambda in 1:len_tune_lambda){
        
        lambda[p] = lambda_list[ilambda]
        out_est = mve(x,K,lambda=lambda,R=Rank.opt,lambda1=0,delta0=10^(-3),tpmax=20,tmax=20)
        #out_est_all[[(p-1)*11+ilambda]]=out_est
        BIC_OUT = myBIC(out_est$L,out_est$mu,sigma=out_est$sigma,Omega=out_est$Omega,sigma.kp=out_est$sigma.kp,omega.kp=out_est$omega.kp,beta=out_est$beta,weight.kp=out_est$weight.kp,x,method=bic_method)
        bic_value_mat_lambda = rbind(bic_value_mat_lambda, c(BIC_OUT, lambda_list[ilambda]))
        rm(out_est)
      }
      bic_value_mat_lambda = bic_value_mat_lambda[-1,]
      bic_value_mat_lambda_all = rbind(bic_value_mat_lambda_all,bic_value_mat_lambda)
      min_BIC_index = which.min(bic_value_mat_lambda[,1])
      lambda.opt[p] = lambda_list[min_BIC_index]
      lambda[p]=0
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
  P$lambda.opt = lambda.opt
  P$time.rank = time.rank
  P$time.lambda = time.lambda
  P$time.lambda1 = time.lambda1
  P$bic_value_mat = bic_value_mat
  
  return(P)
}

myexchangeable = function (p, sigma = 1, rho = 0.1){
  Sigma <- matrix(1, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(i != j ){
        Sigma[i, j] <- sigma^2 * rho      		
      }
      
    }
  }
  return(Sigma)
}

tensor_mean5 <- function(u,p){
  beta1=c(u,u,u,rep(0,p-3))
  beta2=c(rep(0,2),u,u,u,rep(0,p-5))
  beta3=c(rep(0,4),u,u,u,rep(0,p-7))
  beta4=c(rep(0,6),u,u,u,rep(0,p-9))
  #beta5=c(rep(0,4),u,u,u,rep(0,p-7))
  mu10=outer_list(list(beta1,beta1,beta1),3)
  mu20=outer_list(list(beta2,beta2,beta2),3)
  mu30=outer_list(list(beta3,beta3,beta3),3)
  mu40=outer_list(list(beta4,beta4,beta4),3)
  #mu50=outer_list(list(beta5,beta5,beta5),3)
  mu1=mu10+mu20+mu30+mu40#+mu50
  mu2=-1*mu1
  mu3=mu10-mu20+mu30-mu40#+mu50
  mu4=-1*mu3
  list(mu1,mu2,mu3,mu4)
}

K=4
P=3
mp=c(10,10,10)
#totally identity
set.seed(1)
sigma1=myexchangeable(5,1,0.6)
Sigma10=matrix(0,10,10)
Sigma10[1:5,1:5]=Sigma10[6:10,6:10]=sigma1
Sigma0=list(Sigma10,Sigma10,Sigma10)
Sigma=lapply(1:4,function(k)Sigma0)
Sigma.true=lapply(1:K,function(k0) k_list(Sigma[[k0]],3))
omega.true=Sigma

mu.true=mu=tensor_mean5(u=0.80,p=10)


registerDoParallel(cores=25)
sta4 <- foreach(m=1:50) %dopar%{
  cat("simulation",m,"\n")
  t1 = proc.time()
  set.seed(m)
  #data generation and tunning
  DATA1_400 = data_rgen(mu,Sigma,n=400)
  tune1_400 = mytune(x=DATA1_400,K=4,lambda_list=seq(0,0.001,0.0001),tune_method = "sequential")  
  out_est1_400 = mve(DATA1_400,K,lambda=tune1_400$lambda.opt,R=tune1_400$Rank.opt,lambda1=tune1_400$lambda1.opt,delta0=10^(-3),tpmax=20,tmax=20)
  time_ours=(proc.time() - t1)[3]
  L1=out_est1_400$L
  l1=sapply(1:nrow(L1),function(i) which.max(L1[i,]))
  n0=400/4
  l2=rep(1:4,rep(n0,4))
  t2=proc.time()
  x1=t(sapply(DATA1_400,vec))
  cl1=kmeans(x1,4,nstart=100)
  time_kmeans=(proc.time()-t2)[3]
  mu.est1=lapply(1:K,function(k){
    mu.est10=array(cl1$centers[k,],dim=mp)
    as.tensor(mu.est10)
  })
  err1=mycluster_error(l1,l2)
  err2=mycluster_error(cl1$cluster,l2)
  mu.est=out_est1_400$mu
  Sigma.est=out_est1_400$sigma.kp
  omega.est=out_est1_400$Omega
  para_error=myerror(mu.est=mu.est,mu.true=mu.true,Sigma.est=Sigma.est,Sigma.true)
  para_error1=myerror(mu.est=mu.est1,mu.true=mu.true,Sigma.est=Sigma.est,Sigma.true)
  rec=myrecovery(mu.est=mu.est,mu.true,Sigma.est=Sigma.est,Sigma.true,omega.est=omega.est,omega.true)
  c(err1,mu_error=para_error$mu_error,Sigma_error=para_error$Sigma_error,TPR=rec$TPR,FPR=rec$FPR,err2,time_ours,mu_error2=para_error1$mu_error,time_kmeans)
}
save(sta4,file="/nethome/jxz280/tensor/exp11/0.6/sta4.RData")