#data generation
#precision matrix
#identity matrix
sigma_identity <- function(p) diag(p) #p=dimensions
#Chain network
sigma_chain <- function(a,p){
  s_initial=runif(p,0.5,1)
  s=cumsum(s_initial)
  sigma=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      sigma[i,j]=exp(-a*abs(s[i]-s[j]))
    }
  }
  sigma
}
#Power-law network
sigma_powerlaw0 <- function(edge_list,p){
  edge_list1=edge_list[,2:1]
  edge_list2=rbind(edge_list,edge_list1)
  omega=omega_t=sigma=matrix(0,p,p)
  for(n0 in 1:nrow(edge_list2)) omega[edge_list2[n0,1],edge_list2[n0,2]]=runif(1,0.1,0.4)*sample(c(-1,1),1)
  diag(omega)=1
  r_sum=sapply(1:p,function(i) 0.9*(rowSums(abs(omega))-1))
  omega_t0=omega/r_sum
  id0=which(r_sum[,1]==0)
  omega_t0[id0,]=0
  diag(omega_t0)=1
  for(i in 1:p){
    for(j in 1:p){
      omega_t[i,j]=(omega_t0[i,j]+omega_t0[j,i])/2
    }
  }
  for(n0 in 1:nrow(edge_list2)){
    i=edge_list2[n0,1]
    j=edge_list2[n0,2]
    sigma[i,j]=0.9*omega_t[i,j]^(-1)/sqrt(omega_t[i,i]^(-1)*omega_t[j,j]^(-1))
  } 
  diag(sigma)=1
  sigma
}
sigma_powerlaw <- function(edge_list,p){
  edge_list1=edge_list[,2:1]
  edge_list2=rbind(edge_list,edge_list1)
  omega=omega_t=sigma=matrix(0,p,p)
  for(n0 in 1:nrow(edge_list2)) omega[edge_list2[n0,1],edge_list2[n0,2]]=runif(1,0.1,0.4)*sample(c(-1,1),1)
  diag(omega)=1
  r_sum=sapply(1:p,function(i) 0.9*(rowSums(abs(omega))-1))
  omega_t0=omega/r_sum
  id0=which(r_sum[,1]==0)
  omega_t0[id0,]=0
  diag(omega_t0)=1
  for(i in 1:p){
    for(j in 1:p){
      omega_t[i,j]=(omega_t0[i,j]+omega_t0[j,i])/2
    }
  }
  sigma=solve(omega_t)
  sigma
}
edge_list_rgen <- function(n,p,dup=FALSE){
  max_edge=p*(p-1)/2
  edge_list0=matrix(0,max_edge,2)
  edge_list0[,1]=rep(1:(p-1),(p-1):1)
  edge2=lapply(1:(p-1),function(i) (i+1):p )
  edge_list0[,2]=unlist(edge2)
  if(n<0) stop("n should be positive")
  if(n >= max_edge){
    edge_list=edge_list0
  }else{
    id0=sample(1:max_edge,n,replace=FALSE)
    id0_increasing=order(id0,decreasing=FALSE)
    id1=id0[id0_increasing]
    edge_list=edge_list0[id1,]
  }
  return(edge_list)
}
#tensor mean
tensor_mean <- function(u,p){
  beta1=c(u,-u,0.5*u,-0.5*u,rep(0,p-4))
  beta2=c(rep(0,4),u,-u,0.5*u,-0.5*u,rep(0,p-8))
  mu10=outer_list(list(beta1,beta1),2)
  mu20=outer_list(list(beta2,beta2),2)
  mu1=mu10-mu20
  mu2=mu10+mu20
  mu3=-1*mu1
  mu4=-1*mu2
  list(mu1,mu2,mu3,mu4)
}

#data
data_rgen <- function(mu,Sigma,n,out.method="tensor"){
  K=length(mu)
  mp=dim(mu[[1]])
  P=length(mp)
  mp1=sapply(1:length(Sigma[[1]]),function(i) nrow(Sigma[[1]][[i]]) )
  if(sum(mp1!=mp)!=0) stop("the dimensions of mu and Sigma differ")
  sigma.eigen=sapply(1:length(Sigma),function(i)sapply(1:P,function(j)min(eigen(Sigma[[i]][[j]])$values)))
  if(sum(sigma.eigen<=0)>0) stop("Sigma should be all positive")
  nk0=as.integer(n/K)
  nk=rep(nk0,K)
  nk[K]=n-(K-1)*nk0
  m=prod(mp)
  x=mvrnorm(nk[1],mu=vec(mu[[1]]),Sigma=k_list(rev(Sigma[[1]]),P))
  for(k in 2:K) x=rbind(x,mvrnorm(nk[k],mu=vec(mu[[k]]),Sigma=k_list(rev(Sigma[[k]]),P)))
  x0=lapply(1:n,function(i) as.tensor(array(x[i,],dim = mp)) )
  if(out.method=="matrix") return(x)
  else return(x0)
}