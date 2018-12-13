# Testing for true GP; known all parameters;

library(GPfit)
library(MASS)
library(dplyr)
library(purrr)
library(tidyr)
library(reshape2)
library("forcats")
library(mvtnorm)
library(rbenchmark)
library(microbenchmark)
library("coda")
library(ggplot2)
theme_set(theme_bw())
library(MCMCpack)
rm(list=ls())

## read the parameters
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID', 0))
pars=scan(file = 'tsarpars_1119.txt',what = "",skip = TASKID, nlines = 1)%>%as.numeric(.)
#pars=c(40,15,2,0.01,0.01,0.01,0.01,2)
## SImulation dataset
N=pars[1]
V=pars[2]
H=pars[3]
phi_u=pars[4]
phi_x=pars[5]
s_u=pars[6]
s_x=pars[7]
phi_u0=pars[8]
phi_x0=pars[9]
s_u0=pars[10]
s_x0=pars[11]
H_star=pars[12]
a1=pars[13]
a2=pars[14]
set.seed(2)
#c_u=corr_matrix(X = 1:N,beta = log10(k_u),corr = list(type="exponential",power=2))
#c_x=corr_matrix(X = 1:N,beta = log10(k_x),corr = list(type="exponential",power=2))
c_u<-matrix(0,N,N)
#s_u=v_u/(1-phi_u^2)
for (i in 1:N){
  for (j in 1:N){
    c_u[i,j]<-phi_u^(abs(i-j))*s_u
  }}
c_x<-matrix(0,N,N)
#s_x=v_x/(1-phi_x^2)
for (i in 1:N){
  for (j in 1:N){
    c_x[i,j]<-phi_x^(abs(i-j))*s_x
  }}
mu=mvrnorm(n=1,mu=rep(0,N),Sigma = c_u)
X=mvrnorm(n=V*H,mu=rep(0,N),Sigma = c_x)
X_arr=array(data = as.vector(X),dim = c(V,H,N)) 
pi_arr=array(dim=c(V,V,N))
Y_arr=array(dim=c(V,V,N))
for (t in 1:N){
  pi_arr[,,t]=1/(1+exp(-(X_arr[,,t]%*%t(X_arr[,,t])+mu[t])))
  for (i in 2:V){
    for (j in 1:(i-1)){
      Y_arr[i,j,t]=sample(x = c(1,0),size = 1,prob = c(pi_arr[i,j,t],1-pi_arr[i,j,t]))
    }
  }
}
#save(Y_arr,file="Y_arr_sys_sed111.RData")

## prepare
# Naive sampler for PG(1, z)
# (based on the finite approximation of infinite sum)
rpg_naive = function(z, n_terms = 100){
  g = rexp(n_terms, 1)
  out = 1 / (2*(pi^2)) * sum(g / ((1:n_terms - 1/2)^2 + z^2 / (4*(pi^2))))
  return(out)
}


## Inference

N=dim(Y_arr)[3]
V=dim(Y_arr)[1]
Y_arr[,,N]=NA

#a1=a2=2
niter=5000
burnin=1000
set.seed(456)


# MCMC start here

now=proc.time()
#MCMCme=function(N=40,V=15,H_star=10,k_u0=0.05,k_x0=0.05,a1=2,a2=2,niter=5000,burnin=1000,Y_arr=Y_arr){

##prior
c_u0<-matrix(0,N,N)
#s_u0=v_u0/(1-phi_u0^2)
for (i in 1:N){
  for (j in 1:N){
    c_u0[i,j]<-phi_u0^(abs(i-j))*s_u0
  }}
c_x0<-matrix(0,N,N)
#s_x0=v_x0/(1-phi_x0^2)
for (i in 1:N){
  for (j in 1:N){
    c_x0[i,j]<-phi_x0^(abs(i-j))*s_x0
  }}
#c_u0=corr_matrix(X = 1:N,beta = log10(k_u0),corr = list(type="exponential",power=2))
K_mu_inv=chol2inv(chol(c_u0+diag(rep(1e-8,N)))) #inverse
#c_x0=corr_matrix(X = 1:N,beta = log10(k_x0),corr = list(type="exponential",power=2))
K_x_inv=chol2inv(chol(c_x0+diag(rep(1e-8,N)))) #inverse
mu0=mvrnorm(n=1,mu=rep(0,N),Sigma = c_u0)
v1=rgamma(n = 1,shape = a1,rate = 1)
vl=rgamma(n = H_star-1,shape = a2,rate = 1) #(H_star-1)*1 vec
vv=c(v1,vl)
#vv=rep(1,H_star)
tao=accumulate(vv,prod) #in durante code: 1 1 1 1 1 .. as starting value
X_arr0=array(dim=c(V,H_star,N))
for(h in 1:H_star){
  X_arr0[,h,]=mvrnorm(n=V,mu=rep(0,N),Sigma = (1/tao[h])*c_x0)
}

#Store the result: create empty bags

#W_cache=matrix(nrow=V*V*N,ncol=niter)
mu_cache=matrix(nrow=N,ncol=niter)
#X_cache=matrix(nrow=V*H_star*N,ncol=niter)
tao_cache=matrix(nrow=H_star,ncol=niter)
pi_arr_est_cache=matrix(nrow=V*V*N,ncol=niter)
Y_arrN_cache=array(dim=c(V,V,niter))
#preparation:transfer Y_arr to symmetric matrix without diagonal
#for (t in 1:N){
#  Y_arr[,,t][upper.tri(Y_arr[,,t])] = t(Y_arr[,,t])[upper.tri(Y_arr[,,t])]
#}
for (t in 1:N){
  Y_arr[,,t]=xpnd(vech(Y_arr[,,t]),V)
}
pi_arr0=map(1:N,~1/(1+exp(-X_arr0[,,.x]%*%t(X_arr0[,,.x])-mu0[.x])))

for (i in 2:V){
  for (j in 1:(i-1)){
    Y_arr[i,j,N]=sample(x = c(1,0),size = 1,prob = c(pi_arr0[[N]][i,j],1-pi_arr0[[N]][i,j]))
  }
}
Y_arr[,,N]=xpnd(vech(Y_arr[,,N]),V)
##posterior
for (iter in 1:niter){
  
  #sample W:using X_arr0,mu0
  S=map(1:N,~X_arr0[,,.x]%*%t(X_arr0[,,.x])) #S=Xt(X)
  W=array(data = NA,dim=c(V,V,N))
  for (t in 1:N){
    W[,,t][lower.tri(W[,,t])]=map_dbl(S[[t]][lower.tri(S[[t]])]+mu0[t],~rpg_naive(z = .x,n_terms = 100))
  }
  
  #sample mu0:using W
  Sigma_mu0=(apply(X = W,MARGIN = 3,FUN = sum,na.rm=TRUE))%>%diag(.)
  Sigma_mu=chol2inv(chol(Sigma_mu0+K_mu_inv))
  mu_mu=Sigma_mu%*%map_dbl(1:N,~sum(Y_arr[,,.x]-0.5-S[[.x]]*W[,,.x],na.rm = TRUE)) #W-upper NA
  mu0=mvrnorm(1,mu_mu,Sigma_mu)
  mu_cache[,iter]=mu0
  
  #sample X_arr0:using tao,W,X_arr0
  prior_sig_x=diag(tao) %x% K_x_inv
  #for (t in 1:N){
  #  W[,,t][upper.tri(W[,,t])] = t(W[,,t])[upper.tri(W[,,t])]
  #}
  for (t in 1:N){
    W[,,t]=xpnd(vech(W[,,t]),V)
  }
  
  for (v in 1:V){
    X_tilta=matrix(0,nrow = (V-1)*N,ncol = H_star*N)
    w=W[v,,]%>%t(.)%>%as.vector(.)%>%na.omit(.)%>%as.vector()
    Omega=diag(w)
    for (t in 1:N){
      X_tilta[seq(from=t,to=(V-1)*N,by=N),seq(from=t,to=H_star*N,by=N)]=X_arr0[-v,,t]
    }
    Sigma_x0=t(X_tilta)%*%Omega%*%X_tilta+prior_sig_x
    Sigma_x=chol2inv(chol(Sigma_x0+diag(rep(1e-8,N*H_star))))
    y=Y_arr[v,,]%>%t(.)%>%as.vector(.)%>%na.omit(.)%>%as.vector()
    mu_x0=t(X_tilta)%*%(y-0.5-w*rep(mu0,V-1))
    mu_x=Sigma_x%*%mu_x0
    X_v=rmvnorm(n = 1,mean = mu_x,sigma = Sigma_x)
    X_arr0[v,,]=matrix(data = X_v,nrow=H_star,byrow = TRUE)
    
  }
  
  #X_cache[,iter]=X_arr0%>%as.vector()
  
  ##sample tao,v:using Xr_arr0,Xs_arr0,vv
  xKx=map_dbl(.x = 1:H_star,.f = function(l){map_dbl(.x = 1:V,.f = ~ t(X_arr0[.x,l,])%*%K_x_inv%*%X_arr0[.x,l,])%>%sum(.)})
  
  tao_1=c(1,vv[-1])%>%accumulate(.,prod)
  rate1=1+0.5*sum(tao_1*xKx)
  v1=rgamma(n = 1,shape = a1+V*N*H_star/2,rate = rate1)
  rate_l=vector(length = H_star-1)
  for (h in 2:H_star){
    rate_l[h-1]=map_dbl(h:H_star,~prod(vv[1:.x][-h])*xKx[.x])%>%sum(.)*0.5+1
  }
  vl=map_dbl(2:H_star,~rgamma(n = 1,shape = a2+V*N*(H_star-.x+1)/2,rate = rate_l[.x-1]))
  vv=c(v1,vl)
  tao=accumulate(vv,prod)
  tao_cache[,iter]=tao
  
  
  #calculate pi!!!Finally!!!:using X_arr0,mu0,
  pi_arr0=map(1:N,~1/(1+exp(-X_arr0[,,.x]%*%t(X_arr0[,,.x])-mu0[.x])))
  pi_vec=pi_arr0%>%unlist(.)
  pi_arr_est_cache[,iter]=pi_vec

    for (i in 2:V){
      for (j in 1:(i-1)){
        Y_arr[i,j,N]=sample(x = c(1,0),size = 1,prob = c(pi_arr0[[N]][i,j],1-pi_arr0[[N]][i,j]))
      }
    }
  Y_arr[,,N]=xpnd(vech(Y_arr[,,N]),V)
  Y_arrN_cache[,,iter]=Y_arr[,,N]
}
#MCMCres=list(mu_est=mu_cache,pi_est=pi_arr_est_cache)
#save(MCMCres,file="MCMCres_sys_seed456.RData")
#return(MCMCres)
#}

#givemeMCMC=MCMCme(N=40,V=15,H_star=10,k_u0=0.05,k_x0=0.05,a1=2,a2=2,niter=5000,burnin=1000,Y_arr=Y_arr)

future=proc.time()
future-now
save.image(paste0("tsar1119_mis_tau","_ID",TASKID,".RData"),compress = "xz")
