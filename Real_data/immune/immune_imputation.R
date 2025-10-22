#working directory
setwd("/data_immune")
#load data
y=readRDS(file="immune_data.RDS")
p=ncol(y)
n=nrow(y)
ns=c(85, 140, 578, 195)
S=4
######################################################
#initialization
Sigma=(diag(p))*0.5
PCA_est=princomp(y,scores = T)
fact_est=PCA_est 
#shared latent factors and loadings
d=30 # shared latent factors
Lambda=matrix(fact_est$loadings[,1:d], ncol = d, nrow = p)
Lambda
eta=eta_=matrix(NA, ncol = d, nrow = n)
    
for (h in 1:d) {eta[,h]=eta_[,h]=fact_est$scores[h]}
#specific latent factors and loadings
ks=rep(NA, S)
ks=c(3,3,3,3)#n. of specific factors/study 
k=sum(ks)
group=rep(NA,n)
for (s in 1:S) {
    nscumpre=ifelse(s>1,sum(ns[1:(s-1)])+1,1)
    nscum=sum(ns[1:s])
    group[nscumpre:nscum]=s
  }
group

#dummy variables for groups
X=matrix(NA, ncol=S, nrow=n)
X[,1:S]=model.matrix(rep(1,n)~-1+as.factor(group))
tail(X)


specific_loadings=rnorm(p*k, sd=1)
Gamma=matrix(specific_loadings, ncol = k, nrow = p)
phi_=phi=matrix(NA, ncol = k, nrow = n)
for (h in 1:k) {phi[,h]=phi_[,h]=rnorm(n)}

#hyperparameters 
alpha_eta=5
v_eta = c( rbeta(d-1, shape1 = 1, shape2 = alpha_eta), 1)
w_eta = v_eta*c(1,cumprod(1-v_eta[-d]))                     # weights
z_eta = rep(d,d)
alpha_phi=6
v_phi = c( rbeta(k-1, shape1 = 1, shape2 = alpha_phi), 1)
w_phi = v_phi*c(1,cumprod(1-v_phi[-k]))                     # weights
z_phi = rep(k,k)
whichgroup=unique(sapply(1:n ,function (k) (which(X[k,]==1))))
betas=matrix(0, ncol=k, nrow=S, byrow = F)
betas
plogis(betas )
a_lambda=2
b_lambda=2
a_gamma=2
b_gamma=2
a_load=c(rep(a_lambda,d),rep(a_gamma, k))
b_load=c(rep(b_lambda,d),rep(b_gamma, k))

# save all the quantities useful for the MCMC kernel
  state=list(y=y,# response
             Lambda=Lambda,  Lambda_=Lambda, # shared loadings
             eta= (eta),  #sparse and non sparse shared factors
             Gamma=Gamma, #specific loadings
             phi= (phi), phi_= (phi_), #sparse and non sparse specific factors
             Sigma=Sigma,# list of specific covariance matrices
             n=n, ns=ns, X=X, S=S, d=d, k=k,  p=p,
             #prior
             a_sigma=2, b_sigma=2, tau_eta=c(rep(1,d), rep(0, d-d)), 
             tau_phi=c(rep(1,k), rep(0, k-k)), 
             z_eta =z_eta, z_phi= z_phi, 
             w_eta=w_eta, w_phi=w_phi, 
             v_eta=v_eta, v_phi=v_phi, 
             alpha_eta=10, alpha_phi=6, 
             scale_beta=0.1,# equal tO expected n.of active factors
              a_load=a_load, b_load=b_load,
             betas=betas,   
             ps=matrix(rbinom(n*k,1,0.5), ncol=k))
    
    
  
#### arrays for storing results
maxiter=10000
ris_y=array(dim=c(maxiter,dim(state$y)) )
ris_phi=array(dim=c(maxiter,dim(state$phi)) )
ris_eta=array(dim=c(maxiter,dim(state$eta)) )
ris_ps=array(dim=c(maxiter,dim(state$ps)) )
ris_beta=array(dim=c(maxiter,dim(state$betas)) )
ris_eta=array(dim=c(maxiter,dim(state$eta)) )
ris_lambda=array(dim=c(maxiter,dim(state$Lambda)) )
ris_gamma=array(dim=c(maxiter,dim(state$Gamma)) )
ris_tau_eta=array(dim=c(maxiter,length(state$tau_eta)) )
ris_tau_phi=array(dim=c(maxiter,length(state$tau_phi)) )
ris_sigma1=matrix(0,ncol=p,maxiter)
 

y_na=y
# leave out 30 observations, as they were missing.
set.seed(11)
idx_imp=matrix(NA,ncol=2, nrow=30)
idx_imp[,1]=trunc(runif(30,1,n))
idx_imp[,2]=trunc(runif(30,1,p))
idx_imp
 
# function to draw from the predictive distribution, 
# using draws from the posterior (state$...)
conditional_d<- function (i,j) {
  VAR=tcrossprod(t(t(state$Lambda)*state$tau_eta))+
    tcrossprod(t(t(state$Gamma)*(state$tau_phi*state$ps[i,])))+state$Sigma
  V11=VAR[j,j]
  V22=VAR[-j,-j]
  V12=VAR[j,-j]
  V21=VAR[-j,j]
  
  Vc=V11-V12%*%solve(V22)%*%V21
  Mc=mean(y[,j])+ V12%*%solve(V22)%*%(y[i,-j]-apply(y,2,mean)[-j])
  rnorm(1,c(Mc), c(Vc))
}  

#run MCMC with data integration
iter=1
set.seed(1234)
for (iter in iter:maxiter){
  cat(iter)
  for(i in 1:30){
    ii=idx_imp[i,1]
    jj=idx_imp[i,2]
    y_na[idx_imp[i,1],idx_imp[i,2]]=conditional_d(ii,jj)
  }
  
  state$y=y_na
  cat(iter)
  state= Gibbs_Kernel(state)
  ris_y[iter,,]=state$y
  #shrinkage CUSP
  ris_tau_eta[iter,]=state$tau_eta
  ris_tau_phi[iter,]=state$tau_phi
  #factors
  ris_phi[iter,,]=state$phi
  ris_eta[iter,,]=state$eta
  #beta 
  ris_beta[iter,,]=state$betas
  #psi
  ris_ps[iter,,]=state$ps
  #loadings
  ris_lambda[iter,,]=state$Lambda
  ris_gamma[iter,,]=state$Gamma
  #sigma
  ris_sigma1[iter,]=diag(state$Sigma)
  #print iteration and number of active factors
  print(state$tau_eta)
  image(state$ps)
  print(state$tau_phi)
  if(iter%%200==0){
    print(iter)
    #save progress
    save.image("immune_res_imputation.RData")}

  }

# results
#load("immune_res_imputation.RData")

MSE=rep(0,10000)
for(iter in iter:10000) MSE[iter]=  sum((ris_y[iter,,]-y)^2)/30
plot(MSE, type="l")

boxplot((MSE[1:8000]))
boxplot((MSE[10000:8000]))
 

diff=matrix(0, 2000,30)
for(iter in 8001:10000){
for(ii in 1:30) diff[iter-8000,ii]=(y[idx_imp[ii,1],idx_imp[ii,2]]-ris_y[iter,idx_imp[ii,1],idx_imp[ii,2]]) 
}
boxplot((diff))



pred=matrix(0, 2000,30)
for(iter in 8001:10000){
  for(ii in 1:30) pred[iter-8000,ii]=(ris_y[iter,idx_imp[ii,1],idx_imp[ii,2]]) 
}
yobs=matrix(0,1,30)
 
  for(ii in 1:30) yobs[1,ii]=(y[idx_imp[ii,1],idx_imp[ii,2]]) 

# figure 9
par(mfrow=c(1,2))

plot(MSE, type="l", xlab="iteration")
mean(MSE)
boxplot((pred), ylab="prediction")
points(1:30,yobs, col=2,pch=8,lwd=0.5,cex=01.8)
