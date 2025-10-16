#this file reproduces the simulation study with large p  (p=1000)
rm(list=ls())
source("sampler.R")

##########################################################
# seeds used to generate the data (scenario A) for p=1000
seeds=1201:1210
##########################################################
# 'run1sim': this function creates the data,
# initializes the model,
# saves the true parameters (factor, loadings, variance..)  in a .RDS file
# runs MCMC for 10k iterations and
# stores the results for all quantities as a .RData file 
run1sim<-function(setting){ # setting will be an integer from 1 to 10
    set.seed(seeds[setting])
    S=3 # studies
    ns=rep(33,S) #units per study
    n=sum(ns) # total units 
    p=1000 # observed variables
    Sigma <- (diag(p))*0.1 #errors
    eps=matrix(mvtnorm::rmvnorm(n, sigma = Sigma), ncol = p, nrow=n)
    dim(eps)
    #shared latent factors and loadings
    d=100# shared latent factors
    shared_loadings=rnorm(p*d,sd = 1)
    Lambda=matrix(shared_loadings, ncol = d, nrow = p)
    eta=matrix(0, ncol = d, nrow = n)
    for (h in 1:d) {eta[,h]=rnorm(n,1)}
    #specific latent factors and loadings
    ks=rep(NA, S)
    ks=c(5,5,5)#specific factors study  
    k=sum(ks)
    group=rep(NA,n)
    for (s in 1:S) {
      nscumpre=ifelse(s>1,sum(ns[1:(s-1)])+1,1)
      nscum=sum(ns[1:s])
      group[nscumpre:nscum]=s
    }
    #dummy variables group
    #X=matrix(NA, ncol=S+1, nrow=n)
    X=model.matrix(rep(1,n)~-1+as.factor(group))
    dim(X)
    specific_loadings=rnorm(p*k, sd=1)
    Gamma=matrix(specific_loadings, ncol = k, nrow = p)
    #specific factors 
    phi_=phi=matrix(0, ncol = k, nrow = n)
    for (h in 1:k) {phi[,h]=phi_[,h]=rnorm(n,1)}
    
    #sparsify
    for (s in 1:k) {
      phi[,s]=phi[,s]*(group==s)
    }
    phi[,4:k]=0
    for(h in 1:k){
      sdh=sd(phi_[,h])
      Gamma[,h]= Gamma[,h]*(sdh)
      phi_[,h]= phi_[,h]/(sdh)
      phi[,h]= phi[,h]/(sdh)
      
    }
    for(h in 1:d){
      sdh=sd(eta[,h])
      Lambda[,h]= Lambda[,h]*(sdh)
      # Lambda_[,h]= Lambda_[,h]*(sdh)
      eta[,h]= eta[,h]/(sdh)
    }
    
    Lambda_=Lambda 
    Lambda
    Lambda[,4:d]=Lambda_[,4:d]= 0
    head(Lambda_)
    head(phi)
    #observations
    y=t(sapply(1:n, function (i) Lambda%*%eta[i,]+ Gamma%*%phi[i,]))+eps 
    dim(y)
    #hyperparameters 
    alpha_eta=3
    v_eta = c( rbeta(d-1, shape1 = 1, shape2 = alpha_eta), 1)
    w_eta = v_eta*c(1,cumprod(1-v_eta[-d]))                     # weights
    z_eta = rep(d,d)
    alpha_phi=3
    v_phi = c( rbeta(k-1, shape1 =1, shape2 = alpha_phi), 1)
    w_phi = v_phi*c(1,cumprod(1-v_phi[-k]))                     # weights
    z_phi = rep(k,k)
    betas=matrix(0, ncol=k, nrow=S)
    (plogis(betas ))
  

#####################################################
# save all the useful quantities for the MCMC kernel
state=list(y=y,# response
             Lambda=Lambda,  Lambda_=Lambda_, # shared loadings  #sparse and non sparse  
             eta=(eta) , #shared factors
             Gamma=Gamma, #specific loadings
             phi= (phi), phi_= phi_, #sparse and non sparse specific factors
             tau_eta=c(rep(1,d-1),rep(0,1)),
             tau_phi= c(rep(1,k-1),rep(0,1)), 
             z_eta =z_eta, z_phi= z_phi, 
             w_eta=w_eta, w_phi=w_phi, 
             v_eta=v_eta, v_phi=v_phi,
             Sigma=Sigma,# list of specific covariance matrices
             n=n, ns=ns, X=X, S=S, d=d, k=k,  p=p, 
             #prior
             a_lambda=3, b_lambda=1, a_gamma=3, b_gamma=1,
             a_sigma=5, b_sigma=10,  
             alpha_eta=3, alpha_phi=3, # equal to expected n.of active factors 
             betas=betas, 
             ps=matrix(rbinom(n*k,1,0.50), nrow=n)
  )
  
copy_state=state
saveRDS(copy_state, paste("copy_state_p1000", setting, sep=""))
  
state$Lambda_=jitter(state$Lambda_, amount = 1)
state$Lambda=state$Lambda_
state$Gamma=jitter(state$Gamma,amount =1)
state$eta=jitter(state$eta, amount = 1)
state$phi_=jitter(state$phi_, amount = 1)
state$phi=state$phi_
state$Sigma=jitter(state$Sigma, amount = 0.1)
maxiter=10000
  
if(1==1){ris_sigma1=matrix(0,ncol=p,maxiter)
  ris_beta=array(0,dim = c( S, k,maxiter))
  
  ris_tau_eta=matrix(0,ncol=d,maxiter)
  ris_tau_phi=matrix(0,ncol=k,maxiter)
  ris_phi=array(0, dim = c( k, n, maxiter))  
  ris_eta=array(0, dim = c( d, n, maxiter))  
  ris_psi=array(0, dim = c( k, n, maxiter))  
  
  ris_lambda=array(0, dim = c( d,p, maxiter))   
  ris_gamma=array(0, dim = c( k,p, maxiter))   
  }


# Gibbs sampler
par(mar=c(2,2,2,2))
time<-system.time(for (iter in 1:maxiter){
    # if(iter==1) image(state$ps, main=0)
    cat(iter)
    state= Gibbs_Kernel(state)
    #shrinkage tau
    ris_sigma1[iter,]=diag(state$Sigma)
    ris_tau_eta[iter,]=state$tau_eta
    ris_tau_phi[iter,]=state$tau_phi
    #factors
    ris_phi[,,iter]=state$phi 
    
    ris_eta[,,iter]=state$eta 
    
    #beta 
    ris_beta[,,iter]=state$betas 
     
    #theta psi
    ris_psi[,,iter]=state$ps 
    
    #loading
    ris_lambda[,,iter]=state$Lambda 
    
    ris_gamma[,,iter]=state$Gamma 
   
    #print number of active factors
    if(iter%%2==0){ 
      image(state$ps)
      print(state$tau_eta)
      print(state$tau_phi)}
   })

print(time)
image=list(ris_sigma1=ris_sigma1,
           ris_tau_eta=ris_tau_eta,
           ris_tau_phi=ris_tau_phi,
           #factors
           ris_phi=ris_phi,
          ris_eta =ris_eta,
           #beta 
           ris_beta=ris_beta,
           #theta psi
           ris_psi=ris_psi,
           #loadings
           ris_lambda=ris_lambda ,
           ris_gamma= ris_gamma 
  
)
saveRDS(image, paste(setting,"largep1000.RDS",sep = ""))
}
i
# 10 replications 
for (i in 1:10) run1sim(i)
