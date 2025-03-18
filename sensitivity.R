#SUPPLEMENTARY S1: sensitivity analysis
library(Rcpp) 
library(MASS)
library(mvtnorm)
library(MCMCpack)
library(calculus)
library(unbiasedmcmc) 
library(pgdraw)
library(ggplot2) 


### prior for elements of sigma
library(vioplot)
par(mfrow=c(1,2))
set.seed(11)
sample_20_2=rinversegamma(100000,20,2)
sample_5_5=rinversegamma(100000,5,5)


plot(1 ,type = "n", xlim = c(0.005, 2), ylim = range(c(sample_20_2)), 
     xaxt = "n",   xlab = "", ylab = "")
vioplot(sample_20_2,   add = TRUE,   axes = FALSE)
text(1.5,0.3,"a=20, b=2",cex=0.75)
plot(1 ,type = "n", xlim = c(0.005, 2), ylim = range(c(sample_5_5)), 
     xaxt = "n",   xlab = "", ylab = "")
vioplot(sample_5_5,   add = TRUE,   axes = FALSE)
text(1.5,16.42,"a=5, b=5",cex=0.75)

# to perform Sensitivity analysis, change lines 437-446


##########
# setup
standardize=TRUE # always in the sensitivity analysis

#polya gamma data augmentation
cppFunction('
List pgg_m_sigma_(const Eigen::Map<Eigen::MatrixXd>  & omega,
                       const Eigen::Map<Eigen::MatrixXd>  & X,
                       const Eigen::Map<Eigen::MatrixXd>  & invB,
                       const Eigen::Map<Eigen::VectorXd>  & KTkappaplusinvBtimesb){
  int n = X.rows();
  int p = X.cols();
  // The matrix A stores XT Omega X + B^{-1}, that is, Sigma^{-1}
  Eigen::MatrixXd A(p,p);
  for (int j1 = 0; j1 < p; j1 ++){
    for (int j2 = j1; j2 < p; j2 ++){
      A(j1,j2) = invB(j1, j2);
      for (int i = 0; i < n; i++){
        A(j1,j2) = A(j1,j2) + X(i,j1) * X(i,j2) * omega(i);
      }
      A(j2,j1) = A(j1,j2);
    }
  }
  Eigen::LLT<Eigen::MatrixXd> lltofA(A);
  Eigen::MatrixXd lower = lltofA.matrixL();
  Eigen::VectorXd x = lltofA.solve(KTkappaplusinvBtimesb);
  return List::create(Named("m")=x,
                      Named("Sigma_inverse") = A,
                      Named("Cholesky_inverse") = lower,
                      Named("Cholesky") = lower.inverse());
}', depends = "RcppEigen")




# sampler for gaussian data
#--------------------------------------------
# state: state of the chain, list containig
#
# list(n=n, #n. of units, scalar
#      p=p, # n. of observed variables per unit
#      S=S, #n. of groups
#      ns=ns, #n. of units, vector of length S 
#      X=X, # design matrix n x S of dummy variables for groups
#      d=d, k=k, # max number of shared and specific factors
#      y=y, # n x p matrix of responses

# initialization of parameters
#      Lambda=Lambda,  Lambda_=Lambda, # p x d matrix of shared loadings (sparse and non sparse)
#      eta= eta,  n x d matrix of factors
#      Gamma=Gamma, # p x d matrix of specific loadings
#      phi= phi, phi_= phi_, # n x k sparse and non sparse specific factors
#      ps=matrix(rbinom(n*k,1,0.5), ncol=k)# matrix of local activation for phi 
#      Sigma=Sigma,# p x p covariance matrix
#      betas=betas,  matrix of dimension S x k)
#    
#prior hyper parameters
#      scale_beta=0.1 # suggested: 1/n #scale hyperparameter for beta#      
#      a_sigma=2, b_sigma=2,# scalars, InverseGamma hyperparameters for Sigma
#      a_load=a_load, b_load=b_load,#(d+k) vectors of InverseGamma hyperparameters for Lambda and Gamma
#      alpha_eta=10, alpha_phi=6, # hyperparameters of the CUSP process (n. active factors)

#initialization CUSP
#      tau_eta=c(rep(1,d), rep(0, d-d)), # global shrinkage initialization
#      tau_phi=c(rep(1,k), rep(0, k-k)), 
#      z_eta =z_eta, z_phi= z_phi, # vectors of length d and k
#      w_eta=w_eta, w_phi=w_phi, # vectors of length d and k
#      v_eta=v_eta, v_phi=v_phi, # vectors of length d and k)

Gibbs_Kernel=function(state){
  
  if(is.null(state$scale_beta)) state$scale_beta=0.1
  if(is.null(state$X)) state$X=matrix(0.1, ncol=state$d, nrow = state$n)
  if(is.null(state$ps)) state$s=matrix(rbinom(state$n*state$k),0.1, ncol=state$k)
  
  #1.update factors
  I=diag(state$d+state$k)
  invS=diag(1/diag(state$Sigma))
  mean_update=sapply(1:state$n, 
                     function (i) solve(I+rbind(t(state$Lambda),t(state$Gamma)*state$tau_phi*state$ps[i,] )%*%invS%*%
                                          t(rbind(t(state$Lambda), t(state$Gamma)*state$tau_phi*state$ps[i,] )))%*%
                       (rbind(t(state$Lambda),t(state$Gamma)*state$tau_phi*state$ps[i,] ))%*%invS%*%
                       (state$y[i,])) 
  var_update=lapply(1:state$n,
                    function (i) solve(I+rbind(t(state$Lambda),t(state$Gamma)*state$tau_phi*state$ps[i,] )%*%invS%*%
                                         t(rbind(t(state$Lambda), t(state$Gamma)*state$tau_phi*state$ps[i,] ))) )
  
  factors=sapply(1:state$n, function(i) rmvnorm(1,rep(0,state$d+state$k), var_update[[i]]))+(mean_update)
  state$eta=t(factors)[,1:state$d]
  state$phi_=t(factors)[,-c(1:state$d)] #non sparse
  state$phi=t(factors)[,-c(1:state$d)]*t(state$tau_phi*t(state$ps)) #sparse
  
  
  if(standardize==T){
    #standardize
    for(h in 1:state$d){
      sdh=sd(state$eta[,h])
      state$Lambda_[,h]= state$Lambda_[,h]*(sdh)
      state$Lambda[,h]= state$Lambda[,h]*(sdh)
      state$eta[,h]= state$eta[,h]/(sdh)
    }
    for(h in 1:state$k){
      sdh=sd(state$phi_[,h])
      state$Gamma[,h]= state$Gamma[,h]*(sdh)
      state$phi_[,h]= state$phi_[,h]/(sdh)
      state$phi[,h]= state$phi[,h]/(sdh)
    }
  }
  #2.update Sigmas
  Ytil = state$y - tcrossprod(state$eta,state$Lambda)-tcrossprod(state$phi,state$Gamma)
  invsig = rgamma(state$p, state$a_sigma+state$n/2, state$b_sigma+0.5*colSums(Ytil^2))
  state$Sigma = diag(1/invsig)
  
  #3.update betas
  pgg_m_and_sigma <- function(omega, precomputed){
    return(pgg_m_sigma_(omega, precomputed$X, precomputed$invB, precomputed$KTkappaplusinvBtimesb))
  }
  pgg_precomputation <- function(Y, X, b, B){
    invB <- solve(B)
    invBtimesb <- invB %*% (b)
    Ykappa <- matrix(Y - rep(0.5, length(Y)), ncol=1)
    XTkappa <- t(X) %*% Ykappa
    KTkappaplusinvBtimesb <- XTkappa + (invBtimesb)
    return(list(n=nrow(X), p=ncol(X), X=X, Y=Y, b=b, B=B,
                invB=invB, invBtimesb=invBtimesb, KTkappaplusinvBtimesb=KTkappaplusinvBtimesb))
  } 
  pred =  state$X%*%(state$betas)
  logit_phi = plogis(pred)
  ps_= matrix(1, nrow = state$n, ncol = state$k)
  logit_phi0 = logit_phi[which(state$ps==0)]
  p_constant=1 
  which_zero = which(runif(length(logit_phi0))<
                       ((1-logit_phi0)/(1- logit_phi0*p_constant)))
  ps_[ which(state$ps==0)[which_zero] ] = 0
  
  for(h in 1:state$k){
    state$betas=matrix(state$betas, ncol=state$k)
    betas=state$betas[,h]
    #variance prior for betas
    B=diag(length(betas))*(state$n^-1)*state$scale_beta# -0.0095/state$S*(state$n^-1) 
    pgg_precomputed <- pgg_precomputation(ps_[,h], state$X, as.matrix(betas), B)
    pgg_kernel <- function(beta){
      zs <- (pgg_precomputed$X %*% beta)
      w <- pgdraw::pgdraw(1, zs)
      res <- pgg_m_and_sigma(w, pgg_precomputed)
      beta <- unbiasedmcmc:::fast_rmvnorm_chol(1, res$m, res$Cholesky)[1,]
      return(list(beta = beta))
    }
    beta_ <- pgg_kernel(betas)
    state$betas[,h]=beta_$beta
  }
  
  #4. Update precision matrix of loadings
  a_load=c(rep(state$a_lambda, state$d),rep(state$a_gamma, state$k))
  b_load=c(rep(state$b_lambda, state$d),rep(state$b_gamma, state$k))
  loadings=cbind(state$Lambda_, state$Gamma)# non sparse
  Prec = diag(rgamma(state$d+state$k,a_load+0.5*state$p, b_load+0.5*colSums(loadings^2)))
  
  #5. Update the loadings  
  factors=cbind(t(t(state$eta)*(state$tau_eta)), t(t(state$phi))) #sparse
  for (j in 1:state$p){
    var_update=solve((Prec)+t(factors)%*%factors*1/state$Sigma[j,j])
    mean_update=var_update%*%t(factors)%*%(state$y[,j])*1/state$Sigma[j,j]
    Loadings=rmvnorm(1,rep(0,state$d+state$k), var_update)+t(mean_update)
    state$Lambda_[j,]=Loadings[1:state$d]
    state$Gamma[j,]=Loadings[-c(1:state$d)]
  }
  state$Lambda= t(state$tau_eta*t(state$Lambda_)) #sparse
  
  #6. update z
  index(state$Lambda_) = c("j","h")
  index(state$eta) = c("i", "h")
  eta_lam = einstein(state$eta, state$Lambda_, drop = F)  # n x p x k
  mu_eta = tcrossprod( state$eta,state$Lambda)
  mu_phi = tcrossprod( state$phi,state$Gamma)
  mu=mu_eta+mu_phi
  sdy=matrix( rep(sqrt(diag(state$Sigma)),state$n), state$n, state$p, byrow=T)
  for(h in 1:state$d){
    mu_0 = mu - state$tau_eta[h]*eta_lam[,,h]
    mu_1 = mu_0 + eta_lam[,,h]
    f0 = sum(dnorm(state$y, mean= mu_0, sd=sdy, log=T))
    f1 = sum(dnorm(state$y, mean= mu_1, sd=sdy, log=T))
    mf = max(c(f0,f1))
    f0 = f0 - mf
    f1 = f1 - mf
    prob_h = exp( c(rep(f0, h), rep(f1, state$d-h)) +log(state$w_eta))
    if (sum(prob_h)==0){
      prob_h = c(rep(0,state$d-1), 1)
    } else{
      prob_h = prob_h/sum(prob_h)
    }
    state$z_eta[h] = which(rmultinom(n=1, size=1, prob=prob_h)==1)
  }
  
  #7 update tau_eta
  state$tau_eta = rep(1,state$d)
  state$tau_eta[state$z_eta <= seq(1,state$d)]=0
  state$Lambda= t( (state$tau_eta)*t(state$Lambda_))
  
  # 8 --   Update v_eta and w_eta -- #
  for(h in 1:(state$d-1)){
    state$v_eta[h] = rbeta(1, shape1 =1+ sum(state$z_eta==h), 
                           shape2 = state$alpha_eta+sum(state$z_eta>h))
  }
  state$v_eta[state$d] = 1
  state$w_eta = state$v_eta*c(1,cumprod(1-state$v_eta[-state$d]))
  
  mu_eta = tcrossprod( state$eta,state$Lambda)
  ps_phi = state$phi_*state$ps
  index(state$Gamma) = c("j", "h")
  index(ps_phi) = c("i","h")
  phi_ps_gamma= (einstein( (ps_phi),(state$Gamma),drop = F))  # n x p x k
  
  mu=mu_eta+mu_phi
  
  # 6bis update zeta_phi
  for(h in 1:state$k){
    mu_0 = mu - state$tau_phi[h]*phi_ps_gamma[,,h]
    mu_1 = mu_0 + phi_ps_gamma[,,h]
    f0 = sum(dnorm(state$y, mean= mu_0, sd=(sdy), log=T))
    f1 = sum(dnorm(state$y, mean= mu_1, sd=(sdy), log=T))
    mf = max(c(f0,f1))
    f0 = f0 - mf
    f1 = f1 - mf
    prob_h = exp( c(rep(f0, h), rep(f1, state$k-h)) +log(state$w_phi))
    if (sum(prob_h)==0){
      prob_h = c(rep(0,state$k-1), 1)
    } else{
      prob_h = prob_h/sum(prob_h)
    }
    state$z_phi[h] = which(rmultinom(n=1, size=1, prob=prob_h)==1)
  }
  #7bis
  state$tau_phi = rep(1,state$k)
  state$tau_phi[state$z_phi <= seq(1,state$k)]=0
  if(length(state$tau_phi)==0) (state$tau_phi=rep(0, state$k))
  state$phi=state$phi_*t(state$tau_phi*t(state$ps)) 
  
  # 8bis --  Update v_phi and w_phi -- #
  for(h in 1:(state$k-1)){
    state$v_phi[h] = rbeta(1, shape1 = 1+sum(state$z_phi==h),
                           shape2 = state$alpha_phi+sum(state$z_phi>h))
  }
  state$v_phi[state$k] = 1
  state$w_phi = state$v_phi*c(1,cumprod(1-state$v_phi[-state$k]))
  
  # 9 update ps
  mu_eta = tcrossprod(state$eta,state$Lambda)
  pred =  state$X%*%(state$betas)
  logit_phi = plogis(pred)
  mu_phi= tcrossprod(state$phi,state$Gamma)
  mu=mu_eta+mu_phi
  tau_phi =t(t(state$phi_)*state$tau_phi )
  index(tau_phi) = c("i","h")
  index(state$Gamma)=c("j", "h")
  
  phi_tau_gam = einstein(tau_phi,state$Gamma, drop = F)  # n x p x k
  for(h in 1:state$k){
    mu_0 = mu - (phi_tau_gam[,,h])*state$ps[,h]
    mu_1 = mu_0 + phi_tau_gam[,,h]
    f0 = rowSums(dnorm(state$y,   mean= mu_0, sd=(sdy), log=T))
    f1 = rowSums(dnorm(state$y,  mean= mu_1, sd=(sdy), log=T))
    mf =      max(c(f0   ,f1   ))
    f0 = f0 - mf
    f1 = f1 - mf
    lp_phi0 = f0 + log(1-logit_phi[,h]*1)
    lp_phi1 = f1 + log(logit_phi[,h]*1)
    sumlog = apply(cbind(lp_phi0, lp_phi1),1, matrixStats::logSumExp)
    state$ps[,h] =  round( runif(state$n) < exp(lp_phi1-sumlog) )
  }
  state$phi=state$phi_*t(state$tau_phi*t(state$ps)) # sparse
  
  
  #reorder active factors (specific)
  if(sum(state$tau_phi)>0){
    idx_act=which(state$tau_phi==1)
    idx_non_act=c(1:state$k)[-idx_act]
    state$phi=state$phi[,c(idx_act,  idx_non_act)]
    state$phi_=state$phi_[,c(idx_act,  idx_non_act)]
    state$ps=state$ps[,c(idx_act,  idx_non_act)]
    state$tau_phi=state$tau_phi[c(idx_act,  idx_non_act)]
    state$betas=state$betas[,c(idx_act,  idx_non_act)]
    state$Gamma=state$Gamma[,c(idx_act,  idx_non_act)]
  }
  
  #reorder active factors (shared)
  idx_act=which(state$tau_eta==1)
  idx_non_act=c(1:state$d)[-idx_act]
  state$eta=state$eta[,c(idx_act,  idx_non_act)]
  state$tau_eta=state$tau_eta[c(idx_act,  idx_non_act)]
  state$Lambda=state$Lambda[,c(idx_act,  idx_non_act)]
  state$Lambda_=state$Lambda_[,c(idx_act,  idx_non_act)]
  
  
  return(state)
}





######################################################
# data generating process (always scenario A for sensitivity analysis)
setting=1
seeds=c(1234:1243) 
set.seed(seeds[setting])
if(1==1){
  S=3 # studies
  ns=rep(20,S) #units per study
  n=sum(ns) # total unitsB=
  p=10# observed variables
  Sigma <- (diag(p))*0.1 #errors
  eps=matrix(mvtnorm::rmvnorm(n, sigma = Sigma), ncol = p, nrow=n)
  
  #shared latent factors and loadings
  d=6# shared latent factors
  shared_loadings=rnorm(p*d,sd = 1)
  Lambda=matrix(shared_loadings, ncol = d, nrow = p)
  eta=matrix(NA, ncol = d, nrow = n)
  for (h in 1:d) {eta[,h]=rnorm(n,1)}
  #specific latent factors and loadings
  ks=rep(NA, S)
  ks=c(2,2,2)#specific factors study  
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
  phi_=phi=matrix(NA, ncol = k, nrow = n)
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
  Lambda[,4:d]= 0
  head(Lambda)
  head(phi)
  #observations
  y=t(sapply(1:n, function (i) Lambda%*%eta[i,]+ Gamma%*%phi[i,]))+eps 
  
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
}


###################################
# save all the quantities useful for the MCMC kernel


state=list(y=y,# response
           Lambda=Lambda,  Lambda_=Lambda_, # shared loadings  #sparse and non sparse  
           eta=(eta) , #shared factors
           Gamma=Gamma, #specific loadings
           phi= (phi), phi_= phi_, #sparse and non sparse specific factors
           tau_eta=c(rep(1,6),rep(0,d-6)),
           tau_phi= c(rep(1,6),rep(0,k-6)), 
           z_eta =z_eta, z_phi= z_phi, 
           w_eta=w_eta, w_phi=w_phi, 
           v_eta=v_eta, v_phi=v_phi,
           Sigma=Sigma,# list of specific covariance matrices
           n=n, ns=ns, X=X, S=S, d=d, k=k,  p=p, 
           #prior
           a_lambda=1, b_lambda=2, 
           a_gamma=1, b_gamma=2,
           
#sensitivity  for sigma 
            # a_sigma=20, b_sigma=2,  
          a_sigma=5, b_sigma=5,  
            #  a_sigma=5, b_sigma=1,  
#sensitivity  for alpha_eta and   alpha_phi
            # alpha_eta=1, alpha_phi=1,
            # alpha_eta=2, alpha_phi=2,
            # alpha_eta=1, alpha_phi=1,
           alpha_eta=4, alpha_phi=4,
           betas=betas, 
           ps=matrix(rbinom(n*k,1,0.50), nrow=n)
)

copy_state=state
 
state$Lambda_=jitter(state$Lambda_, amount = 1)
state$Lambda=state$Lambda_
state$Gamma=jitter(state$Gamma,amount =1)
state$eta=jitter(state$eta, amount = 1)
state$phi_=jitter(state$phi_, amount = 1)
state$phi=state$phi_
state$Sigma=jitter(state$Sigma, amount = 0.1)
maxiter=10000

if(1==1){ris_sigma1=matrix(0,ncol=p,maxiter)
ris_beta1=matrix(0,ncol=S,maxiter)
ris_beta2=matrix(0,ncol=S,maxiter)
ris_beta3=matrix(0,ncol=S,maxiter)
ris_beta4=matrix(0,ncol=S,maxiter)
ris_beta5=matrix(0,ncol=S,maxiter)
ris_beta6=matrix(0,ncol=S,maxiter)
ris_tau_eta=matrix(0,ncol=d,maxiter)
ris_tau_phi=matrix(0,ncol=k,maxiter)
ris_phi1=matrix(0,ncol=n,maxiter)
ris_phi2=ris_phi3=ris_phi4=ris_phi5=ris_phi6=
  ris_ps1=ris_ps2=ris_ps3=ris_ps4=ris_ps5=ris_ps6=matrix(0,ncol=n,maxiter)
ris_eta1=ris_eta2=ris_eta3=ris_eta4=ris_eta5=ris_eta6=matrix(0,ncol=n,maxiter)
ris_lambda1=ris_lambda2=ris_lambda3=ris_lambda4=ris_lambda5=ris_lambda6=matrix(0,ncol=p,maxiter)
ris_gamma1=ris_gamma2=ris_gamma3=ris_gamma4=ris_gamma5=ris_gamma6=matrix(0,ncol=p,maxiter)
}


# Gibbs sampler
setting
set.seed(setting)#setting 
par(mar=c(2,2,2,2))
time<-system.time(for (iter in 1:maxiter){
  # if(iter==1) image(state$ps, main=0)
  cat(iter)
  state= Gibbs_Kernel(state)
  #shrinkage tau
  ris_tau_eta[iter,]=state$tau_eta
  ris_tau_phi[iter,]=state$tau_phi
  #factors
  ris_phi1[iter,]=state$phi[,1]
  ris_phi2[iter,]=state$phi[,2]
  ris_phi3[iter,]=state$phi[,3]
  ris_phi4[iter,]=state$phi[,4]
  ris_phi5[iter,]=state$phi[,5]
  ris_phi6[iter,]=state$phi[,6]
  ris_eta1[iter,]=state$eta[,1]
  ris_eta2[iter,]=state$eta[,2]
  ris_eta3[iter,]=state$eta[,3]
  ris_eta4[iter,]=state$eta[,4]
  ris_eta5[iter,]=state$eta[,5]
  ris_eta6[iter,]=state$eta[,6]
  #beta 
  ris_beta1[iter,]=state$betas[,1]
  ris_beta2[iter,]=state$betas[,2]
  ris_beta3[iter,]=state$betas[,3]
  ris_beta4[iter,]=state$betas[,4]
  ris_beta5[iter,]=state$betas[,5]
  ris_beta6[iter,]=state$betas[,6]
  
  #theta psi
  ris_ps1[iter,]=state$ps[,1]
  ris_ps2[iter,]=state$ps[,2]
  ris_ps3[iter,]=state$ps[,3]
  ris_ps4[iter,]=state$ps[,4]
  ris_ps5[iter,]=state$ps[,5]
  ris_ps6[iter,]=state$ps[,6]
  #loading
  ris_lambda1[iter,]=state$Lambda[,1]
  ris_lambda2[iter,]=state$Lambda[,2]
  ris_lambda3[iter,]=state$Lambda[,3]
  ris_lambda4[iter,]=state$Lambda[,4]
  ris_lambda5[iter,]=state$Lambda[,5]
  ris_lambda6[iter,]=state$Lambda[,6]
  ris_gamma1[iter,]=state$Gamma[,1] 
  ris_gamma2[iter,]=state$Gamma[,2]
  ris_gamma3[iter,]=state$Gamma[,3]
  ris_gamma4[iter,]=state$Gamma[,4]
  ris_gamma5[iter,]=state$Gamma[,5]
  ris_gamma6[iter,]=state$Gamma[,6]
  ris_sigma1[iter,]=diag(state$Sigma)
  #print iteration and number of active factors
  if(iter%%10==0){ 
    print(state$tau_eta)
    print(mean(diag(state$Sigma)))
    print(state$tau_phi)}
})

time 
 

# convergence assessment: plot the  Loglikelihood 
iter=1
Loglik=rep(0, 10000)
for(iter in 1:10000){
  L=cbind(ris_lambda1[iter,],ris_lambda2[iter,],ris_lambda3[iter,],
          ris_lambda4[iter,],ris_lambda5[iter,],ris_lambda6[iter,])
  L=t(t(L)*ris_tau_eta[iter,1:6])
  LLT=L%*%t(L)
  G=cbind(ris_gamma1[iter,],ris_gamma2[iter,],ris_gamma3[iter,],
          ris_gamma4[iter,],ris_gamma5[iter,],ris_gamma6[iter,])
  G=t(t(G)*ris_tau_phi[iter,1:6])
  PS=(cbind(ris_ps1[iter,],ris_ps2[iter,],ris_ps3[iter,],
            ris_ps4[iter,],ris_ps5[iter,],ris_ps6[iter,])) 
  GGT=lapply(1:n, function (ii ) G%*%diag(PS[ii,1:6])%*%
               t(G))
  S=diag(ris_sigma1[iter,])
  Loglik[iter]=0
  
  for(ii in 1:60){
    VAR=LLT+GGT[[ii]]+S
    Loglik[iter]=Loglik[iter]+dmvnorm(y[ii,],rep(0,p), sigma = VAR, log=T)
  }
  print(iter)
}

par(mar=c(2,2,0.2,0.2)+2)
par(mfrow=c(2,1))
plot(Loglik, cex=.4,  xlab="iteration", ,col=rowSums(ris_tau_eta))
legend("bottomright", col=c(5,4,3), lty=1, 
       legend =  paste("d=",c(5:3)))
plot(Loglik, cex=.4,  xlab="iteration",col=rowSums(ris_tau_phi))
legend("bottomright", col=c(5,4,3), lty=1,
       legend =  paste("k=",c(5:3)))

### need to sligthly adjust the legend and the colors
   