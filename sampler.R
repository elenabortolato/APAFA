#Gibbs sampling for flexible MSFA
library(Rcpp)

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



cppFunction('
arma::mat truncnorm_lg(const arma::mat& y_lower, const arma::mat& y_upper,
                       const arma::mat& mu, const arma::vec& sigma,
                       const arma::mat& u_rand) {
  // Dim of matrix:
  int n = y_lower.n_rows;
  int p = y_lower.n_cols;
  // Storage:
  double val = 0;
  arma::mat z_star(n,p);

  for(int t = 0; t < n; ++t) {
    for(int j = 0; j < p; ++j) {
      // Control
      double uptail1 = (y_lower(t,j) - mu(t,j)) * 1 / sigma(j) > 8;
      double uptail2 = (y_upper(t,j) - mu(t,j)) * 1 / sigma(j) > 8;
      // pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
      // true + false = 1, true + true = 2, false + false = 0
      if((uptail1 + uptail2) == 0){
        // Lower and upper limits, transformed via pnorm:
        double F_lower = R::pnorm(y_lower(t,j), mu(t,j), sigma(j), 1, 0);
        double F_upper = R::pnorm(y_upper(t,j), mu(t,j), sigma(j), 1, 0);
        // replace 0 with 0.000001 and 1 with 0.999999
        if (F_lower == 0) {
          F_lower = 0.000001;
        } else if (F_upper == 1) {
          F_lower = 0.999999;
        }
        if (F_upper == 0) {
          F_upper = 0.000001;
        } else if (F_upper == 1) {
          F_upper = 0.999999;
        }
        // Corresponding sampled value:
        val = R::qnorm(F_lower + u_rand(t,j) * (F_upper - F_lower), mu(t,j), sigma(j), 1, 0);
      }
      else {
        double F_lower = R::pnorm(y_lower(t,j), mu(t,j), sigma(j), 0, 0);
        double F_upper = R::pnorm(y_upper(t,j), mu(t,j), sigma(j), 0, 0);
        // replace 0 with 0.000001 and 1 with 0.999999
        if (F_lower == 0) {
          F_lower = 0.000001;
        } else if (F_upper == 1) {
          F_lower = 0.999999;
        }
        if (F_upper == 0) {
          F_upper = 0.000001;
        } else if (F_upper == 1) {
          F_upper = 0.999999;
        }
        // Corresponding sampled value:
        val = R::qnorm(F_lower + u_rand(t,j) * (F_upper - F_lower), mu(t,j), sigma(j), 0, 0);
      }
      z_star(t,j) = std::min(std::max(y_lower(t,j), val), y_upper(t,j));
    }
  }
  
  return z_star;
}', depends = "RcppArmadillo")

cppFunction('
arma::mat runif_mat(const int& rows, const int& cols, const double& minVal, const double& maxVal) {
  arma::mat M(rows, cols);
  for (int j = 0; j < cols; j++) {
    for (int i = 0; i < rows; i++) {
      M(i, j) = R::runif(minVal, maxVal);
    }
  }
  return M;
}', depends = "RcppArmadillo")



# sampler for non gaussian data
Gibbs_Kernel_non_gauss<-function(state, family="logistic"){
  if(family=="logistic"){
    if(all((unique(state$y))%in%c(0,1))){#changes only at first iteration
      state$obs=state$y
      state$y=matrix(rnorm(prod(dim(state$y))), ncol=ncol(state$y))
    }
    
    mean_y=t(sapply(1:state$n, function (i) state$Lambda%*%state$eta[i,]+ state$Gamma%*%state$phi[i,])) 
    runif_mat=matrix(runif(state$n*state$p),state$n, state$p)
    lb=ifelse(state$obs==0,-15, 0)
    ub=ifelse(state$obs==0,0, 15)
    state$y = truncnorm_lg((lb), (ub), mean_y, matrix(1 , state$p), runif_mat)
  }
  
  if(is.null(state$scale_beta)) state$scale_beta=0.1
  if(is.null(state$X)) state$X=matrix(0.1, ncol=state$d, nrow = state$n)
  if(is.null(state$ps)) state$s=matrix(rbinom(state$n*state$k),0.1, ncol=state$k)
  
  #1.update factors
  I=diag(state$d+state$k)
  invS=diag(state$p)
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
  
  #2.update Sigmas
    state$Sigma = diag(state$p)
  
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
    B=diag(length(betas))*(state$n^-1)*0.1# -0.0095/state$S*(state$n^-1) 
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
  
  # 9 update psi
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

 
