#SUPPLEMENTARY: Resistance to identifiability  issues
######################################################
# data generating process (scenario A) ns=20, S=3, d=k=3
# Specific factor loadings are nearly identical.
# (some noise is added from initial identical columns)
# the resulting correlation between specific factor loadings
# in three scenarios will be: 0.72,0.84,0.95

source("sampler.R")


seeds=c(1234) 
replication=1
mat_all=list(3)
for(replication in 1:1 ){
set.seed(seeds[1])
jittering_amount=c(0.2,0.4,0.6)
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
  specific_loadings=rnorm(p, sd=1)
  
  Gamma=matrix(specific_loadings, ncol = k, nrow = p)
   
  # identical specific factor loadings
  Gamma[1:p,]=Gamma[1:p,]+rnorm(p*d, mean = 0, sd=jittering_amount[replication])
 
  if (replication==1) Gamma_replication1=Gamma
  
  print("correlation")
  par(mar=c(3,3,3,3)+1)
  par(mfrow=c(1,3))
  plot(Gamma[1:10,1:2], xlab=expression(Gamma[1]),ylab=expression(Gamma[2]))
  #abline(0,1,col=2,  ylim=c(-1,1))
  plot(Gamma[1:10,c(1,3)], xlab=expression(Gamma[1]),ylab=expression(Gamma[3]))
  plot(Gamma[1:10,2:3], xlab=expression(Gamma[2]),ylab=expression(Gamma[3]))

  print(mean(cor(Gamma[,1:3])[lower.tri(cor(Gamma[,1:3]))]))
  print((cor(Gamma[,1:3])[lower.tri(cor(Gamma[,1:3]))]))
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
  head(Gamma)
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
           a_sigma=20, b_sigma=2,  
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
  #cat(iter)
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
  if(iter%%1000==0){ 
    print(state$tau_eta)
    print(mean(diag(state$Sigma)))
    print(state$tau_phi)}
})

boxplot(ris_sigma1[8000:10000,])


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
 # print(iter)
}

# convergence assessment with
# number of factors found
par(mar=c(2,2,0.2,0.2)+2)
par(mfrow=c(2,1))
plot(Loglik, cex=.4,  xlab="iteration",col=rowSums(ris_tau_phi))
legend("bottomright", col=c(5,4,3), lty=1,
       legend =  paste("k=",c(5:3)))


# posterior means
mat <- cbind(
  colMeans(ris_ps1[8000:10000, ]),
  colMeans(ris_ps2[8000:10000, ]),
  colMeans(ris_ps3[8000:10000, ])
)

mat_all[[replication]]=mat
}

# Set up layout for image + legend
par(mfrow=c(3,1))
layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1))

# Define color palette (0–1 scale)
cols <- colorRampPalette(c("white",  "red"))(100)
zlim <- c(0, 1)  # proportion range

# Plot the images
corr_est=c(0.95,0.84,0.72)
for (replication in 1:3){
image(
  t(apply(mat_all[[replication]], 2, rev)),        
  col = cols,
  zlim = zlim,
  axes = FALSE,
  xlab=expression(Psi),
  main = paste("correlation", corr_est[replication])
)
axis(1, at = seq(0, 1, length.out = ncol(mat)), labels = 1:ncol(mat))
axis(2, at = seq(0, 1, length.out = nrow(mat)), labels = nrow(mat):1)

# --- Add color legend ---
par(mar = c(5, 2, 4, 4))
image(
  z = matrix(seq(zlim[1], zlim[2], length.out = 100), nrow = 1),
  col = cols,
  axes = FALSE
)
axis(4, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1)
 
}





######################################################
# data generating process (scenario A) ns=30
# Specific factor loadings are nearly identical.

seeds=c(1234) 
replication=1
mat_all=list(1)
for(replication in 1:1 ){
  set.seed(seeds[1])
  jittering_amount=c(0.2,0.4,0.6)
  if(1==1){
    S=3 # studies
    ns=rep(40,S) #units per study
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
    
    
    # use same Gamma as before
    if (replication==1) Gamma=Gamma_replication1
    
    
    print("correlation")
    par(mar=c(3,3,3,3)+1)
    par(mfrow=c(1,3))
    plot(Gamma[1:10,1:2], xlab=expression(Gamma[1]),ylab=expression(Gamma[2]))
    #abline(0,1,col=2,  ylim=c(-1,1))
    plot(Gamma[1:10,c(1,3)], xlab=expression(Gamma[1]),ylab=expression(Gamma[3]))
    plot(Gamma[1:10,2:3], xlab=expression(Gamma[2]),ylab=expression(Gamma[3]))
    
    print(mean(cor(Gamma[,1:3])[lower.tri(cor(Gamma[,1:3]))]))
    print((cor(Gamma[,1:3])[lower.tri(cor(Gamma[,1:3]))]))
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
    head(Gamma)
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
             a_sigma=20, b_sigma=2,  
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
    #cat(iter)
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
    if(iter%%1000==0){ 
      print(state$tau_eta)
      print(mean(diag(state$Sigma)))
      print(state$tau_phi)}
  })
  
  boxplot(ris_sigma1[8000:10000,])
  
  
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
    # print(iter)
  }
  
  # convergence assessment with
  # number of factors found
  par(mar=c(2,2,0.2,0.2)+2)
  par(mfrow=c(2,1))
  plot(Loglik, cex=.4,  xlab="iteration",col=rowSums(ris_tau_phi))
  legend("bottomright", col=c(5,4,3), lty=1,
         legend =  paste("k=",c(5:3)))
  
  
  # posterior means
  mat <- cbind(
    colMeans(ris_ps1[8000:10000, ]),
    colMeans(ris_ps2[8000:10000, ]),
    colMeans(ris_ps3[8000:10000, ])
   
  )
  
  mat_all[[replication]]=mat
}

# Set up layout for image + legend
par(mfrow=c(3,1))
layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1))

# Define color palette (0–1 scale)
cols <- colorRampPalette(c("white",  "red"))(100)
zlim <- c(0, 1)  # proportion range

# Plot the images
corr_est=c(0.95,0.84,0.72)
for (replication in 1:1){
  image(
    t(apply(mat_all[[replication]], 2, rev)),        
    col = cols,
    zlim = zlim,
    axes = FALSE,
    xlab=expression(Psi),
    main = paste("ns=40, correlation", corr_est[replication])
  )
  axis(1, at = seq(0, 1, length.out = ncol(mat)), labels = 1:ncol(mat))
  axis(2, at = seq(0, 1, length.out = nrow(mat)), labels = nrow(mat):1)
  
  # --- Add color legend ---
  par(mar = c(5, 2, 4, 4))
  image(
    z = matrix(seq(zlim[1], zlim[2], length.out = 100), nrow = 1),
    col = cols,
    axes = FALSE
  )
  axis(4, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1)
  
}
