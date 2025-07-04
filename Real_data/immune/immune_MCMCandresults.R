library(MSFA)
library(BiocManager)
library(Rcpp)
library(mvtnorm)
library(MCMCpack)
library(calculus)
library(unbiasedmcmc) 
library(pgdraw)
library(ggplot2) 
#working directory
setwd("C:/Users/Asus/Dropbox/FactorModels/data_immune")
#saveRDS(file="immune_data.RDS",y)  
y=readRDS(file="immune_data.RDS")
n=sum(ns) # total units
p=ncol(dat) # observed variables

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
 
#run MCMC
iter=1
set.seed(1234)
for (iter in iter:maxiter){
  cat(iter)
  state= Gibbs_Kernel(state)
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
    save.image("immune_res.RData")}
}

  
#######################################################################################
#LOAD DATA
memory.limit(100000)
load("immune_res.RData")
 
#active factors
colMeans(ris_tau_eta)
sort.list(abs(colMeans(ris_lambda[,,1])))

plot(colMeans(ris_lambda[,,2]))
order(colMeans(ris_lambda[,,2]))
 
L=colMeans(colMeans(ris_tau_eta[5000:10000, ])*(ris_lambda[5000:10000,, ]))
LL=tcrossprod(L,L)

 
G=colMeans(colMeans(ris_tau_phi[5000:10000, ])*(ris_gamma[5000:10000,, ]))
image(G)
 
 
#reorder specific factors according to activation %
 

GG1=tcrossprod(G[,1],G[,1])
GG2=tcrossprod(G[,2],G[,2])
GG3=tcrossprod(G[,3],G[,3])
GG4=tcrossprod(G[,4],G[,4])
GG5=tcrossprod(G[,5],G[,5])
GG6=tcrossprod(G[,6],G[,6])
GG7=tcrossprod(G[,7],G[,7])
GG8=tcrossprod(G[,8],G[,8])
GG9=tcrossprod(G[,9],G[,9])
  
#### better colors
par(mfrow=c(3,3))
GG1[1,1]=max(abs(GG1))
GG1[p,p]=-max(abs(GG1))
GG2[1,1]=max(abs(GG2))
GG2[p,p]=-max(abs(GG2))
GG3[1,1]=max(abs(GG3))
GG3[p,p]=-max(abs(GG3))
GG4[1,1]=max(abs(GG4))
GG4[p,p]=-max(abs(GG4))
GG5[1,1]=max(abs(GG4))
GG5[p,p]=-max(abs(GG5))
GG6[1,1]=max(abs(GG6))
GG6[p,p]=-max(abs(GG6))
GG7[1,1]=max(abs(GG7))
GG7[p,p]=-max(abs(GG7))
GG8[1,1]=max(abs(GG8))
GG8[p,p]=-max(abs(GG8))
GG9[1,1]=max(abs(GG9))
GG9[p,p]=-max(abs(GG9))
 
 
par(mar=c(4,3,3,3))
# Define the matrix to be visualized 
# Define the color palette with red and blue
custom_colors <- c("red","white", "blue")

# Create a custom color function that maps values to colors
color_function <- colorRampPalette(custom_colors)
colors <- color_function(0.3*length(GG1))
# Plot the image with the custom color palette

# Figure 9: contribution of Gamma Gamma
par(mfrow=c(3,3))
image(GG3, axes=F,xlab=expression(Gamma[1]*Gamma[1]^T),col = colors[sort.list(t(GG3))])
image(GG4, axes=F,xlab=expression(Gamma[2]*Gamma[2]^T),col = colors[sort.list(GG4)])
image(GG5, axes=F,xlab=expression(Gamma[3]*Gamma[3]^T),col = colors[sort.list(GG1)])
image(GG1, axes=F,xlab=expression(Gamma[4]*Gamma[4]^T),col = colors[sort.list(GG2)])
image(GG2, axes=F,xlab=expression(Gamma[5]*Gamma[5]^T),col = colors[sort.list(GG5)])
image(GG6, axes=F,xlab=expression(Gamma[6]*Gamma[6]^T),col = colors[sort.list(GG6)])
image(GG7, axes=F,xlab=expression(Gamma[7]*Gamma[7]^T),col = colors[sort.list(GG7)])
image(GG8, axes=F,xlab=expression(Gamma[8]*Gamma[8]^T),col = colors[sort.list(GG8)])
image(GG9, axes=F,xlab=expression(Gamma[9]*Gamma[9]^T),col = colors[sort.list(GG9)])
 


ps_hat=apply(ris_ps,c(2,3),function(x) mean((x)))
phi_hat=apply(ris_phi[5000:10000,,],c(2,3),function(x) mean((x)))

par(mfrow=c(1,1))
par(mar=c(4,4,4,4))
image((phi_hat[,1:9]), col=colors[sort.list(t(phi_hat))],axes=F, xlab="units",ylab=expression(Phi))
abline(v=cumsum(ns/n))
sapply(1:9,function (k) expression(phi[k]))
#axis(philab, at=x.at)
philab=c("phi[1]","phi[2]","phi[3]","phi[4]","phi[5]","phi[6]","phi[7]","phi[8]","phi[9]")
y=1:9
axis(2, at=y/9*1.091-0.1,labels = 1:9)
abline(v=0)
abline(h=1.06)
abline(h=-0.06)



