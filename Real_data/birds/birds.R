rm(list=ls())

#setwd("/Users/elenabortolato/Desktop/FactorModels/birds")
data.directory = "./birds"
da = read.csv(file.path(data.directory, "data.csv"), stringsAsFactors=TRUE)
 

#Route: location
da$Route = as.factor(da$Route)
dim(da)

# Year: year of repeated observations
summary(da$Year)
unique(da$Route)
dim(da)[1]
#914 observations

length(unique(da$Route))
#number of "studies"=200

#levels(da$Route)=c(1:200)
dim(da)[1]/length(unique(da$Route))
#4 observations per location (number of units per study)
colSums(table(da$Route,da$Year))
# more locations for recent years
 
boxplot(rowSums(da[,c(10:59)])~da$Year)
#total counts increasing
 
# We store candidate environmental covariates
# habitat type and spring temperature
# even if also many other environmental variables could be expected to influence the commmunity.
#XData = data.frame(Route = da$Route, hab=da$Habitat, clim = da$AprMay)
# scale the numeric data
#XData$clim = scale(XData$clim)

# 
colnames(da)[-c(1:9)]
Y = matrix((da[,-c(1:9)])>0, nrow = nrow(da))
Y = apply(Y,MARGIN = 2,FUN = as.numeric)
Y
colnames(Y)=colnames(da)[-c(1:9)]
colnames(Y)
# We next read the datafile containing species traits, and include in the TrData dataframe data on migratory strategy and body mass

#alltraits = read.csv(file.path(data.directory, "traits.csv"), stringsAsFactors = TRUE)
#cbind(as.character(alltraits[,1]), (colnames(da)[-c(1:9)]))
#rownames(alltraits)
#TrData = data.frame(Species = alltraits$Species, Migration=alltraits$Migration, LogMass = log(alltraits$Mass))
# scale the numeric data
#TrData$LogMass = scale(TrData$LogMass)


# X and Traits formulae
#XFormula = ~ hab + poly(clim, degree = 2,raw = TRUE) -1
#TrFormula = ~ Migration + LogMass -1


# we now should find groups of birds defined by the philogenetic tree
#install.packages("ape")
library(ape)

phyloTree = read.tree(file.path(data.directory, "CTree.tre"))
plot(phyloTree)

# from pica pica to prunella modularis (40) we have the same ordine: passeriformes
# dendrocopos_major and dryocopus_martius has order piciformes
#  Cuculiformes
#  gruiformes
#  Charadriiformes
#  Columbiformes
#  tetrao tetrix
# 7 different orders with one of them that has 40 elements: not so beautiful.

# to separate the big group we can consider the family or superfamily
# 1-5 corvidae
# 6-18 Sylvioidea
# 19 Certhioidea
# 20-29 muscicapoidea
# 30 rguloidea
# 31-40 passeroidea
# an alternative is using
# corvidae (5) - Sylvioidea (13) - muscicapoidea (10) - passeroidea (10) as further variables.
# we can decide to ignore or consider also picidae (2) - Scolopacidae (3)

# Finally we identified 4 or 6 group variables (depending if we ignore the smallest groups or not).

# New TrData

corvoidea =  sylvioidea = muscicapoidea = passeroidea = picidae = scolopacidae = rep(0,50) 
corvoidea[which(colnames(Y) %in% phyloTree$tip.label[1:5])] = 1
sylvioidea[which(colnames(Y) %in% phyloTree$tip.label[6:18])] = 1
muscicapoidea[which(colnames(Y) %in% phyloTree$tip.label[20:29])] = 1
passeroidea[which(colnames(Y) %in% phyloTree$tip.label[31:40])] = 1
picidae[which(colnames(Y) %in% phyloTree$tip.label[41:42])] = 1
scolopacidae[which(colnames(Y) %in% phyloTree$tip.label[45:47])] = 1

TrData_taxonomy = data.frame(TrData, 
                           corvoidea, sylvioidea,muscicapoidea, passeroidea, picidae,  scolopacidae)
TrFormula_taxonomy = ~ Migration + LogMass + corvoidea + 
 sylvioidea + muscicapoidea + passeroidea +picidae + scolopacidae -1
TrData_taxonomy[,4]

 
# Define quantities of interest and initialize the algorithm
d=10
k=10
p=ncol(Y)

ns=table(da$Route)
S=length(unique(da$Route)) # studies
S
ns 
n=sum(ns) # total units
n 
 
#shared latent factors and loadings
# shared latent factors
shared_loadings=rnorm(p*d,sd = 1)
Lambda=matrix(shared_loadings, ncol = d, nrow = p)
Lambda
eta=matrix(NA, ncol = d, nrow = n)
for (h in 1:d) {eta[,h]=rnorm(n)}
ks=rep(NA, S)
ks=c(rep(1,S))#specific factors study (important just for the ground truth in simulations)
 
 
ks[cumsum(ks)>k]=0
sum(ks)
group=rep(NA,n)

for (s in 1:S) {
  nscumpre=ifelse(s>1,sum(ns[1:(s-1)])+1,1)
  nscum=sum(ns[1:s])
 group[nscumpre:nscum]=s
}

group
table(group)
#dummy variables group
X=matrix(NA, ncol=S, nrow=n)
X[,1:S]=model.matrix(rep(1,n)~-1+as.factor(group))
 
specific_loadings=rnorm(p*k, sd=1)
Gamma=matrix(specific_loadings, ncol = k, nrow = p)
phi_=phi=matrix(NA, ncol = k, nrow = n)
for (h in 1:k) {phi[,h]=phi_[,h]=rnorm(n)}
for (h in 1:k) {phi[,h]=phi_[,h]*(group[h]==h)}
head(phi)

#hyperparameters 
alpha_eta=5
v_eta = c( rbeta(d-1, shape1 = 1, shape2 = alpha_eta), 1)
w_eta = v_eta*c(1,cumprod(1-v_eta[-d]))                     # weights
z_eta = rep(d,d)
alpha_phi=5
v_phi = c( rbeta(k-1, shape1 = 1, shape2 = alpha_phi), 1)
w_phi = v_phi*c(1,cumprod(1-v_phi[-k]))                     # weights
z_phi = rep(k,k)
whichgroup=unique(sapply(1:n ,function (k) (which(X[k,]==1))))
betas=matrix(0, ncol=k, nrow=S) 
plogis(betas )
 
 
# define the state of the chain
state=list(y=matrix(as.numeric(da[,-c(1:9)]>0), nrow = nrow(da)),# response
           Lambda=Lambda+rnorm(prod(dim(Lambda))),  
                               Lambda_=Lambda+rnorm(prod(dim(Lambda)), sd=0.1), # shared loadings
           eta= (eta),  # shared factors
           Gamma=Gamma+rnorm(prod(dim(Gamma)), sd=0.1), #specific loadings
           phi= (phi) , phi_= (phi_)+rnorm(prod(dim(phi_)), sd=0.1), #sparse and non sparse specific factors
        # list of specific covariance matrices
           n=n, ns=ns, X=X, S=S, d=d, k=k,  p=p,
           #prior
        a_lambda=1, b_lambda=2, a_gamma=1, b_gamma=2,
            tau_eta=c(c(rep(1,10)),c(rep(0,d-10))), 
        tau_phi=c(c(rep(1,10)),c(rep(0,k-10))), 
           z_eta =z_eta, z_phi= z_phi, 
           w_eta=w_eta, w_phi=w_phi, 
           v_eta=v_eta, v_phi=v_phi, 
           alpha_eta=5, alpha_phi=5, # equal tO expected n.of active factors
           betas=betas, 
           ps=matrix(rbinom(n*k,1,0.1),ncol=k)   )
 
copy_state=state
 table(state$y)
#######################################################################################################################
#one iteration of Gibbs sampling
library(Rcpp)
library(mvtnorm)
library(MCMCpack)
library(calculus)
library(unbiasedmcmc) 
library(pgdraw)
library(ggplot2)  


Gaussian=F
maxiter=10000
ris_phi=array(dim=c(maxiter,dim(state$phi)) )
ris_eta=array(dim=c(maxiter,dim(state$eta)) )
ris_th=array(dim=c(maxiter,dim(state$th)) )
ris_ps=array(dim=c(maxiter,dim(state$ps)) )
 
ris_beta=array(dim=c(maxiter,dim(state$betas)) )
ris_eta=array(dim=c(maxiter,dim(state$eta)) )
ris_lambda=array(dim=c(maxiter,dim(state$Lambda)) )
ris_gamma=array(dim=c(maxiter,dim(state$Gamma)) )
ris_tau_eta=array(dim=c(maxiter,length(state$tau_eta)) )
ris_tau_phi=array(dim=c(maxiter,length(state$tau_phi)) )
set.seed(111)
for (iter in 1:10000){
  cat(iter)
  state= Gibbs_Kernel_non_gauss(state)
 #shrinkage tau
  ris_tau_eta[iter,]=state$tau_eta
  ris_tau_phi[iter,]=state$tau_phi
  #factors
  ris_phi[iter,,]=state$phi
  ris_eta[iter,,]=state$eta
  
  #beta 
  ris_beta[iter,,]=state$betas
   
  
  #theta and psi
  
  ris_ps[iter,,]=state$ps
   
  #loadings
  ris_lambda[iter,,]=state$Lambda
  ris_gamma[iter,,]=state$Gamma
   
  
  #print iteration and number of active factors
  print(state$tau_eta)
  print(state$tau_phi)
 
}
