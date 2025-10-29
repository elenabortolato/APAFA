"Bird species occurence dataset"
---

Bortolato, E. and Canale, A.
October 2025



```{r, echo=FALSE}
library(dplyr)
#source("../../sampler.R")
```
 

## Description

This script reproduces the analysis of the *birds* dataset (Section 4.1). The dataset records co-occurrence patterns of \(p=50\) bird species observed in Finland over nine years (2006–2014) across $S=200$ locations.

We fit a factor model to study species occurrences using only the *group* information (observation site), as specified in the APAFA framework.

Subsequently, we examine the relationships between the available environmental covariates and the latent factors inferred by the model *a posteriori*.


## Data preprocessing 

Let $Y$ be the $n\times p$ matrix containing  the counts (presence/absence) of species in the S=200 location (group) over time.

 

```{r}
data = read.csv( ("data/data.csv"), stringsAsFactors=TRUE)
Y = data %>%
  select(-c(1:9)) %>%
  as.matrix()
```

We next load the information on the species traits:

```{r}
traits <- read.csv("traits.csv") %>% 
  mutate(Migration = as.factor(Migration)) %>% 
  mutate(Species = as.factor(Species))
TrFormula_taxonomy = ~ Migration + LogMass + corvoidea +
  sylvioidea + muscicapoidea + passeroidea + picidae + scolopacidae - 1
```
 
## Initializations 


We define the latent dimensions and group-specific sample sizes
```{r}
d = 10
k = 10
p = ncol(Y)

ns = table(data$Route) # units in each study
S = length(unique(data$Route)) # studies
n = sum(ns) # total units in all studies
```
 
 
   
Initialization of shared loadings and factors

```{r}
shared_loadings = rnorm(p * d, sd = 1)
Lambda = matrix(shared_loadings, ncol = d, nrow = p)
eta = matrix(NA, ncol = d, nrow = n)
for (h in 1:d) {
  eta[, h] = rnorm(n)
}
```

Initialization of specific loadings and factors and related quantities

```{r}
ks = c(rep(1, S)) #specific factors study (important just for the ground truth in simulations)
ks[cumsum(ks) > k] = 0 # set some factors to "inactive"

#assign group labels (location labels)
group = rep(NA, n)
for (s in 1:S) {
  nscumpre = ifelse(s > 1, sum(ns[1:(s - 1)]) + 1, 1)
  nscum = sum(ns[1:s])
  group[nscumpre:nscum] = s
}
table(group)
# convert into dummy variables the group labels
X = matrix(NA, ncol = S, nrow = n)
X[, 1:S] = model.matrix(rep(1, n) ~ -1 + as.factor(group))

# initialize specific loadings and factors
specific_loadings = rnorm(p * k, sd = 1)
Gamma = matrix(specific_loadings, ncol = k, nrow = p)
phi_ = phi = matrix(NA, ncol = k, nrow = n)
for (h in 1:k) {
  phi[, h] = phi_[, h] = rnorm(n)
}
for (i in 1:n){
for (h in 1:k) {
  phi[i, h] = phi_[i, h] * (group[i] == h) # set 0s blocks of different groups
}
}
```

## Set prior hyperparameters
 

```{r}
alpha_eta = 5 # number of active shared factors a priori
v_eta = c(rbeta(d - 1, shape1 = 1, shape2 = alpha_eta), 1)
w_eta = v_eta * c(1, cumprod(1 - v_eta[-d]))                     # weights
z_eta = rep(d, d)
alpha_phi = 5 # number of active specific factors a priori
v_phi = c(rbeta(k - 1, shape1 = 1, shape2 = alpha_phi), 1)
w_phi = v_phi * c(1, cumprod(1 - v_phi[-k]))                     # weights
z_phi = rep(k, k)
whichgroup = unique(sapply(1:n , function (k)
  (which(X[k, ] == 1))))
betas = matrix(0, ncol = k, nrow = S)
```

## Define the state of the MCMC chain

```{r}
state = list(
  y = Y,
  # response
  Lambda = Lambda + rnorm(prod(dim(Lambda))),
  Lambda_ = Lambda + rnorm(prod(dim(Lambda)), sd =0.1),
  # shared loadings
  eta = (eta),
  # shared factors
  Gamma = Gamma + rnorm(prod(dim(Gamma)), sd = 0.1),
  #specific loadings
  phi = (phi) ,
  phi_ = (phi_) + rnorm(prod(dim(phi_)), sd = 0.1),
  #sparse and non sparse specific factors
  # list of specific covariance matrices
  n = n,
  ns = ns,
  X = X,
  S = S,
  d = d,
  k = k,
  p = p,
  #prior
  a_lambda = 1,
  b_lambda = 2,
  a_gamma = 1,
  b_gamma = 2,
  tau_eta = c(c(rep(1, 10)), c(rep(0, d - 10))),
  tau_phi = c(c(rep(1, 10)), c(rep(0, k - 10))),
  z_eta = z_eta,
  z_phi = z_phi,
  w_eta = w_eta,
  w_phi = w_phi,
  v_eta = v_eta,
  v_phi = v_phi,
  alpha_eta = 5,
  alpha_phi = 5,
  # equal to expected n.of active factors
  betas = betas,
  ps = matrix(rbinom(n * k, 1, 0.1), ncol = k)
)
copy_state = state
```

## Run the Gibbs sampler

In the following we run the Gibbs sampler. Note that `run` is set to `FALSE` to avoid running the Gibbs sampler (approximately 3 or 4 hours depending on your computer). To re-run the analysis change `run = FALSE` into `run = TRUE`.  

```{r}
run = FALSE
Gaussian = F
maxiter = 1
ris_phi = array(dim = c(maxiter, dim(state$phi)))
ris_eta = array(dim = c(maxiter, dim(state$eta)))
ris_th = array(dim = c(maxiter, dim(state$th)))
ris_ps = array(dim = c(maxiter, dim(state$ps)))

ris_beta = array(dim = c(maxiter, dim(state$betas)))
ris_eta = array(dim = c(maxiter, dim(state$eta)))
ris_lambda = array(dim = c(maxiter, dim(state$Lambda)))
ris_gamma = array(dim = c(maxiter, dim(state$Gamma)))
ris_tau_eta = array(dim = c(maxiter, length(state$tau_eta)))
ris_tau_phi = array(dim = c(maxiter, length(state$tau_phi)))

set.seed(111)
maxiter=10000
iter=1
if (run==TRUE){
for (iter in 1:maxiter) {
  cat(iter)
  state = Gibbs_Kernel_non_gauss(state)
  #shrinkage tau
  ris_tau_eta[iter, ] = state$tau_eta
  ris_tau_phi[iter, ] = state$tau_phi
  #factors
  ris_phi[iter, , ] = state$phi
  ris_eta[iter, , ] = state$eta
 
  #beta
  ris_beta[iter, , ] = state$betas

  #local activation psi
  ris_ps[iter, , ] = state$ps
 
  #loadings
  ris_lambda[iter, , ] = state$Lambda
  ris_gamma[iter, , ] = state$Gamma
 
  #print iteration and number of active factors
  if(iter%%100==0){
    print(state$tau_eta)
    print(state$tau_phi)
    }
}
}
if (run==TRUE) save.image("birds_results.RData")
```

# Load the results

```{r}
load("birds_results.RData")
```


# Plots

The first specific factor is associated to the type of habitat


``` r
par(mfrow=c(1,1))
#specific factors vs type of habitat
levels(da[,6])=c("Broadleaved","Conifer","Open","Urban","Wetlands")
for (i in 6:6) plot(colMeans(ris_phi[5000:10000,,1])~da[,i],
                    xlab=colnames(da)[i],ylab=expression(varphi[1]))
```

![](birdsmd_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->


We show the posterior mean of the factors, ordering units by latitude


```{r}
# Define the color palette with red and blue
custom_colors <- c("red","white", "blue")

# Create a custom color function that maps values to colors
color_function <- colorRampPalette(custom_colors)
colors <- color_function(0.3*p^2)
  
image(t(colMeans((ris_phi[5000:10000,,1:5])*
                    ris_ps[5000:10000,,1:5])),
ylab="latitude", col =colors, axes=F , xlab = expression(Phi))
axis(1,0.25*c(1:5)-0.25, 1:5)
```
 


![](birdsmd_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Posterior mean ordering by latitude and habitat

``` r
oor=order(da$Habitat)
cuts=table(da$Habitat[oor])/sum(table(da$Habitat[oor]))
image(t(colMeans((ris_phi[5000:10000,oor,1:5])*
                   ris_ps[5000:10000,oor,1:5])),
      ylab="environment", col =colors, axes=F , xlab = expression(Phi))
abline(h=cumsum(cuts)*1.01-0., lwd=2.985)
axis(1,0.25*c(1:5)-0.25, 1:5)
text(.5,cumsum(cuts)-0.012*c(4,8,2,4,1),c("Broadleaved", "Conifer", "Open", "Urban","Wetlands"))
```

![](birdsmd_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

We also display  lower and upper  bounds of credible intervals (10% and 90%) for specific factors ordering, again, by latitude. Lower bounds (10%)

```{r}
#function to compute 10% quantiles
colQuant10=function (m)apply(m,c(2,3),function(x) quantile(x,0.1))
image(t(colQuant10((ris_phi[5000:10000,,1:5])*
                   ris_ps[5000:10000,,1:5]))
      ,  main="q=0.1",
      ylab="latitude", col =colors, axes=F , xlab = expression(Phi))
axis(1,0.25*c(1:5)-0.25, 1:5)
```
 
![](birdsmd_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->


Upper bounds (90%) 

```{r}
#function to compute 90% quantiles
colQuant90=function (m) apply(m,c(2,3), function(x) quantile(x,0.9))
image(t(colQuant90((ris_phi[5000:10000,,1:5])*
                     ris_ps[5000:10000,,1:5])),
      main="q=0.9",
      ylab="latitude", col =colors, axes=F , xlab = expression(Phi))
axis(1,0.25*c(1:5)-0.25, 1:5)
```

![](birdsmd_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

We produce the same plots, first for the 10% quantiles,  ordering the specific factors by habitat

```{r}
image(t(colQuant10((ris_phi[5000:10000,oor,1:5])*
                     ris_ps[5000:10000,oor,1:5]))
      ,  main="q=0.1",
      ylab="environment", col =colors, axes=F , xlab = expression(Phi))
abline(h=cumsum(cuts)*1.01-0., lwd=2.985)
axis(1,0.25*c(1:5)-0.25, 1:5)
text(.5,cumsum(cuts)-0.012*c(4,8,2,4,1),c("Broadleaved", "Conifer", "Open", "Urban","Wetlands"))
```
 

![](birdsmd_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

then, for the 90% quantiles,  ordering the specific factors by habitat.

```{r}
image(t(colQuant90((ris_phi[5000:10000,oor,1:5])*
                     ris_ps[5000:10000,oor,1:5])),
      main="q=0.9",
      ylab="environment", col =colors, axes=F , xlab = expression(Phi))
abline(h=cumsum(cuts)*1.01-0., lwd=2.985)
axis(1,0.25*c(1:5)-0.25, 1:5)
text(.5,cumsum(cuts)-0.012*c(4,8,2,4,1),c("Broadleaved", "Conifer", "Open", "Urban","Wetlands"))
```
![](birdsmd_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->


