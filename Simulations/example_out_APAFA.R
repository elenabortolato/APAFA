# =========================================================
# Example file to compute performance metrics for the
# simulation studies (Section 3 of the paper)
# =========================================================

# -------------------------------
# 0. Setup
# -------------------------------

# Option to enable plotting
plot <- 0
par(mfrow = c(1,1))  # reset plotting area

# Initialize matrix to store activation patterns across simulations
partition <- matrix(0, ncol = 3, nrow = 45)  # 45 variables, 3 specific factors

# -------------------------------
# 1. Loop over all simulation settings
# -------------------------------
for(setting in 1:10){
  
  # Load the simulation results for the current setting
  name <- paste0("scen1_", setting, ".RData")
  load(name)
  
  # -------------------------------
  # 1a. Define a helper function to compute MAP (maximum a posteriori)
  # For discrete variables: mode
  # For continuous variables: density peak
  # -------------------------------
  MAP <- function(val){
    if(length(unique(val)) < 50){
      t <- table(val)
      as.numeric(names(t))[which.max(t)]
    } else {
      dd <- density(val)
      dd$x[which.max(dd$y)]
    }
  }
  
  # -------------------------------
  # 1b. Compute mean number of active factors
  # ris_tau_eta: latent factor indicators for shared factors
  # ris_tau_phi: latent factor indicators for specific factors
  # -------------------------------
  NFACTOR <- cbind(
    mean(apply(ris_tau_eta, 1, sum)),
    mean(apply(ris_tau_phi, 1, sum))
  )
  
  # -------------------------------
  # 1c. Define helper functions
  # -------------------------------
  # Function to compute trace of a matrix
  tr <- function(A) sum(diag(A))
  
  # Replace NaN values if necessary
  nanzero <- function(a, b = 0){
    if(is.nan(b)) return(a) else return(b)
  }
  
  burnin <- 8000  # burn-in for MCMC samples
  
  # Compute posterior mean of latent factor activity probabilities
  tau_ph <- apply(ris_tau_phi, 2, mean)
  tau_et <- apply(ris_tau_eta, 2, mean)
  
  # Boolean indicator matrix of active specific factors
  copy_state$ps <- copy_state$phi != 0
  
  # -------------------------------
  # 1d. Align columns to avoid label switching
  # -------------------------------
  ps1 <- apply(ris_ps1[-c(1:burnin), ], 2, mean)
  ps2 <- apply(ris_ps2[-c(1:burnin), ], 2, mean)
  ps3 <- apply(ris_ps3[-c(1:burnin), ], 2, mean)
  ps4 <- apply(ris_ps4[-c(1:burnin), ], 2, mean)
  psest <- cbind(ps1, ps2, ps3, ps4)
  
  # Compute distance between estimated and true specific factors
  tab <- sapply(1:length(copy_state$tau_phi), function(x) {
    sapply(1:ncol(psest), function(k) sum(abs(copy_state$ps[,x] - psest[,k])))
  })
  
  # Avoid matching with inactive factors
  tab[, copy_state$tau_phi == 0] <- copy_state$n * 2  # true factor = 0
  tab[tau_ph == 0, ] <- copy_state$n * 2  # estimated factor = 0
  
  # -------------------------------
  # 1e. Match estimated factors to true factors (best alignment)
  # -------------------------------
  ordine <- order(apply(tab, 2, min))
  taken <- NA
  similar_ps <- rep(NA, ncol(copy_state$ps))
  
  for(i in 1:ncol(tab)){
    found <- FALSE
    which <- 1
    while(!found){
      similar_ps[ordine[i]] <- order(tab[, ordine[i]])[which]
      if(similar_ps[[ordine[i]]] %in% taken){
        which <- which + 1
      } else {
        taken <- c(taken, similar_ps[ordine[i]])
        found <- TRUE
      }
    }
  }
  
  # Reorder estimated factors
  psest <- psest[, similar_ps]
  
  # -------------------------------
  # 1f. Compute posterior mean of Lambda and Gamma matrices
  # -------------------------------
  lambda1 <- apply(ris_lambda1[-c(1:burnin), ], 2, mean)
  lambda2 <- apply(ris_lambda2[-c(1:burnin), ], 2, mean)
  lambda3 <- apply(ris_lambda3[-c(1:burnin), ], 2, mean)
  
  gamma1 <- apply(ris_gamma1[-c(1:burnin), ], 2, mean)
  gamma2 <- apply(ris_gamma2[-c(1:burnin), ], 2, mean)
  gamma3 <- apply(ris_gamma3[-c(1:burnin), ], 2, mean)
  gamma4 <- apply(ris_gamma4[-c(1:burnin), ], 2, mean)
  
  Lam <- cbind(lambda1, lambda2, lambda3)
  Gam <- cbind(gamma1, gamma2, gamma3, gamma4)
  Gam <- Gam[, similar_ps]  # reorder columns to match true factors
  
  # -------------------------------
  # 1g. Posterior mean of residual covariances
  # -------------------------------
  Sig1 <- diag(apply(ris_sigma1[-c(1:burnin), ], 2, mean))
  Sig2 <- diag(apply(ris_sigma1[-c(1:burnin), ], 2, mean))
  Sig3 <- diag(apply(ris_sigma1[-c(1:burnin), ], 2, mean))
  
  SIG <- list(Sig1, Sig2, Sig3)
  
  # Reorder SIG according to factor matching
  SIGOR <- SIG
  for(k in 1:length(SIGOR)) SIGOR[[k]] <- SIG[[similar_ps[k]]]
  Sig1 <- SIGOR[[1]]; Sig2 <- SIGOR[[2]]; Sig3 <- SIGOR[[3]]
  
  # -------------------------------
  # 1h. Scale specific factors with tau_ph
  # -------------------------------
  ps <- psest
  ps[is.na(ps)] <- 0
  ps <- t(t(ps) * tau_ph)
  
  # -------------------------------
  # 1i. Compute true covariance components
  # -------------------------------
  LL <- copy_state$Lambda %*% t(copy_state$Lambda)
  dd <- diag(4) * 0
  gr1_ <- dd; gr2_ <- dd; gr3_ <- dd
  gr1_[1,1] <- 1; gr2_[2,2] <- 1; gr3_[3,3] <- 1
  
  GG1 <- copy_state$Gamma %*% gr1_ %*% t(copy_state$Gamma)
  GG2 <- copy_state$Gamma %*% gr2_ %*% t(copy_state$Gamma)
  GG3 <- copy_state$Gamma %*% gr3_ %*% t(copy_state$Gamma)
  
  LLGG1 <- LL + GG1
  LLGG2 <- LL + GG2
  LLGG3 <- LL + GG3
  
  true1 <- LLGG1 + copy_state$Sigmas[[1]]
  true2 <- LLGG2 + copy_state$Sigmas[[2]]
  true3 <- LLGG3 + copy_state$Sigmas[[3]]
  
  # -------------------------------
  # 1j. Compute performance metrics (RV coefficients)
  # -------------------------------
  gamphi <- lapply(1:n, function(i) Gam %*% diag(ps[i, ]) %*% t(Gam))
  lamet <- lapply(1:n, function(i) Lam %*% diag(th[i, ]) %*% t(Lam))
  
  # Overall similarity between estimated and true shared factors
  norm_eta <- mean(sapply(1:n, function(i) tr(t(lamet[[i]]) %*% LL) /
                             sqrt(tr(t(lamet[[i]]) %*% lamet[[i]]) * tr(t(LL) %*% LL))))
  
  # Similarity for each specific factor group
  norm1 <- mean(sapply(1:10, function(i) tr(t(lamet[[i]] + gamphi[[i]] + Sig1) %*% true1) /
                          sqrt(tr(t(lamet[[i]] + gamphi[[i]] + Sig1) %*% (lamet[[i]] + gamphi[[i]] + Sig1)) *
                               tr(t(true1) %*% true1))))
  
  norm2 <- mean(sapply(11:20, function(i) tr(t(lamet[[i]] + gamphi[[i]] + Sig2) %*% true2) /
                           sqrt(tr(t(lamet[[i]] + gamphi[[i]] + Sig2) %*% (lamet[[i]] + gamphi[[i]] + Sig2)) *
                                tr(t(true2) %*% true2))))
  
  norm3 <- mean(sapply(21:30, function(i) tr(t(lamet[[i]] + gamphi[[i]] + Sig3) %*% true3) /
                           sqrt(tr(t(lamet[[i]] + gamphi[[i]] + Sig3) %*% (lamet[[i]] + gamphi[[i]] + Sig3)) *
                                tr(t(true3) %*% true3))))
  
  # Similarity for gamma components only
  

  normGG1=mean(sapply(1:10,function(i) tr(gamphi[[i]]%*%GG1)/
                                                 sqrt(tr(t(gamphi[[i]])%*%gamphi[[i]])%*%tr(t(GG1)%*%GG1))), na.rm = T)
  normGG2=mean( sapply(11:20,function(i) tr(gamphi[[i]]%*%GG2)/
                                                 sqrt(tr(t(gamphi[[i]])%*%gamphi[[i]])%*%tr(t(GG2)%*%GG2))), na.rm = T)
  normGG3=mean(sapply(21:30,function(i)tr(gamphi[[i]]%*%GG3)/
                                                 sqrt(tr(t(gamphi[[i]])%*%gamphi[[i]])%*%tr(t(GG3)%*%GG3))), na.rm = T)

  
  # number of active factors
  NFACTOR=cbind(NFACTOR)
  # RV coefficient
  TRACE=cbind(rbind(apply((cbind(norm1, norm2 ,norm3 ,normGG1 , normGG2,normGG3)),2,mean)),
              rbind(apply((cbind(norm_eta)),2,mean)))
  
  # Prepare the ROC CURVE plot
  partition=partition+((ps[,copy_state$tau_phi==1]))
  auc=pROC::auc(c(c(copy_state$ps[,copy_state$tau_phi==1])!=0),c(ps[,copy_state$tau_phi==1]))
  err=1-mean(abs(c(c(copy_state$ps[,copy_state$tau_phi==1])!=0)-c(ps[,copy_state$tau_phi==1])))
  
  auc=c(auc,err) # AUC and relative error
  pROC::roc
  rocc=pROC::roc(c(c(copy_state$phi[,copy_state$tau_phi==1])!=0),c(ps[,copy_state$tau_phi==1]))
  
  #plot roc curve
  if(setting==1){
   plot(rocc)}
  else{
    plot(rocc,add=T)
  }
  #plot activation pattern of specific factors
  image((ps[,copy_state$tau_phi==1]))

  NFACTOR_TRACE=cbind(NFACTOR,TRACE)
  NFACTOR_TRACE 
  nameNFACTOR_TRACE=paste(paste(paste("scen1_",setting, sep = ""),"NFACTOR_TRACE",sep=""),".RDS", sep="")
  nameNFACTOR_TRACE
  saveRDS(object = NFACTOR_TRACE, nameNFACTOR_TRACE)  
  
  nameNMATRIX=paste(paste(paste("scen1_",setting, sep = ""),"MATRIX",sep=""),".RDS", sep="")
    saveRDS(object = auc, nameNMATRIX)  
}
partition
parition=partition/10
par(mfrow=c(1,1))
image(partition)
# table
iqr_<- function (x) return(rbind(quantile(x, c(0.5,0.25, 0.75))))



#### average results over 10 datasets: produce the metric performance
#### of table 1
setting=1
nameNFACTOR_TRACE=paste(paste(paste("scen1_",setting, sep = ""),"NFACTOR_TRACE",sep=""),".RDS", sep="")
FT=readRDS(nameNFACTOR_TRACE, nameNFACTOR_TRACE)
FT_=FT
nameNMATRIX=paste(paste(paste("scen1_",setting, sep = ""),"MATRIX",sep=""),".RDS", sep="")
 
AUC=cbind(rep(NA,10),rep(NA,10))
mat=readRDS(nameNMATRIX)
AUC[1,]=mat[1:2]
  
# prepare activation pattern figure 
for(setting in 2:10){
  nameNFACTOR_TRACE=paste(paste(paste("scen1_",setting, sep = ""),"NFACTOR_TRACE",sep=""),".RDS", sep="")
  nameNFACTOR_TRACE
  FT=FT+readRDS(nameNFACTOR_TRACE)
  FT_=rbind(FT_,readRDS(nameNFACTOR_TRACE))
  nameNMATRIX=paste(paste(paste("scen1_",setting, sep = ""),"MATRIX",sep=""),".RDS", sep="")
  mat=readRDS(nameNMATRIX)
  AUC[setting,]=mat[1:2]
}

FT=FT/10
FT_=apply(FT_, 2, iqr_)

apply(AUC,2,iqr_)
FT=FT
FT_=FT_[,-c(6:8)]

xtable::xtable(FT)
xtable::xtable(FT_)
xtable::xtable(mat)

par(mfrow=c(1,1))
image(mat, xlab = "", ylab = "",oldstyle = F)

