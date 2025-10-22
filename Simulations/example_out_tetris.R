# ===============================
# Simulation Analysis for TETRIS
# ===============================

# Clear workspace
rm(list = ls())

# Set working directory
setwd("Simulations/large")

# Initialize parameters
plot <- 0
NFACTOR_ <- matrix(NA, nrow = 10, ncol = 2)  # Store number of factors (shared / specific)
NFACTOR <- matrix(NA, nrow = 10, ncol = 2)
p <- 6  # Assuming 'p' is number of variables/features

# -------------------------------
# 1. Compute average number of factors for each setting
# -------------------------------
for (setting in c(1, 2, 4:10)) {
  cat("Processing setting:", setting, "\n")
  
  # Load simulation results
  name <- paste0("scen2TETRIS2_", setting, ".RData")
  load(name)
  
  # Compute number of shared (==3) and total factors
  NFACTOR_[setting, 1] <- mean(unlist(lapply(1:2000, function(x) sum(colSums(run$A[[x]]) == 3))))
  NFACTOR_[setting, 2] <- mean(unlist(lapply(1:2000, function(x) ncol(run$A[[x]])))) - NFACTOR_[setting, 1]
}

NFACTOR_

# -------------------------------
# 2. Compute estimated covariance matrices and performance measures
# -------------------------------
for (setting in c(1, 2, 4:10)) {
  cat("Processing covariance estimation for setting:", setting, "\n")
  
  # Load simulation results
  name <- paste0("scen2TETRIS2_", setting, ".RData")
  load(name)
  
  out <- run
  
  # Helper function: maximum a posteriori (MAP) for discrete or continuous
  MAP <- function(val) {
    if (length(unique(val)) < 50) {
      t <- table(val)
      as.numeric(names(t))[which.max(t)]
    } else {
      dd <- density(val)
      dd$x[which.max(dd$y)]
    }
  }
  
  # -------------------------------
  # 2a. Factor estimation
  # -------------------------------
  A <- choose.A(run, alpha_IBP = 6, S = 3)
  run_fixed <- tetris(data_y, alpha = 6, beta = 2, fixed = TRUE, A_fixed = A)
  Lambda <- getLambda(run_fixed, A)
  
  # -------------------------------
  # 2b. Compute Psi averages for each component
  # -------------------------------
  Sig_list <- lapply(out$Psi, function(sig) {
    sig_mat <- matrix(unlist(sig), ncol = p, byrow = FALSE)
    diag(apply(sig_mat, 2, mean))
  })
  
  Sig1 <- Sig_list[[1]]
  Sig2 <- Sig_list[[2]]
  Sig3 <- Sig_list[[3]]
  
  # -------------------------------
  # 2c. Function to compute trace
  # -------------------------------
  tr <- function(A) sum(diag(A))
  
  # -------------------------------
  # 2d. Compute estimated covariance matrices for each latent factor group
  # -------------------------------
  compute_LLGG <- function(A_list, Lambda_list, group_idx) {
    temp <- lapply(A_list, function(x) which(x[group_idx, ] == 1))
    temp <- lapply(seq_along(temp), function(x) Lambda_list[[x]][, temp[[x]]] %*% t(Lambda_list[[x]][, temp[[x]]]))
    
    # Average across simulations
    est <- matrix(0, ncol = p, nrow = p)
    for (j in 1:p) {
      for (h in 1:p) {
        est[j, h] <- mean(sapply(1:length(temp), function(x) temp[[x]][j, h]))
      }
    }
    return(est)
  }
  
  LLGG1est <- compute_LLGG(out$A, out$Lambda, 1)
  LLGG2est <- compute_LLGG(out$A, out$Lambda, 2)
  LLGG3est <- compute_LLGG(out$A, out$Lambda, 3)
  
  # -------------------------------
  # 2e. True covariance matrices from copy_state
  # -------------------------------
  LL <- copy_state$Lambda %*% diag(c(1, 1, 1, 0)) %*% t(copy_state$Lambda)
  
  # Group matrices
  make_group_matrix <- function(gr_idx) {
    dd <- diag(p) * 0
    dd[gr_idx, gr_idx] <- 1
    copy_state$Gamma %*% dd %*% t(copy_state$Gamma)
  }
  
  GG1 <- make_group_matrix(1)
  GG2 <- make_group_matrix(2)
  GG3 <- make_group_matrix(3)
  
  true1 <- LL + GG1 + copy_state$Sigma
  true2 <- LL + GG2 + copy_state$Sigma
  true3 <- LL + GG3 + copy_state$Sigma
  
  # -------------------------------
  # 2f. Performance measures
  # -------------------------------
  norm_eta <- tr(tcrossprod(copy_state$Lambda) %*% tcrossprod(Lambda)) /
    sqrt(tr(tcrossprod(copy_state$Lambda) %*% tcrossprod(copy_state$Lambda)) *
           tr(tcrossprod(Lambda) %*% tcrossprod(Lambda)))
  
  norm1 <- tr(t(LLGG1est + Sig1) %*% true1) / sqrt(tr(t(LLGG1est + Sig1) %*% (LLGG1est + Sig1)) *
                                                   tr(t(true1) %*% true1))
  
  norm2 <- tr(t(LLGG2est + Sig2) %*% true2) / sqrt(tr(t(LLGG2est + Sig2) %*% (LLGG2est + Sig2)) *
                                                   tr(t(true2) %*% true2))
  
  norm3 <- tr(t(LLGG3est + Sig3) %*% true3) / sqrt(tr(t(LLGG3est + Sig3) %*% (LLGG3est + Sig3)) *
                                                   tr(t(true3) %*% true3))
  
  # Store results
  TRACE <- cbind(matrix(c(mean(c(norm1, norm2, norm3)), mean(norm_eta)), nrow = 1))
  NFACTOR_TRACE <- cbind(NFACTOR_[setting, 1], NFACTOR_[setting, 2], TRACE)
  
  # Save results
  nameNFACTOR_TRACE <- paste0("scen2TETRIS2_", setting, "NFACTOR_TRACE.RDS")
  saveRDS(object = NFACTOR_TRACE, file = nameNFACTOR_TRACE)
}

# -------------------------------
# 3. Function to compute median and IQR
# -------------------------------
iqr_ <- function(x) {
  rbind(quantile(x, c(0.5, 0.25, 0.75)))
}

# -------------------------------
# 4. Average results across 10 datasets
# -------------------------------
setting <- 1
FT <- readRDS(paste0("scen1TETRIS1_", setting, "NFACTOR_TRACE.RDS"))
FT_ <- FT

for (setting in 2:10) {
  tmp <- readRDS(paste0("scen1TETRIS1_", setting, "NFACTOR_TRACE.RDS"))
  FT <- FT + tmp
  FT_ <- rbind(FT_, tmp)
}

FT <- FT / 10  # Average over datasets
FT_iqr <- apply(FT_[, -1], 2, iqr_)

# -------------------------------
# 5. Output tables
# -------------------------------
library(xtable)
xtable(FT)
xtable(FT_iqr[, -c(4:6)])
xtable(apply(NFACTOR_, 2, iqr_))

# -------------------------------
# 6. Plot covariance image (example)
# -------------------------------
par(mfrow = c(1, 1))
image(mat, xlab = "", ylab = "", oldstyle = FALSE)
