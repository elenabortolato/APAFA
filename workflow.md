## Example of usage: simulated data

Load the file sampler.R that contains the main Gibbs sampler method and
all library dependencies

    # Define the path to the sampler script
    #sampler_file <- "sampler.R"

    # Check if the file exists before sourcing
    #if (file.exists(sampler_file)) {
      # Load the sampler that contains the main Gibbs sampler method
    #  source(sampler_file)
    #  message("Sampler successfully loaded.")
    #} else {
    #  stop(paste("File not found:", sampler_file))
    #}
      
Generate synthetic data

    #-------------------------------------------------------
    # Example data generation for APAFA (shared + specific factors)
    #-------------------------------------------------------

    # --- Reproducibility ---
    set.seed(12345)

    # --- Study and data setup ---
    S  <- 3              # Number of studies (groups)
    ns <- rep(20, S)     # Sample size per study
    n  <- sum(ns)        # Total number of observations
    p  <- 10             # Number of observed variables (dimensions)

    # --- Residual covariance matrix ---
    Sigma <- diag(p) * 0.1                # Idiosyncratic (error) variance
    eps   <- mvtnorm::rmvnorm(n, sigma=Sigma)  # Random noise (n x p matrix)

    # --- Shared latent factors ---
    d <- 6                                # Number of shared latent factors
    shared_loadings <- rnorm(p * d, sd=1) # Random loadings for shared structure
    Lambda <- matrix(shared_loadings, nrow=p, ncol=d)
    eta    <- matrix(NA, nrow=n, ncol=d)  # Shared latent factor matrix (n x d)

    # Generate shared factors (standard normal)
    for (h in 1:d) {
      eta[, h] <- rnorm(n, mean=1)
    }

    # --- Study-specific latent factors ---
    ks <- rep(2, S)                       # Two specific factors per study
    k  <- sum(ks)                         # Total number of specific factors across studies

    # Group membership index for each observation
    group <- rep(NA, n)
    for (s in 1:S) {
      nscumpre <- ifelse(s > 1, sum(ns[1:(s-1)]) + 1, 1)
      nscum    <- sum(ns[1:s])
      group[nscumpre:nscum] <- s
    }

    # --- Dummy variables for group indicators ---
    X <- model.matrix(rep(1, n) ~ -1 + as.factor(group))
    dim(X)

    ## [1] 60  3

    # --- Specific factor loadings ---
    specific_loadings <- rnorm(p * k, sd=1)
    Gamma <- matrix(specific_loadings, nrow=p, ncol=k)

    # --- Initialize specific factors ---
    phi_ <- phi <- matrix(NA, nrow=n, ncol=k)
    for (h in 1:k) {
      phi[, h]  <- phi_[, h] <- rnorm(n, mean=1)
    }

    # --- Sparsify specific factors ---
    # Keep factors active only for their corresponding study
    for (s in 1:k) {
      phi[, s] <- phi[, s] * (group == s)
    }
    # Set inactive (higher-indexed) specific factors to zero
    phi[, 4:k] <- 0

    # --- Scale specific factors and loadings ---
    for (h in 1:k) {
      sdh <- sd(phi_[, h])
      Gamma[, h] <- Gamma[, h] * sdh     # Adjust loadings to match factor scale
      phi_[, h]  <- phi_[, h] / sdh
      phi[, h]   <- phi[, h] / sdh
    }

    # --- Scale shared factors and loadings ---
    for (h in 1:d) {
      sdh <- sd(eta[, h])
      Lambda[, h] <- Lambda[, h] * sdh
      eta[, h] <- eta[, h] / sdh
    }

    # --- Enforce sparsity on shared structure (optional) ---
    Lambda_ <- Lambda
    Lambda[, 4:d] <- 0  # Only first 3 shared factors active

    # --- Inspect structure ---
    head(Lambda)

    ##            [,1]        [,2]        [,3] [,4] [,5] [,6]
    ## [1,] -1.8398113 -0.48559055  1.29003918    0    0    0
    ## [2,]  0.2378274  0.24249345  0.18294380    0    0    0
    ## [3,] -0.5225391 -1.47957933 -0.04995755    0    0    0
    ## [4,] -0.7445867  0.97934627 -0.56959783    0    0    0
    ## [5,] -0.1494594  0.07015995  0.63896298    0    0    0
    ## [6,] -1.4857867  0.85071198 -0.05741765    0    0    0

    head(phi)

    ##            [,1] [,2] [,3] [,4] [,5] [,6]
    ## [1,] 0.15995409    0    0    0    0    0
    ## [2,] 1.36763630    0    0    0    0    0
    ## [3,] 0.06742471    0    0    0    0    0
    ## [4,] 2.38846642    0    0    0    0    0
    ## [5,] 0.12766600    0    0    0    0    0
    ## [6,] 0.58229441    0    0    0    0    0

    # --- Generate observed data ---
    # y_i = Λ * η_i + Γ * φ_i + ε_i
    y <- t(sapply(1:n, function(i)
      Lambda %*% eta[i, ] + Gamma %*% phi[i, ])) + eps

    # --- Hyperparameters for stick-breaking priors ---
    # These control the expected number of active factors

    # Shared factors
    alpha_eta <- 3
    v_eta <- c(rbeta(d - 1, shape1=1, shape2=alpha_eta), 1)
    w_eta <- v_eta * c(1, cumprod(1 - v_eta[-d]))   # Stick-breaking weights
    z_eta <- rep(d, d)

    # Study-specific factors
    alpha_phi <- 3
    v_phi <- c(rbeta(k - 1, shape1=1, shape2=alpha_phi), 1)
    w_phi <- v_phi * c(1, cumprod(1 - v_phi[-k]))
    z_phi <- rep(k, k)

    # --- Regression parameters for activation (logistic link) ---
    betas <- matrix(0, nrow=S, ncol=k)
    plogis(betas)  # Activation probabilities

    ##      [,1] [,2] [,3] [,4] [,5] [,6]
    ## [1,]  0.5  0.5  0.5  0.5  0.5  0.5
    ## [2,]  0.5  0.5  0.5  0.5  0.5  0.5
    ## [3,]  0.5  0.5  0.5  0.5  0.5  0.5

    #-------------------------------------------------------
    # Initialize state list for MCMC kernel
    #-------------------------------------------------------

    state <- list(
      # -----------------------------
      # Observed data
      # -----------------------------
      y       = y,         # Observed responses (n x p matrix)
      
      # -----------------------------
      # Shared latent factor structure
      # -----------------------------
      Lambda  = Lambda,    # Shared factor loadings (sparse, active factors)
      Lambda_ = Lambda_,   # Non-sparse copy of Lambda for reference
      eta     = eta,       # Shared latent factors (n x d)

      # -----------------------------
      # Study-specific latent factor structure
      # -----------------------------
      Gamma   = Gamma,     # Specific factor loadings (p x k)
      phi     = phi,       # Sparse specific factors (n x k)
      phi_    = phi_,      # Non-sparse specific factors for reference

      # -----------------------------
      # Sparsity / activation indicators
      # -----------------------------
      tau_eta = c(rep(1, 6), rep(0, d - 6)),   # Indicator for active shared factors
      tau_phi = c(rep(1, 6), rep(0, k - 6)),   # Indicator for active specific factors
      z_eta   = z_eta,                          # Number of active shared factors (stick-breaking)
      z_phi   = z_phi,                          # Number of active specific factors
      w_eta   = w_eta,                          # Stick-breaking weights (shared)
      w_phi   = w_phi,                          # Stick-breaking weights (specific)
      v_eta   = v_eta,                          # Beta variables for stick-breaking (shared)
      v_phi   = v_phi,                          # Beta variables for stick-breaking (specific)

      # -----------------------------
      # Error / covariance structure
      # -----------------------------
      Sigma   = Sigma,     # Residual covariance matrix (p x p)

      # -----------------------------
      # Dimensions and study info
      # -----------------------------
      n       = n,         # Total number of observations
      ns      = ns,        # Number of observations per study
      X       = X,         # Dummy variable matrix for group indicators (n x S)
      S       = S,         # Number of studies
      d       = d,         # Number of shared latent factors
      k       = k,         # Number of specific latent factors
      p       = p,         # Number of observed variables

      # -----------------------------
      # Prior hyperparameters
      # -----------------------------
      a_lambda  = 1,       # Shape for Gamma prior on shared loadings
      b_lambda  = 2,       # Rate for Gamma prior on shared loadings
      a_gamma   = 1,       # Shape for Gamma prior on specific loadings
      b_gamma   = 2,       # Rate for Gamma prior on specific loadings
      a_sigma   = 20,      # Shape for Inverse-Gamma prior on residual variance
      b_sigma   = 2,       # Rate for Inverse-Gamma prior
      alpha_eta = 3,       # Stick-breaking hyperparameter for shared factors
      alpha_phi = 3,       # Stick-breaking hyperparameter for specific factors

      # -----------------------------
      # Factor activation parameters
      # -----------------------------
      betas = betas,       # Logistic regression parameters for factor activation
      ps    = matrix(rbinom(n * k, 1, 0.50), nrow=n) # Binary indicators of factor activity per observation
    )

    copy_state = state
    #save the state before initialization (containing the true parameters)
    #saveRDS(copy_state, paste("copy_state", setting))

Now we add Uniform\[-amount, amount\] noise to initialize the model.

    state$Lambda_ = jitter(state$Lambda_, amount = 1)
    state$Lambda = state$Lambda_
    state$Gamma = jitter(state$Gamma, amount = 1)
    state$eta = jitter(state$eta, amount = 1)
    state$phi_ = jitter(state$phi_, amount = 1)
    state$phi = state$phi_
    state$Sigma = jitter(state$Sigma, amount = 0.1)

Define the number of Gibbs sampler iterations and pre - allocate
matrices to store results of MCMC

    maxiter = 5000

    ris_sigma1 = matrix(0, ncol = p, maxiter) # will contain the entries of the diagonal matrix Sigma

    ris_beta1 = matrix(0, ncol = S, maxiter) # k vectors of coefficient to model the activation of specific factors among groups
    ris_beta2 = matrix(0, ncol = S, maxiter)
    ris_beta3 = matrix(0, ncol = S, maxiter)
    ris_beta4 = matrix(0, ncol = S, maxiter)
    ris_beta5 = matrix(0, ncol = S, maxiter)
    ris_beta6 = matrix(0, ncol = S, maxiter)

    ris_tau_eta = matrix(0, ncol = d, maxiter) # d activation indicators
    ris_tau_phi = matrix(0, ncol = k, maxiter) # k activation indicators

    #specific factors
    ris_phi1 = ris_phi2 = ris_phi3 = ris_phi4 = ris_phi5 = ris_phi6 =
    ris_ps1 = ris_ps2 = ris_ps3 = ris_ps4 = ris_ps5 = ris_ps6 =
    matrix(0, ncol = n, maxiter)
    #common factors
    ris_eta1 = ris_eta2 = ris_eta3 = ris_eta4 = ris_eta5 = ris_eta6 = matrix(0, ncol =
    n, maxiter)

    #loadings
    ris_lambda1 = ris_lambda2 = ris_lambda3 = ris_lambda4 = ris_lambda5 =
    ris_lambda6 = matrix(0, ncol = p, maxiter)
    ris_gamma1 = ris_gamma2 = ris_gamma3 = ris_gamma4 = ris_gamma5 = ris_gamma6 =
    matrix(0, ncol = p, maxiter)

## Run the sampler

    # Gibbs sampler
    time <- system.time(for (iter in 1:maxiter) {
    state = Gibbs_Kernel(state)
    # if(iter==1) image(state$ps, main=0) # plot activation pattern of local factors (psi)
    ris_tau_eta[iter, ] = state$tau_eta
    ris_tau_phi[iter, ] = state$tau_phi
    #factors
    ris_phi1[iter, ] = state$phi[, 1]
    ris_phi2[iter, ] = state$phi[, 2]
    ris_phi3[iter, ] = state$phi[, 3]
    ris_phi4[iter, ] = state$phi[, 4]
    ris_phi5[iter, ] = state$phi[, 5]
    ris_phi6[iter, ] = state$phi[, 6]
    ris_eta1[iter, ] = state$eta[, 1]
    ris_eta2[iter, ] = state$eta[, 2]
    ris_eta3[iter, ] = state$eta[, 3]
    ris_eta4[iter, ] = state$eta[, 4]
    ris_eta5[iter, ] = state$eta[, 5]
    ris_eta6[iter, ] = state$eta[, 6]
    #beta
    ris_beta1[iter, ] = state$betas[, 1]
    ris_beta2[iter, ] = state$betas[, 2]
    ris_beta3[iter, ] = state$betas[, 3]
    ris_beta4[iter, ] = state$betas[, 4]
    ris_beta5[iter, ] = state$betas[, 5]
    ris_beta6[iter, ] = state$betas[, 6]

    #theta psi
    ris_ps1[iter, ] = state$ps[, 1]
    ris_ps2[iter, ] = state$ps[, 2]
    ris_ps3[iter, ] = state$ps[, 3]
    ris_ps4[iter, ] = state$ps[, 4]
    ris_ps5[iter, ] = state$ps[, 5]
    ris_ps6[iter, ] = state$ps[, 6]
    #loading
    ris_lambda1[iter, ] = state$Lambda[, 1]
    ris_lambda2[iter, ] = state$Lambda[, 2]
    ris_lambda3[iter, ] = state$Lambda[, 3]
    ris_lambda4[iter, ] = state$Lambda[, 4]
    ris_lambda5[iter, ] = state$Lambda[, 5]
    ris_lambda6[iter, ] = state$Lambda[, 6]
    ris_gamma1[iter, ] = state$Gamma[, 1]
    ris_gamma2[iter, ] = state$Gamma[, 2]
    ris_gamma3[iter, ] = state$Gamma[, 3]
    ris_gamma4[iter, ] = state$Gamma[, 4]
    ris_gamma5[iter, ] = state$Gamma[, 5]
    ris_gamma6[iter, ] = state$Gamma[, 6]
    ris_sigma1[iter, ] = diag(state$Sigma)
    #print iteration, number of active factors, average idyosincratic error variance
    if (iter %% 100 == 0) {
    cat(iter)
    cat("\n")
    #shrinkage tau
    cat("tau[eta]")
    print(state$tau_eta)
    cat("tau[phi]")
    print(state$tau_phi)
    cat("sigma^2")
    print(mean(diag(state$Sigma)))
    }
    })

    ## 100
    ## tau[eta][1] 1 1 1 1 1 0
    ## tau[phi][1] 1 1 1 1 1 0
    ## sigma^2[1] 0.1028437
    ## 200
    ## tau[eta][1] 1 1 1 1 1 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.0880105
    ....
    ## 4900
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.09978126
    ## 5000
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1053641

    time

    ##    user  system elapsed 
    ## 113.556   2.846 119.173

## Plots: convergence assessments and uncertainty quantification based on posterior

    # convergence assessment: Compute running loglikelihood
    iter = 1
    Loglik = rep(0, maxiter)
    for (iter in 1:maxiter) {
      L = cbind(
        ris_lambda1[iter, ],
        ris_lambda2[iter, ],
        ris_lambda3[iter, ],
        ris_lambda4[iter, ],
        ris_lambda5[iter, ],
        ris_lambda6[iter, ]
      )
      L = t(t(L) * ris_tau_eta[iter, 1:6])
      LLT = L %*% t(L)
      G = cbind(ris_gamma1[iter, ],
                ris_gamma2[iter, ],
                ris_gamma3[iter, ],
                ris_gamma4[iter, ],
                ris_gamma5[iter, ],
                ris_gamma6[iter, ])
      G = t(t(G) * ris_tau_phi[iter, 1:6])
      PS = (cbind(ris_ps1[iter, ], ris_ps2[iter, ], ris_ps3[iter, ], ris_ps4[iter, ], ris_ps5[iter, ], ris_ps6[iter, ]))
      GGT = lapply(1:n, function (ii)
        G %*% diag(PS[ii, 1:6]) %*%
          t(G))
      S = diag(ris_sigma1[iter, ])
      Loglik[iter] = 0
      
      for (ii in 1:state$n) {
        VAR = LLT + GGT[[ii]] + S
        Loglik[iter] = Loglik[iter] + dmvnorm(y[ii, ], rep(0, p), sigma = VAR, log =
                                                T)
      }
      #if (iter%%100==0) print(iter)
    }

    
    

    # Plot the loglikelihood
    par(mar = c(2, 2, 0.2, 0.2) + 2)
    par(mfrow = c(2, 1))
    plot(Loglik,pch="|",
         cex = 1,
         xlab = "iteration",main="LLik and number of shared factors", 
         ,
         col = rowSums(ris_tau_eta))
    legend(
      "bottomright",
      col = c(5, 4, 3),
      lty = 1,
      legend =  paste("d=", c(5:3))
    )
    plot(Loglik,
         cex = 1,pch="|",
         xlab = "iteration",main="LLik and number of specific factors",
         col = rowSums(ris_tau_phi))
    legend(
      "bottomright",
      col = c(5, 4, 3),
      lty = 1,
      legend =  paste("k=", c(5:3))
    )

![](workflow_files/figure-markdown_strict/unnamed-chunk-9-1.png)

Posterior distribution of the entries of the covariance matrix compared
with true values (red marks)

    ### uncertainty quantification
    iter0 = round(maxiter*0.8)
    len=abs(iter0-maxiter)
    VARR = VARRTRUE = matrix(0, len, p^2)
    iter=1
    for (iter in 1:len) {
    iter_ = iter + iter0
    L = cbind(
    ris_lambda1[iter+iter0, ],
    ris_lambda2[iter+iter0, ],
    ris_lambda3[iter+iter0, ],
    ris_lambda4[iter+iter0, ],
    ris_lambda5[iter+iter0, ],
    ris_lambda6[iter+iter0, ]
    )
    L = t(t(L) * ris_tau_eta[iter_, 1:d])
    LLT = L %*% t(L)

    S = diag(ris_sigma1[iter_, ])
    VAR = LLT + S
    VARTRUE = copy_state$Lambda %*% t(copy_state$Lambda) + copy_state$Sigma
    VARR[iter, ] =  c(VAR)
    VARRTRUE[iter, ] =  (VARTRUE)

    print(iter)
    }

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    # ....
    ## [1] 999
    ## [1] 1000

    or = order(colMeans(VARR))
    par(mar = c(4, 4, 4, 4))

    boxplot(VARR[, or] ,
    outline = T,
    ylab = expression(lambda[ij] + sigma[ii]))
    points(1:(p^2), VARTRUE[or], col = 2, cex = 1.5, pch="x")

![](workflow_files/figure-markdown_strict/unnamed-chunk-10-1.png)
