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
    ## 300
    ## tau[eta][1] 1 1 1 1 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.0931883
    ## 400
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.09666213
    ## 500
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.09272252
    ## 600
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.1061118
    ## 700
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.09270123
    ## 800
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.09364668
    ## 900
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.1037738
    ## 1000
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.1039276
    ## 1100
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.1034119
    ## 1200
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.09718746
    ## 1300
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.1040795
    ## 1400
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.1051697
    ## 1500
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 1 0 0
    ## sigma^2[1] 0.1018728
    ## 1600
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1069601
    ## 1700
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.09598138
    ## 1800
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1006346
    ## 1900
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1036444
    ## 2000
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.09346617
    ## 2100
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.102946
    ## 2200
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1125706
    ## 2300
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1028534
    ## 2400
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.09795223
    ## 2500
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1145823
    ## 2600
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1071655
    ## 2700
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.100833
    ## 2800
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1023072
    ## 2900
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.09668988
    ## 3000
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1000898
    ## 3100
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1057075
    ## 3200
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1038873
    ## 3300
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1121723
    ## 3400
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.115452
    ## 3500
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1063707
    ## 3600
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1062338
    ## 3700
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1088476
    ## 3800
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1078105
    ## 3900
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.09301477
    ## 4000
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1003133
    ## 4100
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.09958615
    ## 4200
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.099363
    ## 4300
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1180885
    ## 4400
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1026818
    ## 4500
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1000509
    ## 4600
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1018197
    ## 4700
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1036277
    ## 4800
    ## tau[eta][1] 1 1 1 0 0 0
    ## tau[phi][1] 1 1 1 0 0 0
    ## sigma^2[1] 0.1121025
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
      if (iter%%100==0) print(iter)
    }

    ## [1] 100
    ## [1] 200
    ## [1] 300
    ## [1] 400
    ## [1] 500
    ## [1] 600
    ## [1] 700
    ## [1] 800
    ## [1] 900
    ## [1] 1000
    ## [1] 1100
    ## [1] 1200
    ## [1] 1300
    ## [1] 1400
    ## [1] 1500
    ## [1] 1600
    ## [1] 1700
    ## [1] 1800
    ## [1] 1900
    ## [1] 2000
    ## [1] 2100
    ## [1] 2200
    ## [1] 2300
    ## [1] 2400
    ## [1] 2500
    ## [1] 2600
    ## [1] 2700
    ## [1] 2800
    ## [1] 2900
    ## [1] 3000
    ## [1] 3100
    ## [1] 3200
    ## [1] 3300
    ## [1] 3400
    ## [1] 3500
    ## [1] 3600
    ## [1] 3700
    ## [1] 3800
    ## [1] 3900
    ## [1] 4000
    ## [1] 4100
    ## [1] 4200
    ## [1] 4300
    ## [1] 4400
    ## [1] 4500
    ## [1] 4600
    ## [1] 4700
    ## [1] 4800
    ## [1] 4900
    ## [1] 5000

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
    ## [1] 9
    ## [1] 10
    ## [1] 11
    ## [1] 12
    ## [1] 13
    ## [1] 14
    ## [1] 15
    ## [1] 16
    ## [1] 17
    ## [1] 18
    ## [1] 19
    ## [1] 20
    ## [1] 21
    ## [1] 22
    ## [1] 23
    ## [1] 24
    ## [1] 25
    ## [1] 26
    ## [1] 27
    ## [1] 28
    ## [1] 29
    ## [1] 30
    ## [1] 31
    ## [1] 32
    ## [1] 33
    ## [1] 34
    ## [1] 35
    ## [1] 36
    ## [1] 37
    ## [1] 38
    ## [1] 39
    ## [1] 40
    ## [1] 41
    ## [1] 42
    ## [1] 43
    ## [1] 44
    ## [1] 45
    ## [1] 46
    ## [1] 47
    ## [1] 48
    ## [1] 49
    ## [1] 50
    ## [1] 51
    ## [1] 52
    ## [1] 53
    ## [1] 54
    ## [1] 55
    ## [1] 56
    ## [1] 57
    ## [1] 58
    ## [1] 59
    ## [1] 60
    ## [1] 61
    ## [1] 62
    ## [1] 63
    ## [1] 64
    ## [1] 65
    ## [1] 66
    ## [1] 67
    ## [1] 68
    ## [1] 69
    ## [1] 70
    ## [1] 71
    ## [1] 72
    ## [1] 73
    ## [1] 74
    ## [1] 75
    ## [1] 76
    ## [1] 77
    ## [1] 78
    ## [1] 79
    ## [1] 80
    ## [1] 81
    ## [1] 82
    ## [1] 83
    ## [1] 84
    ## [1] 85
    ## [1] 86
    ## [1] 87
    ## [1] 88
    ## [1] 89
    ## [1] 90
    ## [1] 91
    ## [1] 92
    ## [1] 93
    ## [1] 94
    ## [1] 95
    ## [1] 96
    ## [1] 97
    ## [1] 98
    ## [1] 99
    ## [1] 100
    ## [1] 101
    ## [1] 102
    ## [1] 103
    ## [1] 104
    ## [1] 105
    ## [1] 106
    ## [1] 107
    ## [1] 108
    ## [1] 109
    ## [1] 110
    ## [1] 111
    ## [1] 112
    ## [1] 113
    ## [1] 114
    ## [1] 115
    ## [1] 116
    ## [1] 117
    ## [1] 118
    ## [1] 119
    ## [1] 120
    ## [1] 121
    ## [1] 122
    ## [1] 123
    ## [1] 124
    ## [1] 125
    ## [1] 126
    ## [1] 127
    ## [1] 128
    ## [1] 129
    ## [1] 130
    ## [1] 131
    ## [1] 132
    ## [1] 133
    ## [1] 134
    ## [1] 135
    ## [1] 136
    ## [1] 137
    ## [1] 138
    ## [1] 139
    ## [1] 140
    ## [1] 141
    ## [1] 142
    ## [1] 143
    ## [1] 144
    ## [1] 145
    ## [1] 146
    ## [1] 147
    ## [1] 148
    ## [1] 149
    ## [1] 150
    ## [1] 151
    ## [1] 152
    ## [1] 153
    ## [1] 154
    ## [1] 155
    ## [1] 156
    ## [1] 157
    ## [1] 158
    ## [1] 159
    ## [1] 160
    ## [1] 161
    ## [1] 162
    ## [1] 163
    ## [1] 164
    ## [1] 165
    ## [1] 166
    ## [1] 167
    ## [1] 168
    ## [1] 169
    ## [1] 170
    ## [1] 171
    ## [1] 172
    ## [1] 173
    ## [1] 174
    ## [1] 175
    ## [1] 176
    ## [1] 177
    ## [1] 178
    ## [1] 179
    ## [1] 180
    ## [1] 181
    ## [1] 182
    ## [1] 183
    ## [1] 184
    ## [1] 185
    ## [1] 186
    ## [1] 187
    ## [1] 188
    ## [1] 189
    ## [1] 190
    ## [1] 191
    ## [1] 192
    ## [1] 193
    ## [1] 194
    ## [1] 195
    ## [1] 196
    ## [1] 197
    ## [1] 198
    ## [1] 199
    ## [1] 200
    ## [1] 201
    ## [1] 202
    ## [1] 203
    ## [1] 204
    ## [1] 205
    ## [1] 206
    ## [1] 207
    ## [1] 208
    ## [1] 209
    ## [1] 210
    ## [1] 211
    ## [1] 212
    ## [1] 213
    ## [1] 214
    ## [1] 215
    ## [1] 216
    ## [1] 217
    ## [1] 218
    ## [1] 219
    ## [1] 220
    ## [1] 221
    ## [1] 222
    ## [1] 223
    ## [1] 224
    ## [1] 225
    ## [1] 226
    ## [1] 227
    ## [1] 228
    ## [1] 229
    ## [1] 230
    ## [1] 231
    ## [1] 232
    ## [1] 233
    ## [1] 234
    ## [1] 235
    ## [1] 236
    ## [1] 237
    ## [1] 238
    ## [1] 239
    ## [1] 240
    ## [1] 241
    ## [1] 242
    ## [1] 243
    ## [1] 244
    ## [1] 245
    ## [1] 246
    ## [1] 247
    ## [1] 248
    ## [1] 249
    ## [1] 250
    ## [1] 251
    ## [1] 252
    ## [1] 253
    ## [1] 254
    ## [1] 255
    ## [1] 256
    ## [1] 257
    ## [1] 258
    ## [1] 259
    ## [1] 260
    ## [1] 261
    ## [1] 262
    ## [1] 263
    ## [1] 264
    ## [1] 265
    ## [1] 266
    ## [1] 267
    ## [1] 268
    ## [1] 269
    ## [1] 270
    ## [1] 271
    ## [1] 272
    ## [1] 273
    ## [1] 274
    ## [1] 275
    ## [1] 276
    ## [1] 277
    ## [1] 278
    ## [1] 279
    ## [1] 280
    ## [1] 281
    ## [1] 282
    ## [1] 283
    ## [1] 284
    ## [1] 285
    ## [1] 286
    ## [1] 287
    ## [1] 288
    ## [1] 289
    ## [1] 290
    ## [1] 291
    ## [1] 292
    ## [1] 293
    ## [1] 294
    ## [1] 295
    ## [1] 296
    ## [1] 297
    ## [1] 298
    ## [1] 299
    ## [1] 300
    ## [1] 301
    ## [1] 302
    ## [1] 303
    ## [1] 304
    ## [1] 305
    ## [1] 306
    ## [1] 307
    ## [1] 308
    ## [1] 309
    ## [1] 310
    ## [1] 311
    ## [1] 312
    ## [1] 313
    ## [1] 314
    ## [1] 315
    ## [1] 316
    ## [1] 317
    ## [1] 318
    ## [1] 319
    ## [1] 320
    ## [1] 321
    ## [1] 322
    ## [1] 323
    ## [1] 324
    ## [1] 325
    ## [1] 326
    ## [1] 327
    ## [1] 328
    ## [1] 329
    ## [1] 330
    ## [1] 331
    ## [1] 332
    ## [1] 333
    ## [1] 334
    ## [1] 335
    ## [1] 336
    ## [1] 337
    ## [1] 338
    ## [1] 339
    ## [1] 340
    ## [1] 341
    ## [1] 342
    ## [1] 343
    ## [1] 344
    ## [1] 345
    ## [1] 346
    ## [1] 347
    ## [1] 348
    ## [1] 349
    ## [1] 350
    ## [1] 351
    ## [1] 352
    ## [1] 353
    ## [1] 354
    ## [1] 355
    ## [1] 356
    ## [1] 357
    ## [1] 358
    ## [1] 359
    ## [1] 360
    ## [1] 361
    ## [1] 362
    ## [1] 363
    ## [1] 364
    ## [1] 365
    ## [1] 366
    ## [1] 367
    ## [1] 368
    ## [1] 369
    ## [1] 370
    ## [1] 371
    ## [1] 372
    ## [1] 373
    ## [1] 374
    ## [1] 375
    ## [1] 376
    ## [1] 377
    ## [1] 378
    ## [1] 379
    ## [1] 380
    ## [1] 381
    ## [1] 382
    ## [1] 383
    ## [1] 384
    ## [1] 385
    ## [1] 386
    ## [1] 387
    ## [1] 388
    ## [1] 389
    ## [1] 390
    ## [1] 391
    ## [1] 392
    ## [1] 393
    ## [1] 394
    ## [1] 395
    ## [1] 396
    ## [1] 397
    ## [1] 398
    ## [1] 399
    ## [1] 400
    ## [1] 401
    ## [1] 402
    ## [1] 403
    ## [1] 404
    ## [1] 405
    ## [1] 406
    ## [1] 407
    ## [1] 408
    ## [1] 409
    ## [1] 410
    ## [1] 411
    ## [1] 412
    ## [1] 413
    ## [1] 414
    ## [1] 415
    ## [1] 416
    ## [1] 417
    ## [1] 418
    ## [1] 419
    ## [1] 420
    ## [1] 421
    ## [1] 422
    ## [1] 423
    ## [1] 424
    ## [1] 425
    ## [1] 426
    ## [1] 427
    ## [1] 428
    ## [1] 429
    ## [1] 430
    ## [1] 431
    ## [1] 432
    ## [1] 433
    ## [1] 434
    ## [1] 435
    ## [1] 436
    ## [1] 437
    ## [1] 438
    ## [1] 439
    ## [1] 440
    ## [1] 441
    ## [1] 442
    ## [1] 443
    ## [1] 444
    ## [1] 445
    ## [1] 446
    ## [1] 447
    ## [1] 448
    ## [1] 449
    ## [1] 450
    ## [1] 451
    ## [1] 452
    ## [1] 453
    ## [1] 454
    ## [1] 455
    ## [1] 456
    ## [1] 457
    ## [1] 458
    ## [1] 459
    ## [1] 460
    ## [1] 461
    ## [1] 462
    ## [1] 463
    ## [1] 464
    ## [1] 465
    ## [1] 466
    ## [1] 467
    ## [1] 468
    ## [1] 469
    ## [1] 470
    ## [1] 471
    ## [1] 472
    ## [1] 473
    ## [1] 474
    ## [1] 475
    ## [1] 476
    ## [1] 477
    ## [1] 478
    ## [1] 479
    ## [1] 480
    ## [1] 481
    ## [1] 482
    ## [1] 483
    ## [1] 484
    ## [1] 485
    ## [1] 486
    ## [1] 487
    ## [1] 488
    ## [1] 489
    ## [1] 490
    ## [1] 491
    ## [1] 492
    ## [1] 493
    ## [1] 494
    ## [1] 495
    ## [1] 496
    ## [1] 497
    ## [1] 498
    ## [1] 499
    ## [1] 500
    ## [1] 501
    ## [1] 502
    ## [1] 503
    ## [1] 504
    ## [1] 505
    ## [1] 506
    ## [1] 507
    ## [1] 508
    ## [1] 509
    ## [1] 510
    ## [1] 511
    ## [1] 512
    ## [1] 513
    ## [1] 514
    ## [1] 515
    ## [1] 516
    ## [1] 517
    ## [1] 518
    ## [1] 519
    ## [1] 520
    ## [1] 521
    ## [1] 522
    ## [1] 523
    ## [1] 524
    ## [1] 525
    ## [1] 526
    ## [1] 527
    ## [1] 528
    ## [1] 529
    ## [1] 530
    ## [1] 531
    ## [1] 532
    ## [1] 533
    ## [1] 534
    ## [1] 535
    ## [1] 536
    ## [1] 537
    ## [1] 538
    ## [1] 539
    ## [1] 540
    ## [1] 541
    ## [1] 542
    ## [1] 543
    ## [1] 544
    ## [1] 545
    ## [1] 546
    ## [1] 547
    ## [1] 548
    ## [1] 549
    ## [1] 550
    ## [1] 551
    ## [1] 552
    ## [1] 553
    ## [1] 554
    ## [1] 555
    ## [1] 556
    ## [1] 557
    ## [1] 558
    ## [1] 559
    ## [1] 560
    ## [1] 561
    ## [1] 562
    ## [1] 563
    ## [1] 564
    ## [1] 565
    ## [1] 566
    ## [1] 567
    ## [1] 568
    ## [1] 569
    ## [1] 570
    ## [1] 571
    ## [1] 572
    ## [1] 573
    ## [1] 574
    ## [1] 575
    ## [1] 576
    ## [1] 577
    ## [1] 578
    ## [1] 579
    ## [1] 580
    ## [1] 581
    ## [1] 582
    ## [1] 583
    ## [1] 584
    ## [1] 585
    ## [1] 586
    ## [1] 587
    ## [1] 588
    ## [1] 589
    ## [1] 590
    ## [1] 591
    ## [1] 592
    ## [1] 593
    ## [1] 594
    ## [1] 595
    ## [1] 596
    ## [1] 597
    ## [1] 598
    ## [1] 599
    ## [1] 600
    ## [1] 601
    ## [1] 602
    ## [1] 603
    ## [1] 604
    ## [1] 605
    ## [1] 606
    ## [1] 607
    ## [1] 608
    ## [1] 609
    ## [1] 610
    ## [1] 611
    ## [1] 612
    ## [1] 613
    ## [1] 614
    ## [1] 615
    ## [1] 616
    ## [1] 617
    ## [1] 618
    ## [1] 619
    ## [1] 620
    ## [1] 621
    ## [1] 622
    ## [1] 623
    ## [1] 624
    ## [1] 625
    ## [1] 626
    ## [1] 627
    ## [1] 628
    ## [1] 629
    ## [1] 630
    ## [1] 631
    ## [1] 632
    ## [1] 633
    ## [1] 634
    ## [1] 635
    ## [1] 636
    ## [1] 637
    ## [1] 638
    ## [1] 639
    ## [1] 640
    ## [1] 641
    ## [1] 642
    ## [1] 643
    ## [1] 644
    ## [1] 645
    ## [1] 646
    ## [1] 647
    ## [1] 648
    ## [1] 649
    ## [1] 650
    ## [1] 651
    ## [1] 652
    ## [1] 653
    ## [1] 654
    ## [1] 655
    ## [1] 656
    ## [1] 657
    ## [1] 658
    ## [1] 659
    ## [1] 660
    ## [1] 661
    ## [1] 662
    ## [1] 663
    ## [1] 664
    ## [1] 665
    ## [1] 666
    ## [1] 667
    ## [1] 668
    ## [1] 669
    ## [1] 670
    ## [1] 671
    ## [1] 672
    ## [1] 673
    ## [1] 674
    ## [1] 675
    ## [1] 676
    ## [1] 677
    ## [1] 678
    ## [1] 679
    ## [1] 680
    ## [1] 681
    ## [1] 682
    ## [1] 683
    ## [1] 684
    ## [1] 685
    ## [1] 686
    ## [1] 687
    ## [1] 688
    ## [1] 689
    ## [1] 690
    ## [1] 691
    ## [1] 692
    ## [1] 693
    ## [1] 694
    ## [1] 695
    ## [1] 696
    ## [1] 697
    ## [1] 698
    ## [1] 699
    ## [1] 700
    ## [1] 701
    ## [1] 702
    ## [1] 703
    ## [1] 704
    ## [1] 705
    ## [1] 706
    ## [1] 707
    ## [1] 708
    ## [1] 709
    ## [1] 710
    ## [1] 711
    ## [1] 712
    ## [1] 713
    ## [1] 714
    ## [1] 715
    ## [1] 716
    ## [1] 717
    ## [1] 718
    ## [1] 719
    ## [1] 720
    ## [1] 721
    ## [1] 722
    ## [1] 723
    ## [1] 724
    ## [1] 725
    ## [1] 726
    ## [1] 727
    ## [1] 728
    ## [1] 729
    ## [1] 730
    ## [1] 731
    ## [1] 732
    ## [1] 733
    ## [1] 734
    ## [1] 735
    ## [1] 736
    ## [1] 737
    ## [1] 738
    ## [1] 739
    ## [1] 740
    ## [1] 741
    ## [1] 742
    ## [1] 743
    ## [1] 744
    ## [1] 745
    ## [1] 746
    ## [1] 747
    ## [1] 748
    ## [1] 749
    ## [1] 750
    ## [1] 751
    ## [1] 752
    ## [1] 753
    ## [1] 754
    ## [1] 755
    ## [1] 756
    ## [1] 757
    ## [1] 758
    ## [1] 759
    ## [1] 760
    ## [1] 761
    ## [1] 762
    ## [1] 763
    ## [1] 764
    ## [1] 765
    ## [1] 766
    ## [1] 767
    ## [1] 768
    ## [1] 769
    ## [1] 770
    ## [1] 771
    ## [1] 772
    ## [1] 773
    ## [1] 774
    ## [1] 775
    ## [1] 776
    ## [1] 777
    ## [1] 778
    ## [1] 779
    ## [1] 780
    ## [1] 781
    ## [1] 782
    ## [1] 783
    ## [1] 784
    ## [1] 785
    ## [1] 786
    ## [1] 787
    ## [1] 788
    ## [1] 789
    ## [1] 790
    ## [1] 791
    ## [1] 792
    ## [1] 793
    ## [1] 794
    ## [1] 795
    ## [1] 796
    ## [1] 797
    ## [1] 798
    ## [1] 799
    ## [1] 800
    ## [1] 801
    ## [1] 802
    ## [1] 803
    ## [1] 804
    ## [1] 805
    ## [1] 806
    ## [1] 807
    ## [1] 808
    ## [1] 809
    ## [1] 810
    ## [1] 811
    ## [1] 812
    ## [1] 813
    ## [1] 814
    ## [1] 815
    ## [1] 816
    ## [1] 817
    ## [1] 818
    ## [1] 819
    ## [1] 820
    ## [1] 821
    ## [1] 822
    ## [1] 823
    ## [1] 824
    ## [1] 825
    ## [1] 826
    ## [1] 827
    ## [1] 828
    ## [1] 829
    ## [1] 830
    ## [1] 831
    ## [1] 832
    ## [1] 833
    ## [1] 834
    ## [1] 835
    ## [1] 836
    ## [1] 837
    ## [1] 838
    ## [1] 839
    ## [1] 840
    ## [1] 841
    ## [1] 842
    ## [1] 843
    ## [1] 844
    ## [1] 845
    ## [1] 846
    ## [1] 847
    ## [1] 848
    ## [1] 849
    ## [1] 850
    ## [1] 851
    ## [1] 852
    ## [1] 853
    ## [1] 854
    ## [1] 855
    ## [1] 856
    ## [1] 857
    ## [1] 858
    ## [1] 859
    ## [1] 860
    ## [1] 861
    ## [1] 862
    ## [1] 863
    ## [1] 864
    ## [1] 865
    ## [1] 866
    ## [1] 867
    ## [1] 868
    ## [1] 869
    ## [1] 870
    ## [1] 871
    ## [1] 872
    ## [1] 873
    ## [1] 874
    ## [1] 875
    ## [1] 876
    ## [1] 877
    ## [1] 878
    ## [1] 879
    ## [1] 880
    ## [1] 881
    ## [1] 882
    ## [1] 883
    ## [1] 884
    ## [1] 885
    ## [1] 886
    ## [1] 887
    ## [1] 888
    ## [1] 889
    ## [1] 890
    ## [1] 891
    ## [1] 892
    ## [1] 893
    ## [1] 894
    ## [1] 895
    ## [1] 896
    ## [1] 897
    ## [1] 898
    ## [1] 899
    ## [1] 900
    ## [1] 901
    ## [1] 902
    ## [1] 903
    ## [1] 904
    ## [1] 905
    ## [1] 906
    ## [1] 907
    ## [1] 908
    ## [1] 909
    ## [1] 910
    ## [1] 911
    ## [1] 912
    ## [1] 913
    ## [1] 914
    ## [1] 915
    ## [1] 916
    ## [1] 917
    ## [1] 918
    ## [1] 919
    ## [1] 920
    ## [1] 921
    ## [1] 922
    ## [1] 923
    ## [1] 924
    ## [1] 925
    ## [1] 926
    ## [1] 927
    ## [1] 928
    ## [1] 929
    ## [1] 930
    ## [1] 931
    ## [1] 932
    ## [1] 933
    ## [1] 934
    ## [1] 935
    ## [1] 936
    ## [1] 937
    ## [1] 938
    ## [1] 939
    ## [1] 940
    ## [1] 941
    ## [1] 942
    ## [1] 943
    ## [1] 944
    ## [1] 945
    ## [1] 946
    ## [1] 947
    ## [1] 948
    ## [1] 949
    ## [1] 950
    ## [1] 951
    ## [1] 952
    ## [1] 953
    ## [1] 954
    ## [1] 955
    ## [1] 956
    ## [1] 957
    ## [1] 958
    ## [1] 959
    ## [1] 960
    ## [1] 961
    ## [1] 962
    ## [1] 963
    ## [1] 964
    ## [1] 965
    ## [1] 966
    ## [1] 967
    ## [1] 968
    ## [1] 969
    ## [1] 970
    ## [1] 971
    ## [1] 972
    ## [1] 973
    ## [1] 974
    ## [1] 975
    ## [1] 976
    ## [1] 977
    ## [1] 978
    ## [1] 979
    ## [1] 980
    ## [1] 981
    ## [1] 982
    ## [1] 983
    ## [1] 984
    ## [1] 985
    ## [1] 986
    ## [1] 987
    ## [1] 988
    ## [1] 989
    ## [1] 990
    ## [1] 991
    ## [1] 992
    ## [1] 993
    ## [1] 994
    ## [1] 995
    ## [1] 996
    ## [1] 997
    ## [1] 998
    ## [1] 999
    ## [1] 1000

    or = order(colMeans(VARR))
    par(mar = c(4, 4, 4, 4))

    boxplot(VARR[, or] ,
    outline = T,
    ylab = expression(lambda[ij] + sigma[ii]))
    points(1:(p^2), VARTRUE[or], col = 2, cex = 1.5, pch="x")

![](workflow_files/figure-markdown_strict/unnamed-chunk-10-1.png)
