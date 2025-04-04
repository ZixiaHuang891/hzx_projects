#######
####### BIOST 571 Final Project Simulations

#######
#######
####### Here, we run 1,000 simulations under each of five strategies:
  ### 1. Naive approach: Using error-prone data
  ### 2. IPW analysis on SRS-selected validation sample
  ### 3. IPW analysis on optimally-selected validation sample
  ### 4. IPW analysis with generalized raking on SRS-selected validation sample
  ### 5. IPW analysis with generalized raking on optimally-selected validation sample

#######
#######
####### Setup

#### Load libraries
library(Matrix)
library(dplyr)
library(geepack)
library(MIIPW)
library(survey)
library(optimall)
library(MASS)

#### Helper Functions

## Function to compute Equation 7 from paper
compute_generalized_raking_weights <- function(h_star, cluster_ids, R, pi_k) {
  # Convert h_star into a data frame with clear names
  h_star_df <- as.data.frame(h_star)
  # h_star_df <- h_star_df[, 2, drop = FALSE]  # Keep only the second column
  names(h_star_df) <- "h_star_X"
  
  # Create the full dataset
  data <- data.frame(cluster = cluster_ids, R = R, h_star_df)
  
  # Phase 1 totals (population) - Include intercept!
  phase1_totals <- c(`(Intercept)` = nrow(h_star_df), h_star_X = sum(h_star_df$h_star_X))
  
  # Or try removing intercept
  phase1_totals_noint <- c(h_star_X = sum(h_star_df$h_star_X))
  
  # Phase 2 data (sample)
  phase2_data <- subset(data, R == 1)
  
  # Aggregate by cluster (needed for correct dimension matching)
  phase2_cluster_data <- aggregate(. ~ cluster, data = phase2_data[, -2], sum)
  
  phase2_cluster_data$pi_k <- pi_k[match(phase2_cluster_data$cluster, names(pi_k))]
  
  
  # Design object
  numeric_clusters <- phase2_cluster_data$cluster
  phase2_cluster_data$cluster <- as.factor(phase2_cluster_data$cluster)
  design <- svydesign(id = ~cluster, data = phase2_cluster_data, probs = ~pi_k)
  phase2_cluster_data$cluster <- numeric_clusters
  
  # Calibration step (raking)
  rake_design <- tryCatch({
    calibrate(design, 
              formula = ~ h_star_X, 
              population = phase1_totals, 
              calfun = "raking")
  }, error = function(e) {
    # If it fails, return NULL to indicate failure
    return(NULL)
  })
  
  # If rake_design is NULL (error in first step), run the second calibration
  if (is.null(rake_design)) {
    rake_design <- calibrate(design, 
                             formula = ~ h_star_X - 1,  # Exclude intercept with -1
                             population = phase1_totals_noint, 
                             calfun = "raking")
  }
  
  # Calibrated weights
  calibrated_weights <- weights(rake_design)
  names(calibrated_weights) <- phase2_cluster_data$cluster
  
  return(calibrated_weights)
}

# New inf_fun GEE
inf_fun_gee <- function(fit, cluster_ids, rho = 0.1) {
  # Extract model matrix and residuals
  dm <- model.matrix(fit)  # Design matrix (N x p)
  eta <- dm %*% coef(fit)  # Linear predictor
  mu <- fit$fitted.values  # Fitted probabilities
  res <- fit$y - mu        # Response residuals
  
  # Compute diagonal weight matrix W (variance of Y | X)
  W_diag <- mu * (1 - mu)     
  W <- diag(as.vector(W_diag))  # (N x N) weight matrix
  
  # Compute D = diag(mu * (1 - mu)) %*% dm
  D <- W %*% dm  # (N x p)
  
  # Initialize working correlation matrix V
  N <- length(mu)
  V <- matrix(0, nrow = N, ncol = N)
  clusters <- unique(cluster_ids)
  
  for (cl in clusters) {
    idx <- which(cluster_ids == cl)
    n_cl <- length(idx)
    
    # Working correlation matrix R (exchangeable structure)
    R_cl <- matrix(rho, n_cl, n_cl) + diag(1 - rho, n_cl)
    
    # Compute V_cl = sqrt(W) %*% R_cl %*% sqrt(W)
    W_cl <- diag(W_diag[idx])  # Weight submatrix for cluster
    V_cl <- sqrt(W_cl) %*% R_cl %*% sqrt(W_cl)
    
    V[idx, idx] <- V_cl  # Insert into full V matrix
  }
  
  # Compute Fisher information matrix using V
  I_hat <- t(D) %*% solve(V) %*% D  # (p x p)
  
  # Solve for the inverse of I_hat
  I_hat_inv <- solve(I_hat)  # (p x p)
  
  # Compute score function U
  U <- t(D) %*% solve(V) %*% res  # (p x 1)
  
  # Influence function
  infl <- D %*% I_hat_inv %*% U  # (N x p) * (p x 1) = (N x 1)
  
  return(infl)
}


#######
#######
####### Strategy 1: Naive approach
one_sim_naive_SRS <- function(K, Nk_range, beta0, beta1, beta2, beta3, gamma, 
                           sampling_frac = 1/4, seed){
  set.seed(seed)
  ### Generate Data
  # Simulation parameters
  K <- K            # Number of clusters
  
  # True beta values
  beta0 <- beta0
  beta1 <- beta1      # Parameter of interest
  beta2 <- beta2
  beta3 <- beta3
  
  gamma <- gamma  # Effect of Vk on X*ki
  
  # Generate categorical cluster-level covariate
  cluster_cov <- sample(c("A", "B"), prob = c(0.6, 0.4), K, replace=TRUE)
  cluster_cov <- sort(cluster_cov)
  # Assign cluster sizes
  Nk_A_range <- Nk_range 
  Nk_B_range <- Nk_range 
  
  # Count the number of clusters in each stratum
  num_A <- sum(cluster_cov == "A")
  num_B <- sum(cluster_cov == "B")
  
  # Assign cluster sizes accordingly
  Nk_A <- sample(Nk_A_range, num_A, replace=TRUE)  # Sizes for A clusters
  Nk_B <- sample(Nk_B_range, num_B, replace=TRUE)  # Sizes for B clusters
  Nk <- c(Nk_A, Nk_B)
  
  N <- sum(Nk)
  Vk <- rep(cluster_cov, times=Nk)
  
  # Generate cluster-level random effect Rk
  Rk <- rnorm(K, mean=0, sd = 0.3)
  cluster_ids <- rep(1:K, times=Nk)
  Rk_full <- rep(Rk, times=Nk)
  
  # Generate individual-level covariates
  covs <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(1, 0.15, 0.15, 1), nrow=2))
  Xki <- covs[,1]
  Wki <- covs[,2]
  
  # Compute true linear predictor and mu
  eta <- beta0 + beta1*Xki + beta2*Wki + beta3*ifelse(Vk == "B", 1, 0) + Rk_full
  mu <- exp(eta) / (1 + exp(eta))
  
  # True binary outcome Yki
  Yki <- rbinom(N, 1, mu)
  
  # Generate error-prone observed data
  Xki_star <- Xki + gamma*ifelse(Vk == "B", 1,0) + 0.3*Wki +  rnorm(N, mean=0, sd=0.5)
  Yki_star <- ifelse(Yki == 1, rbinom(N, 1, 0.9), rbinom(N, 1, 0.05))
  
  # Phase 1 data (Y*, X*, W, Vk)
  phase1_data <- data.frame(Y_star=Yki_star, X_star=Xki_star, W=Wki, V=Vk, cluster=cluster_ids)
  
  naive_model <- geeglm(Y_star ~ X_star + W + V, family=binomial, data=phase1_data,
                        id = cluster, corstr = "exchangeable")
  final_results <- summary(naive_model)
  
  return(c(final_results$coefficients["X_star", "Estimate"],
    final_results$coefficients["X_star", "Std.err"]))
}

# run 1000 simulations and store results in a dataframe
num_cores <- detectCores() - 2
cl <- makeCluster(num_cores)

# Export necessary variables and functions to the cluster
clusterExport(cl, c("one_sim_naive_SRS", "inf_fun_gee",
                    "compute_generalized_raking_weights"))

# Run the simulations in parallel
sim_results_naive_list <- parLapply(cl, 1:1000, function(i) {
  library(survey)
  library(MASS)
  library(geepack)
  one_sim_naive_SRS(120, 5:7, -2, 0.7, 0.1, 0.8, 0.5, 1/3, seed = i)
})

# Combine the results into a single data frame
sim_results_naive <- as.data.frame(do.call(rbind, sim_results_naive_list))

# Stop the cluster
stopCluster(cl)
# save(sim_results_naive, file = "/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/SecondYearCoursework/BIOST571/Final Project/sim_results_naive_SRS.rds")

mean(sim_results_naive[,1])
mean(sim_results_naive[,2])
var(sim_results_naive[,1])

# Coverage
sim_results_naive$lower <- sim_results_naive[,1] - 1.96*sim_results_naive[,2]
sim_results_naive$upper <- sim_results_naive[,1] + 1.96*sim_results_naive[,2]
sim_results_naive$cover <- ifelse(0.7 >= sim_results_naive$lower & 0.7 <= sim_results_naive$upper, 1, 0)
sum(sim_results_naive$cover)/1000


#######
#######
####### Strategy 2: IPW analysis on SRS-selected validation sample
one_sim_IPW_SRS <- function(K, Nk_range, beta0, beta1, beta2, beta3, gamma, 
                           sampling_frac = 1/4, seed){
  set.seed(seed)
  ### Generate Data
  # Simulation parameters
  K <- K            # Number of clusters
  
  # True beta values
  beta0 <- beta0
  beta1 <- beta1      # Parameter of interest
  beta2 <- beta2
  beta3 <- beta3
  
  gamma <- gamma  # Effect of Vk on X*ki
  
  # Generate categorical cluster-level covariate
  cluster_cov <- sample(c("A", "B"), prob = c(0.6, 0.4), K, replace=TRUE)
  cluster_cov <- sort(cluster_cov)
  # Assign cluster sizes
  Nk_A_range <- Nk_range 
  Nk_B_range <- Nk_range 
  
  # Count the number of clusters in each stratum
  num_A <- sum(cluster_cov == "A")
  num_B <- sum(cluster_cov == "B")
  
  # Assign cluster sizes accordingly
  Nk_A <- sample(Nk_A_range, num_A, replace=TRUE)  # Sizes for A clusters
  Nk_B <- sample(Nk_B_range, num_B, replace=TRUE)  # Sizes for B clusters
  Nk <- c(Nk_A, Nk_B)
  
  N <- sum(Nk)
  Vk <- rep(cluster_cov, times=Nk)
  
  # Generate cluster-level random effect Rk
  Rk <- rnorm(K, mean=0, sd = 0.3)
  cluster_ids <- rep(1:K, times=Nk)
  Rk_full <- rep(Rk, times=Nk)
  
  # Generate individual-level covariates
  covs <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(1, 0.15, 0.15, 1), nrow=2))
  Xki <- covs[,1]
  Wki <- covs[,2]
  
  # Compute true linear predictor and mu
  eta <- beta0 + beta1*Xki + beta2*Wki + beta3*ifelse(Vk == "B", 1, 0) + Rk_full
  mu <- exp(eta) / (1 + exp(eta))
  
  # True binary outcome Yki
  Yki <- rbinom(N, 1, mu)
  
  # Generate error-prone observed data
  Xki_star <- Xki + gamma*ifelse(Vk == "B", 1,0) + 0.3*Wki +  rnorm(N, mean=0, sd=0.5)
  Yki_star <- ifelse(Yki == 1, rbinom(N, 1, 0.9), rbinom(N, 1, 0.05))
  
  # Phase 1 data (Y*, X*, W, Vk)
  phase1_data <- data.frame(Y_star=Yki_star, X_star=Xki_star, W=Wki, V=Vk, cluster=cluster_ids)
  
  #######
  ####### Phase 2
  
  ## Sample Phase 2
  # Select clusters for phase 2 validation (SRS)
  phase2_clusters <- sample(unique(cluster_ids), size=round(K*sampling_frac), replace=FALSE)
  phase1_data$R <- ifelse(phase1_data$cluster %in% phase2_clusters, 1, 0)
  
  # True data collected in phase 2
  phase2_data <- subset(phase1_data, R==1)
  phase2_data$Y <- Yki[phase1_data$R==1]
  phase2_data$X <- Xki[phase1_data$R==1]
  
  # Phase-2 inclusion probabilities pi_k (simplified equal probs)
  pi_k <- rep(length(phase2_clusters)/K, K)
  names(pi_k) <- 1:K
  
  # Final GEE estimation using true data (phase 2) and IPW weights
  gee_IPW_Survey <- function(data, weights){
    data$pik <- weights[as.character(data$cluster)]
    data$pik_inv <- 1 / weights[as.character(data$cluster)]
    pik_inv <- 1 / 1 / weights[as.character(data$cluster)]
    data$cluster <- as.factor(data$cluster)
    mydesign <- svydesign(ids = ~cluster, probs = ~pik, data = data)
    model <- svyglm(Y ~ X + W + V, family=quasibinomial, design = mydesign)
    return(summary(model))
  }
  
  #final_results2 <- gee_generalized_raking_SIPW(phase2_data, calibrated_weights)$beta
  final_results3 <- gee_IPW_Survey(phase2_data, pi_k)
  
  return(c(#final_results2["X",1],
    final_results3$coefficients["X", "Estimate"],
    final_results3$coefficients["X", "Std. Error"], length(cluster_ids)))
}

# run 1000 simulations and store results in a dataframe
num_cores <- detectCores() - 2
cl <- makeCluster(num_cores)

# Export necessary variables and functions to the cluster
clusterExport(cl, c("one_sim_IPW_SRS", "inf_fun_gee",
                    "compute_generalized_raking_weights"))

# Run the simulations in parallel
sim_results_list_IPW_SRS <- parLapply(cl, 1:1000, function(i) {
  library(survey)
  library(MASS)
  library(geepack)
  one_sim_IPW_SRS(120, 5:7, -2, 0.7, 0.1, 0.8, 0.5, 1/3, seed = i)
})

# Combine the results into a single data frame
sim_results_IPW_SRS <- as.data.frame(do.call(rbind, sim_results_list_IPW_SRS))

# Stop the cluster
stopCluster(cl)
save(sim_results_IPW_SRS, file = "/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/SecondYearCoursework/BIOST571/Final Project/sim_results_IPW_SRS.rds")

mean(sim_results_IPW_SRS[,1])
mean(sim_results_IPW_SRS[,2])
var(sim_results_IPW_SRS[,1])

# Coverage
sim_results_IPW_SRS$lower <- sim_results_IPW_SRS[,1] - 1.96*sim_results_IPW_SRS[,2]
sim_results_IPW_SRS$upper <- sim_results_IPW_SRS[,1] + 1.96*sim_results_IPW_SRS[,2]
sim_results_IPW_SRS$cover <- ifelse(0.7 >= sim_results_IPW_SRS$lower & 0.7 <= sim_results_IPW_SRS$upper, 1, 0)
sum(sim_results_IPW_SRS$cover)/1000

#######
#######
####### Strategy 3: IPW analysis on optimally-selected validation sample
one_sim_IPW_OALL <- function(K, Nk_range, beta0, beta1, beta2, beta3, gamma,
                            sampling_frac, seed){
  set.seed(seed)
  ### Generate Data
  # Simulation parameters
  K <- K            # Number of clusters
  
  # True beta values
  beta0 <- beta0
  beta1 <- beta1      # Parameter of interest
  beta2 <- beta2
  beta3 <- beta3
  
  gamma <- gamma  # Effect of Vk on X*ki
  
  # Generate categorical cluster-level covariate
  cluster_cov <- sample(c("A", "B"), prob = c(0.6, 0.4), K, replace=TRUE)
  cluster_cov <- sort(cluster_cov)
  # Assign cluster sizes
  Nk_A_range <- Nk_range
  Nk_B_range <- Nk_range 
  
  # Count the number of clusters in each stratum
  num_A <- sum(cluster_cov == "A")
  num_B <- sum(cluster_cov == "B")
  
  # Assign cluster sizes accordingly
  Nk_A <- sample(Nk_A_range, num_A, replace=TRUE)  # Sizes for A clusters
  Nk_B <- sample(Nk_B_range, num_B, replace=TRUE)  # Sizes for B clusters
  Nk <- c(Nk_A, Nk_B)
  
  N <- sum(Nk)
  Vk <- rep(cluster_cov, times=Nk)
  
  # Generate cluster-level random effect Rk
  Rk <- rnorm(K, mean=0, sd = 0.5)
  cluster_ids <- rep(1:K, times=Nk)
  Rk_full <- rep(Rk, times=Nk)
  
  # Generate individual-level covariates
  covs <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(1, 0.15, 0.15, 1), nrow=2))
  Xki <- covs[,1]
  Wki <- covs[,2]
  
  # Compute true linear predictor and mu
  eta <- beta0 + beta1*Xki + beta2*Wki + beta3*ifelse(Vk == "B", 1, 0) + Rk_full
  mu <- exp(eta) / (1 + exp(eta))
  
  # True binary outcome Yki
  Yki <- rbinom(N, 1, mu)
  
  # Generate error-prone observed data
  Xki_star <- Xki + gamma*ifelse(Vk == "B", 1,0) + 0.3*Wki +  rnorm(N, mean=0, sd=0.5)
  Yki_star <- ifelse(Yki == 1, rbinom(N, 1, 0.9), rbinom(N, 1, 0.05))
  
  #var(h_star[phase1_data$strata == "high"])
  #var(h_star[phase1_data$strata == "low"])
  
  # Generate strata
  test <- dplyr::group_by(data.frame(cluster = cluster_ids, Yki_star, Xki_star,
                                     Vk), cluster) %>%
    dplyr::summarize(mean_Y = mean(Yki_star), Vk = Vk[1]) 
  
  quantileA <- quantile(test$mean_Y[test$Vk == "A"], 0.4)
  quantileB <- quantile(test$mean_Y[test$Vk == "B"], 0.35)
  
  test <- test %>%
    dplyr::mutate(strata = ifelse(Vk == "A" & mean_Y > quantileA, "high.A",
                                  ifelse(Vk == "A" & mean_Y <= quantileA, "low.A",
                                         ifelse(Vk == "B" & mean_Y > quantileB, "high.B",
                                                "low.B"))))
  
  # Phase 1 data (Y*, X*, W, Vk)
  phase1_data <- data.frame(Y_star=Yki_star, X_star=Xki_star, W=Wki, V=Vk, cluster=cluster_ids) %>%
    dplyr::left_join(dplyr::select(test, cluster, strata), by = "cluster")
  
  #######
  ####### Phase 2
  
  # Compute influence functions for optimum allocation
  # Calculate influence functions h_star based on phase 1 estimates
  
  # Initial GEE estimation with error-prone data (phase 1)
  gee_phase1 <- geeglm(Y_star ~ X_star + W + V, family=binomial, data=phase1_data, 
                       id=cluster, corstr="exchangeable")
  
  beta_phase1 <- coef(gee_phase1)
  X_phase1 <- model.matrix(~ X_star + W + V, data=phase1_data)
  # h_star <- compute_h_ki_logistic(phase1_data$Y_star, X_phase1, beta_phase1, phase1_data$cluster)
  h_star <- inf_fun_gee(gee_phase1, phase1_data$cluster)
  
  ######
  ###### Now, using h_star, I compute everything needed for the optimum allocation
  compute_A_j <- function(h_star, cluster_ids, I_j) {
    A_j <- 0
    for (k in I_j) {
      indices_k <- which(cluster_ids == k)
      h_k <- h_star[indices_k]
      A_j <- A_j + sum(outer(h_k, h_k, "*")) # Compute sum of all pairwise products within cluster
    }
    return(A_j)
  }
  
  compute_B_j <- function(h_star, cluster_ids, I_j) {
    B_j <- 0
    cluster_list <- unique(cluster_ids)
    
    for (k in I_j) {
      indices_k <- which(cluster_ids == k)
      h_k <- h_star[indices_k]
      
      for (k_prime in I_j[I_j != k]) { # Only consider k' != k
        indices_k_prime <- which(cluster_ids == k_prime)
        h_k_prime <- h_star[indices_k_prime]
        
        B_j <- B_j + sum(outer(h_k, h_k_prime, "*")) # Compute sum of pairwise products between clusters
      }
    }
    return(B_j)
  }
  
  compute_S_j <- function(h_star, cluster_ids, strata_df) {
    # Get unique strata
    unique_strata <- unique(strata_df$strata)
    
    # Initialize result storage
    results <- data.frame(stratum = unique_strata, S_j = NA)
    
    for (j in unique_strata) {
      # Get clusters in this stratum
      I_j <- strata_df$cluster[strata_df$strata == j]
      K_j <- length(I_j)  # Number of clusters in this stratum
      
      # Compute A_j and B_j
      A_j <- compute_A_j(h_star, cluster_ids, I_j)
      B_j <- compute_B_j(h_star, cluster_ids, I_j)
      
      # Compute S_j
      S_j <- A_j - B_j / (K_j - 1)
      
      # Store result
      results$S_j[results$stratum == j] <- S_j
    }
    
    return(results)
  }
  # Now apply it
  strata_df <- dplyr::select(test, cluster, strata)
  S_j_results <- compute_S_j(h_star, cluster_ids, strata_df)
  
  Nj <- strata_df %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(n = n())
  myNj <- Nj$n
  
  S_j_results <- dplyr::arrange(S_j_results, stratum)
  S_j_results$Nj <- myNj
  
  # Optimal allocation
  cluster_level_df <- strata_df
  cluster_level_df$cluster <- as.character(cluster_level_df$cluster)
  cluster_level_df <- as.data.frame(cluster_level_df)
  
  allocation <- optimall::optimum_allocation(S_j_results, strata = "stratum", 
                                             sd_h = "S_j",
                                             N_h = "Nj", 
                                             nsample = round(K*sampling_frac),
                                             ndigits = 10,
                                             method = "WrightII")
  
  ### TEMP
  # S_j_results_temp <- S_j_results
  # S_j_results_temp$S_j_new <- c(3,3,1, 1)
  # allocation <- optimall::optimum_allocation(S_j_results_temp, strata = "stratum", 
  #                                            sd_h = "S_j_new",
  #                                            N_h = "Nj", 
  #                                            nsample = round(K*sampling_frac),
  #                                            ndigits = 10,
  #                                            method = "WrightII")
  
  cluster_level_df <- sample_strata(cluster_level_df, design_data = allocation,
                                    strata = "strata", id = "cluster",
                                    design_strata = "strata",
                                    n_allocated = "stratum_size")
  phase2_clusters <- unique(dplyr::filter(cluster_level_df, sample_indicator == 1)$id)
  phase1_data$R <- ifelse(phase1_data$cluster %in% phase2_clusters, 1, 0)
  
  # Inclusion probs
  allocation$pi_k <- allocation$stratum_size/allocation$npop
  phase1_data$pi_k <- dplyr::left_join(phase1_data, 
                                       dplyr::select(allocation,
                                                     "strata", "pi_k"), 
                                       by = c("strata" = "strata"))$pi_k
  
  # True data collected in phase 2
  phase2_data <- subset(phase1_data, R==1)
  phase2_data$Y <- Yki[phase1_data$R==1]
  phase2_data$X <- Xki[phase1_data$R==1]
  phase2_data$W <- Wki[phase1_data$R==1]
  phase2_data$V <- Vk[phase1_data$R==1]
  
  # Get vector of inclusion probs
  pi_k_df <- unique(dplyr::select(phase2_data, cluster, pi_k))
  pi_k <- pi_k_df$pi_k
  names(pi_k) <- pi_k_df$cluster
  
  # Final GEE estimation using true data (phase 2) and IPW weights
  gee_IPW_Survey <- function(data, weights){
    data$pik <- weights[as.character(data$cluster)]
    data$pik_inv <- 1 / weights[as.character(data$cluster)]
    pik_inv <- 1 / 1 / weights[as.character(data$cluster)]
    data$cluster <- as.factor(data$cluster)
    mydesign <- svydesign(ids = ~cluster, probs = ~pik, data = data)
    model <- svyglm(Y ~ X + W + V, family=quasibinomial, design = mydesign)
    return(summary(model))
  }

  #final_results2 <- gee_generalized_raking_SIPW(phase2_data, calibrated_weights)$beta
  final_results3 <- gee_IPW_Survey(phase2_data, pi_k)
  
  return(c(#final_results2["X",1],
    final_results3$coefficients["X", "Estimate"],
    final_results3$coefficients["X", "Std. Error"], length(cluster_ids)))
}


# run 1000 simulations and store results in a dataframe
num_cores <- detectCores() - 2
cl <- makeCluster(num_cores)

# Export necessary variables and functions to the cluster
clusterExport(cl, c("one_sim_IPW_OALL", "inf_fun_gee"))

# Run the simulations in parallel
sim_results_IPW_list <- parLapply(cl, 1:1000, function(i) {
  library(survey)
  library(MASS)
  library(geepack)
  library(dplyr)
  library(optimall)
  one_sim_IPW_OALL(120, 5:7, -2, 0.7, 0.1, 0.8, 0.5, 1/3, seed = i)
})

# Combine the results into a single data frame
sim_results_IPW <- as.data.frame(do.call(rbind, sim_results_IPW_list))

# Stop the cluster
stopCluster(cl)
# save(sim_results_IPW, file = "/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/SecondYearCoursework/BIOST571/Final Project/sim_results_IPW_OALL.rds")

mean(sim_results_IPW[,1])
mean(sim_results_IPW[,2])
var(sim_results_IPW[,1])

# Coverage
sim_results_IPW$lower <- sim_results_IPW[,1] - 1.96*sim_results_IPW[,2]
sim_results_IPW$upper <- sim_results_IPW[,1] + 1.96*sim_results_IPW[,2]
sim_results_IPW$cover <- ifelse(0.7 >= sim_results_IPW$lower & 0.7 <= sim_results_IPW$upper, 1, 0)
sum(sim_results_IPW$cover)/1000

#######
#######
####### Strategy 3: IPW analysis with generalized raking on SRS-selected validation sample
one_sim_GR_SRS <- function(K, Nk_range, beta0, beta1, beta2, beta3, gamma, 
                           sampling_frac = 1/4, seed){
  set.seed(seed)
  ### Generate Data
  # Simulation parameters
  K <- K            # Number of clusters
  
  # True beta values
  beta0 <- beta0
  beta1 <- beta1      # Parameter of interest
  beta2 <- beta2
  beta3 <- beta3
  
  gamma <- gamma  # Effect of Vk on X*ki
  
  # Generate categorical cluster-level covariate
  cluster_cov <- sample(c("A", "B"), prob = c(0.6, 0.4), K, replace=TRUE)
  cluster_cov <- sort(cluster_cov)
  # Assign cluster sizes
  Nk_A_range <- Nk_range 
  Nk_B_range <- Nk_range 
  
  # Count the number of clusters in each stratum
  num_A <- sum(cluster_cov == "A")
  num_B <- sum(cluster_cov == "B")
  
  # Assign cluster sizes accordingly
  Nk_A <- sample(Nk_A_range, num_A, replace=TRUE)  # Sizes for A clusters
  Nk_B <- sample(Nk_B_range, num_B, replace=TRUE)  # Sizes for B clusters
  Nk <- c(Nk_A, Nk_B)
  
  N <- sum(Nk)
  Vk <- rep(cluster_cov, times=Nk)
  
  # Generate cluster-level random effect Rk
  Rk <- rnorm(K, mean=0, sd = 0.3)
  cluster_ids <- rep(1:K, times=Nk)
  Rk_full <- rep(Rk, times=Nk)
  
  # Generate individual-level covariates
  covs <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(1, 0.15, 0.15, 1), nrow=2))
  Xki <- covs[,1]
  Wki <- covs[,2]
  
  # Compute true linear predictor and mu
  eta <- beta0 + beta1*Xki + beta2*Wki + beta3*ifelse(Vk == "B", 1, 0) + Rk_full
  mu <- exp(eta) / (1 + exp(eta))
  
  # True binary outcome Yki
  Yki <- rbinom(N, 1, mu)
  
  # Generate error-prone observed data
  Xki_star <- Xki + gamma*ifelse(Vk == "B", 1,0) + 0.3*Wki +  rnorm(N, mean=0, sd=0.5)
  Yki_star <- ifelse(Yki == 1, rbinom(N, 1, 0.9), rbinom(N, 1, 0.05))
  
  # Phase 1 data (Y*, X*, W, Vk)
  phase1_data <- data.frame(Y_star=Yki_star, X_star=Xki_star, W=Wki, V=Vk, cluster=cluster_ids)
  
  #######
  ####### Phase 2
  
  ## Sample Phase 2
  # Select clusters for phase 2 validation (Neyman allocation approximation here)
  phase2_clusters <- sample(unique(cluster_ids), size=round(K*sampling_frac), replace=FALSE)
  phase1_data$R <- ifelse(phase1_data$cluster %in% phase2_clusters, 1, 0)
  
  # True data collected in phase 2
  phase2_data <- subset(phase1_data, R==1)
  phase2_data$Y <- Yki[phase1_data$R==1]
  phase2_data$X <- Xki[phase1_data$R==1]
  
  # Initial GEE estimation with error-prone data (phase 1)
  gee_phase1 <- geeglm(Y_star ~ X_star + W + V, family=binomial, 
                       data=phase1_data,
                       id=cluster, corstr="exchangeable")
  
  # Calculate influence functions h_star based on phase 1 estimates
  beta_phase1 <- coef(gee_phase1)
  X_phase1 <- model.matrix(~ X_star + W + V, data=phase1_data)
  # h_star <- compute_h_ki_logistic(phase1_data$Y_star, X_phase1, beta_phase1, phase1_data$cluster)
  h_star <- inf_fun_gee(gee_phase1, phase1_data$cluster)
  
  # Phase-2 inclusion probabilities pi_k (simplified equal probs)
  pi_k <- rep(length(phase2_clusters)/K, K)
  names(pi_k) <- 1:K
  
  # Generalized raking weights
  calibrated_weights <- compute_generalized_raking_weights(h_star, cluster_ids, phase1_data$R, pi_k)
  
  # Final GEE estimation using true data (phase 2) and calibrated weights
  gee_generalized_raking <- function(data, calibrated_weights){
    data$gk_inv <- 1 / calibrated_weights[as.character(data$cluster)]
    model <- suppressWarnings(geeglm(Y ~ X + W + V, family=binomial, data=data, id=cluster,
                                     weights=gk_inv, corstr="exchangeable"))
    return(summary(model))
  }
  
  # gee_generalized_raking_SIPW <- function(data, calibrated_weights){
  #   data$gk_inv <- calibrated_weights[as.character(data$cluster)]
  #   data <- data %>%
  #     dplyr::group_by(cluster) %>%
  #     dplyr::mutate(myvisit = row_number()) %>%
  #     dplyr::ungroup()
  #   gk_inv <- 1 / calibrated_weights[as.character(data$cluster)]
  #   model <- SIPW(Y ~ X + W + V, family='binomial', data=data, id='cluster', visit = 'myvisit',
  #                   weights='gk_inv', corstr="exchangeable")
  #   return(model)
  # }
  
  gee_generalized_raking_Survey <- function(data, calibrated_weights){
    data$gk_inv <- 1 / calibrated_weights[as.character(data$cluster)]
    gk_inv <- 1 / calibrated_weights[as.character(data$cluster)]
    data$cluster <- as.factor(data$cluster)
    mydesign <- svydesign(ids = ~cluster, probs = ~gk_inv, data = data,
                          calibration.formula = ~1)
    model <- svyglm(Y ~ X + W + V, family=quasibinomial, design = mydesign,
                    weights=gk_inv)
    return(summary(model))
  }
  
  # final_results <- gee_generalized_raking(phase2_data, calibrated_weights)
  # final_results2 <- gee_generalized_raking_SIPW(phase2_data, calibrated_weights)$beta
  final_results3 <- gee_generalized_raking_Survey(phase2_data, calibrated_weights)
  
  return(c(#final_results$coefficients["X", "Estimate"],
           final_results3$coefficients["X", "Estimate"],
           final_results3$coefficients["X", "Std. Error"],
           length(phase2_clusters)))
}

# run 1000 simulations and store results in a dataframe
num_cores <- detectCores() - 2
cl <- makeCluster(num_cores)

# Export necessary variables and functions to the cluster
clusterExport(cl, c("one_sim_GR_SRS", "inf_fun_gee",
                    "compute_generalized_raking_weights"))

# Run the simulations in parallel
sim_results_list <- parLapply(cl, 1:1000, function(i) {
  library(survey)
  library(MASS)
  library(geepack)
  one_sim_GR_SRS(120, 5:7, -2, 0.7, 0.1, 0.8, 0.5, 1/3, seed = i)
})

# Combine the results into a single data frame
sim_results <- as.data.frame(do.call(rbind, sim_results_list))

# Stop the cluster
stopCluster(cl)
sim_results_GR_SRS <- sim_results
# save(sim_results_GR_SRS, file = "/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/SecondYearCoursework/BIOST571/Final Project/sim_results_GR_SRS.rds")

mean(sim_results[,1])
mean(sim_results[,2])
var(sim_results[,1])

# Coverage
sim_results$lower <- sim_results[,1] - 1.96*sim_results[,2]
sim_results$upper <- sim_results[,1] + 1.96*sim_results[,2]
sim_results$cover <- ifelse(0.7 >= sim_results$lower & 0.7 <= sim_results$upper, 1, 0)
sum(sim_results$cover)/1000

#######
#######
####### Strategy 4: IPW analysis with generalized raking on optimally-selected validation sample

one_sim_GR_OALL <- function(K, Nk_range, beta0, beta1, beta2, beta3, gamma,
                            sampling_frac, seed){
  set.seed(seed)
  ### Generate Data
  # Simulation parameters
  K <- K            # Number of clusters
  
  # True beta values
  beta0 <- beta0
  beta1 <- beta1      # Parameter of interest
  beta2 <- beta2
  beta3 <- beta3
  
  gamma <- gamma  # Effect of Vk on X*ki
  
  # Generate categorical cluster-level covariate
  cluster_cov <- sample(c("A", "B"), prob = c(0.6, 0.4), K, replace=TRUE)
  cluster_cov <- sort(cluster_cov)
  # Assign cluster sizes
  Nk_A_range <- Nk_range
  Nk_B_range <- Nk_range 
  
  # Count the number of clusters in each stratum
  num_A <- sum(cluster_cov == "A")
  num_B <- sum(cluster_cov == "B")
  
  # Assign cluster sizes accordingly
  Nk_A <- sample(Nk_A_range, num_A, replace=TRUE)  # Sizes for A clusters
  Nk_B <- sample(Nk_B_range, num_B, replace=TRUE)  # Sizes for B clusters
  Nk <- c(Nk_A, Nk_B)
  
  N <- sum(Nk)
  Vk <- rep(cluster_cov, times=Nk)
  
  # Generate cluster-level random effect Rk
  Rk <- rnorm(K, mean=0, sd = 0.5)
  cluster_ids <- rep(1:K, times=Nk)
  Rk_full <- rep(Rk, times=Nk)
  
  # Generate individual-level covariates
  covs <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(1, 0.15, 0.15, 1), nrow=2))
  Xki <- covs[,1]
  Wki <- covs[,2]
  
  # Compute true linear predictor and mu
  eta <- beta0 + beta1*Xki + beta2*Wki + beta3*ifelse(Vk == "B", 1, 0) + Rk_full
  mu <- exp(eta) / (1 + exp(eta))
  
  # True binary outcome Yki
  Yki <- rbinom(N, 1, mu)
  
  # Generate error-prone observed data
  Xki_star <- Xki + gamma*ifelse(Vk == "B", 1,0) + 0.3*Wki +  rnorm(N, mean=0, sd=0.5)
  Yki_star <- ifelse(Yki == 1, rbinom(N, 1, 0.9), rbinom(N, 1, 0.05))
  
  #var(h_star[phase1_data$strata == "high"])
  #var(h_star[phase1_data$strata == "low"])
  
  # Generate strata
  test <- dplyr::group_by(data.frame(cluster = cluster_ids, Yki_star, Xki_star,
                                     Vk), cluster) %>%
    dplyr::summarize(mean_Y = mean(Yki_star), Vk = Vk[1]) 
  
  quantileA <- quantile(test$mean_Y[test$Vk == "A"], 0.4)
  quantileB <- quantile(test$mean_Y[test$Vk == "B"], 0.35)
  
  test <- test %>%
    dplyr::mutate(strata = ifelse(Vk == "A" & mean_Y > quantileA, "high.A",
                                  ifelse(Vk == "A" & mean_Y <= quantileA, "low.A",
                                         ifelse(Vk == "B" & mean_Y > quantileB, "high.B",
                                                "low.B"))))
  
  # Phase 1 data (Y*, X*, W, Vk)
  phase1_data <- data.frame(Y_star=Yki_star, X_star=Xki_star, W=Wki, V=Vk, cluster=cluster_ids) %>%
    dplyr::left_join(dplyr::select(test, cluster, strata), by = "cluster")
  
  #######
  ####### Phase 2
  
  # Compute influence functions for optimum allocation
  # Calculate influence functions h_star based on phase 1 estimates
  
  # Initial GEE estimation with error-prone data (phase 1)
  gee_phase1 <- geeglm(Y_star ~ X_star + W + V, family=binomial, data=phase1_data, 
                       id=cluster, corstr="exchangeable")
  
  beta_phase1 <- coef(gee_phase1)
  X_phase1 <- model.matrix(~ X_star + W + V, data=phase1_data)
  # h_star <- compute_h_ki_logistic(phase1_data$Y_star, X_phase1, beta_phase1, phase1_data$cluster)
  h_star <- inf_fun_gee(gee_phase1, phase1_data$cluster)
  
  ######
  ###### Now, using h_star, I compute everything needed for the optimum allocation
  compute_A_j <- function(h_star, cluster_ids, I_j) {
    A_j <- 0
    for (k in I_j) {
      indices_k <- which(cluster_ids == k)
      h_k <- h_star[indices_k]
      A_j <- A_j + sum(outer(h_k, h_k, "*")) # Compute sum of all pairwise products within cluster
    }
    return(A_j)
  }
  
  compute_B_j <- function(h_star, cluster_ids, I_j) {
    B_j <- 0
    cluster_list <- unique(cluster_ids)
    
    for (k in I_j) {
      indices_k <- which(cluster_ids == k)
      h_k <- h_star[indices_k]
      
      for (k_prime in I_j[I_j != k]) { # Only consider k' != k
        indices_k_prime <- which(cluster_ids == k_prime)
        h_k_prime <- h_star[indices_k_prime]
        
        B_j <- B_j + sum(outer(h_k, h_k_prime, "*")) # Compute sum of pairwise products between clusters
      }
    }
    return(B_j)
  }
  
  compute_S_j <- function(h_star, cluster_ids, strata_df) {
    # Get unique strata
    unique_strata <- unique(strata_df$strata)
    
    # Initialize result storage
    results <- data.frame(stratum = unique_strata, S_j = NA)
    
    for (j in unique_strata) {
      # Get clusters in this stratum
      I_j <- strata_df$cluster[strata_df$strata == j]
      K_j <- length(I_j)  # Number of clusters in this stratum
      
      # Compute A_j and B_j
      A_j <- compute_A_j(h_star, cluster_ids, I_j)
      B_j <- compute_B_j(h_star, cluster_ids, I_j)
      
      # Compute S_j
      S_j <- A_j - B_j / (K_j - 1)
      
      # Store result
      results$S_j[results$stratum == j] <- S_j
    }
    
    return(results)
  }
  # Now apply it
  strata_df <- dplyr::select(test, cluster, strata)
  S_j_results <- compute_S_j(h_star, cluster_ids, strata_df)
  
  Nj <- strata_df %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(n = n())
  myNj <- Nj$n
  
  S_j_results <- dplyr::arrange(S_j_results, stratum)
  S_j_results$Nj <- myNj
  
  # Optimal allocation
  cluster_level_df <- strata_df
  cluster_level_df$cluster <- as.character(cluster_level_df$cluster)
  cluster_level_df <- as.data.frame(cluster_level_df)
  
  allocation <- optimall::optimum_allocation(S_j_results, strata = "stratum", 
                                             sd_h = "S_j",
                                             N_h = "Nj", 
                                             nsample = round(K*sampling_frac),
                                             ndigits = 10,
                                             method = "WrightII")
  
  cluster_level_df <- sample_strata(cluster_level_df, design_data = allocation,
                                    strata = "strata", id = "cluster",
                                    design_strata = "strata",
                                    n_allocated = "stratum_size")
  phase2_clusters <- unique(dplyr::filter(cluster_level_df, sample_indicator == 1)$id)
  phase1_data$R <- ifelse(phase1_data$cluster %in% phase2_clusters, 1, 0)
  
  # Inclusion probs
  allocation$pi_k <- allocation$stratum_size/allocation$npop
  phase1_data$pi_k <- dplyr::left_join(phase1_data, 
                                       dplyr::select(allocation,
                                                     "strata", "pi_k"), 
                                       by = c("strata" = "strata"))$pi_k
  
  # True data collected in phase 2
  phase2_data <- subset(phase1_data, R==1)
  phase2_data$Y <- Yki[phase1_data$R==1]
  phase2_data$X <- Xki[phase1_data$R==1]
  phase2_data$W <- Wki[phase1_data$R==1]
  phase2_data$V <- Vk[phase1_data$R==1]
  
  # Get vector of inclusion probs
  pi_k_df <- unique(dplyr::select(phase2_data, cluster, pi_k))
  pi_k <- pi_k_df$pi_k
  names(pi_k) <- pi_k_df$cluster
  
  # Generalized raking weights
  calibrated_weights <- compute_generalized_raking_weights(h_star, cluster_ids, phase1_data$R, pi_k)
  
  # Final GEE estimation using true data (phase 2) and calibrated weights
  gee_generalized_raking <- function(data, calibrated_weights){
    data$gk_inv <- 1 / calibrated_weights[as.character(data$cluster)]
    model <- suppressWarnings(geeglm(Y ~ X + W + V, family=binomial, data=data, id=cluster,
                                     weights=gk_inv, corstr="exchangeable"))
    return(summary(model))
  }
  
  # gee_generalized_raking_SIPW <- function(data, calibrated_weights){
  #   data$gk_inv <- calibrated_weights[as.character(data$cluster)]
  #   data <- data %>%
  #     dplyr::group_by(cluster) %>%
  #     dplyr::mutate(myvisit = row_number()) %>%
  #     dplyr::ungroup()
  #   gk_inv <- 1 / calibrated_weights[as.character(data$cluster)]
  #   model <- SIPW(Y ~ X + W + V, family='binomial', data=data, id='cluster', visit = 'myvisit',
  #                   weights='gk_inv', corstr="exchangeable")
  #   return(model)
  # }
  
  gee_generalized_raking_Survey <- function(data, calibrated_weights){
    data$gk_inv <- 1 / calibrated_weights[as.character(data$cluster)]
    gk_inv <- 1 / calibrated_weights[as.character(data$cluster)]
    data$cluster <- as.factor(data$cluster)
    mydesign <- svydesign(ids = ~cluster, probs = ~gk_inv, data = data,
                          calibration.formula = ~1)
    model <- svyglm(Y ~ X + W + V, family=quasibinomial, design = mydesign,
                    weights=gk_inv)
    return(summary(model))
  }
  
  final_results3 <- gee_generalized_raking_Survey(phase2_data, calibrated_weights)
  
  return(c(#final_results2["X",1],
           final_results3$coefficients["X", "Estimate"],
           final_results3$coefficients["X", "Std. Error"], length(cluster_ids)))
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

# Export necessary variables and functions to the cluster
clusterExport(cl, c("one_sim_GR_OALL", "inf_fun_gee",
                    "compute_generalized_raking_weights"))

# Run the simulations in parallel
sim_results_list2 <- parLapply(cl, 1:1000, function(i) {
  library(survey)
  library(MASS)
  library(geepack)
  library(dplyr)
  library(optimall)
  one_sim_GR_OALL(120, 5:7, -2, 0.7, 0.1, 0.8, 0.5, 1/3, seed = i)
})

# Stop the cluster
stopCluster(cl)

# Combine the results into a single data frame
sim_results2 <- as.data.frame(do.call(rbind, sim_results_list2))

mean(sim_results2[,1])
mean(sim_results2[,1])
var(sim_results2[,1])

# Coverage
sim_results2$lower <- sim_results2[,1] - 1.96*sim_results2[,2]
sim_results2$upper <- sim_results2[,1] + 1.96*sim_results2[,2]
sim_results2$cover <- ifelse(0.7 >= sim_results2$lower & 0.7 <= sim_results2$upper, 1, 0)
sum(sim_results2$cover)/1000

# Save
sim_results_GR_OALL <- sim_results2
#save(sim_results_GR_OALL, file = "/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/SecondYearCoursework/BIOST571/Final Project/sim_results_GR_OALL.rds")

