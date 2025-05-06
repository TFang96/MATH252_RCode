# Load Required Libraries
library(mlbench)
library(ggplot2)
library(dbscan)
library(cluster)
library(mclust)
library(reshape2)
library(readr)
library(fpc)
library(aricode)

# Load Datasets (all six)
data("Glass")
glass <- Glass
wine <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv", sep = ";")
magic_data <- read_csv("C:/Users/Maggie/Downloads/Spring2025/MATH252/project/magic04.data")
beans <- read.csv("C:/Users/Maggie/Downloads/Spring2025/MATH252/project/bean.csv", header = TRUE, stringsAsFactors = FALSE)
seeds <- read.table("C:/Users/Maggie/Downloads/Spring2025/MATH252/project/seeds_dataset.txt", header = FALSE)
wholesale_data <- read.csv("C:/Users/Maggie/Downloads/Spring2025/MATH252/project/wholesalecustomers.csv", header = TRUE, stringsAsFactors = FALSE)

# Mahalanobis Distance-Based Outlier Detection 
calculate_mahalanobis_outliers <- function(data_frame) {
  numeric_data <- data_frame[, sapply(data_frame, is.numeric)]
  n_dims <- ncol(numeric_data)
  scaled_data <- scale(numeric_data)
  
  # Use robust error handling for covariance matrix calculation
  cov_matrix <- tryCatch({
    cov(scaled_data)
  }, error = function(e) {
    diag(diag(cov(scaled_data))) + 0.001 * diag(n_dims)
  })
  
  # Calculate Mahalanobis distance and outlier flags
  mahal_dist <- mahalanobis(scaled_data, center = colMeans(scaled_data), cov = cov_matrix)
  mahal_outliers <- mahal_dist > qchisq(0.95, df = n_dims)
  
  # Return outlier flags, distances, and percentage
  list(
    mahal_outliers = mahal_outliers,
    mahal_outlier_percent = sum(mahal_outliers) / length(mahal_outliers) * 100,
    mahal_distances = mahal_dist
  )
}

# Prepare Datasets and Generate Outlier Plots
datasets_list <- list(
  Glass = glass,
  Wine = wine,
  Magic = magic_data,
  Beans = beans,
  Seeds = seeds,
  Wholesale = wholesale_data
)

processed_datasets <- list()

for(dataset_name in names(datasets_list)) {
  dataset <- datasets_list[[dataset_name]]
  numeric_data <- dataset[, sapply(dataset, is.numeric)]
  outlier_metrics <- calculate_mahalanobis_outliers(numeric_data)
  
  processed_datasets[[dataset_name]] <- list(
    data = dataset,
    numeric_data = numeric_data,
    outlier_metrics = outlier_metrics
  )
  
  # Summary output
  cat("\nDataset:", dataset_name)
  cat("\nDimensions:", nrow(numeric_data), "x", ncol(numeric_data))
  cat("\nMahalanobis outliers (", ncol(numeric_data), " dimensions): ", 
      round(outlier_metrics$mahal_outlier_percent, 2), "%\n")
  
  # Correlation heatmap
  cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")
  cor_melted <- melt(cor_matrix)
  heatmap_plot <- ggplot(data = cor_melted, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
    coord_fixed() +
    ggtitle(paste(dataset_name, "- Correlation Heatmap"))
  print(heatmap_plot)
  
  # Pairwise scatterplots with outliers highlighted
  par(mfrow = c(1, 1))
  max_vars <- min(6, ncol(numeric_data))
  selected_vars <- colnames(numeric_data)[1:max_vars]
  cat("Variables selected for pairs plot:", paste(selected_vars, collapse = ", "), "\n")
  pairs(numeric_data[, selected_vars, drop = FALSE],
        col = ifelse(outlier_metrics$mahal_outliers, "red", "blue"),
        pch = 19, cex = 0.5,
        gap = 0.5,
        cex.labels = 1.2,
        main = paste(dataset_name, "Dataset: Pairwise Relationships with Outliers Highlighted"))
}

# Data Preparation for Clustering 
prepare_for_clustering <- function(dataset_info) {
  x <- scale(dataset_info$numeric_data)
  x[is.na(x)] <- 0  # Replace NA with 0
  
  # Determine labels if present, otherwise derive them
  if("Type" %in% colnames(dataset_info$data)) {
    labels <- as.numeric(dataset_info$data$Type)
  } else if("Class" %in% colnames(dataset_info$data)) {
    labels <- as.numeric(as.factor(dataset_info$data$Class))
  } else {
    categorical_cols <- names(dataset_info$data)[sapply(dataset_info$data, is.factor) | sapply(dataset_info$data, is.character)]
    if(length(categorical_cols) > 0) {
      labels <- as.numeric(as.factor(dataset_info$data[[categorical_cols[1]]]))
      cat("Using column", categorical_cols[1], "as class labels\n")
    } else {
      set.seed(123)
      k_value <- min(max(3, ncol(x) / 2), 10)
      labels <- kmeans(x, centers = k_value)$cluster
      cat("Note: Using k-means clustering with k =", k_value, "to generate labels\n")
    }
  }
  
  # Regenerate labels if only one unique label
  if(length(unique(labels)) < 2) {
    set.seed(123)
    labels <- kmeans(x, centers = 3)$cluster
    cat("Warning: Had to regenerate labels using k-means with k=3\n")
  }
  
  list(
    x = x,
    labels = labels,
    outlier_density = dataset_info$outlier_metrics$mahal_outlier_percent
  )
}

datasets_for_clustering <- lapply(processed_datasets, prepare_for_clustering)

# Implementation of a custom contaminated normal mixture model
fit_contaminated_normal <- function(x, k) {
  n <- nrow(x)
  d <- ncol(x)
  
  cat("Using custom contaminated normal mixture model implementation\n")
  
  # Try to fit the model
  tryCatch({
    # For large datasets, take a random subsample for model fitting
    if(n > 5000) {
      cat("Large dataset detected, using a subsample of 5000 observations for model fitting\n")
      set.seed(123)
      subsample_idx <- sample(1:n, 5000)
      x_sample <- x[subsample_idx, ]
      n_sample <- 5000
    } else {
      x_sample <- x
      n_sample <- n
      subsample_idx <- 1:n
    }
    
    # Set random seed for reproducibility
    set.seed(123)
    
    # Step 1: Initialize with k-means
    cat("Initializing with k-means clustering\n")
    km <- kmeans(x_sample, centers = k, nstart = 10)
    
    # Step 2: Prepare data structures for EM algorithm
    # Initial mixing proportions
    pi_g <- km$size / n_sample
    
    # Initial mean vectors
    mu_g <- km$centers
    
    # Initial covariance matrices (one per cluster)
    sigma_g <- array(0, dim = c(d, d, k))
    for (g in 1:k) {
      cluster_points <- x_sample[km$cluster == g, , drop = FALSE]
      if (nrow(cluster_points) > 1) {
        sigma_g[, , g] <- cov(cluster_points)
      } else {
        sigma_g[, , g] <- diag(d)
      }
      
      # Ensure covariance matrix is positive definite
      eigen_values <- eigen(sigma_g[, , g])$values
      if (any(eigen_values <= 1e-10)) {
        sigma_g[, , g] <- sigma_g[, , g] + diag(0.01, d)
      }
    }
    
    # Initial contamination parameters
    # alpha_g: weight of good component in mixture (0.5-0.99)
    # eta_g: inflation factor for covariance in bad component (>1)
    alpha_g <- rep(0.95, k)  # Start with 95% of data in good component
    eta_g <- rep(3.0, k)     # Bad component has 3x covariance
    
    # Step 3: EM algorithm for contaminated normal mixture
    max_iter <- 50  # Reduced for large datasets
    epsilon <- 1e-5  # Relaxed convergence threshold
    log_lik_old <- -Inf
    
    # Data structures for posterior probabilities
    z_ig <- matrix(0, nrow = n_sample, ncol = k)      # Cluster memberships
    v_ig <- matrix(0, nrow = n_sample, ncol = k)      # Component (good/bad) memberships
    
    # For numerical stability in log calculations
    tiny <- .Machine$double.eps * 100
    
    for (iter in 1:max_iter) {
      # E-step: Calculate posterior probabilities
      
      # First, calculate densities for each observation in each cluster and component
      dens_good <- matrix(0, nrow = n_sample, ncol = k)  # Density in good component
      dens_bad <- matrix(0, nrow = n_sample, ncol = k)   # Density in bad component
      
      for (g in 1:k) {
        # Calculate determinants and inverses once per iteration
        sigma_g_inv <- try(solve(sigma_g[, , g]), silent = TRUE)
        if (inherits(sigma_g_inv, "try-error")) {
          # If matrix is singular, add a small regularization
          sigma_g[, , g] <- sigma_g[, , g] + diag(0.01, d)
          sigma_g_inv <- solve(sigma_g[, , g])
        }
        
        sigma_g_det <- det(sigma_g[, , g])
        if (sigma_g_det <= 0) {
          # If determinant is zero or negative, regularize
          sigma_g[, , g] <- sigma_g[, , g] + diag(0.01, d)
          sigma_g_det <- det(sigma_g[, , g])
          sigma_g_inv <- solve(sigma_g[, , g])
        }
        
        # Determinant of the bad component's covariance
        sigma_g_bad_det <- sigma_g_det * (eta_g[g]^d)
        
        # Calculate Mahalanobis distances
        for (i in 1:n_sample) {
          # Center the observation
          centered <- x_sample[i, ] - mu_g[g, ]
          
          # Mahalanobis distance (squared)
          mahal_dist <- sum((centered %*% sigma_g_inv) * centered)
          
          # Density in good component
          dens_good[i, g] <- (2 * pi)^(-d/2) * sigma_g_det^(-0.5) * 
            exp(-0.5 * mahal_dist)
          
          # Density in bad component (inflated covariance)
          mahal_dist_bad <- mahal_dist / eta_g[g]
          dens_bad[i, g] <- (2 * pi)^(-d/2) * sigma_g_bad_det^(-0.5) * 
            exp(-0.5 * mahal_dist_bad)
        }
      }
      
      # Calculate overall posterior probabilities
      for (i in 1:n_sample) {
        # Denominator for z_ig (overall probability across all clusters and components)
        denom_z <- 0
        for (g in 1:k) {
          denom_z <- denom_z + pi_g[g] * (alpha_g[g] * dens_good[i, g] + 
                                            (1 - alpha_g[g]) * dens_bad[i, g])
        }
        
        # Calculate z_ig (cluster membership)
        for (g in 1:k) {
          z_ig[i, g] <- pi_g[g] * (alpha_g[g] * dens_good[i, g] + 
                                     (1 - alpha_g[g]) * dens_bad[i, g]) / 
            (denom_z + tiny)
        }
        
        # Calculate v_ig (component membership within cluster)
        for (g in 1:k) {
          numer_v <- alpha_g[g] * dens_good[i, g]
          denom_v <- alpha_g[g] * dens_good[i, g] + (1 - alpha_g[g]) * dens_bad[i, g]
          v_ig[i, g] <- numer_v / (denom_v + tiny)
        }
      }
      
      # Calculate log-likelihood
      log_lik <- 0
      for (i in 1:n_sample) {
        sum_i <- 0
        for (g in 1:k) {
          sum_i <- sum_i + pi_g[g] * (alpha_g[g] * dens_good[i, g] + 
                                        (1 - alpha_g[g]) * dens_bad[i, g])
        }
        log_lik <- log_lik + log(sum_i + tiny)
      }
      
      # Check for convergence
      if (abs(log_lik - log_lik_old) < epsilon) {
        cat("EM algorithm converged after", iter, "iterations\n")
        break
      }
      log_lik_old <- log_lik
      
      # M-step: Update parameters
      
      # Update mixing proportions
      n_g <- colSums(z_ig)
      pi_g <- n_g / n_sample
      
      # Update means and covariances
      for (g in 1:k) {
        # Weighted observations for this component
        z_v_ig <- z_ig[, g] * v_ig[, g]
        sum_z_v <- sum(z_v_ig)
        
        if (sum_z_v > tiny) {
          # Update mean
          mu_g[g, ] <- colSums(z_v_ig * x_sample) / sum_z_v
          
          # Update covariance
          sigma_g[, , g] <- matrix(0, d, d)
          for (i in 1:n_sample) {
            centered <- x_sample[i, ] - mu_g[g, ]
            sigma_g[, , g] <- sigma_g[, , g] + z_v_ig[i] * 
              (centered %*% t(centered))
          }
          sigma_g[, , g] <- sigma_g[, , g] / sum_z_v
          
          # Ensure covariance matrix is positive definite
          eigen_values <- eigen(sigma_g[, , g])$values
          if (any(eigen_values <= 1e-10)) {
            sigma_g[, , g] <- sigma_g[, , g] + diag(0.01, d)
          }
        }
        
        # Update alpha (weight of good component)
        alpha_g[g] <- sum(z_ig[, g] * v_ig[, g]) / n_g[g]
        
        # Bound alpha to avoid degenerate solutions
        alpha_g[g] <- max(min(alpha_g[g], 0.99), 0.5)
      }
      
      # No direct update for eta in standard EM
      # Typically requires numerical optimization
      # For simplicity, we keep eta fixed
    }
    
    if (iter == max_iter) {
      cat("EM algorithm reached maximum iterations without convergence\n")
    }
    
    # Step 4: If we used a subsample, extend results to full dataset
    if (n > 5000) {
      # For subsample, get cluster assignments and component memberships
      sample_classifications <- apply(z_ig, 1, which.max)
      sample_v_ig <- v_ig
      
      # Initialize full dataset structures
      classifications <- rep(0, n)
      v_ig_full <- matrix(0, nrow = n, ncol = k)
      z_ig_full <- matrix(0, nrow = n, ncol = k)
      
      # First, assign subsample points
      classifications[subsample_idx] <- sample_classifications
      v_ig_full[subsample_idx, ] <- sample_v_ig
      z_ig_full[subsample_idx, ] <- z_ig
      
      # For non-sampled points, calculate posterior based on fitted model
      non_sample_idx <- setdiff(1:n, subsample_idx)
      
      # Initialize matrices for densities
      dens_good_full <- matrix(0, nrow = length(non_sample_idx), ncol = k)
      dens_bad_full <- matrix(0, nrow = length(non_sample_idx), ncol = k)
      
      for (g in 1:k) {
        # Get inverse of covariance matrix
        sigma_g_inv <- solve(sigma_g[, , g])
        sigma_g_det <- det(sigma_g[, , g])
        sigma_g_bad_det <- sigma_g_det * (eta_g[g]^d)
        
        # Calculate densities for each non-sampled point
        for (i_idx in 1:length(non_sample_idx)) {
          i <- non_sample_idx[i_idx]
          centered <- x[i, ] - mu_g[g, ]
          
          # Mahalanobis distance (squared)
          mahal_dist <- sum((centered %*% sigma_g_inv) * centered)
          
          # Density in good component
          dens_good_full[i_idx, g] <- (2 * pi)^(-d/2) * sigma_g_det^(-0.5) * 
            exp(-0.5 * mahal_dist)
          
          # Density in bad component (inflated covariance)
          mahal_dist_bad <- mahal_dist / eta_g[g]
          dens_bad_full[i_idx, g] <- (2 * pi)^(-d/2) * sigma_g_bad_det^(-0.5) * 
            exp(-0.5 * mahal_dist_bad)
        }
      }
      
      # Calculate posterior probabilities for non-sampled points
      for (i_idx in 1:length(non_sample_idx)) {
        i <- non_sample_idx[i_idx]
        
        # Denominator for z_ig
        denom_z <- 0
        for (g in 1:k) {
          denom_z <- denom_z + pi_g[g] * (alpha_g[g] * dens_good_full[i_idx, g] + 
                                            (1 - alpha_g[g]) * dens_bad_full[i_idx, g])
        }
        
        # Calculate z_ig (cluster membership)
        for (g in 1:k) {
          z_ig_full[i, g] <- pi_g[g] * (alpha_g[g] * dens_good_full[i_idx, g] + 
                                          (1 - alpha_g[g]) * dens_bad_full[i_idx, g]) / 
            (denom_z + tiny)
        }
        
        # Calculate v_ig (component membership within cluster)
        for (g in 1:k) {
          numer_v <- alpha_g[g] * dens_good_full[i_idx, g]
          denom_v <- alpha_g[g] * dens_good_full[i_idx, g] + (1 - alpha_g[g]) * dens_bad_full[i_idx, g]
          v_ig_full[i, g] <- numer_v / (denom_v + tiny)
        }
      }
      
      # Assign cluster based on maximum posterior probability
      for (i in non_sample_idx) {
        classifications[i] <- which.max(z_ig_full[i, ])
      }
      
      # Use full dataset matrices
      z_ig <- z_ig_full
      v_ig <- v_ig_full
      
    } else {
      # For regular dataset, just take cluster with highest posterior
      classifications <- apply(z_ig, 1, which.max)
    }
    
    # Identify outliers based on component membership probabilities
    # If probability of being in the "good" component is low, it's an outlier
    outliers <- rep(FALSE, n)
    for (i in 1:n) {
      g <- classifications[i]
      if (v_ig[i, g] < 0.5) {  # More likely to be in bad component than good
        outliers[i] <- TRUE
      }
    }
    
    # Calculate cluster centers (means)
    centers <- mu_g
    
    # Return results
    return(list(
      classification = classifications,
      posterior = z_ig,                         # Posterior cluster probabilities
      component_posterior = v_ig,               # Posterior component probabilities
      mahalanobis_distances = rep(0, n),        # Not explicitly calculated
      outliers = outliers,
      n_outliers = sum(outliers),
      outlier_percent = sum(outliers) / n * 100,
      centers = centers,
      covariances = sigma_g,
      mixing_props = pi_g,
      alpha = alpha_g,                         # Weight of good component
      eta = eta_g,                             # Inflation factor for bad component
      log_likelihood = log_lik,
      success = TRUE
    ))
  }, error = function(e) {
    cat("Error in custom contaminated normal mixture model:", e$message, "\n")
    cat("Falling back to mclust model\n")
    
    # Use a subsample for large datasets with mclust too
    if(n > 5000) {
      cat("Large dataset detected, using a subsample of 5000 observations for mclust\n")
      set.seed(123)
      subsample_idx <- sample(1:n, 5000)
      x_sample <- x[subsample_idx, ]
    } else {
      x_sample <- x
      subsample_idx <- 1:n
    }
    
    # Fallback to mclust if custom implementation fails
    library(mclust)
    mclust_result <- Mclust(x_sample, G = k)
    
    # For subsample, extend results to full dataset
    if(n > 5000) {
      # Initialize classification and posterior for full dataset
      classification <- rep(0, n)
      
      # First, assign subsample points
      classification[subsample_idx] <- mclust_result$classification
      
      # For non-sampled points, find nearest cluster center
      non_sample_idx <- setdiff(1:n, subsample_idx)
      
      # Calculate distances to cluster centers
      centers <- t(mclust_result$parameters$mean)
      for(i in non_sample_idx) {
        distances <- apply(centers, 1, function(center) {
          sum((x[i, ] - center)^2)  # Euclidean distance squared
        })
        classification[i] <- which.min(distances)
      }
      
      # Create simple posterior matrix for full dataset
      posterior <- matrix(0, nrow = n, ncol = k)
      for(i in 1:n) {
        posterior[i, classification[i]] <- 1
      }
    } else {
      classification <- mclust_result$classification
      posterior <- mclust_result$z
    }
    
    # Calculate outliers based on likelihood
    # Points with low likelihood are potential outliers
    if(n <= 5000) {
      density_values <- mclust_result$z %*% mclust_result$parameters$pro
      outlier_threshold <- quantile(density_values, 0.05)  # Bottom 5%
      outliers <- density_values < outlier_threshold
    } else {
      # For large datasets, use simple distance-based outlier detection
      outliers <- rep(FALSE, n)
      for(i in 1:n) {
        # Calculate distance to assigned cluster center
        g <- classification[i]
        if(g <= k) {  # Make sure cluster index is valid
          center <- centers[g, ]
          dist_squared <- sum((x[i, ] - center)^2)
          
          # Mark as outlier if distance is large
          # This is a simple approximation when we can't calculate full likelihoods
          outliers[i] <- dist_squared > median(dist_squared) * 3
        }
      }
    }
    
    return(list(
      classification = classification,
      posterior = posterior,
      mahalanobis_distances = rep(0, n),
      outliers = outliers,
      n_outliers = sum(outliers),
      outlier_percent = sum(outliers) / n * 100,
      centers = if(n <= 5000) t(mclust_result$parameters$mean) else centers,
      success = TRUE,
      model_type = "mclust_fallback"
    ))
  })
}

# Clustering Benchmark Function 
benchmark_clust <- function(x, true_lbls) {
  # Ensure all required packages are loaded
  required_packages <- c("dbscan", "cluster", "aricode", "fpc")
  
  for(pkg in required_packages) {
    if(!require(pkg, character.only = TRUE)) {
      cat("Installing package:", pkg, "\n")
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  # Check specific functions
  if(!exists("ARI", mode = "function")) {
    cat("Loading ARI from aricode package\n")
    if(require("aricode")) {
      ARI <- aricode::ARI
    } else {
      stop("Cannot load ARI function. Please install the aricode package.")
    }
  }
  
  if(!exists("NMI", mode = "function")) {
    cat("Loading NMI from aricode package\n")
    if(require("aricode")) {
      NMI <- aricode::NMI
    } else {
      stop("Cannot load NMI function. Please install the aricode package.")
    }
  }
  
  res <- list()
  start_time <- Sys.time()
  
  # For large datasets, use a simpler parameter search to save time
  n <- nrow(x)
  if(n > 5000) {
    cat("Large dataset detected, using simplified parameter search for DBSCAN/OPTICS\n")
    eps_values <- seq(0.3, 1.0, by = 0.2)
    min_pts_values <- c(5, 10)
  } else {
    eps_values <- seq(0.3, 1.5, by = 0.1)
    min_pts_values <- c(3, 4, 5, 10)
  }
  
  # DBSCAN Tuning and Evaluation
  best_eps <- 0.5
  best_ari <- -1
  
  for (min_pts in min_pts_values) {
    for (eps in eps_values) {
      db_temp <- dbscan::dbscan(x, eps = eps, minPts = min_pts)
      if (length(unique(db_temp$cluster)) > 1) {
        ari_temp <- aricode::ARI(db_temp$cluster, true_lbls)  # Use with namespace
        if (ari_temp > best_ari) {
          best_ari <- ari_temp
          best_eps <- eps
          best_min_pts <- min_pts
        }
      }
    }
  }
  
  db <- dbscan::dbscan(x, eps = best_eps, minPts = best_min_pts)
  db_time <- difftime(Sys.time(), start_time, units = "secs")
  db_clusters <- db$cluster
  valid_db <- db_clusters > 0
  
  # Calculate silhouette only for smaller datasets
  if(n <= 5000 && length(unique(db_clusters[valid_db])) > 1) {
    avg_sil_db <- mean(cluster::silhouette(db_clusters[valid_db], dist(x[valid_db, ]))[, 3])
  } else {
    avg_sil_db <- NA
  }
  
  res$DBSCAN_ARI <- aricode::ARI(db_clusters, true_lbls)
  res$DBSCAN_NMI <- aricode::NMI(true_lbls, db_clusters)
  res$DBSCAN_params <- paste("eps =", best_eps, ", minPts =", best_min_pts)
  res$DBSCAN_time <- as.numeric(db_time)
  res$DBSCAN_n_clusters <- length(unique(db_clusters[db_clusters > 0]))
  res$DBSCAN_noise <- sum(db_clusters == 0) / length(db_clusters) * 100
  res$DBSCAN_silhouette <- avg_sil_db
  
  # OPTICS Clustering (only for datasets smaller than 10000 rows)
  if(n <= 10000) {
    start_time <- Sys.time()
    opt <- dbscan::optics(x, eps = 10, minPts = best_min_pts)
    
    if(n > 5000) {
      eps_cl_values <- seq(0.2, 1.0, by = 0.2)
    } else {
      eps_cl_values <- seq(0.1, 2.0, by = 0.1)
    }
    
    best_eps_cl <- 0.5
    best_ari <- -1
    for (eps_cl in eps_cl_values) {
      cl_opt_temp <- dbscan::extractDBSCAN(opt, eps_cl = eps_cl)$cluster
      if (length(unique(cl_opt_temp)) > 1) {
        ari_temp <- aricode::ARI(cl_opt_temp, true_lbls)
        if (ari_temp > best_ari) {
          best_ari <- ari_temp
          best_eps_cl <- eps_cl
        }
      }
    }
    cl_opt <- dbscan::extractDBSCAN(opt, eps_cl = best_eps_cl)$cluster
    opt_time <- difftime(Sys.time(), start_time, units = "secs")
    valid_opt <- cl_opt > 0
    
    # Calculate silhouette only for smaller datasets
    if(n <= 5000 && length(unique(cl_opt[valid_opt])) > 1) {
      avg_sil_opt <- mean(cluster::silhouette(cl_opt[valid_opt], dist(x[valid_opt, ]))[, 3])
    } else {
      avg_sil_opt <- NA
    }
    
    res$OPTICS_ARI <- aricode::ARI(cl_opt, true_lbls)
    res$OPTICS_NMI <- aricode::NMI(true_lbls, cl_opt)
    res$OPTICS_params <- paste("minPts =", best_min_pts, ", eps_cl =", best_eps_cl)
    res$OPTICS_time <- as.numeric(opt_time)
    res$OPTICS_n_clusters <- length(unique(cl_opt[cl_opt > 0]))
    res$OPTICS_noise <- sum(cl_opt == 0) / length(cl_opt) * 100
    res$OPTICS_silhouette <- avg_sil_opt
  } else {
    # Skip OPTICS for very large datasets
    cat("Skipping OPTICS for large dataset (n >", n, ")\n")
    res$OPTICS_ARI <- NA
    res$OPTICS_NMI <- NA
    res$OPTICS_params <- NA
    res$OPTICS_time <- NA
    res$OPTICS_n_clusters <- NA
    res$OPTICS_noise <- NA
    res$OPTICS_silhouette <- NA
  }
  
  # PAM Clustering (only for datasets smaller than 10000 rows)
  if(n <= 10000) {
    start_time <- Sys.time()
    k <- length(unique(true_lbls))
    pam_result <- cluster::pam(x, k, pamonce = TRUE)  # Using faster implementation
    pam_time <- difftime(Sys.time(), start_time, units = "secs")
    res$PAM_ARI <- aricode::ARI(pam_result$clustering, true_lbls)
    res$PAM_NMI <- aricode::NMI(true_lbls, pam_result$clustering)
    res$PAM_params <- paste("k =", k)
    res$PAM_time <- as.numeric(pam_time)
    res$PAM_silhouette <- mean(pam_result$silinfo$avg.width)
  } else {
    # Skip PAM for very large datasets
    cat("Skipping PAM for large dataset (n >", n, ")\n")
    res$PAM_ARI <- NA
    res$PAM_NMI <- NA
    res$PAM_params <- NA
    res$PAM_time <- NA
    res$PAM_silhouette <- NA
  }
  
  # CN Clustering using custom implementation
  start_time <- Sys.time()
  k <- length(unique(true_lbls))
  cn_model <- fit_contaminated_normal(x, k)
  cn_clusters <- cn_model$classification
  cn_time <- difftime(Sys.time(), start_time, units = "secs")
  
  # Calculate silhouette only for smaller datasets
  if(n <= 5000 && length(unique(cn_clusters)) > 1) {
    avg_sil_cn <- mean(cluster::silhouette(cn_clusters, dist(x))[, 3])
  } else {
    avg_sil_cn <- NA
  }
  
  res$CN_ARI <- aricode::ARI(cn_clusters, true_lbls)
  res$CN_NMI <- aricode::NMI(true_lbls, cn_clusters)
  res$CN_params <- paste("k =", k)
  res$CN_time <- as.numeric(cn_time)
  res$CN_n_clusters <- length(unique(cn_clusters))
  res$CN_BIC <- NA  # BIC not directly provided by our implementation
  res$CN_silhouette <- avg_sil_cn
  
  return(res)
}

# Run Benchmarking Across all Datasets
results <- lapply(datasets_for_clustering, function(d) {
  cat("\nBenchmarking dataset with outlier density:", round(d$outlier_density, 2), "%\n")
  benchmark_clust(d$x, d$labels)
})

# Compile Results into a Data Frame
results_df <- do.call(rbind, lapply(names(results), function(dataset_name) {
  result <- results[[dataset_name]]
  df <- data.frame(Dataset = dataset_name,
                   Outlier_Density = datasets_for_clustering[[dataset_name]]$outlier_density,
                   stringsAsFactors = FALSE)
  for (name in names(result)) {
    df[[name]] <- result[[name]]
  }
  df
}))

print(results_df)

# Prepare and Print Comparison Tables
comparison_scores <- data.frame(
  Dataset = results_df$Dataset,
  DBSCAN_ARI = results_df$DBSCAN_ARI,
  OPTICS_ARI = results_df$OPTICS_ARI,
  PAM_ARI = results_df$PAM_ARI,
  CN_ARI = results_df$CN_ARI,
  DBSCAN_NMI = results_df$DBSCAN_NMI,
  OPTICS_NMI = results_df$OPTICS_NMI,
  PAM_NMI = results_df$PAM_NMI,
  CN_NMI = results_df$CN_NMI,
  DBSCAN_Silhouette = results_df$DBSCAN_silhouette,
  OPTICS_Silhouette = results_df$OPTICS_silhouette,
  PAM_Silhouette = results_df$PAM_silhouette,
  CN_Silhouette = results_df$CN_silhouette
)

cat("Adjusted Rand Index (ARI), Normalized Mutual Information (NMI), and Silhouette Width Comparison:\n")
print(comparison_scores)

# Summary of Best Performing Algorithms by Metric
cat("\n\nClustering Algorithm Summary\n")
cat("Dataset outlier characteristics:\n")
for(ds in names(datasets_for_clustering)) {
  cat(paste0("- ", ds, ": ", round(datasets_for_clustering[[ds]]$outlier_density, 2), "% outlier density\n"))
}

# Best algorithm by ARI
cat("\nBest performing algorithm by dataset (based on ARI):\n")
ari_scores <- comparison_scores
for(i in 1:nrow(ari_scores)) {
  algo_scores <- as.numeric(ari_scores[i, 2:5])
  algo_names <- colnames(ari_scores)[2:5]
  valid_indices <- which(!is.na(algo_scores))
  if(length(valid_indices) > 0) {
    best_algo <- algo_names[which.max(algo_scores[valid_indices])]
    best_score <- max(algo_scores[valid_indices], na.rm = TRUE)
    cat(paste0("- ", ari_scores$Dataset[i], ": ", best_algo, 
               " (ARI = ", round(best_score, 3), ")\n"))
  } else {
    cat(paste0("- ", ari_scores$Dataset[i], ": No valid ARI scores available\n"))
  }
}

# Best algorithm by NMI
cat("\nBest performing algorithm by dataset (based on NMI):\n")
nmi_scores <- comparison_scores
for(i in 1:nrow(nmi_scores)) {
  algo_scores <- as.numeric(ari_scores[i, 6:9])
  algo_names <- colnames(nmi_scores)[6:9]
  valid_indices <- which(!is.na(algo_scores))
  if(length(valid_indices) > 0) {
    best_algo <- algo_names[which.max(algo_scores[valid_indices])]
    best_score <- max(algo_scores[valid_indices], na.rm = TRUE)
    cat(paste0("- ", nmi_scores$Dataset[i], ": ", best_algo, 
               " (NMI = ", round(best_score, 3), ")\n"))
  } else {
    cat(paste0("- ", nmi_scores$Dataset[i], ": No valid NMI scores available\n"))
  }
}

# Best algorithm by Silhouette Score
cat("\nBest performing algorithm by dataset (based on Silhouette Score):\n")
sil_scores <- comparison_scores
for(i in 1:nrow(sil_scores)) {
  algo_scores <- as.numeric(sil_scores[i, 10:13])
  algo_names <- colnames(sil_scores)[10:13]
  valid_indices <- which(!is.na(algo_scores))
  if(length(valid_indices) > 0) {
    best_algo <- algo_names[which.max(algo_scores[valid_indices])]
    best_score <- max(algo_scores[valid_indices], na.rm = TRUE)
    cat(paste0("- ", sil_scores$Dataset[i], ": ", best_algo, 
               " (Silhouette = ", round(best_score, 3), ")\n"))
  } else {
    cat(paste0("- ", sil_scores$Dataset[i], ": No valid Silhouette scores available\n"))
  }
}

# Generate Summary Tables
# Create a table of ARI scores
ari_table <- data.frame(
  Dataset = results_df$Dataset,
  Outlier_Density = round(results_df$Outlier_Density, 2),
  DBSCAN = round(results_df$DBSCAN_ARI, 3),
  OPTICS = round(results_df$OPTICS_ARI, 3),
  PAM = round(results_df$PAM_ARI, 3),
  CN = round(results_df$CN_ARI, 3)
)

# Create a table of NMI scores
nmi_table <- data.frame(
  Dataset = results_df$Dataset,
  Outlier_Density = round(results_df$Outlier_Density, 2),
  DBSCAN = round(results_df$DBSCAN_NMI, 3),
  OPTICS = round(results_df$OPTICS_NMI, 3),
  PAM = round(results_df$PAM_NMI, 3),
  CN = round(results_df$CN_NMI, 3)
)

# Create a table of Silhouette scores
sil_table <- data.frame(
  Dataset = results_df$Dataset,
  Outlier_Density = round(results_df$Outlier_Density, 2),
  DBSCAN = round(results_df$DBSCAN_silhouette, 3),
  OPTICS = round(results_df$OPTICS_silhouette, 3),
  PAM = round(results_df$PAM_silhouette, 3),
  CN = round(results_df$CN_silhouette, 3)
)

cat("\n\nARI Scores by Dataset and Algorithm:\n")
print(ari_table, row.names = FALSE)

cat("\n\nNMI Scores by Dataset and Algorithm:\n")
print(nmi_table, row.names = FALSE)

cat("\n\nSilhouette Scores by Dataset and Algorithm:\n")
print(sil_table, row.names = FALSE)

# Create visualizations if ggplot2 is available
if(require(ggplot2)) {
  # Convert results to long format for plotting
  ari_long <- reshape2::melt(ari_table, 
                             id.vars = c("Dataset", "Outlier_Density"),
                             variable.name = "Algorithm", 
                             value.name = "ARI")
  
  # Plot ARI scores by algorithm and dataset
  ari_plot <- ggplot(ari_long, aes(x = Dataset, y = ARI, fill = Algorithm)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_minimal() +
    labs(title = "ARI Scores by Dataset and Algorithm",
         subtitle = "Higher values indicate better clustering performance",
         y = "Adjusted Rand Index") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(ari_plot)
  
  # Create a plot comparing performance against outlier density
  outlier_plot <- ggplot(ari_long, aes(x = Outlier_Density, y = ARI, color = Algorithm)) +
    geom_point(size = 3) +
    geom_smooth(method = "loess", se = FALSE) +
    theme_minimal() +
    labs(title = "Algorithm Performance vs. Outlier Density",
         x = "Outlier Density (%)",
         y = "Adjusted Rand Index")
  
  print(outlier_plot)
}

