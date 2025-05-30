library(mlbench)
library(ggplot2)
library(dbscan)
library(cluster)
library(mclust)       
library(reshape2)    
library(gridExtra)
library(MASS)       # For mvrnorm
library(readr) # For read_csv

# Loading datasets - using actual standard datasets when possible

# Glass dataset
data("Glass")
glass <- Glass

# Wine dataset 
wine <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv", sep = ";")

# Magic dataset
magic_data <- read_csv("C:/Users/tf245/Documents/GitHub/MATH252_RCode/Project/magic+gamma+telescope/magic04.data")

# Beans dataset 
beans <- read.csv("C:/Users/tf245/Documents/GitHub/MATH252_RCode/Project/bean.csv", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

# Seeds dataset
seeds <- read.table("seeds_dataset.txt", 
                    header = FALSE)

# Wholesale dataset
wholesale_data <- read.csv("Wholesale customers data.csv", 
                           header = TRUE, 
                           stringsAsFactors = FALSE)

# Function to calculate Mahalanobis distance outliers
calculate_mahalanobis_outliers <- function(data_frame) {
  # Only for numeric columns
  numeric_data <- data_frame[, sapply(data_frame, is.numeric)]
  n_dims <- ncol(numeric_data)
  
  scaled_data <- scale(numeric_data)
  
  cov_matrix <- tryCatch({
    cov(scaled_data)
  }, error = function(e) {
    diag(diag(cov(scaled_data))) + 0.001 * diag(n_dims)
  })
  
  mahal_dist <- mahalanobis(scaled_data, 
                            center = colMeans(scaled_data), 
                            cov = cov_matrix)
  
  mahal_outliers <- mahal_dist > qchisq(0.95, df = n_dims)
  mahal_outlier_percent <- sum(mahal_outliers) / length(mahal_outliers) * 100
  
  metrics <- list()
  metrics$mahal_outliers <- mahal_outliers
  metrics$mahal_outlier_percent <- mahal_outlier_percent
  metrics$mahal_distances <- mahal_dist
  
  return(metrics)
}

# Process each dataset
datasets_list <- list(
  Glass = glass,
  Wine = wine,
  Magic = magic_data,
  Beans = beans,
  Seeds = seeds,
  Wholesale = wholesale_data
)

# Pre-process datasets and calculate outlier metrics
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
  
  cat("\nDataset:", dataset_name)
  cat("\nDimensions:", nrow(numeric_data), "x", ncol(numeric_data))
  cat("\nMahalanobis outliers (", ncol(numeric_data), " dimensions): ", 
      round(outlier_metrics$mahal_outlier_percent, 2), "%\n")
  
  cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")
  cor_melted <- melt(cor_matrix)
  
  heatmap_plot <- ggplot(data = cor_melted, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1)) +
    coord_fixed() +
    ggtitle(paste(dataset_name, "- Correlation Heatmap"))
  
  print(heatmap_plot)

  par(mfrow = c(1, 1))
  max_vars <- min(6, ncol(numeric_data))
  selected_vars <- colnames(numeric_data)[1:max_vars]
  
  cat("Variables selected for pairs plot:", paste(selected_vars, collapse=", "), "\n")
  
  pairs(numeric_data[, selected_vars, drop = FALSE], 
        col = ifelse(outlier_metrics$mahal_outliers, "red", "blue"),
        pch = 19, cex = 0.5,
        gap = 0.5,
        cex.labels = 1.2,
        main = paste(dataset_name, "Dataset: Pairwise Relationships with Outliers Highlighted"))
}



# Prepare for clustering ------------------------------------------
# Function to prepare data for clustering
prepare_for_clustering <- function(dataset_info) {
  # Standardize all numeric features
  x <- scale(dataset_info$numeric_data)
  
  # Handle potential NA values from scaling
  x[is.na(x)] <- 0
  
  # Extract or create labels (if available)
  if("Type" %in% colnames(dataset_info$data)) {
    labels <- as.numeric(dataset_info$data$Type)
  } else if("Class" %in% colnames(dataset_info$data)) {
    labels <- as.numeric(as.factor(dataset_info$data$Class))
  } else {
    # Try to find any categorical column that could serve as labels
    categorical_cols <- names(dataset_info$data)[sapply(dataset_info$data, is.factor) | 
                                                   sapply(dataset_info$data, is.character)]
    
    if(length(categorical_cols) > 0) {
      # Use the first categorical column as labels
      labels <- as.numeric(as.factor(dataset_info$data[[categorical_cols[1]]]))
      cat("Using column", categorical_cols[1], "as class labels\n")
    } else {
      # Create dummy labels if no class information is available
      set.seed(123)  # For reproducibility
      # Use a larger number of clusters for potentially more complex datasets
      k_value <- min(max(3, ncol(x) / 2), 10)  # Between 3 and 10 based on dimensionality
      labels <- kmeans(x, centers = k_value)$cluster
      cat("Note: Using k-means clustering with k=", k_value, "to generate labels\n")
    }
  }
  
  # Make sure we have reasonable number of clusters (at least 2)
  if(length(unique(labels)) < 2) {
    set.seed(123)
    labels <- kmeans(x, centers = 3)$cluster
    cat("Warning: Had to regenerate labels using k-means with k=3\n")
  }
  
  return(list(
    x = x, 
    labels = labels, 
    outlier_density = dataset_info$outlier_metrics$mahal_outlier_percent
  ))
}

# Prepare all datasets for clustering
datasets_for_clustering <- lapply(processed_datasets, prepare_for_clustering)

# Benchmark function -------------------------------------
benchmark_clust <- function(x, true_lbls) {
  res <- list()
  
  # Timer for performance measurement
  start_time <- Sys.time()
  
  # DBSCAN with parameter tuning
  # Try a range of eps values to find optimal
  eps_values <- seq(0.3, 1.5, by = 0.1)  # Wider range of eps values
  best_eps <- 0.5
  best_ari <- -1
  
  # Try different minPts values for more robust clustering
  min_pts_values <- c(3, 4, 5, 10)
  
  for (min_pts in min_pts_values) {
    for (eps in eps_values) {
      db_temp <- dbscan(x, eps = eps, minPts = min_pts)
      # Skip if all points are noise (cluster = 0)
      if (length(unique(db_temp$cluster)) > 1) {
        ari_temp <- adjustedRandIndex(db_temp$cluster, true_lbls)
        if (ari_temp > best_ari) {
          best_ari <- ari_temp
          best_eps <- eps
          best_min_pts <- min_pts
        }
      }
    }
  }
  
  # Use best eps and minPts for final DBSCAN
  db <- dbscan(x, eps = best_eps, minPts = best_min_pts)
  db_time <- difftime(Sys.time(), start_time, units = "secs")
  
  res$DBSCAN_ARI <- adjustedRandIndex(db$cluster, true_lbls)
  res$DBSCAN_params <- paste("eps =", best_eps, ", minPts =", best_min_pts)
  res$DBSCAN_time <- as.numeric(db_time)
  res$DBSCAN_n_clusters <- length(unique(db$cluster[db$cluster > 0]))
  res$DBSCAN_noise <- sum(db$cluster == 0) / length(db$cluster) * 100  # % noise points
  
  # OPTICS with parameter tuning
  start_time <- Sys.time()
  # Use a larger eps_bound for OPTICS to ensure more points are processed
  opt <- optics(x, eps = 10, minPts = best_min_pts)
  
  # Try different extraction parameters with a wider range
  eps_cl_values <- seq(0.1, 2.0, by = 0.1)
  best_eps_cl <- 0.5  # Default
  best_ari <- -1
  
  for (eps_cl in eps_cl_values) {
    cl_opt_temp <- extractDBSCAN(opt, eps_cl = eps_cl)$cluster
    # Skip if all points are noise
    if (length(unique(cl_opt_temp)) > 1) {
      ari_temp <- adjustedRandIndex(cl_opt_temp, true_lbls)
      if (ari_temp > best_ari) {
        best_ari <- ari_temp
        best_eps_cl <- eps_cl
      }
    }
  }
  
  # If we couldn't find a good extraction parameter, try Xi method
  if (best_ari == -1) {
    xi_values <- seq(0.01, 0.2, by = 0.01)
    for (xi in xi_values) {
      cl_opt_temp <- extractXi(opt, xi = xi)$cluster
      if (length(unique(cl_opt_temp)) > 1) {
        ari_temp <- adjustedRandIndex(cl_opt_temp, true_lbls)
        if (ari_temp > best_ari) {
          best_ari <- ari_temp
          best_xi <- xi
          # Set extraction method to Xi
          extraction_method <- "Xi"
        }
      }
    }
    # Use the Xi extraction if it worked better
    if (exists("extraction_method") && extraction_method == "Xi") {
      cl_opt <- extractXi(opt, xi = best_xi)$cluster
      extraction_params <- paste("xi =", best_xi)
    } else {
      # Fall back to DBSCAN extraction
      cl_opt <- extractDBSCAN(opt, eps_cl = best_eps_cl)$cluster
      extraction_params <- paste("eps_cl =", best_eps_cl)
    }
  } else {
    # Use DBSCAN extraction
    cl_opt <- extractDBSCAN(opt, eps_cl = best_eps_cl)$cluster
    extraction_params <- paste("eps_cl =", best_eps_cl)
  }
  
  opt_time <- difftime(Sys.time(), start_time, units = "secs")
  
  res$OPTICS_ARI <- adjustedRandIndex(cl_opt, true_lbls)
  
  # Update the params string based on the extraction method used
  if(exists("extraction_method") && extraction_method == "Xi") {
    res$OPTICS_params <- paste("minPts =", best_min_pts, ",", extraction_params)
  } else {
    res$OPTICS_params <- paste("minPts =", best_min_pts, ",", extraction_params)
  }
  
  res$OPTICS_time <- as.numeric(opt_time)
  res$OPTICS_n_clusters <- length(unique(cl_opt[cl_opt > 0]))
  res$OPTICS_noise <- sum(cl_opt == 0) / length(cl_opt) * 100  # % noise points
  
  # PAM (k-medoids)
  start_time <- Sys.time()
  k <- length(unique(true_lbls))
  pam_result <- pam(x, k)
  pam_time <- difftime(Sys.time(), start_time, units = "secs")
  
  res$PAM_ARI <- adjustedRandIndex(pam_result$clustering, true_lbls)
  res$PAM_params <- paste("k =", k)
  res$PAM_time <- as.numeric(pam_time)
  res$PAM_n_clusters <- k
  res$PAM_silhouette <- mean(pam_result$silinfo$avg.width)  # Average silhouette width
  
  # Contaminated Normal Distribution approach using mclust
  start_time <- Sys.time()
  k <- length(unique(true_lbls))
  
  # Function to implement contaminated normal model using mclust
  fit_contaminated_normal <- function(data, k) {
    tryCatch({
      # VVV: allows variable volumes, variable shapes, variable orientations
      # This model is better for handling outliers
      model1 <- Mclust(data, G = k, modelNames = "VVV")
      
      # Get model with different covariance structure - higher chance of modeling contamination
      model2 <- Mclust(data, G = k, modelNames = "VVI")
      
      # Choose the model with better BIC
      if (!is.null(model1) && !is.null(model2)) {
        if (model1$bic > model2$bic) {
          model <- model1
        } else {
          model <- model2
        }
      } else if (!is.null(model1)) {
        model <- model1
      } else if (!is.null(model2)) {
        model <- model2
      } else {
        # If both failed, try a more basic model
        model <- Mclust(data, G = k)
      }
      
      # If we have a model, get the cluster assignments
      if (!is.null(model)) {
        # Calculate Mahalanobis distance for each observation to identify outliers
        mahal_distances <- numeric(nrow(data))
        
        # For each observation, calculate its Mahalanobis distance to assigned cluster center
        for (i in 1:nrow(data)) {
          cluster_i <- model$classification[i]
          mean_i <- model$parameters$mean[, cluster_i]
          
          # Get covariance matrix for this cluster (handle different model structures)
          if (model$modelName %in% c("VVV", "VVI", "VVE")) {
            # Different covariance for each cluster
            if (length(dim(model$parameters$variance$sigma)) == 3) {
              cov_i <- model$parameters$variance$sigma[,, cluster_i]
            } else {
              # Handle diagonal case
              cov_i <- diag(model$parameters$variance$sigma[, cluster_i])
            }
          } else {
            # Same covariance for all clusters or other structure
            if (is.list(model$parameters$variance)) {
              if ("sigma" %in% names(model$parameters$variance)) {
                if (is.matrix(model$parameters$variance$sigma)) {
                  cov_i <- model$parameters$variance$sigma
                } else if (is.array(model$parameters$variance$sigma) && 
                           length(dim(model$parameters$variance$sigma)) == 3) {
                  cov_i <- model$parameters$variance$sigma[,, cluster_i]
                } else {
                  # Diagonal case
                  cov_i <- diag(model$parameters$variance$sigma)
                }
              } else {
                # Default to identity matrix if structure unclear
                cov_i <- diag(ncol(data))
              }
            } else {
              # Default to identity matrix
              cov_i <- diag(ncol(data))
            }
          }
          
          # Make sure covariance matrix is invertible
          if (any(is.na(cov_i)) || det(cov_i) <= 1e-10) {
            cov_i <- diag(ncol(data))
          }
          
          # Calculate Mahalanobis distance
          data_i <- as.numeric(data[i,])
          mahal_distances[i] <- mahalanobis(matrix(data_i, ncol = length(data_i)), 
                                            mean_i, cov_i)
        }
        
        # Identify outliers using chi-square distribution threshold
        # Points with distance greater than 95th percentile of chi-square are outliers
        outliers <- mahal_distances > qchisq(0.95, df = ncol(data))
        outlier_pct <- sum(outliers) / length(outliers) * 100
        
        # Create result
        result <- list(
          cluster = model$classification,
          outliers = outliers,
          outlier_pct = outlier_pct,
          modelName = model$modelName,
          BIC = model$bic,
          parameters = list(mean = model$parameters$mean, 
                            variance = model$parameters$variance),
          model = model
        )
        
        return(result)
      } else {
        # If all models failed, use kmeans as a fallback
        km <- kmeans(data, centers = k, nstart = 10)
        result <- list(
          cluster = km$cluster,
          outliers = rep(FALSE, nrow(data)),
          outlier_pct = 0,
          modelName = "Fallback_KMeans",
          BIC = NA,
          parameters = NULL,
          model = km
        )
        return(result)
      }
    }, error = function(e) {
      warning(paste("Contaminated normal model failed:", e$message))
      # Fallback to kmeans if all else fails
      km <- kmeans(data, centers = k, nstart = 10)
      result <- list(
        cluster = km$cluster,
        outliers = rep(FALSE, nrow(data)),
        outlier_pct = 0,
        modelName = "Fallback_KMeans",
        BIC = NA,
        parameters = NULL,
        model = km
      )
      return(result)
    })
  }
  
  # Fit the contaminated normal model
  cn_model <- fit_contaminated_normal(x, k)
  cn_time <- difftime(Sys.time(), start_time, units = "secs")
  
  # Extract results
  res$CN_ARI <- adjustedRandIndex(cn_model$cluster, true_lbls)
  res$CN_params <- paste("k =", k, ", model =", cn_model$modelName)
  res$CN_time <- as.numeric(cn_time)
  res$CN_n_clusters <- length(unique(cn_model$cluster))
  res$CN_BIC <- cn_model$BIC
  res$CN_outlier_pct <- cn_model$outlier_pct
  
  return(res)
}

# Run benchmarks & display enhanced results -----------------------
results <- lapply(datasets_for_clustering, function(d) {
  cat("\nBenchmarking dataset with outlier density:", 
      round(d$outlier_density, 2), "%\n")
  benchmark_clust(d$x, d$labels)
})

# Process results into a comprehensive comparison table
results_df <- do.call(rbind, lapply(names(results), function(dataset_name) {
  result <- results[[dataset_name]]
  
  # Create a default structure with NA values for missing fields
  df <- data.frame(
    Dataset = dataset_name,
    Outlier_Density = datasets_for_clustering[[dataset_name]]$outlier_density,
    
    DBSCAN_ARI = NA,
    DBSCAN_Params = NA,
    DBSCAN_Time = NA,
    DBSCAN_Clusters = NA,
    DBSCAN_Noise_Pct = NA,
    
    OPTICS_ARI = NA,
    OPTICS_Params = NA,
    OPTICS_Time = NA,
    OPTICS_Clusters = NA,
    OPTICS_Noise_Pct = NA,
    
    PAM_ARI = NA,
    PAM_Params = NA,
    PAM_Time = NA,
    PAM_Silhouette = NA,
    
    CN_ARI = NA,
    CN_Params = NA,
    CN_Time = NA,
    CN_BIC = NA,
    CN_Outlier_Pct = NA,
    
    stringsAsFactors = FALSE
  )
  
  # Only fill in the fields that exist in the result
  if (!is.null(result$DBSCAN_ARI)) df$DBSCAN_ARI <- result$DBSCAN_ARI
  if (!is.null(result$DBSCAN_params)) df$DBSCAN_Params <- result$DBSCAN_params
  if (!is.null(result$DBSCAN_time)) df$DBSCAN_Time <- result$DBSCAN_time
  if (!is.null(result$DBSCAN_n_clusters)) df$DBSCAN_Clusters <- result$DBSCAN_n_clusters
  if (!is.null(result$DBSCAN_noise)) df$DBSCAN_Noise_Pct <- result$DBSCAN_noise
  
  if (!is.null(result$OPTICS_ARI)) df$OPTICS_ARI <- result$OPTICS_ARI
  if (!is.null(result$OPTICS_params)) df$OPTICS_Params <- result$OPTICS_params
  if (!is.null(result$OPTICS_time)) df$OPTICS_Time <- result$OPTICS_time
  if (!is.null(result$OPTICS_n_clusters)) df$OPTICS_Clusters <- result$OPTICS_n_clusters
  if (!is.null(result$OPTICS_noise)) df$OPTICS_Noise_Pct <- result$OPTICS_noise
  
  if (!is.null(result$PAM_ARI)) df$PAM_ARI <- result$PAM_ARI
  if (!is.null(result$PAM_params)) df$PAM_Params <- result$PAM_params
  if (!is.null(result$PAM_time)) df$PAM_Time <- result$PAM_time
  if (!is.null(result$PAM_silhouette)) df$PAM_Silhouette <- result$PAM_silhouette
  
  if (!is.null(result$CN_ARI)) df$CN_ARI <- result$CN_ARI
  if (!is.null(result$CN_params)) df$CN_Params <- result$CN_params
  if (!is.null(result$CN_time)) df$CN_Time <- result$CN_time
  if (!is.null(result$CN_BIC)) df$CN_BIC <- result$CN_BIC
  if (!is.null(result$CN_outlier_pct)) df$CN_Outlier_Pct <- result$CN_outlier_pct
  
  return(df)
}))

print(results_df)

# Extract just the ARI scores for comparison
ari_scores <- data.frame(
  Dataset = results_df$Dataset,
  DBSCAN = results_df$DBSCAN_ARI,
  OPTICS = results_df$OPTICS_ARI,
  PAM = results_df$PAM_ARI,
  CN = results_df$CN_ARI
)

# Print ARI score comparison
print("Adjusted Rand Index (ARI) Comparison:")
print(ari_scores)

# Summary
cat("\n\n========== CLUSTERING ALGORITHM COMPARATIVE ANALYSIS ==========\n")
cat("Dataset outlier characteristics:\n")
for(ds in names(datasets_for_clustering)) {
  cat(paste0("- ", ds, ": ", round(datasets_for_clustering[[ds]]$outlier_density, 2), 
             "% outlier density\n"))
}

cat("\nBest performing algorithm by dataset (based on ARI):\n")
for(i in 1:nrow(ari_scores)) {
  algo_scores <- as.numeric(ari_scores[i, 2:5])
  algo_names <- colnames(ari_scores)[2:5]
  
  # Handle NA values
  valid_indices <- which(!is.na(algo_scores))
  if(length(valid_indices) > 0) {
    valid_scores <- algo_scores[valid_indices]
    valid_names <- algo_names[valid_indices]
    
    best_algo <- valid_names[which.max(valid_scores)]
    best_score <- max(valid_scores, na.rm = TRUE)
    
    cat(paste0("- ", ari_scores$Dataset[i], ": ", best_algo, 
               " (ARI = ", round(best_score, 3), ")\n"))
  } else {
    cat(paste0("- ", ari_scores$Dataset[i], ": No valid ARI scores available\n"))
  }
}
