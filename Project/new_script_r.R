library(mlbench)
library(ggplot2)
library(dbscan)
library(cluster)
library(mclust)       
library(reshape2)    
library(gridExtra)

# Loading datasets

# =========================
# Load and preprocess MAGIC dataset
# =========================
# Read and clean the data
# Load and clean the MAGIC data
col_names_magic <- c(
  "fLength", "fWidth", "fSize", "fConc", "fConc1", 
  "fAsym", "fM3Long", "fM3Trans", "fAlpha", "fDist", "Class"
)
magic_data <- read_csv("C:/Users/tf245/Documents/GitHub/MATH252_RCode/Project/magic+gamma+telescope/magic04.data",
                       col_names = col_names_magic)

magic_data <- na.omit(magic_data)  # Remove any incomplete rows

magic_labels <- as.numeric(as.factor(magic_data$Class))

magic_features <- scale(magic_data[, 1:10])

# =========================
# Load the wine dataset
# =========================
data(wine)
wine <- wine

# Function to calculate outlier metrics with proper dimensionality
calculate_outlier_metrics <- function(data_frame) {
  # Only for numeric columns
  numeric_data <- data_frame[, sapply(data_frame, is.numeric)]
  n_dims <- ncol(numeric_data)  # Get the full dimensionality
  
  metrics <- list()
  
  # IQR-based outlier count per feature
  outlier_counts <- sapply(numeric_data, function(x) {
    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower_bound <- q1 - 1.5 * iqr
    upper_bound <- q3 + 1.5 * iqr
    sum(x < lower_bound | x > upper_bound, na.rm = TRUE)
  })
  
  # Percentage of outliers per feature
  outlier_percentages <- outlier_counts / nrow(numeric_data) * 100
  
  # Overall outlier density (average percentage across features)
  overall_outlier_density <- mean(outlier_percentages)
  
  # Feature-wise skewness to understand distribution
  skewness_values <- sapply(numeric_data, function(x) {
    n <- length(x)
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    sum((x - m)^3, na.rm = TRUE) / (n * s^3)
  })
  
  # Calculate Mahalanobis distance using full dimensionality
  scaled_data <- scale(numeric_data)
  mahal_dist <- mahalanobis(scaled_data, 
                            center = colMeans(scaled_data), 
                            cov = cov(scaled_data))
  
  # Identify outliers using chi-square with proper degrees of freedom
  mahal_outliers <- mahal_dist > qchisq(0.95, df = n_dims)  # Using full dimensionality
  mahal_outlier_percent <- sum(mahal_outliers) / length(mahal_outliers) * 100
  
  metrics$outlier_counts <- outlier_counts
  metrics$outlier_percentages <- outlier_percentages
  metrics$overall_outlier_density <- overall_outlier_density
  metrics$skewness <- skewness_values
  metrics$mahal_outliers <- mahal_outliers
  metrics$mahal_outlier_percent <- mahal_outlier_percent
  metrics$mahal_distances <- mahal_dist
  
  return(metrics)
}

# Preprocess the magic dataset to ensure all columns are properly handled
# First, make sure we're only working with numeric columns for analysis
magic_numeric <- magic_data[, sapply(magic_data, is.numeric)]
cat("Numeric variables in magic dataset:", paste(names(magic_numeric), collapse=", "), "\n")

# Calculate outlier metrics for both datasets
wine_data_outlier_metrics <- calculate_outlier_metrics(wine)
magic_outlier_metrics <- calculate_outlier_metrics(magic_numeric)

# Visualization for outlier patterns
# wine dataset outlier visualization
wine_data_num_all <- wine[, sapply(wine, is.numeric)]
n_dims_wine <- ncol(wine_data_num_all)

# Boxplots for outlier visualization
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
for (nm in names(wine_data_num_all)) {
  boxplot(wine_data_num_all[[nm]], main = paste0(nm, "\n(", 
                                             round(wine_data_outlier_metrics$outlier_percentages[nm], 1), "% outliers)"), 
          horizontal = TRUE, col = "lightblue")
}

# Scatter plot with outliers using full-dimensional Mahalanobis distance
wine_data_2 <- as.data.frame(scale(wine_data_num_all))[, 1:2]
colnames(wine_data_2) <- c("X1", "X2")

# Use the full-dimensional Mahalanobis outliers
wine_data_2$outlier <- wine_data_outlier_metrics$mahal_outliers

# Scatter plot
p1 <- ggplot(wine_data_2, aes(x = X1, y = X2, color = outlier)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("blue", "red"), 
                     labels = c("Regular", "Outlier")) +
  labs(
    title = "Wine Dataset: Feature Distribution with Outliers",
    subtitle = paste("Mahalanobis outliers (", n_dims_wine, " dimensions): ", 
                     round(wine_data_outlier_metrics$mahal_outlier_percent, 2), "%"),
    x = names(wine_data_num_all)[1],
    y = names(wine_data_num_all)[2],
    color = "Point Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
p1

# magic dataset outlier visualization
magic_num_all <- magic_numeric  # Use the preprocessed numeric-only data
n_dims_magic <- ncol(magic_num_all)
n_magic_vars <- ncol(magic_num_all)
rows <- ceiling(sqrt(n_magic_vars))
cols <- ceiling(n_magic_vars / rows)

# Boxplots for outlier visualization
par(mfrow = c(rows, cols), mar = c(4, 4, 2, 1))
for (nm in names(magic_num_all)) {
  boxplot(magic_num_all[[nm]], main = paste0(nm, "\n(", 
                                             round(magic_outlier_metrics$outlier_percentages[nm], 1), "% outliers)"), 
          horizontal = TRUE, col = "lightgreen")
}

# Scatter plot with outlier highlighting for first two features
magic_2 <- as.data.frame(scale(magic_num_all))[, 1:2]
colnames(magic_2) <- c("X1", "X2")

# Use the full-dimensional Mahalanobis outliers
magic_2$outlier <- magic_outlier_metrics$mahal_outliers

# Create enhanced scatter plot
p2 <- ggplot(magic_2, aes(x = X1, y = X2, color = outlier)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("darkgreen", "red"), 
                     labels = c("Regular", "Outlier")) +
  labs(
    title = "Dry magic Dataset: Feature Distribution with Outliers",
    subtitle = paste("Mahalanobis outliers (", n_dims_magic, " dimensions): ", 
                     round(magic_outlier_metrics$mahal_outlier_percent, 2), "%"),
    x = names(magic_num_all)[1],
    y = names(magic_num_all)[2],
    color = "Point Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
p2

# Print comparative outlier plots
grid.arrange(p1, p2, ncol = 2)

# Add pairs plots for better visualization (as suggested by professor)
# For bc_data dataset
par(mfrow = c(1, 1))  # Reset the plotting environment
if(ncol(wine_data_num_all) <= 8) {
  # Color points by outlier status
  pairs(wine_data_num_all, col = ifelse(bc_data_outlier_metrics$mahal_outliers, "red", "blue"),
        pch = 19, cex = 0.7,
        main = "Wine Dataset: Pairwise Relationships with Outliers Highlighted")
} else {
  # If too many columns, select a subset
  cols_subset <- names(wine_data_num_all)[1:min(8, ncol(wine_data_num_all))]
  pairs(wine_data_num_all[, cols_subset], 
        col = ifelse(wine_data_outlier_metrics$mahal_outliers, "red", "blue"),
        pch = 19, cex = 0.7,
        main = "Wine Dataset: Pairwise Relationships (Subset) with Outliers Highlighted")
}

# For magic dataset 
par(mfrow = c(1, 1))

# Select first 6 variables for better visualization (can be changed)
magic_vars_for_pairs <- names(magic_num_all)[1:min(6, ncol(magic_num_all))]
cat("Magic variables selected for pairs plot:", paste(magic_vars_for_pairs, collapse=", "), "\n")

# Explicitly create the subset data frame
magic_subset <- magic_num_all[, bean_vars_for_pairs, drop = FALSE]

# Create the pairs plot with specific settings
pairs(magic_subset, 
      col = ifelse(magic_outlier_metrics$mahal_outliers, "red", "darkgreen"),
      pch = 19, cex = 0.5,  # Smaller points for clarity
      gap = 0.5,            # More space between panels
      cex.labels = 1.2,     # Larger variable names
      main = "magic Dataset: Pairwise Relationships with Outliers Highlighted")

# Prepare for clustering ------------------------------------------
# Ensure features and labels come from the same cleaned data

# Clean and align breast cancer data
wine_clean <- na.omit(wine)
wine_data_x <- scale(wine_clean[, sapply(wine_clean, is.numeric)])
wine_data_labels <- as.numeric(wine_clean$Type)

# Clean and align magic data
magic_clean <- na.omit(magic_data)
magic_x <- scale(magic_clean[, 1:10])
magic_labels <- as.numeric(as.factor(magic_clean$Class))

# Final dataset list
datasets <- list(
  wine = list(x = wine_data_x, labels = wine_data_labels, 
                 outlier_density = wine_data_outlier_metrics$overall_outlier_density),
  magic = list(x = magic_x, labels = magic_labels, 
               outlier_density = magic_outlier_metrics$overall_outlier_density)
)

# Benchmark function -------------------------------------
benchmark_clust <- function(x, true_lbls) {
  library(dbscan)
  library(cluster)
  library(mclust)
  res <- list()
  
  # Timer for performance measurement
  start_time <- Sys.time()
  
  # DBSCAN with parameter tuning
  # Try a range of eps values to find optimal
  eps_values <- seq(0.3, 0.8, by = 0.1)
  best_eps <- 0.5
  best_ari <- -1
  
  for (eps in eps_values) {
    min_pts <- 5  # Default minPts
    db_temp <- dbscan(x, eps = eps, minPts = min_pts)
    # Skip if all points are noise (cluster = 0)
    if (length(unique(db_temp$cluster)) > 1) {
      # Debug check to catch vector length mismatch
      if (length(db_temp$cluster) != length(true_lbls)) {
        stop(paste(
          "Vector length mismatch:",
          "db_temp$cluster =", length(db_temp$cluster),
          "true_lbls =", length(true_lbls)
        ))
      }
      ari_temp <- adjustedRandIndex(db_temp$cluster, true_lbls)
      if (ari_temp > best_ari) {
        best_ari <- ari_temp
        best_eps <- eps
      }
    }
  }
  
  # Use best eps for final DBSCAN
  db <- dbscan(x, eps = best_eps, minPts = 5)
  db_time <- difftime(Sys.time(), start_time, units = "secs")
  
  res$DBSCAN_ARI <- adjustedRandIndex(db$cluster, true_lbls)
  res$DBSCAN_params <- paste("eps =", best_eps, ", minPts = 5")
  res$DBSCAN_time <- as.numeric(db_time)
  res$DBSCAN_n_clusters <- length(unique(db$cluster[db$cluster > 0]))
  res$DBSCAN_noise <- sum(db$cluster == 0) / length(db$cluster) * 100  # % noise points
  
  # OPTICS with parameter tuning
  start_time <- Sys.time()
  opt <- optics(x, eps = best_eps, minPts = 5)
  
  # Try different extraction parameters
  eps_cl_values <- seq(0.3, 0.8, by = 0.1)
  best_eps_cl <- best_eps  # Default
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
  
  cl_opt <- extractDBSCAN(opt, eps_cl = best_eps_cl)$cluster
  opt_time <- difftime(Sys.time(), start_time, units = "secs")
  
  res$OPTICS_ARI <- adjustedRandIndex(cl_opt, true_lbls)
  res$OPTICS_params <- paste("eps =", best_eps, ", eps_cl =", best_eps_cl)
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
results <- lapply(datasets, function(d) {
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
    Outlier_Density = datasets[[dataset_name]]$outlier_density,
    
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

# summary
cat("\n\n========== CLUSTERING ALGORITHM COMPARATIVE ANALYSIS ==========\n")
cat("Dataset outlier characteristics:\n")
for(ds in names(datasets)) {
  cat(paste0("- ", ds, ": ", round(datasets[[ds]]$outlier_density, 2), 
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
