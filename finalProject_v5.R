library(mlbench)
library(ggplot2)
library(dbscan)
library(cluster)
library(ContaminatedMixt)
library(mclust)
library(reshape2)
library(readr)
library(fpc)
library(aricode)

# Load datasets
data("Glass")
glass <- Glass
wine <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv", sep = ";")
magic_data <- read_csv("magic04.data")
beans <- read.csv("bean.csv", header = TRUE, stringsAsFactors = FALSE)
seeds <- read.table("seeds_dataset.txt", header = FALSE)
wholesale_data <- read.csv("wholesalecustomers.csv", header = TRUE, stringsAsFactors = FALSE)

# Mahalanobis outlier detection
calculate_mahalanobis_outliers <- function(data_frame) {
  numeric_data <- data_frame[, sapply(data_frame, is.numeric)]
  n_dims <- ncol(numeric_data)
  scaled_data <- scale(numeric_data)
  cov_matrix <- tryCatch({
    cov(scaled_data)
  }, error = function(e) {
    diag(diag(cov(scaled_data))) + 0.001 * diag(n_dims)
  })
  mahal_dist <- mahalanobis(scaled_data, center = colMeans(scaled_data), cov = cov_matrix)
  mahal_outliers <- mahal_dist > qchisq(0.95, df = n_dims)
  list(
    mahal_outliers = mahal_outliers,
    mahal_outlier_percent = sum(mahal_outliers) / length(mahal_outliers) * 100,
    mahal_distances = mahal_dist
  )
}

# Organize datasets
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
  
  cat("\nDataset:", dataset_name)
  cat("\nDimensions:", nrow(numeric_data), "x", ncol(numeric_data))
  cat("\nMahalanobis outliers (", ncol(numeric_data), " dimensions): ", 
      round(outlier_metrics$mahal_outlier_percent, 2), "%\n")
  
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

# Prepare datasets for clustering
prepare_for_clustering <- function(dataset_info) {
  x <- scale(dataset_info$numeric_data)
  x[is.na(x)] <- 0
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

# Clustering benchmark
benchmark_clust <- function(x, true_lbls) {
  res <- list()
  start_time <- Sys.time()
  eps_values <- seq(0.3, 1.5, by = 0.1)
  best_eps <- 0.5
  best_ari <- -1
  min_pts_values <- c(3, 4, 5, 10)
  
  for (min_pts in min_pts_values) {
    for (eps in eps_values) {
      db_temp <- dbscan::dbscan(x, eps = eps, minPts = min_pts)
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
  
  db <- dbscan::dbscan(x, eps = best_eps, minPts = best_min_pts)
  db_time <- difftime(Sys.time(), start_time, units = "secs")
  db_clusters <- db$cluster
  valid_db <- db_clusters > 0
  avg_sil_db <- if (length(unique(db_clusters[valid_db])) > 1)
    mean(silhouette(db_clusters[valid_db], dist(x[valid_db, ]))[, 3]) else NA
  
  res$DBSCAN_ARI <- adjustedRandIndex(db_clusters, true_lbls)
  res$DBSCAN_NMI <- NMI(true_lbls, db_clusters)
  res$DBSCAN_params <- paste("eps =", best_eps, ", minPts =", best_min_pts)
  res$DBSCAN_time <- as.numeric(db_time)
  res$DBSCAN_n_clusters <- length(unique(db_clusters[db_clusters > 0]))
  res$DBSCAN_noise <- sum(db_clusters == 0) / length(db_clusters) * 100
  res$DBSCAN_silhouette <- avg_sil_db
  
  # OPTICS
  start_time <- Sys.time()
  opt <- optics(x, eps = 10, minPts = best_min_pts)
  eps_cl_values <- seq(0.1, 2.0, by = 0.1)
  best_eps_cl <- 0.5
  best_ari <- -1
  for (eps_cl in eps_cl_values) {
    cl_opt_temp <- extractDBSCAN(opt, eps_cl = eps_cl)$cluster
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
  valid_opt <- cl_opt > 0
  avg_sil_opt <- if (length(unique(cl_opt[valid_opt])) > 1)
    mean(silhouette(cl_opt[valid_opt], dist(x[valid_opt, ]))[, 3]) else NA
  
  res$OPTICS_ARI <- adjustedRandIndex(cl_opt, true_lbls)
  res$OPTICS_NMI <- NMI(true_lbls, cl_opt)
  res$OPTICS_params <- paste("minPts =", best_min_pts, ", eps_cl =", best_eps_cl)
  res$OPTICS_time <- as.numeric(opt_time)
  res$OPTICS_n_clusters <- length(unique(cl_opt[cl_opt > 0]))
  res$OPTICS_noise <- sum(cl_opt == 0) / length(cl_opt) * 100
  res$OPTICS_silhouette <- avg_sil_opt
  
  # PAM
  start_time <- Sys.time()
  k <- length(unique(true_lbls))
  pam_result <- pam(x, k)
  pam_time <- difftime(Sys.time(), start_time, units = "secs")
  res$PAM_ARI <- adjustedRandIndex(pam_result$clustering, true_lbls)
  res$PAM_NMI <- NMI(true_lbls, pam_result$clustering)
  res$PAM_params <- paste("k =", k)
  res$PAM_time <- as.numeric(pam_time)
  res$PAM_silhouette <- mean(pam_result$silinfo$avg.width)
  
  # CN (Gaussian Mixture)
  start_time <- Sys.time()
  pca_result <- prcomp(x, scale. = TRUE)
  x_reduced <- pca_result$x[, 1:min(5, ncol(pca_result$x))]
  model <- tryCatch({
    CNmixt(x_reduced, G = k, initialization = "kmeans", alphafix = TRUE, model = "EEE")
  }, error = function(e) {
    cat("CNmixt failed: ", e$message, "\n")
    return(NULL)
  })
  
  # Proceed only if the object is valid and classification is complete
  if (!is.null(model) &&
      "classification" %in% names(model) &&
      length(model$classification) == nrow(x)) {
    
    cn_clusters <- model$classification
    # Continue with ARI/NMI etc.
    
  } else {
    cat("CNmixt did not return a valid classification. Falling back to Mclust.\n")
    model <- Mclust(x, G = k)
    cn_clusters <- model$classification
    # Use fallback model for metrics
  }
  cn_time <- difftime(Sys.time(), start_time, units = "secs")
  avg_sil_cn <- if (length(unique(cn_clusters)) > 1)
    mean(silhouette(cn_clusters, dist(x))[, 3]) else NA
  res$CN_ARI <- adjustedRandIndex(cn_clusters, true_lbls)
  res$CN_NMI <- NMI(true_lbls, cn_clusters)
  res$CN_params <- paste("k =", k, ", model =", model$modelName)
  res$CN_time <- as.numeric(cn_time)
  res$CN_n_clusters <- length(unique(cn_clusters))
  res$CN_BIC <- model$bic
  res$CN_outlier_pct <- NA
  res$CN_silhouette <- avg_sil_cn
  
  return(res)
}

# Run benchmarks
results <- lapply(datasets_for_clustering, function(d) {
  cat("\nBenchmarking dataset with outlier density:", round(d$outlier_density, 2), "%\n")
  benchmark_clust(d$x, d$labels)
})

# Results table
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

# ARI & Silhouette comparison
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

cat("\n\nClustering Algorithm Summary\n")
cat("Dataset outlier characteristics:\n")
for(ds in names(datasets_for_clustering)) {
  cat(paste0("- ", ds, ": ", round(datasets_for_clustering[[ds]]$outlier_density, 2), "% outlier density\n"))
}

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
    cat(paste0("- ", nmi_scores$Dataset[i], ": No valid ARI scores available\n"))
  }
}

cat("\nBest performing algorithm by dataset (based on Silhouette Score):\n")
sil_scores <- comparison_scores
for(i in 1:nrow(sil_scores)) {
  algo_scores <- as.numeric(sil_scores[i, 10:13])  # Columns for silhouette scores
  algo_names <- colnames(sil_scores)[10:13]        # Matching algorithm names
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

