library(readxl) ## this is needed to read excel files
library(readr) ## this is needed to read the csv
library(dbscan) ## this is needed to run dbscan
library(cluster) ## this is needed to run pam
library(ContaminatedMixt) # needed to run contaminated normal
library(NbClust) # needed to check the number of clusters a dataset has
real_estate_data <- read_excel("C:/Users/tf245/Documents/GitHub/MATH252_RCode/Project/real+estate+valuation+data+set/Real estate valuation data set.xlsx") # read in real estate data

# Read the .names file as lines of text
lines_magic_names <- readLines("C:/Users/tf245/Documents/GitHub/MATH252_RCode/Project/magic+gamma+telescope/magic04.names")

# Extract lines that start with a number and a period
name_lines_magic <- grep("^\\d+\\.", lines_magic_names, value = TRUE)

# Use regex to pull the actual variable names
col_names_magic <- sub("^\\d+\\.\\s*([a-zA-Z0-9_]+):.*", "\\1", name_lines_magic)

#read the data file
magic_data <- read_csv("C:/Users/tf245/Documents/GitHub/MATH252_RCode/Project/magic+gamma+telescope/magic04.data", col_names = col_names_magic)
magic_data_numeric <- magic_data[, 1:10]

#check for the number of clusters
#NbClust(real_estate_data, distance = "euclidean", min.nc = 2, max.nc = 15, method = "kmeans")
#NbClust(magic_data_numeric, distance = "euclidean", min.nc = 2, max.nc = 15, method = "kmeans")

## cluster real estate data ##
real_estate_scaled <- scale(real_estate_data) # scale the data
real_estate_dbscan <- dbscan(real_estate_scaled, eps = 2.75, minPts = 8) #dbscan
real_estate_optics <- optics(real_estate_scaled, minPts = 8) #optics
real_estate_pam <- pam(real_estate_data, k = 3)

## cluster the magic data ##
magic_data_scaled <- scale(magic_data)
magic_dbscan <- dbscan(magic_data_scaled, eps = 2.75, minPts = 8)
magic_optics <- optics(magic_scaled, minPts = 8)
magic_pam <- pam(magic_data, k = 2)


## CN Mixture

fit_best_CN_model <- function(data, G_range = 2:4) {
  tryCatch({
    best_model <- NULL
    best_bic <- -Inf
    
    for (g in G_range) {
      model_names <- c("EEE", "VVV", "EEV", "VEV")
      
      for (model_name in model_names) {
        fit <- tryCatch({
          fit_CNmixt(data, G = g, model = model_name)
        }, error = function(e) NULL)
        
        if (!is.null(fit)) {
          if (!is.null(fit$bic) && fit$bic > best_bic) {
            best_bic <- fit$bic
            best_model <- fit
          }
        }
      }
    }
    
    if (is.null(best_model)) {
      warning("All CNmixt models failed. Falling back to kmeans.")
      km <- kmeans(data, centers = mean(G_range), nstart = 10)
      return(list(
        cluster = km$cluster,
        outliers = rep(FALSE, nrow(data)),
        outlier_pct = 0,
        modelName = "Fallback_KMeans",
        BIC = NA,
        parameters = NULL,
        model = km
      ))
    }
    
    cluster <- best_model$classification
    outliers <- best_model$outliers
    outlier_pct <- sum(outliers) / length(outliers) * 100
    
    return(list(
      cluster = cluster,
      outliers = outliers,
      outlier_pct = outlier_pct,
      modelName = best_model$modelName,
      BIC = best_model$bic,
      G = best_model$G,  # Save best G
      parameters = best_model$parameters,
      model = best_model
    ))
    
  }, error = function(e) {
    warning(paste("CNmixt fitting failed:", e$message))
    km <- kmeans(data, centers = mean(G_range), nstart = 10)
    return(list(
      cluster = km$cluster,
      outliers = rep(FALSE, nrow(data)),
      outlier_pct = 0,
      modelName = "Fallback_KMeans",
      BIC = NA,
      parameters = NULL,
      model = km
    ))
  })
}

#CN for real estate
best_fit_cn_res_estate <- fit_best_CN_model(real_estate_data)

#CN for magic
best_fit_cn_res_magic <- fit_best_CN_model(magic_data_numeric)
