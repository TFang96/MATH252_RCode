# Name: Ziheng (Tony) Fang
# MATH 252
# Spring 2025
# Assignment 3
library(MASS) # needed to simulate clusters
library(mclust) # needed to fit data on a gmm
library(MixGHD) # needed for ARI function
library(MixtureMissing) #needed for select_mixture function
#### First work on the cluster that fits well with a GMM. 
#generate mean and sigma for the two clusters
mu1 <- c(0, 0, 0) #mean for cluster 1
mu2 <- c(5, 5, 5) #mean for cluster 2
sigma <- diag(3) #sigma for both clusters

cluster1 <- mvrnorm(n = 100, mu = mu1, Sigma = sigma) #simulate cluster 1
cluster2 <- mvrnorm(n = 100, mu = mu2, Sigma = sigma) #simualte cluster 2

good_normal <- rbind(cluster1, cluster2) #combine the clusters into 1 data set

labels_good <- c(rep(1, 100), rep(2, 100)) #cluster assignments/labels for data points in the good normal cluster

#A GMM cannot capture data that is skewed or has heavy outliers. 

mu1 <- c(0, 0, 0) #mean for cluster 1
mu2 <- c(5, 5, 5) #mean for cluster 2

set.seed(123) #control randomness in generation

cluster1_skewed <- cbind(rexp(100, rate = 1), rnorm(100, mean = 0, sd = 1), rnorm(100, mean = 0, sd = 1)) #simulate a skewed cluster

cluster2_t <- cbind(rt(100, df = 1) * 0.2 + 5, rt(100, df = 1) + 5, rt(100, df = 1) + 5) #simulate a t-distribution w/heavy tails

bad_normal <- rbind(cluster1_skewed, cluster2_t) #combine the two clusters

labels_bad <- c(rep(1, 100), rep(2, 100)) #labels for all the data points in the bad dataset

# Add some extreme outliers
outliers <- matrix(runif(30 * 3, min = -5000, max = 5000), ncol = 3) #generate the outliers
bad_skewed_outliers <- rbind(bad_normal, outliers) # add outliers to data set
labels_bad_outliers <- c(labels_bad, rep(0, 30))  # label 0 = noise/outlier

### Fit good data ###
good_fitted <- Mclust(good_normal, G = 2, mnames = "vvv") #fit the good data using a GMM
bad_fitted <- Mclust(bad_skewed_outliers, G = 2, mnames = "VVV") #fit the bad data (skewed and heavy outlier) using a GMM

cat("ARI for good data", ARI(labels_good, good_fitted$classification), "\n") #print the ARI for the dataset with no outliers and skewness
cat("ARI for dataset that is skewed and has outliers", ARI(labels_bad_outliers, bad_fitted$classification), "\n") #print the ARI for the dataset that is skewed and has outliers.

## select best models to fit the data ##
best_choice_normal <- select_mixture(good_normal, G = 2, model = c("N", "t", "St", "CN"), criterion = "BIC", init_method = "kmeans")
best_choice_skewed_and_outliers <- select_mixture(bad_skewed_outliers, G = 2, model = c("N", "t", "St", "CN"), criterion = "BIC", init_method = "kmeans")

