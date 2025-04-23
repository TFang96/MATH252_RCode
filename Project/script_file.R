library(readxl) ## this is needed to read excel files
library(readr) ## this is needed to read the csv
library(dbscan) ## this is needed to run dbscan
library(cluster) ## this is needed to run pam
library(ContaminatedMixt) # needed to run contaminated normal
real_estate_data <- read_excel("C:/Users/tf245/Documents/GitHub/MATH252/Project/real+estate+valuation+data+set/Real estate valuation data set.xlsx") # read in real estate data

# Read the .names file as lines of text
lines_magic_names <- readLines("C:/Users/tf245/Documents/GitHub/MATH252/Project/magic+gamma+telescope/magic04.names")

# Extract lines that start with a number and a period
name_lines_magic <- grep("^\\d+\\.", lines_magic_names, value = TRUE)

# Use regex to pull the actual variable names
col_names_magic <- sub("^\\d+\\.\\s*([a-zA-Z0-9_]+):.*", "\\1", name_lines_magic)

#read the data file
magic_data <- read_csv("C:/Users/tf245/Documents/GitHub/MATH252/Project/magic+gamma+telescope/magic04.data", col_names = col_names_magic)

## Work on the real estate data ##
real_estate_scaled <- scale(real_estate_data) # scale the data
real_estate_dbscan <- dbscan(real_estate_scaled, eps = 2.75, minPts = 8) #dbscan
real_estate_optics <- optics(real_estate_scaled, minPts = 8) #optics
real_estate_pam <- pam(real_estate_data, k = 2)

