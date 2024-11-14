##Script to run factor analysis on significant predictors, selected from EN!
library(ggplot2)
library(reshape2)
library(corrplot)
library(ggdendro)

## Data prep ##

RNGkind(kind = "Mersenne-Twister")
set.seed(1234)

setwd("~/Desktop/NeuroPathPredict/V1/Model_V1/Input_data")
df_in<- read.csv("XY_MFG_V1_right_std_1021.csv")

# Retrieve the names from the text file
list_SP_30 <- readLines("EN1_SigPred_30.txt")

# Subset the data frame to include only the significant predictors
df_list <- df_in[, colnames(df_in) %in% list_SP_30]

# Determine the number of factors to extract for FA
num_factors <- 6

# Perform Factor Analysis
factor_analysis_result <- factanal(df_list, factors = num_factors, rotation = "varimax", start = NULL)

# Print FA results
print(factor_analysis_result)

# Plot the loadings
loadings <- as.data.frame(factor_analysis_result$loadings[, 1:num_factors])
loadings$Variable <- rownames(loadings)

# Melt the loadings data frame for ggplot2
loadings_melted <- melt(loadings, id.vars = "Variable")

# Plot the factor loadings
ggplot(loadings_melted, aes(x = Variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Factor Loadings",
       x = "Variables",
       y = "Loadings") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# Plot the covariance matrix
cov_matrix <- cov(df_list)
corrplot(cov_matrix, method = "color", type = "upper", order = "hclust", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45, 
         title = "Covariance Matrix", mar = c(0,0,1,0), number.cex = 0.4, tl.cex = 0.6)

# Plot the dendrogram
dist_matrix <- dist(t(loadings[, 1:num_factors]))
hc <- hclust(dist_matrix)
dendro_data <- dendro_data(hc)
ggdendrogram(dendro_data, rotate = TRUE, size = 2) +
  labs(title = "Dendrogram of Factor Loadings") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Calculate eigenvalues for the scree plot
cor_matrix <- cor(df_list)
eigenvalues <- eigen(cor_matrix)$values
scree_data <- data.frame(Factor = 1:length(eigenvalues), Eigenvalue = eigenvalues)

# Plot the scree plot
ggplot(scree_data, aes(x = Factor, y = Eigenvalue)) +
  geom_point() +
  geom_line() +
  labs(title = "Scree Plot", x = "Factor", y = "Eigenvalue") +
  theme_minimal()

# Cluster variables based on factor loadings
threshold <- 0.6  # Define a threshold for high loadings
clusters <- apply(loadings[, 1:num_factors], 1, function(x) {
  factor <- which.max(abs(x))
  if (abs(x[factor]) >= threshold) {
    return(factor)
  } else {
    return(NA)
  }
})

# Create a data frame with variable names and their corresponding clusters
cluster_df <- data.frame(
  Variable = rownames(loadings),
  Cluster = clusters
)

# Remove variables that do not meet the threshold
cluster_df <- cluster_df[!is.na(cluster_df$Cluster), ]

# Print the clustered variables
print(cluster_df)

# Calculate and plot covariance matrices for each cluster
for (i in 1:4) {
  cluster_vars <- cluster_df$Variable[cluster_df$Cluster == i]
  if (length(cluster_vars) > 1) {
    cluster_data <- df_list[, cluster_vars]
    cov_matrix <- cov(cluster_data)
    corrplot(cov_matrix, method = "color", type = "upper", order = "hclust", 
             addCoef.col = "black", tl.col = "black", tl.srt = 45, 
             title = paste("Covariance Matrix for Cluster", i), mar = c(0,0,1,0), 
             number.cex = 0.5, tl.cex = 0.7)
  } else {
    print(paste("Cluster", i, "has less than 2 variables, skipping covariance matrix."))
  }
}
