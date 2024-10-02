# Load necessary libraries
library(readr)       
library(dplyr)       
library(vegan)        
library(ggplot2)     

# Load the data
data <- read_csv('data.csv')  # Read the CSV file into a data frame

# Select relevant columns for PCA
morpho_columns <- c('Morpho1', 'Morpho2', 'Morpho3')  # Define the columns containing morphotype data
abundance_data <- data %>% select(all_of(morpho_columns))  # Select only the morphotype columns from the data

# Hellinger transformation
abundance_data_hellinger <- decostand(abundance_data, method = "hellinger")  # Apply Hellinger transformation to the data

# Perform PCA
pca_result <- rda(abundance_data_hellinger)  # Perform PCA on the Hellinger-transformed data

# Extract PCA scores
pca_scores <- scores(pca_result, display = "sites")  # Extract the PCA scores for the samples

# Convert to data frame for easier plotting
pca_scores_df <- as.data.frame(pca_scores)  # Convert the PCA scores to a data frame

# Add sample names for easier identification
pca_scores_df$sample_name <- data$sample_name  # Add the sample names to the PCA scores data frame

# Plot PCA results
ggplot(pca_scores_df, aes(x = PC1, y = PC2, label = sample_name)) +  # Create a ggplot object with PC1 and PC2 as axes
  geom_point() +  # Add points to the plot
  geom_text(vjust = -0.5, hjust = 0.5) +  # Add text labels for the sample names
  theme_minimal() +  # Use a minimal theme for the plot
  labs(title = "PCA of Hellinger-Transformed Morphotype Data",  # Add a title to the plot
       x = "Principal Component 1",  # Label the x-axis
       y = "Principal Component 2")  # Label the y-axis