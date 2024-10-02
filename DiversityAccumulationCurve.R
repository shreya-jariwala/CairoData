library(dplyr)
library(tidyr)
library(vegan)

data <- read.csv("CairoData_clean.csv")

# Select relevant columns: 'sample_name' and 'Morpho3'
morphotype_data <- data %>%
  select(sample_name, Morpho3) %>%
  
  # Filter out rows where 'Morpho2' is NA (missing values)
  filter(!is.na(Morpho3)) %>%
  
  # Group the data by 'sample_name' and 'Morpho3'
  group_by(sample_name, Morpho3) %>%
  
  # Summarize the data by counting the number of occurrences of each 'Morpho2' in each 'sample_name'
  summarise(Count = n(), .groups = 'drop') %>%
  
  # Pivot the data to a wider format where each unique 'Morpho2' becomes a column
  # Fill missing values with 0 (since some 'Morpho2' types might not be present in some samples)
  pivot_wider(names_from = Morpho3, values_from = Count, values_fill = list(Count = 0))

# Convert to matrix
morphotype_matrix <- as.matrix(morphotype_data[,-1])
rownames(morphotype_matrix) <- morphotype_data$sample_name

# Calculate the minimum sample size in the dataset
min_sample_size <- min(rowSums(morphotype_matrix))

# Ensure the minimum sample size is at least 20
sample <- max(20, min_sample_size)

# Filter out samples with very low counts
min_count_threshold <- 10  # Adjust this threshold as needed
filtered_morphotype_matrix <- morphotype_matrix[rowSums(morphotype_matrix) >= min_count_threshold, ]

# Calculate rarefied species richness for each sample
rarefied_richness <- rarefy(filtered_morphotype_matrix, sample = sample)

# Print the rarefied richness values
print(rarefied_richness)

# Generate rarefaction curves
rare_curves <- rarecurve(filtered_morphotype_matrix, step = 5, sample = sample, lwd = 2, cex = 0.6)

# Add a title to the plot
title(main = "Rarefaction Curves for Each Sample")
