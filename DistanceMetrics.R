library(tidyverse)
library(vegan)

#Raw abundances (morphotype counts) -- use Bray-Curtis Distance

data_count <- as.data.frame(read_csv("Morpho2_counts_byrock.csv")) #import data

rocknames <- as.vector(data_count$rock_name) #make vector of rock names

data_count <- data_count %>% #sort by rock names for consistency
  arrange(rocknames)

rownames(data_count) <- data_count[,1] #set rock names to row names and delete first column
data_count <- data_count[,-1] 

data_count <- as.matrix(data_count) #format as matrix for input into vegan

dist_bray <- vegdist(data_count, method = "bray")

##visualize distances with NMDS -- not really that useful but could be interesting if you get ggplot to label the rock names
nmds <- metaMDS(dist_bray)
scores(nmds) %>%    
  as_tibble(rownames = "Rock") %>%
  ggplot(aes(x=NMDS1,y=NMDS2)) +
  geom_point()

dist_bray_df <- as.data.frame(as.matrix(dist_bray))
colnames(dist_bray_df) <- rownames(dist_bray_df)
write.csv(dist_bray_df, file = "bray_curtis.csv", row.names = TRUE)

##----------------------------------------------#

#Proportions (morphotype frequencies) -- use Hellinger Distance

data_freq <- as.data.frame(read_csv("Morpho2_freqs_byrock.csv")) #import data

rocknames <- as.vector(data_freq$rock_name) #make vector of rock names

data_freq <- data_freq %>% #sort by rock names for consistency
  arrange(rocknames)

rownames(data_freq) <- data_freq[,1] #set rock names to row names and delete first column
data_freq <- data_freq[,-1] 

data_freq <- as.matrix(data_freq) #format as matrix for input into vegan

dist_hel<- vegdist(data_freq, method = "hellinger")

dist_hel_df <- as.data.frame(as.matrix(dist_hel))
colnames(dist_hel_df) <- rownames(dist_hel_df)
write.csv(dist_hel_df, file = "hellinger.csv", row.names = TRUE)

########################################################
#Code for making pairwise stratigraphic matrix:

data_strat <- as_tibble(read.csv("dataStratFiltered.csv")) #read in the wide dataframe that was filtered for rocks with strat. info
data_strat <- data_strat %>%
  arrange(rock_name)
#Make table with every rock name and its corresponding stratigraphic height.
rock_heights <- data_strat %>% 
  group_by(rock_name) %>%
  distinct(strat_height)

rocksnames_strat <- as.vector(rock_heights$rock_name)

rock_heights2 <- as.matrix(rock_heights)
rock_heights2 <- as.numeric(rock_heights2[,2])

strat_pairwise <- vegdist(rock_heights2, method = "euclidean")
strat_pairwise <- as.matrix(strat_pairwise)
rownames(strat_pairwise) <- rocksnames_strat
colnames(strat_pairwise) <- rocksnames_strat
######################################
#Code similar to first section (make an ecological distance matrix) but filtered for samples w/ stratigraphic info
data_count_strat <- as.data.frame(read_csv("Morpho2_counts_byrock_Strat.csv")) #import data

rocknames_strat <- as.vector(data_count_strat$rock_name) #make vector of rock names

data_count_strat <- data_count_strat %>% #sort by rock names for consistency
  arrange(rocknames_strat)

rownames(data_count_strat) <- data_count_strat[,1] #set rock names to row names and delete first column
data_count_strat <- data_count_strat[,-1] 

data_count_strat <- as.matrix(data_count_strat) #format as matrix for input into vegan

dist_bray_strat <- vegdist(data_count_strat, method = "bray")

dist_bray_df_strat <- as.data.frame(as.matrix(dist_bray_strat))
colnames(dist_bray_df_strat) <- rownames(dist_bray_df_strat)
write.csv(dist_bray_df_strat, file = "bray_curtis_StratFiltered.csv", row.names = TRUE)

dist_bray_strat <- as.matrix(dist_bray_strat)
rownames(dist_bray_strat) <- rocksnames_strat
colnames(dist_bray_strat) <- rocksnames_strat

#same for Hellinger:
data_freq_strat <- as.data.frame(read_csv("Morpho2_freqs_byrock_Strat.csv")) #import data
rocknames_strat <- as.vector(data_freq_strat$rock_name) #make vector of rock names

data_freq_strat <- data_freq_strat %>% #sort by rock names for consistency
  arrange(rocknames_strat)

rownames(data_freq_strat) <- data_freq_strat[,1] #set rock names to row names and delete first column
data_freq_strat <- data_freq_strat[,-1] 

data_freq_strat <- as.matrix(data_freq_strat) #format as matrix for input into vegan

dist_hel_strat <- vegdist(data_freq_strat, method = "hellinger")

dist_hel_df_strat <- as.data.frame(as.matrix(dist_hel_strat))
colnames(dist_hel_df_strat) <- rownames(dist_hel_df_strat)
write.csv(dist_hel_df_strat, file = "hellinger_StratFiltered.csv", row.names = TRUE)

dist_hel_strat <- as.matrix(dist_hel_strat)
rownames(dist_hel_strat) <- rocksnames_strat
colnames(dist_hel_strat) <- rocksnames_strat

#####MANTEL TEST FOR STRATIGRAPHIC HEIGHTS###################

mantel_strat_count <- mantel(dist_bray_strat, strat_pairwise, method = "pearson", permutations = 999)
print(mantel_strat_count)

mantel_strat_freq <- mantel(dist_hel_strat, strat_pairwise, method = "pearson", permutations = 999)
print(mantel_strat_freq)

