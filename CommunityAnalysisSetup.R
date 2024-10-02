###-----------------------------------------------------------------------------------------------------###
##          Purpose of this script: format data into a table with rows for each rock sample,             ##
##   and columns for proportional counts of each Morpho2 category, to plug into statistical analyses.    ##
###-----------------------------------------------------------------------------------------------------###

library(tidyverse)

data <- as_tibble(read.csv("CairoData_clean.csv")) #read in product of 'dfSplitMorphoColumns.R'

#make sure R is interpreting the relevant categorical variable columns as factor types:
data <- data %>% 
  mutate(Morpho2=as.factor(Morpho2)) %>%
  mutate(rock_name=as.factor(rock_name))

morphotypes <- levels(data$Morpho2) #make a vector of morphotype names, which will form the relevant columns in the "wide" table format we're going to create

data_wide <- data %>% pivot_wider(names_from = Morpho2, values_from = no_specimens, values_fn = list) #breaks up Morpho2 column into separate columns for each category

##Make our final tables: group specimens from the same rock sample, and summarize each morphotype column...
#as raw counts for rock:
community.counts <- data_wide %>% 
                      group_by(rock_name) %>%
                      summarise(across(morphotypes, ~ sum(unlist(.x))))
#and as proportions of total specimens from each rock:
community.freqs <- data_wide %>% 
                      group_by(rock_name) %>%
                      summarise(across(morphotypes, ~ sum(unlist(.x))/n()))

#Output tables as .csv files:
write_csv(community.counts, "Morpho2_counts_byrock.csv")
write_csv(community.freqs, "Morpho2_freqs_byrock.csv")

##Make count and frequency tables filtering out rocks without strat. measurements:
has.strat <- data_wide %>%
  filter(is.na(strat_height)==FALSE)

write_csv(has.strat, "dataStratFiltered.csv")

#as raw counts for rock:
community.counts.strat <- has.strat %>% 
  group_by(rock_name) %>%
  summarise(across(morphotypes, ~ sum(unlist(.x))))

#and as proportions of total specimens from each rock:
community.freqs.strat <- has.strat %>% 
  group_by(rock_name) %>%
  summarise(across(morphotypes, ~ sum(unlist(.x))/n()))

#Output tables as .csv files:
write_csv(community.counts.strat, "Morpho2_counts_byrock_Strat.csv")
write_csv(community.freqs.strat, "Morpho2_freqs_byrock_Strat.csv")
