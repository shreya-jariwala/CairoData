library(tidyverse)
library(stringr)

df <- read.csv("CairoData_May29_2024.csv")

df <- df %>% 
  rename(Morpho3=Morphotype) #rename Morphotype column

df <- df %>% 
  mutate(Morpho2=str_trunc(Morpho3, 3, side = "right", ellipsis = "")) #new column with only 2 broader categories

df <- df %>% 
  mutate(Morpho1=str_trunc(Morpho3, 2, side = "right", ellipsis = "")) #new column with only broadest category

write_csv(df, "CairoData_clean.csv")

