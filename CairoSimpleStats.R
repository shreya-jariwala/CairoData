library(tidyverse)
library(stringr)

data <- as_tibble(read.csv("CairoData_clean.csv")) #read in product of 'dfSplitMorphoColumns.R'

#set categorical variables as factor types:
data <- data %>% 
  mutate(Morpho1=as.factor(Morpho1)) %>%
  mutate(Morpho2=as.factor(Morpho2)) %>%
  mutate(Morpho3=as.factor(Morpho3)) %>%
  mutate(Preservation=as.factor(Preservation)) %>% 
  mutate(microsite_name=as.factor(microsite_name)) %>% 
  mutate(rock_name=as.factor(rock_name)) %>% 
  mutate(sample_name=as.factor(sample_name))

hist(data$strat_height) #how are specimens distributed stratigraphically?

Morpho2hist <- ggplot(data, aes(x=fct_infreq(Morpho2))) + geom_bar() #Distribution of morphotypes
Morpho2hist

Morpho2hist_by_site <- Morpho2hist + facet_wrap(~microsite_name, scales = 'free') #Distribution of morphotypes split by microsite
Morpho2hist_by_site

Morpho2hist_by_rock <- Morpho2hist + facet_wrap(~rock_name, scales = 'free') #Distribution of morphotypes split by rock sample
Morpho2hist_by_rock

#Counts and proportions of different fossil types:
data <- data %>% 
  mutate(Is.Arthropod = (Morpho1=="Ap")|(Morpho1=="Tr")|(Morpho1=="Sp")|(Morpho1=="Va")) %>% #arthropod
  mutate(Is.Concho = (Morpho2=="Va1")|(Morpho2=="Va2")|(Morpho2=="Va3")|(Morpho2=="Va4")|(Morpho2=="Va5")| Morpho3=="Va") %>% #clam shrimp
  mutate(Is.Cairocaris = (Morpho2=="ApB")|(Morpho3=="ApT2")|(Morpho3=="ApT3")|(Morpho3=="ApT4")|(Morpho3=="ApT7")|(Morpho3=="TrA1")|(Morpho3=="TrA2")|(Morpho3=="TrS2")|(Morpho3=="TrP1")|(Morpho3=="SpS1")) %>% #Cairocaris
  mutate(Is.TerrAr = (Morpho3=="ApT5")|(Morpho2=="ApW")|(Morpho3=="ApT6")|(Morpho3=="SpT3")) %>% #terrestrial arthropod 
  mutate(Is.Euryp = (Morpho3=="ApT8")|(Morpho3=="ApT13")|(Morpho3=="ApT14")) #eurypterid 
#prepare columns for a table summarizing the fossil type statistics:
categories <- c("Arthropod", "Clam Shrimp", "Cairocaris", "Eurypterid", "TerrArthr", "TerrPlants", "EggSacs", "Spore-likes", "Frisbeekites")
counts <- c(sum(data$Is.Arthropod), sum(data$Is.Concho), sum(data$Is.Cairocaris), sum(data$Is.Euryp), sum(data$Is.TerrAr), sum(data$Morpho2=="BrF",data$Morpho2=="PhS",data$Morpho2=="PhI"), sum(data$Morpho3=="OrS1"), sum(data$Morpho2=="OrR"), sum(data$Morpho3=="OrF1"))
proportions <- counts/as.numeric(count(data))

#create tables summarizing overall and arthropod-specific stats:
FossilTypes <- tibble(Category = categories, Count = counts, Proportion = proportions)
ArthropodTypes <- tibble(Category = categories[2:5], Count = counts[2:5], Proportion = counts[2:5]/counts[1])
ArthropodTypes <- ArthropodTypes %>%
  add_row(Category = "Other", Count = counts[1]-sum(counts[2:5]), Proportion = (counts[1]-sum(counts[2:5]))/counts[1])

#make a  pie chart showing arthropod proportions:
ArthropodTypes <- ArthropodTypes %>% 
  arrange(desc(Category)) %>%
  mutate(percent = Proportion *100) %>%
  mutate(ypos = cumsum(percent)- 0.5*percent )
ArthroBreakdown <- ggplot(ArthropodTypes, aes(x="",y=Count, fill=Category)) + 
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0) + theme_void() + theme(legend.position = "none") +
  geom_text(aes(x=1.7, label=Category),position = position_stack(vjust = 0.4),color="black",size=4.5) + scale_fill_brewer(palette = "Set1")
ArthroBreakdown

#make the above summaries for each individual site (CSA, CSC, CSD, CSJ):
arthropods.by.site <- data %>% group_by(microsite_name) %>% summarise(Count = sum(Is.Arthropod))
conchos.by.site <- data %>% group_by(microsite_name) %>% summarise(Count = sum(Is.Concho))
cairocaris.by.site <- data %>% group_by(microsite_name) %>% summarise(Count = sum(Is.Cairocaris))
euryp.by.site <- data %>% group_by(microsite_name) %>% summarise(Count = sum(Is.Euryp))
terrar.by.site <- data %>% group_by(microsite_name) %>% summarise(Count = sum(Is.TerrAr))
terrpl.by.site <- data %>% group_by(microsite_name) %>% summarise(Count = sum(Morpho2=="BrF", Morpho2=="PhS", Morpho2=="PhI"))
eggsac.by.site <- data %>% group_by(microsite_name) %>% summarise(Count = sum(Morpho3=="OrS1"))
spore.by.site <- data %>% group_by(microsite_name) %>% summarise(Count = sum(Morpho2=="OrR"))
frisbee.by.site <- data %>% group_by(microsite_name) %>% summarise(Count = sum(Morpho3=="OrF1"))

site_names <- levels(data$microsite_name)
FossilTypes.by.site <- tibble(Site=factor(), Category=factor(), Count=numeric(), Proportion=numeric())
ArthropodTypes.by.site <- tibble(Site=factor(), Category=factor(), Count=numeric(), Proportion=numeric())
for (i in 1:length(site_names)) {
  counts.tmp <- c(arthropods.by.site$Count[i], conchos.by.site$Count[i], cairocaris.by.site$Count[i], euryp.by.site$Count[i], terrar.by.site$Count[i], terrpl.by.site$Count[i], eggsac.by.site$Count[i], spore.by.site$Count[i], frisbee.by.site$Count[i])
  
  temp.list <- tibble(Category = categories, Count = counts.tmp, Proportion = counts.tmp/as.numeric(sum(data$microsite_name==site_names[i], na.rm = TRUE)))
  temp.list <- temp.list %>% mutate(Site=site_names[i], .before=Category)
  FossilTypes.by.site <- FossilTypes.by.site %>% add_row(temp.list)
  #name <- paste("FossilTypes", site_names[i], sep = "_")
  #assign(name, temp.list)
  
  temp.list <- tibble(Category = categories[2:5], Count = counts.tmp[2:5], Proportion = counts.tmp[2:5]/counts.tmp[1])
  temp.list <- temp.list %>% mutate(Site=site_names[i], .before=Category)
  ArthropodTypes.by.site <- ArthropodTypes.by.site %>% add_row(temp.list)
  ArthropodTypes.by.site <- ArthropodTypes.by.site %>%
    add_row(Site = site_names[i], Category = "Other", Count = counts.tmp[1]-sum(counts.tmp[2:5]), Proportion = (counts.tmp[1]-sum(counts.tmp[2:5]))/counts.tmp[1])
  #name <- paste("ArthropodTypes", site_names[i], sep = "_")
  #assign(name, temp.list)
  
  #rm(name)
  rm(temp.list)
}

ArthropodTypes.by.site <- ArthropodTypes.by.site %>% 
  arrange(desc(Category)) %>%
  mutate(percent = Proportion *100) %>%
  mutate(ypos = cumsum(percent)- 0.5*percent )
ArthroBreakdown.by.site <- ggplot(ArthropodTypes.by.site, aes(x="",y=Count, fill=Category)) + 
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0) + theme_void() + theme(legend.position = "none") +
  geom_text(aes(x=1.7, label=Category),position = position_stack(vjust = 0.4),color="black",size=4.5) + scale_fill_brewer(palette = "Set1") +
  facet_wrap(~Site, scales = 'free')
ArthroBreakdown.by.site

#make the above summaries for individual rock samples:
arthropods.by.rock <- data %>% group_by(rock_name) %>% summarise(Count = sum(Is.Arthropod))
conchos.by.rock <- data %>% group_by(rock_name) %>% summarise(Count = sum(Is.Concho))
cairocaris.by.rock <- data %>% group_by(rock_name) %>% summarise(Count = sum(Is.Cairocaris))
euryp.by.rock <- data %>% group_by(rock_name) %>% summarise(Count = sum(Is.Euryp))
terrar.by.rock <- data %>% group_by(rock_name) %>% summarise(Count = sum(Is.TerrAr))
terrpl.by.rock <- data %>% group_by(rock_name) %>% summarise(Count = sum(Morpho2=="BrF", Morpho2=="PhS", Morpho2=="PhI"))
eggsac.by.rock <- data %>% group_by(rock_name) %>% summarise(Count = sum(Morpho3=="OrS1"))
spore.by.rock <- data %>% group_by(rock_name) %>% summarise(Count = sum(Morpho2=="OrR"))
frisbee.by.rock<- data %>% group_by(rock_name) %>% summarise(Count = sum(Morpho3=="OrF1"))

rock_names <- levels(data$rock_name)
FossilTypes.by.rock <- tibble(Rock=factor(), Category=factor(), Count=numeric(), Proportion=numeric())
ArthropodTypes.by.rock <- tibble(Rock=factor(), Category=factor(), Count=numeric(), Proportion=numeric())
for (i in 1:length(rock_names)) {
  counts.tmp <- c(arthropods.by.rock$Count[i], conchos.by.rock$Count[i], cairocaris.by.rock$Count[i], euryp.by.rock$Count[i], terrar.by.rock$Count[i], terrpl.by.rock$Count[i], eggsac.by.rock$Count[i], spore.by.rock$Count[i], frisbee.by.rock$Count[i])
  
  temp.list <- tibble(Category = categories, Count = counts.tmp, Proportion = counts.tmp/as.numeric(sum(data$rock_name==rock_names[i], na.rm = TRUE)))
  temp.list <- temp.list %>% mutate(Rock=rock_names[i], .before=Category)
  FossilTypes.by.rock <- FossilTypes.by.rock %>% add_row(temp.list)
  
  temp.list <- tibble(Category = categories[2:5], Count = counts.tmp[2:5], Proportion = counts.tmp[2:5]/counts.tmp[1])
  temp.list <- temp.list %>% mutate(Rock=rock_names[i], .before=Category)
  ArthropodTypes.by.rock <- ArthropodTypes.by.rock %>% add_row(temp.list)
  ArthropodTypes.by.rock <- ArthropodTypes.by.rock %>%
    add_row(Rock=rock_names[i], Category = "Other", Count = counts.tmp[1]-sum(counts.tmp[2:5]), Proportion = (counts.tmp[1]-sum(counts.tmp[2:5]))/counts.tmp[1])

  rm(temp.list)
}

ArthropodTypes.by.rock <- ArthropodTypes.by.rock %>% 
  arrange(desc(Category)) %>%
  mutate(percent = Proportion *100) %>%
  mutate(ypos = cumsum(percent)- 0.5*percent )
ArthroBreakdown.by.rock <- ggplot(ArthropodTypes.by.rock, aes(x="",y=Count, fill=Category)) + 
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0) + theme_void() + 
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~Rock, scales = 'free')
ArthroBreakdown.by.rock

#Make bins for strat. heights:
data <- data %>% 
  mutate(Strat.Bin = cut(data$strat_height, 7, ordered_result = TRUE), .before = strat_height)

#make summaries for stratigraphic bins:
arthropods.by.strat <- data %>% group_by(Strat.Bin) %>% summarise(Count = sum(Is.Arthropod))
conchos.by.strat <- data %>% group_by(Strat.Bin) %>% summarise(Count = sum(Is.Concho))
cairocaris.by.strat <- data %>% group_by(Strat.Bin) %>% summarise(Count = sum(Is.Cairocaris))
euryp.by.strat <- data %>% group_by(Strat.Bin) %>% summarise(Count = sum(Is.Euryp))
terrar.by.strat <- data %>% group_by(Strat.Bin) %>% summarise(Count = sum(Is.TerrAr))
terrpl.by.strat <- data %>% group_by(Strat.Bin) %>% summarise(Count = sum(Morpho2=="BrF", Morpho2=="PhS", Morpho2=="PhI"))
eggsac.by.strat <- data %>% group_by(Strat.Bin) %>% summarise(Count = sum(Morpho3=="OrS1"))
spore.by.strat <- data %>% group_by(Strat.Bin) %>% summarise(Count = sum(Morpho2=="OrR"))
frisbee.by.strat <- data %>% group_by(Strat.Bin) %>% summarise(Count = sum(Morpho3=="OrF1"))

strat_bins <- levels(data$Strat.Bin)
FossilTypes.by.strat <- tibble(Height=factor(), Category=factor(), Count=numeric(), Proportion=numeric())
ArthropodTypes.by.strat <- tibble(Height=factor(), Category=factor(), Count=numeric(), Proportion=numeric())
for (i in 1:(length(strat_bins)-1)) {
  counts.tmp <- c(arthropods.by.strat$Count[i], conchos.by.strat$Count[i], cairocaris.by.strat$Count[i], euryp.by.strat$Count[i], terrar.by.strat$Count[i], 
                  terrpl.by.strat$Count[i], eggsac.by.strat$Count[i], spore.by.strat$Count[i], frisbee.by.strat$Count[i])
  
  temp.list <- tibble(Category = categories, Count = counts.tmp, Proportion = counts.tmp/as.numeric(sum(data$rock_name==strat_bins[i], na.rm = TRUE)))
  temp.list <- temp.list %>% mutate(Height=strat_bins[i], .before=Category)
  FossilTypes.by.strat <- FossilTypes.by.strat %>% add_row(temp.list)
  
  temp.list <- tibble(Category = categories[2:5], Count = counts.tmp[2:5], Proportion = counts.tmp[2:5]/counts.tmp[1])
  temp.list <- temp.list %>% mutate(Height=strat_bins[i], .before=Category)
  ArthropodTypes.by.strat <- ArthropodTypes.by.strat %>% add_row(temp.list)
  ArthropodTypes.by.strat <- ArthropodTypes.by.strat %>%
    add_row(Height=strat_bins[i], Category = "Other", Count = counts.tmp[1]-sum(counts.tmp[2:5]), Proportion = (counts.tmp[1]-sum(counts.tmp[2:5]))/counts.tmp[1])
  
  rm(temp.list)
}

ArthropodTypes.by.strat <- ArthropodTypes.by.strat %>% 
  arrange(desc(Height)) %>%
  mutate(percent = Proportion *100) %>%
  mutate(ypos = cumsum(percent)- 0.5*percent )
ArthroBreakdown.by.strat <- ggplot(ArthropodTypes.by.strat, aes(x="",y=Count, fill=Category)) + 
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0) + theme_void() + 
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~Height, scales = 'free')
ArthroBreakdown.by.strat

#make subset table for samples from Susan's block of layers ("CLST") with additional 'layer' column
data_CLST <- data %>%
  filter(str_sub(sample_name, 1, 4)=="CLST")
data_CLST <- data_CLST %>%
  mutate(layer = as.numeric(str_sub(sample_name, 6)), .before = strat_height)
