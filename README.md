# as3

#INTRODUCTION 
#Carnivora constitutes a significant order within Mammalia. Carnivora is one of the few mammalian orders that include terrestrial and aquatic species and appears regularly on all continents. It also has one of the largest size ranges of any mammalian order, from 35 g to 5000 kg (Nyakatura & Bininda-Emonds, 2012). The order Carnivora has two clades: Feliformia (cat-like carnivorans) and Caniformia (dog-like carnivorans). The size and diet of carnivores vary widely, from little stoats to bears and from meat-eating like a cheetah to vegetable-eating like a panda. Carnivores also have a broad range of habitat tolerance. They may be found in all terrestrial and aquatic environments, from the domestic cat (Felis catus) to the great white shark (Carcharodon carcharias), from the African wild dog (Lycaon pictus) to the polar bear (Ursus maritimus). Even though the order Carnivora originated in North America roughly 42 million years ago from members of the family Miacidae, its species presently inhabit a variety of continents, latitudes, and longitudes. The question is whether America still has the highest abundance of carnivores among the five continents of the Americas, Africa, Asia, Europe, and Oceania.
#A study was conducted in 2005 on the evolution of eight Carnivora families that have migrated across the Northern Hemisphere to determine if continental invasions also lead to explosive diversification dynamics, using a Bayesian approach to estimate speciation and extinction rates from a large dataset of fossil occurrences. The radiation of arriving lineages from North America to Eurasia correlates with the loss of existing lineages or stages of climatic change, indicating that ecological variations across continents may be responsible for the different rates of diversification. The fossil record also implies that North America and Eurasia are involved in faunal exchanges (Pires et al., 2015). Due to the migration of carnivores across the three continents, which one of the two Eurasian continents will have a species composition most comparable to that of the continent of origin (America)?

#choose taxonomic group and input the Carnivora data 
library(readr)
Carnivora <- read_tsv("https://www.boldsystems.org/index.php/API_Public/combined?taxon=Carnivora&format=tsv")
write_tsv(Carnivora, "Carnivora_BOLD_data.tsv")
Carnivora <- read_tsv("Carnivora_BOLD_data.tsv")

#Part1: Does America continent have the highest abundance of carnivores among the five continents of the Americas, Africa, Asia, Europe, and Oceania?
install.packages("tidyverse")
library(tidyverse)
#since the table includes 80 columns, only main columns including sampleid, family_name, spp_name, species_name, country, bin_uri are extract
data <- data.frame(sampleid=Carnivora$sampleid,family_name=Carnivora$family_name, spp_name=Carnivora$species_name, bin_uri=Carnivora$bin_uri, country=Carnivora$country)
view(data)

#country and bin_uri are main objective for this test -> remove all rows that contains NA values in these 2 rows
install.packages("tidyr")
library(tidyr)
data.filtered <- data %>% 
  drop_na(bin_uri, country)
View(data.filtered)

#countrycode package converts "country" from data.filtered into continents 
install.packages("countrycode") 
library("countrycode") 
#continent is created based on data.filtered df created earlier, the command will find the country names and convert them into compatible continent. 
continent <- countrycode(data.filtered$country, origin = "country.name" , destination = "continent")

#return a dataframe for further analysis
df.continent <- as.data.frame(continent)
view(df.continent)


#calculate the length of country column to check if it is reasonable to remove some unrelated data
length(data.filtered$country) #there are 2381 rows in total

#the testing continents include Americas, Africa, Asia, Europe, and Oceania. Even though Antarctica is also a continent, it is more general than the 5 listed because it includes both Europe and Oceania continents.
#using sum to check the occurrence of "Antarctica" in "country" column 
sum(data.filtered$country == "Antarctica") #there are only 3 samples from Antartica -> the elements are acceptable be removed.

#other excluded rows: Exception - Culture, Exception - Zoological Park, North Pacific Ocean - the reason to be excluded is because they do not show the country -> the package cannot classify them into appropriate continents 
sum(data.filtered$country == "Exception - Culture") #there are only 7 samples from Exception-Culture -> the elements are acceptable to be removed 
sum(data.filtered$country == "Exception - Zoological Park") #there are only 9 samples from Exception-Zoological Park -> the elements are acceptable to be removed 

#similar to Antarctia, North Pacific Ocean includes two continents, Asia and Oceania -> the samples origin is not specific enough  
sum(data.filtered$country == "North Pacific Ocean") #there are only 1 samples from North Pacific Ocean -> the elements are acceptable to be removed 

#append column using tible -> add continent column next to existing data.filtered table to create a complete data set 
data.filtered.update <- data.filtered %>% 
  add_column(df.continent)
view(data.filtered.update)

#use table() function to find the continent that has the highest abundance of species 
freq_table <- data.frame(table(data.filtered.update$continent))
view(freq_table)

#using barplot to visualize frequency table
barplot(Freq ~ Var1, freq_table, col=c("cyan3"), xlab="Continent", ylab="Abundance", main = "Figure 1. Barplot of Carnivora abundance among the five main continents") 

#############################################################
#Part2: Species accumulation curve and data filtering 
library(vegan)

BINs.by.continent <- data.filtered.update %>%
  group_by(continent, bin_uri) %>%    #extract only "continent" and "bin_uri" columns 
  count(bin_uri)  #count the number of individuals in each BIN/species
view(BINs.by.continent) 

#remove all rows that have NA values in continent because further analysis require numeric values only
BINs.by.continent.update <- BINs.by.continent %>% 
  filter(!is.na(continent))

#converting the data to a table of species variables and continents, each observations represents the number of individuals of a species in a continent. 
BINs.spread.by.continent <- pivot_wider(data = BINs.by.continent.update, names_from = bin_uri, values_from = n)
view(BINs.spread.by.continent)

#NA value means that species is not available in that continent. Further analysis requires numeric values -> NA converts to 0 
BINs.spread.by.continent[is.na(BINs.spread.by.continent)] <- 0

#remove rowname header, convert the column names into row names for further analysis because it 'x' must be numeric
BINs.spread.by.continent.update <- BINs.spread.by.continent %>%
  remove_rownames %>%
  column_to_rownames(var = "continent")

#use specaccum to find the number of species at continent sites 
BINs.spread.by.continent.update.accum <- specaccum(BINs.spread.by.continent.update) #one error, Oceania element does not contain sufficient number of observations; detailed explain in discussion section
view(BINs.spread.by.continent.update)
 
#plot the accumulation curve to comprehend the BIN composition of survey plots and predict species abundance.
plot(BINs.spread.by.continent.update.accum, xlab = "Continents Sampled", ylab = "BIN Richness", main = "Site-Based Accumulation Curve, with Continents as Sites")

#############################################################
#Part 3: Europe-Americas correlation vs Asia-Americas correlation

#filter all observations that have the continent equal "Americas", "Europe" and "Asia" for the last test
BINs.spread.by.continent.subset <- BINs.spread.by.continent %>% 
filter(continent == "Europe" | continent == "Americas" | continent == "Asia")


#There are some species that's not available in these 3 continents -> update the data by checking if there is at least one observation available for each variable
BINs.spread.by.continent.subset <- BINs.spread.by.continent.subset[, colSums(BINs.spread.by.continent.subset != 0) > 0] %>% 
remove_rownames() %>% 
column_to_rownames(var="continent")

install.packages("FactoMineR")
install.packages("factoextra")
library(FactoMineR)
library(factoextra)

#aim to test for correlation between continents -> transpose df to new tibble 
t_BINs.spread.by.continent.subset <- t(BINs.spread.by.continent.subset)

#perform Principal Component Analysis to understand the correlations between continents
res.pca <- PCA(t_BINs.spread.by.continent.subset, graph =FALSE)

#print res.pca
print(res.pca)

#create a pca-biplot using res.pca data 
fviz_pca_biplot(res.pca, label="var") + theme_minimal() + labs(title="PCA-Biplot of Americas, Europe, and Asia continents", x="PC1", y="PC2")

# RESULT AND DISCUSSION
# Given that North America is the origin of Carnivora order, one of the predictions for this experiment is that carnivores will be most abundant there. With its long history and extensive exposure, the order expects to be more able to adapt to the environment of the American continent than others. In addition, a study by Pires et al. revealed a continental faunal of carnivores between Eurasia and North America (2015). Therefore, it is reasonable to expect North America, Europe, and Asia to contain the highest proportion of carnivores. The Americas, Europe, and Asia have the highest abundance of carnivores, respectively, among the five continents shown in Figure 1. The predictions were accurate.
# The accumulation curve shown in Figure 2 aims to comprehend the BIN composition of survey plots and predict species abundance. For several reasons, the accumulation curve does not provide the optimal outcome. First, the data set has a small size. It is acceptable if certain species are unavailable on a specific continent. However, even if they are available on one continent, the numbers of variables only contain from one to five observations per species. Additionally, Oceania contains just 101 observations, of which 100 correspond to the first variable, BOLD:AAA1542. The data was insufficient. Hence the algorithm was unable to display the Oceania data on the accumulation curve and provide a meaningful prediction for this experiment.
# According to the first two tests, America and Europe have the highest abundance, so it is reasonable to assume that these two variables will have a stronger correlation than America and Asia. The correlations between continents were explored using a PCA test. The findings demonstrate that the Americas and Asia are strongly correlated, but Europe has a lower correlation with the other variables. This finding contrasts with the expectations, indicating that even though one continent has higher individual numbers, it does not necessarily have a higher probability of possessing the same BINs as another continent. Even if the data set is insufficient for a prediction test, it is adequate to answer all questions listed. A few requirements must be applied to a data set before it can qualify for a larger project. First, the number of observations after filtering must be raised to at least 30 per species on a single continent to be considered statistical significance. Second, the order Carnivora includes 296 species categorized into 16 families (Hassanin et al., 2021), yet the data set only includes 143 species. It implies that not only must the number of observations be raised, but also the number of variables. A larger project could perhaps compare the BIN composition of each continent now with that of 50 years ago to see whether climate change and global warming have a role in the changes and, if so, to determine if the changes are statistically significant.

# REFERENCES
# Hassanin, A., Veron, G., Ropiquet, A., Jansen van Vuuren, B., Lécu, A., Goodman, S. M., Haider, J., & Nguyen, T. T. (2021, February 16). Evolutionary history of Carnivora (Mammalia, Laurasiatheria) inferred from mitochondrial genomes. PLOS ONE, 16(2), e0240770. https://doi.org/10.1371/journal.pone.0240770
# Noonan, M. J., Newman, C., Buesching, C. D., & Macdonald, D. W. (2015, October 13). Evolution and function of fossoriality in the Carnivora: implications for group-living. Frontiers in Ecology and Evolution, 3. https://doi.org/10.3389/fevo.2015.00116
# Nyakatura, K., & Bininda-Emonds, O. R. (2012, February 27). Updating the evolutionary history of Carnivora (Mammalia): a new species-level supertree complete with divergence time estimates. BMC Biology, 10(1). https://doi.org/10.1186/1741-7007-10-12
# Pires, M. M., Silvestro, D., & Quental, T. B. (2015, October 22). Continental faunal exchange and the asymmetrical radiation of carnivores. Proceedings of the Royal Society B: Biological Sciences, 282(1817), 20151952. https://doi.org/10.1098/rspb.2015.1952

