# as3

#choose taxonomic group and input the Carnivora data 

library(readr)

Carnivora <- read_tsv("https://www.boldsystems.org/index.php/API_Public/combined?taxon=Carnivora&format=tsv")

write_tsv(Carnivora, "Carnivora_BOLD_data.tsv")

Carnivora <- read_tsv("Carnivora_BOLD_data.tsv")

#Part1: Does America continent have the highest abundance of carnivores among the five continents of the Americas, Africa, Asia, Europe, and Oceania?

#install.packages("tidyverse")

library(tidyverse)

#since the table includes 80 columns, only main columns including sampleid, family_name, spp_name, species_name, country, bin_uri are extract

data <- data.frame(sampleid=Carnivora$sampleid,family_name=Carnivora$family_name, spp_name=Carnivora$species_name, bin_uri=Carnivora$bin_uri, country=Carnivora$country)

view(data)

#country and bin_uri are main objective for this test -> remove all rows that contains NA values in these 2 rows

#install.packages("tidyr")

library(tidyr)

data.filtered <- data %>% 
  drop_na(bin_uri, country)

#countrycode package converts "country" from data.filtered into continents 

install.packages("countrycode") 

library("countrycode") 

#continent is created based on data.filtered df created earlier, the command will find the country names and convert them into compatible continent. 

continent <- countrycode(data.filtered$country, origin = "country.name" , destination = "continent")

#return a dataframe for further analysis

df.continent <- as.data.frame(continent)

#calculate the length of country column to check if it is reasonable to remove some unrelated data

length(data.filtered$country) #there are 2381 rows in total

#the testing continents include Americas, Africa, Asia, Europe, and Oceania. Even though Antarctica is also a continent, it is more general than the 5 listed because it includes both Europe and Oceania continents.

#using sum to check the occurrence of "Antarctica" in "country" column 

sum(data.filtered$country == "Antarctica") #there are only 3 samples from Antartica -> the elements are acceptable be removed.

#other excluded rows: Exception - Culture, Exception - Zoological Park, North Pacific Ocean - the reason to be excluded is because they do not show the 

country -> the package cannot classify them into appropriate continents 

sum(data.filtered$country == "Exception - Culture") #there are only 7 samples from Exception-Culture -> the elements are acceptable to be removed 

sum(data.filtered$country == "Exception - Zoological Park") #there are only 9 samples from Exception-Zoological Park -> the elements are acceptable to be removed 

#similar to Antarctia, North Pacific Ocean includes two continents, Asia and Oceania -> the samples origin is not specific enough  

sum(data.filtered$country == "North Pacific Ocean") #there are only 1 samples from North Pacific Ocean -> the elements are acceptable to be removed 

#append column using tible -> add continent column next to existing data.filtered table to create a complete data set 

data.filtered.update <- data.filtered %>% 
        add_column(df.continent)

#use table() function to find the continent that has the highest abundance of species 

freq_table <- data.frame(table(data.filtered.update$continent))

#using barplot to visualize frequency table

barplot(Freq ~ Var1, freq_table, col=c("cyan3"), xlab="Continent", ylab="Abundance", main = "Figure 1. Barplot of Carnivora abundance among the five main continents") 

#############################################################

#Part2: Species accumulation curve and data filtering 

library(vegan)

BINs.by.continent <- data.filtered.update %>%
  group_by(continent, bin_uri) %>%    #extract only "continent" and "bin_uri" columns 
  count(bin_uri)  #count the number of individuals in each BIN/species


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

#install.packages("FactoMineR")

#install.packages("factoextra")

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

#REFERENCES

#Hassanin, A., Veron, G., Ropiquet, A., Jansen van Vuuren, B., Lécu, A., Goodman, S. M., Haider, J., & Nguyen, T. T. (2021, February 16). Evolutionary history of Carnivora (Mammalia, Laurasiatheria) inferred from mitochondrial genomes. PLOS ONE, 16(2), e0240770. https://doi.org/10.1371/journal.pone.0240770

#Noonan, M. J., Newman, C., Buesching, C. D., & Macdonald, D. W. (2015, October 13). Evolution and function of fossoriality in the Carnivora: implications for group-living. Frontiers in Ecology and Evolution, 3. https://doi.org/10.3389/fevo.2015.00116

#Nyakatura, K., & Bininda-Emonds, O. R. (2012, February 27). Updating the evolutionary history of Carnivora (Mammalia): a new species-level supertree complete with divergence time estimates. BMC Biology, 10(1). https://doi.org/10.1186/1741-7007-10-12

#Pires, M. M., Silvestro, D., & Quental, T. B. (2015, October 22). Continental faunal exchange and the asymmetrical radiation of carnivores. Proceedings of the Royal Society B: Biological Sciences, 282(1817), 20151952. https://doi.org/10.1098/rspb.2015.1952

