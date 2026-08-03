# Per trap type statistics, per Zone statistics

library(pracma)
library(tidyverse)
library(lubridate)

folderDataLocal = "Data"
folderOutput = "Outputs"

# load data
totDFmod <- readRDS(file = paste0(folderDataLocal, "/totDFmod_ElDorado_Sepulveda.rds"))

# Whole data ----

# trap type
traps = unique(totDFmod$TrapType)
species = unique(totDFmod$Species)
species = species[-which(is.na(species))]

#df correspondence trap-species
correspondenceTrapDF = data.frame(species = rep(species, times = length(traps)),
                                  trapType = rep(traps, each = length(species)),
                                  match = F)

## loop----
# histogram of species per trap 

for(trapi in traps){
  
  cat(trapi, "\n")
  
  histDF <- totDFmod %>%
    filter(!is.na(Species)) %>% 
    filter(TrapType == trapi) %>%
    mutate(GenusSpecies = paste(Genus, Species)) %>%
    group_by(GenusSpecies, Species) %>%
    summarise(totAbundance = sum(AvgAbundance)) %>%
    ungroup() %>%
    mutate(perc = paste0(round(100*totAbundance/sum(totAbundance), 3), "%"))
  
  speciesi = unique(histDF$Species)
  
  for(speciesii in speciesi){
    correspondenceTrapDF$match[which(correspondenceTrapDF$trapType == trapi & correspondenceTrapDF$species %in% speciesii)] = T
  }
  
  ggplot(histDF, aes(x = totAbundance , y = GenusSpecies, label = perc))+
    xlim(c(0, 1.05*max(histDF$totAbundance)))+
    geom_col(stat = "identity")+ 
    geom_text(hjust = -0.1,    # nudge above top of bar
              size = 3)+
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = "gray90"))
  
  ggsave(filename = paste0(folderOutput, "/A - Histogram of species for ", trapi,".png"),
         device = "png", width = 7, height = 5)

}

ggplot(correspondenceTrapDF, aes(x = trapType, y = species, fill = match))+
  geom_tile()
  

ggsave(filename = paste0(folderOutput, "/A - Correspondence table.png"),
       device = "png", width = 7, height = 5)

# Only selected traps and periods----

totDFsel <- readRDS(file = paste0(folderDataLocal, "/totDFsel_ElDorado_Sepulveda.rds"))

# trap type
traps = unique(totDFsel$TrapType)
species = unique(totDFsel$Species)
species = species[-which(is.na(species))]

#df correspondence trap-species
correspondenceTrapDF = data.frame(species = rep(species, times = length(traps)),
                                  trapType = rep(traps, each = length(species)),
                                  match = F)

## loop----
# histogram of species per trap 

for(trapi in traps){
  
  cat(trapi, "\n")
  
  histDF <- totDFsel %>%
    filter(!is.na(Species)) %>% 
    filter(TrapType == trapi) %>%
    mutate(GenusSpecies = paste(Genus, Species)) %>%
    group_by(GenusSpecies, Species) %>%
    summarise(totAbundance = sum(AvgAbundance)) %>%
    ungroup() %>%
    mutate(perc = paste0(round(100*totAbundance/sum(totAbundance), 3), "%"))
  
  speciesi = unique(histDF$Species)
  
  for(speciesii in speciesi){
    correspondenceTrapDF$match[which(correspondenceTrapDF$trapType == trapi & correspondenceTrapDF$species %in% speciesii)] = T
  }
  
  ggplot(histDF, aes(x = totAbundance , y = GenusSpecies, label = perc))+
    xlim(c(0, 1.05*max(histDF$totAbundance)))+
    geom_col(stat = "identity")+ 
    geom_text(hjust = -0.1,    # nudge above top of bar
              size = 3)+
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = "gray90"))
  
  ggsave(filename = paste0(folderOutput, "/D - Histogram of species for ", trapi," after selection.png"),
         device = "png", width = 7, height = 5)
  
}


ggplot(correspondenceTrapDF, aes(x = trapType, y = species, fill = match))+
  geom_tile()


ggsave(filename = paste0(folderOutput, "/D - Correspondence table after selection.png"),
       device = "png", width = 7, height = 5)