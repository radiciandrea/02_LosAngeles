# Per trap type statistics, per Zone statistics

library(pracma)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)

folderDataLocal = "Data"
folderOutput = "Outputs"

# load data
totDFmod <- readRDS(file = paste0(folderDataLocal, "/totDFmod_ElDorado_Sepulveda.rds"))

# Whole data ----

totDFmod <- totDFmod %>%
  filter(!is.na(Species)) %>%
  mutate(genusSpecies = paste(Genus, Species))

# trap type
traps = unique(totDFmod$TrapType)
genusSpecies = unique(totDFmod$genusSpecies)

#surveillanceEffort per trap
trapDFmod <- totDFmod %>%
  group_by(TrapType) %>%
  summarize(Sites = length(unique(siteCode)),
            Weeks = length(unique(collectionWeek)),
            Effort = sum(totNightTraps))%>%
  ungroup() 

gsites = ggplot(trapDFmod, aes(x = TrapType, y = Sites))+
  geom_col() +
  labs(y = "Spatial records (sites)", x = "Traps" ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  ggtitle("")

gweeks = ggplot(trapDFmod, aes(x = TrapType, y = Weeks))+
  geom_col() +
  labs(y = "Temporal records (weeks)", x = "Traps" ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  ggtitle("Measures of trapping effort")

geffort = ggplot(trapDFmod, aes(x = TrapType, y = Effort))+
  geom_col() +
  labs(y = "Effort (total trap-nights)", x = "Traps" ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  ggtitle("")

gtot = gsites + gweeks + geffort 

ggsave(plot = gtot, filename = paste0(folderOutput, "/A - Trapping effort.png"),
       device = "png", width = 7, height = 5)

# whole DF

histTotDF <- totDFmod %>%
  mutate(genusSpecies = paste(Genus, Species)) %>%
  group_by(genusSpecies) %>%
  summarise(totAbundance = sum(AvgAbundance)) %>%
  ungroup() %>%
  mutate(perc = totAbundance/sum(totAbundance)) 

#df correspondence trap-genusSpecies
correspondenceTrapDF = data.frame(genusSpecies = rep(genusSpecies, times = length(traps)),
                                  trapType = rep(traps, each = length(genusSpecies)),
                                  match = "d) absent")

## loop----
# histogram of genusSpecies per trap 

for(trapi in traps){
  
  cat(trapi, "\n")
  
  histDF <- totDFmod %>%
    filter(TrapType == trapi) %>%
    group_by(genusSpecies) %>%
    summarise(totAbundance = sum(AvgAbundance)) %>%
    ungroup() %>%
    mutate(perc = totAbundance/sum(totAbundance)) %>%
    mutate(percLab = paste0(round(perc, 3), "%"))
  
  genusSpeciesi = unique(histDF$genusSpecies)
  
  ggplot(histDF, aes(x = totAbundance , y = genusSpecies, label = percLab))+
    xlim(c(0, 1.05*max(histDF$totAbundance)))+
    geom_col(stat = "identity")+ 
    geom_text(hjust = -0.1,    # nudge above top of bar
              size = 3)+
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = "gray90"))
  
  ggsave(filename = paste0(folderOutput, "/A - Histogram of species for ", trapi,".png"),
         device = "png", width = 7, height = 5)
  
  for(genusSpeciesii in genusSpeciesi){
    
    percTrapii = histDF %>% filter(genusSpecies == genusSpeciesii) %>% pull(perc)
    perc = histTotDF %>% filter(genusSpecies == genusSpeciesii) %>% pull(perc)
    
    if (percTrapii > 10*perc){
      match = "a) more present than average"
    } else if(percTrapii > 0.1*perc) {
      match = "b) present as average"
    } else if(percTrapii > 0.) {
      match = "c) less present than average"
    } else {
      match = "d) absent"
    }
    
    correspondenceTrapDF$match[which(correspondenceTrapDF$trapType == trapi &
                                       correspondenceTrapDF$genusSpecies == genusSpeciesii)] =  match
  }

}

ggplot(correspondenceTrapDF, aes(x = trapType, y = genusSpecies, fill = match))+
  geom_tile() +
  scale_fill_discrete()+
  ggtitle("Order of magnitude of detection")
  

ggsave(filename = paste0(folderOutput, "/A - Correspondence table.png"),
       device = "png", width = 7, height = 5)

# Only selected traps and periods----

totDFmod <- readRDS(file = paste0(folderDataLocal, "/totDFsel_ElDorado_Sepulveda.rds"))

totDFmod <- totDFmod %>%
  filter(!is.na(Species)) %>%
  mutate(genusSpecies = paste(Genus, Species))

# # trap type
# traps = unique(totDFmod$TrapType)
# genusSpecies = unique(totDFmod$genusSpecies)

#surveillanceEffort per trap
trapDFmod <- totDFmod %>%
  group_by(TrapType) %>%
  summarize(Sites = length(unique(siteCode)),
            Weeks = length(unique(collectionWeek)),
            Effort = sum(totNightTraps))%>%
  ungroup() 

gsites = ggplot(trapDFmod, aes(x = TrapType, y = Sites))+
  geom_col() +
  labs(y = "Spatial records (sites)", x = "Traps" ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  ggtitle("")

gweeks = ggplot(trapDFmod, aes(x = TrapType, y = Weeks))+
  geom_col() +
  labs(y = "Temporal records (weeks)", x = "Traps" ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  ggtitle("Measures of trapping effort")

geffort = ggplot(trapDFmod, aes(x = TrapType, y = Effort))+
  geom_col() +
  labs(y = "Effort (total trap-nights)", x = "Traps" ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  ggtitle("")

gtot = gsites + gweeks + geffort 

ggsave(plot = gtot, filename = paste0(folderOutput, "/E - Trapping effort after selection.png"),
       device = "png", width = 7, height = 5)

# whole DF

histTotDF <- totDFmod %>%
  mutate(genusSpecies = paste(Genus, Species)) %>%
  group_by(genusSpecies) %>%
  summarise(totAbundance = sum(AvgAbundance)) %>%
  ungroup() %>%
  mutate(perc = totAbundance/sum(totAbundance)) 

#df correspondence trap-genusSpecies
correspondenceTrapDF = data.frame(genusSpecies = rep(genusSpecies, times = length(traps)),
                                  trapType = rep(traps, each = length(genusSpecies)),
                                  match = "d) absent")

## loop----
# histogram of genusSpecies per trap 

for(trapi in traps){
  
  cat(trapi, "\n")
  
  histDF <- totDFmod %>%
    filter(TrapType == trapi) %>%
    group_by(genusSpecies) %>%
    summarise(totAbundance = sum(AvgAbundance)) %>%
    ungroup() %>%
    mutate(perc = totAbundance/sum(totAbundance)) %>%
    mutate(percLab = paste0(round(perc, 3), "%"))
  
  genusSpeciesi = unique(histDF$genusSpecies)
  
  ggplot(histDF, aes(x = totAbundance , y = genusSpecies, label = percLab))+
    xlim(c(0, 1.05*max(histDF$totAbundance)))+
    geom_col(stat = "identity")+ 
    geom_text(hjust = -0.1,    # nudge above top of bar
              size = 3)+
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = "gray90"))
  
  ggsave(filename = paste0(folderOutput, "/E - Histogram of species for ", trapi," after selection.png"),
         device = "png", width = 7, height = 5)
  
  for(genusSpeciesii in genusSpeciesi){
    
    percTrapii = histDF %>% filter(genusSpecies == genusSpeciesii) %>% pull(perc)
    perc = histTotDF %>% filter(genusSpecies == genusSpeciesii) %>% pull(perc)
    
    if (percTrapii > 10*perc){
      match = "a) more present than average"
    } else if(percTrapii > 0.1*perc) {
      match = "b) present as average"
    } else if(percTrapii > 0.) {
      match = "c) less present than average"
    } else {
      match = "d) absent"
    }
    
    correspondenceTrapDF$match[which(correspondenceTrapDF$trapType == trapi &
                                       correspondenceTrapDF$genusSpecies == genusSpeciesii)] =  match
  }
  
}

ggplot(correspondenceTrapDF, aes(x = trapType, y = genusSpecies, fill = match))+
  geom_tile() +
  scale_fill_discrete()+
  ggtitle("Order of magnitude of detection")


ggsave(filename = paste0(folderOutput, "/E - Correspondence table after selection.png"),
       device = "png", width = 7, height = 5)
