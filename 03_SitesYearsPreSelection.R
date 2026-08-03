# Select only recent data.

library(pracma)
library(tidyverse)
library(lubridate)

folderDataLocal = "Data"
folderOutput = "Outputs"

# load data
totDFmod <- readRDS(file = paste0(folderDataLocal, "/totDFmod_ElDorado_Sepulveda.rds"))
siteWeeksDF <- readRDS(file = paste0(folderDataLocal, "/siteWeeksDF_ElDorado_Sepulveda.rds"))

## Traps and year selection --

DFsel1 <- siteWeeksDF %>%
  mutate(year = as.numeric(substr(week, 1, 4)))%>%
  group_by(year)%>%
  summarise(nsites = n()/52)%>%
  ungroup()

ggplot(DFsel1, aes(x=year, y=nsites))+
  geom_bar(stat = "identity")

# let's consider 203 to 2022
yearSel = 2003:2022

# also: activity period

DFsel2 <- siteWeeksDF %>%
  group_by(site)%>%
  summarise(nweeks = n())%>%
  ungroup()

ggplot(DFsel2, aes(x=site, y=nweeks, label = nweeks))+
  geom_bar(stat = "identity")+
  geom_text(vjust = -0.2,    # nudge above top of bar
            size = 3)

# let's consider at least 20 

siteSel = DFsel2 %>% filter(nweeks > 20) %>% pull(site)

siteWeeksDFsel <- siteWeeksDF %>%
  filter(site %in% siteSel) %>%
  filter(datesLabels >= as.Date(paste0(min(yearSel), "-01-01"))) %>%
  filter(datesLabels <= as.Date(paste0(max(yearSel), "-12-31")))

totDFsel <- totDFmod %>%
  mutate(year = as.numeric(substr(collectionWeek, 1, 4)))%>%
  filter(year %in% yearSel) %>%
  filter(siteCode %in% siteSel)

# histogram of species
histDF <- totDFsel %>%
  filter(!is.na(Species)) %>%
  mutate(GenusSpecies = paste(Genus, Species)) %>%
  group_by(GenusSpecies, Species) %>%
  summarise(totAbundance = sum(AvgAbundance)) %>%
  ungroup() %>%
  mutate(perc = paste0(round(100*totAbundance/sum(totAbundance), 3), "%"))

ggplot(histDF, aes(x = totAbundance , y = GenusSpecies, label = perc))+
  xlim(c(0, 1.05*max(histDF$totAbundance)))+
  geom_col(stat = "identity")+ 
  geom_text(hjust = -0.1,    # nudge above top of bar
            size = 3)+
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/D - Histogram after selection.png"), device = "png", width = 7, height = 5)

# plot to save for all species

for(si in species){
  
  gsName = histDF %>% filter(Species == si) %>% pull(GenusSpecies)
  gsPerc = histDF %>% filter(Species == si) %>% pull(perc)
  
  ggplot(data = siteWeeksDFsel, aes(x = datesLabels, y = site, fill = .data[[si]]))+
    geom_tile()+
    scale_fill_viridis_c(option = "H", name="N/trap/night")+
    geom_tile(data = siteWeeksDFsel %>% filter(.data[[si]] == 0), aes(x = datesLabels, y = site), fill = "gray90")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 5.5),
          panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = "gray90"))+
    facet_wrap(.~area, scales = "free_y", space = "free_y")+
    ggtitle(paste0("Detection of ", gsName, " - ", gsPerc, " of the total sampled species"))
  
  ggsave(filename = paste0(folderOutput, "/E - ", gsName, " after selection.png"), device = "png", width = 14, height = 7)
  
}

# save

saveRDS(totDFsel, file = paste0(folderDataLocal, "/totDFsel_ElDorado_Sepulveda.rds"))
saveRDS(siteWeeksDFsel, file = paste0(folderDataLocal, "/siteWeeksDFsel_ElDorado_Sepulveda.rds"))
