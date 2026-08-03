# Move to "format long"

library(pracma)
library(tidyverse)
library(lubridate)
library(ISOweek)

folderDataLocal = "Data"
folderOutput = "Outputs"

# load data
totDF <- readRDS(file = paste0(folderDataLocal, "/totDF_ElDorado_Sepulveda.rds"))
totDFmod <- readRDS(file = paste0(folderDataLocal, "/totDFmod_ElDorado_Sepulveda.rds"))

# Variables definition ----

species = unique(totDFmod$Species)
species = species[-which(is.na(species))]
sites = unique(totDFmod$siteCode)
#traps = unique(totDFmod$TrapType)

# let's plot the actovity period for each trap (in weeksRep)
yearStart = substr(min(totDFmod$collectionWeek), 1, 4)
yearEnd = substr(max(totDFmod$collectionWeek), 1, 4)

seriesLength = as.numeric(as.Date(paste0(yearEnd, "-12-31")) - as.Date(paste0(yearStart, "-01-01")))
progressiveFirstDayOfTheWeek = seq(1, seriesLength, 7)

dates = as.Date(progressiveFirstDayOfTheWeek, format = "%Y-%m-%d", origin = paste0(yearStart, "-01-01"))
weeks = ISOweek(dates)

# Dataframe of abundances ----

siteWeeksTempng = data.frame(site = as.factor(rep(sites, times = length(weeks))),
                             datesLabels = rep(dates, each = length(sites)),
                             week = rep(weeks, each = length(sites)),
                             active = 0,
                             area = NA)

speciesTempng = data.frame(matrix(data = 0, nrow = nrow(siteWeeksTempng), ncol = length(species)))
names(speciesTempng) = species

siteWeeksDF = cbind(siteWeeksTempng, speciesTempng)

nc = ncol(siteWeeksTempng)

## Filling ----
# not the fastyest way, but: 

tic()
for(i in 1:nrow(totDFmod)){
  
  wi = totDFmod$collectionWeek[i]
  sci = totDFmod$siteCode[i]
  ai = totDFmod$Area[i]
  si = totDFmod$Species[i]
  
  #which point of the DataFrame
  r = which(siteWeeksDF$week == wi & siteWeeksDF$site == sci)
  
  siteWeeksDF$active[r] = 1
  siteWeeksDF$area[r] = ai
  
  # species
  if(!is.na(si)){
    siteWeeksDF[r, nc + which(species == si)] = totDFmod$AvgAbundance[i]
  }

}
toc() #80-90 sec

# histogram of species
histDF <- totDFmod %>%
  filter(!is.na(Species)) %>%
  mutate(GenusSpecies = paste(Genus, Species)) %>%
  group_by(GenusSpecies, Species) %>%
  summarise(totAbundance = sum(AvgAbundance)) %>%
  ungroup() %>%
  mutate(perc = paste0(round(100*totAbundance/sum(totAbundance), 1), "%"))

ggplot(histDF, aes(x = totAbundance , y = GenusSpecies, label = perc))+
  xlim(c(0, 1.05*max(histDF$totAbundance)))+
  geom_col(stat = "identity")+ 
  geom_text(hjust = -0.1,    # nudge above top of bar
            size = 3)+
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/A - Histogram.png"), device = "png", width = 7, height = 5)

# well... quinquefasciatus is the winner

# Image preprocessing

siteWeeksDF <- siteWeeksDF %>%
  filter(!is.na(area))

saveRDS(siteWeeksDF, file = paste0(folderDataLocal, "/siteWeeksDF_ElDorado_Sepulveda.rds"))

siteWeeksDF <- readRDS(file = paste0(folderDataLocal, "/siteWeeksDF_ElDorado_Sepulveda.rds"))

# Plot ----

## Surveillance ----

ggplot(data = siteWeeksDF, aes(x = datesLabels, y = site, fill = active))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5.5),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Active surveillance")

ggsave(filename = paste0(folderOutput, "/B - ActiveSurveillance.png"), device = "png", width = 14, height = 7)

## Ae. aegypti----

ggplot(data = siteWeeksDF, aes(x = datesLabels, y = site, fill = aegypti))+
  geom_tile()+
  scale_fill_viridis_c(option = "H")+
  geom_tile(data = siteWeeksDF %>% filter(aegypti == 0), aes(x = datesLabels, y = site), fill = "gray90")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5.5),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Detection of Ae. aegypti")

## C. quinquefasciatus----

ggplot(data = siteWeeksDF, aes(x = datesLabels, y = site, fill = quinquefasciatus))+
  geom_tile()+
  scale_fill_viridis_c(option = "H")+
  geom_tile(data = siteWeeksDF %>% filter(quinquefasciatus == 0), aes(x = datesLabels, y = site), fill = "gray90")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5.5),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Detection of C. quinquefasciatus")

# plot to save for all species

for(si in species){
  
  gsName = histDF %>% filter(Species == si) %>% pull(GenusSpecies)
  gsPerc = histDF %>% filter(Species == si) %>% pull(perc)
  
  ggplot(data = siteWeeksDF, aes(x = datesLabels, y = site, fill = .data[[si]]))+
    geom_tile()+
    scale_fill_viridis_c(option = "H", name="N/trap/night")+
    geom_tile(data = siteWeeksDF %>% filter(.data[[si]] == 0), aes(x = datesLabels, y = site), fill = "gray90")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 5.5),
          panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = "gray90"))+
    facet_wrap(.~area, scales = "free_y", space = "free_y")+
    ggtitle(paste0("Detection of ", gsName, " - ", gsPerc, " of the total sampled species"))
  
    ggsave(filename = paste0(folderOutput, "/C - ", gsName, ".png"), device = "png", width = 14, height = 7)
  
}

## Plot of 1 trajectory ----

# let's consider 1 trap (site) and 1 species for plotting

# species = Culex quinquefasciatus
# site = 2800

siteMaxCq = siteWeeksDF %>%
  select(c("site", "quinquefasciatus")) %>%
  group_by(site)%>%
  summarise(sumAvg = sum(quinquefasciatus)) %>%
  ungroup() %>%
  filter(sumAvg == max(sumAvg)) %>%
  pull(site)

ggplot(data = siteWeeksDF %>%
         filter(site == siteMaxCq),
       aes(x = datesLabels,
           y = quinquefasciatus)) +
  geom_point()

hist(totDFmod$AvgAbundance)

# Surveillance season ----

collectionDF <- totDFmod %>%
  mutate(year = as.numeric(substr(collectionWeek, 1,4))) %>%
  select(c("collectionWeek", "year")) %>%
  unique() %>%
  group_by(year) %>%
  summarize(nWeeks = n()) %>%
  ungroup() %>%
  mutate(Area = "all", type = "all year around")
  
collectionAreaDF <- totDFmod %>%
  mutate(year = as.numeric(substr(collectionWeek, 1,4))) %>%
  select(c("collectionWeek", "year", "Area")) %>%
  unique() %>%
  group_by(year, Area) %>%
  summarize(nWeeks = n()) %>%
  ungroup() %>%
  mutate(type = "all year around")

# eary surveillance = before 15st of may, before week 11

earlyCollectionDF <- totDFmod %>%
  mutate(year = as.numeric(substr(collectionWeek, 1, 4))) %>%
  mutate(collectionWeek = as.numeric(substr(collectionWeek, 7,8))) %>%
  filter(collectionWeek < 11 ) %>%
  select(c("collectionWeek", "year")) %>%
  unique() %>%
  group_by(year) %>%
  summarize(nWeeks = n()) %>%
  ungroup() %>%
  mutate(Area = "all", type = "early (may)")

earlyCollectionAreaDF <- totDFmod %>%
  mutate(year = as.numeric(substr(collectionWeek, 1,4))) %>%
  mutate(collectionWeek = as.numeric(substr(collectionWeek, 7,8))) %>%
  filter(collectionWeek < 11 ) %>%
  select(c("collectionWeek", "year", "Area")) %>%
  unique() %>%
  group_by(year, Area) %>%
  summarize(nWeeks = n()) %>%
  ungroup() %>%
  mutate(type = "early (may)")

# late surveillance = after 15st of october, after week 41

lateCollectionDF <- totDFmod %>%
  mutate(year = as.numeric(substr(collectionWeek, 1, 4))) %>%
  mutate(collectionWeek = as.numeric(substr(collectionWeek, 7,8))) %>%
  filter(collectionWeek > 41 ) %>%
  select(c("collectionWeek", "year")) %>%
  unique() %>%
  group_by(year) %>%
  summarize(nWeeks = n()) %>%
  ungroup() %>%
  mutate(Area = "all", type = "late (october)")

lateCollectionAreaDF <- totDFmod %>%
  mutate(year = as.numeric(substr(collectionWeek, 1,4))) %>%
  mutate(collectionWeek = as.numeric(substr(collectionWeek, 7,8))) %>%
  filter(collectionWeek > 41 ) %>%
  select(c("collectionWeek", "year", "Area")) %>%
  unique() %>%
  group_by(year, Area) %>%
  summarize(nWeeks = n()) %>%
  ungroup() %>%
  mutate(type = "late (october)")

totCollectionDF = rbind(collectionDF, collectionAreaDF,
                     earlyCollectionDF, earlyCollectionAreaDF,
                     lateCollectionDF, lateCollectionAreaDF)

ggplot(data = totCollectionDF, aes(x = year, y = nWeeks, color = Area))+
  geom_point()+
  geom_line()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5.5),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~type, scales = "free_y", nrow = 3)+
  ggtitle("Surveilled weeks per year")

ggsave(filename = paste0(folderOutput, "/B - ActiveWeekSurveillance.png"), device = "png", width = 14, height = 7)
