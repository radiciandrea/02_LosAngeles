# Move to "format long"

library(pracma)
library(tidyverse)
library(lubridate)
library(ISOweek)

folderDataLocal = "Data"

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

siteWeeksTempDF = data.frame(site = as.factor(rep(sites, times = length(weeks))),
                     week = rep(weeks, each = length(sites)),
                     active = 0,
                     area = NA)

speciesTempDF = data.frame(matrix(data = 0, nrow = nrow(siteWeeksDF), ncol = length(species)))
names(speciesTempDF) = species

siteWeeksDF = cbind(siteWeeksTempDF, speciesTempDF)

nc = ncol(siteWeeksTempDF)

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
toc() #90 sec

# histogram of species
histDF <- totDFmod %>%
  group_by(Species) %>%
  summarise(tot = sum(AvgAbundance)) %>%
  ungroup()

ggplot(histDF, aes(x = tot , y = Species))+
  geom_col(stat = "identity")

# well... quinquefasciatus is the winner

# Image preprocessing

siteWeeksDF <- siteWeeksDF %>%
  filter(!is.na(area))

# Plot ----

## Surveillance ----

ggplot(data = siteWeeksDF, aes(x = week, y = site, fill = active))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Active surveillance")

## Ae. aegypti----

ggplot(data = siteWeeksDF, aes(x = week, y = site, fill = aegypti))+
  geom_tile()+
  scale_fill_viridis_c(option = "A", direction = -1)+
  geom_tile(data = siteWeeksDF %>% filter(aegypti == 0), aes(x = week, y = site), fill = "gray90")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Detection of Ae. aegypti")

## C. quinquefasciatus----

ggplot(data = siteWeeksDF, aes(x = week, y = site, fill = quinquefasciatus))+
  geom_tile()+
  scale_fill_viridis_c(option = "A", direction = -1)+
  geom_tile(data = siteWeeksDF %>% filter(quinquefasciatus == 0), aes(x = week, y = site), fill = "gray90")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Detection of C. quinquefasciatus")

## Plot of 1 trajectory ----

# let's consider 1 trap (site) and 1 species for plotting

# species = Culex quinquefasciatus
# site = 2800

ggplot(data = totDFmod %>%
         filter(siteCode == 2655,
                Species == "quinquefasciatus"),
       aes(x = collectionWeek,
           y = AvgAbundance)) +
  geom_point()

hist(totDFmod$AvgAbundance)

saveRDS(siteWeeksDF, file = paste0(folderDataLocal, "/siteWeeksDF_ElDorado_Sepulveda.rds"))
