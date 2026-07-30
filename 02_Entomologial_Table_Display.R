# Surveillance: compute surveillance degradation (delay)

library(pracma)
library(tidyverse)
library(lubridate)
library(ISOweek)

folderDataLocal = "Data"

# load data
totDFmod <- readRDS(file = paste0(folderDataLocal, "/totDFmod_ElDorado_Sepulveda.rds"))

# Variables definition ----

species = unique(totDFmod$species)
sites = unique(totDFmod$site_code)
traps = unique(totDFmod$trap_type)

# let's plot the actovity period for each trap (in weeks)
year_start = year(min(totDFmod$collection_date))
year_end = year(max(totDFmod$collection_date))

years = rep(year_start:year_end, each = 52)
weeks = rep(1:52, times = length(year_start:year_end))

# Dataframe of presences ----

sitesDF = data.frame(site = as.factor(rep(sites, times = length(weeks))),
                     year = rep(years, each = length(sites)),
                     week = rep(weeks, each = length(sites)),
                     active = 0,
                     area = NA,
                     quinquefasciatus = 0,
                     tarsalis = 0,
                     stigmatosoma = 0,
                     aegypti = 0)

## Filling ----
# not the fastyest way, but: 

tic()
for(i in 1:nrow(totDFmod)){
  mi = week(totDFmod$collection_date[i])
  yi = year(totDFmod$collection_date[i])
  si = totDFmod$site_code[i]
  ai = totDFmod$area[i]
  quinquefasciatusi = (totDFmod$species[i] == "Culex quinquefasciatus")*(totDFmod$num_count[i])
  tarsalisi = (totDFmod$species[i] == "Culex tarsalis")*(totDFmod$num_count[i])
  stigmatosomai = (totDFmod$species[i] == "Culex stigmatosoma")*(totDFmod$num_count[i])
  aegyptii = (totDFmod$species[i] == "Aedes aegypti")*(totDFmod$num_count[i])
  
  r = which(sitesDF$year == yi & sitesDF$week == mi & sitesDF$site == si)
  
  sitesDF$active[r] = 1
  sitesDF$area[r] = ai
  
  # species
  sitesDF$quinquefasciatus[r] = sitesDF$quinquefasciatus[r] + quinquefasciatusi
  sitesDF$tarsalis[r] = sitesDF$tarsalis[r] + tarsalisi
  sitesDF$stigmatosoma[r] =  sitesDF$stigmatosoma[r] + stigmatosomai
  sitesDF$aegypti[r]= sitesDF$aegypti[r] + aegyptii
}
toc() #9 sec

sitesDF$progYear = :ISOweek2date(
  sprintf("%d-W%02d-1", sitesDF$year, sitesDF$week)
)

# histogram of species
histDF <- totDFmod %>%
  group_by(species) %>%
  summarise(tot = sum(num_count)) %>%
  ungroup()

ggplot(histDF, aes(x = species, y = tot))+
  geom_col(stat = "identity")

# well... quinquefasciatus is the winner

# Image preprocessing

sitesDF <- sitesDF %>%
  filter(!is.na(area))

# Plot ----

## Surveillance ----

ggplot(data = sitesDF, aes(x = progYear, y = site, fill = active))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Active surveillance")

## Ae. aegypti----

ggplot(data = sitesDF, aes(x = progYear, y = site, fill = aegypti))+
  geom_tile()+
  scale_fill_viridis_c(option = "A", direction = -1)+
  geom_tile(data = sitesDF %>% filter(aegypti == 0), aes(x = progYear, y = site), fill = "gray90")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Detection of Ae. aegypti")

## C. quinquefasciatus----

ggplot(data = sitesDF, aes(x = progYear, y = site, fill = quinquefasciatus))+
  geom_tile()+
  scale_fill_viridis_c(option = "A", direction = -1)+
  geom_tile(data = sitesDF %>% filter(quinquefasciatus == 0), aes(x = progYear, y = site), fill = "gray90")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Detection of C. quinquefasciatus")

## C. stigmatosoma----

ggplot(data = sitesDF, aes(x = progYear, y = site, fill = stigmatosoma))+
  geom_tile()+
  scale_fill_viridis_c(option = "A", direction = -1)+
  geom_tile(data = sitesDF %>% filter(stigmatosoma == 0), aes(x = progYear, y = site), fill = "gray90")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Detection of C. stigmatosoma")

## C. tarsalis----

ggplot(data = sitesDF, aes(x = progYear, y = site, fill = tarsalis))+
  geom_tile()+
  scale_fill_viridis_c(option = "A", direction = -1)+
  geom_tile(data = sitesDF %>% filter(tarsalis == 0), aes(x = progYear, y = site), fill = "gray90")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Detection of C. tarsalis")

## Plot of 1 trajectory ----

# let's consider 1 trap (site) and 1 species for plotting

# species = Culex quinquefasciatus
# site = 2800

ggplot(data = totDFmod %>%
         filter(site_code == 2655,
                species == "Culex quinquefasciatus"),
       aes(x = collection_date,
           y = num_count)) +
  geom_point()

# what does it mean "50"? 

hist(totDFmod$num_count)

saveRDS(sitesDF, file = paste0(folderDataLocal, "/sitesDF_ElDorado_Sepulveda.rds"))
