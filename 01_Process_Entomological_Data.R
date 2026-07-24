# preprocess entomological observation files

library(xlsx) #4.5
library(pracma)
library(tidyverse)

# datafodler

folderData = "C:/Users/2024ar003/Desktop/Alcuni file permanenti/Post_doc_biodivecity/Dati/US_Los_Angeles"
folderDataLocal = "Data"

dir.create(folderDataLocal)
# Preprocessing ----

# #  save data as RDS (done once for all)
# tic()
# ElDoradoDF = read.xlsx(file = paste0(folderData, "/El Dorado Park abundance and arbovirus data.xlsx"),
#                        sheetName = "Sheet2")
# SepulvedaDF = read.xlsx(file = paste0(folderData, "/Sepulveda Basin abundance and arbovirus data.xlsx"),
#                         sheetName = "Sheet2")
# toc()
# 
# # read.xlsx is super slow: 
# 
# # the two tables have the same items (rbind is ok)
# 
# ElDoradoDF$area = "El Dorado"
# ElDoradoDF$type = "green area"
# SepulvedaDF$area = "Sepulveda"
# SepulvedaDF$type = "green area"
# 
# #rbind
# 
# totDF = rbind(ElDoradoDF, SepulvedaDF)
# 
# saveRDS(totDF, file = paste0(folderDataLocal, "/Abundance_ElDorado_Sepulveda.rds"))

totDF <- readRDS(file = paste0(folderDataLocal, "/Abundance_ElDorado_Sepulveda.rds"))

# Variables definition ----

species = unique(totDF$species)
sites = unique(totDF$site_code)
traps = unique(totDF$trap_type)

# let's plot the actovity period for each trap (in months)
year_start = year(min(totDF$collection_date))
year_end = year(max(totDF$collection_date))

years = rep(year_start:year_end, each = 12)
months = rep(1:12, times = length(year_start:year_end))

# Dataframe of presences ----

sitesDF = data.frame(site = as.factor(rep(sites, times = length(months))),
                     year = rep(years, each = length(sites)),
                     month = rep(months, each = length(sites)),
                     active = 0,
                     area = NA,
                     quinquefasciatus = 0,
                     tarsalis = 0,
                     stigmatosoma = 0,
                     aegypti = 0)

## Filling ----
# not the fastyest way, but: 

for(i in 1:nrow(totDF)){
  mi = month(totDF$collection_date[i])
  yi = year(totDF$collection_date[i])
  si = totDF$site_code[i]
  ai = totDF$area[i]
  quinquefasciatusi = (totDF$species[i] == "Culex quinquefasciatus")
  tarsalisi = (totDF$species[i] == "Culex tarsalis")
  stigmatosomai = (totDF$species[i] == "Culex stigmatosoma")
  aegyptii = (totDF$species[i] == "Aedes aegypti")
  
  r = which(sitesDF$year == yi & sitesDF$month == mi & sitesDF$site == si)
  
  sitesDF$active[r] = 1
  sitesDF$area[r] = ai
  sitesDF$quinquefasciatus[r] = pmax(sitesDF$quinquefasciatus[r], quinquefasciatusi)
  sitesDF$tarsalis[r] = pmax(sitesDF$tarsalis[r],tarsalisi)
  sitesDF$stigmatosoma[r] = pmax(sitesDF$stigmatosoma[r], stigmatosomai)
  sitesDF$aegypti[r]= pmax(sitesDF$aegypti[r], aegyptii)
}

# Image preprocessing

sitesDF <- sitesDF %>%
  filter(!is.na(area))

sitesDF$progYear = round(sitesDF$month/12 + sitesDF$year,2)

# Plot ----

## Surveillance ----

ggplot(data = sitesDF, aes(x = progYear, y = site, fill = active))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")+
  facet_wrap(.~area, scales = "free_y")+
  ggtitle("Active surveillance")

## Ae. aegypti----

ggplot(data = sitesDF, aes(x = progYear, y = site, fill = aegypti))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")+
  facet_wrap(.~area, scales = "free_y")+
  ggtitle("Detection of Ae. aegypti")

## C. quinquefasciatus----
 
ggplot(data = sitesDF, aes(x = progYear, y = site, fill = quinquefasciatus))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")+
  facet_wrap(.~area, scales = "free_y")+
  ggtitle("Detection of C. quinquefasciatus")

## Plot of 1 trajectory ----

# let's consider 1 trap (site) and 1 species for plotting

# species = Culex quinquefasciatus
# site = 2800

ggplot(data = totDF %>%
         filter(site_code == 2655,
                species == "Culex quinquefasciatus"),
       aes(x = collection_date,
           y = num_count)) +
  geom_point()

# what does it mean "50"? 

hist(totDF$num_count)



# can I compute diversity indices by month (how can I do that? Max abundance is limited to 50...)

# Detection delay

correct_date_detection = (sitesDF %>%
  filter(aegypti > 0) %>%
  filter(progYear == min(progYear)) %>%
  pull(progYear))[1]

# delay with 3 site less
date_detection = c()

tic()
for(i in 1:length(sites)){
  s1 = sites[i]
  sites1 = sites[-i]
  for(j in 1:length(sites1)){
    s2 = sites1[j]
    sites2 = sites1[-j]
    for(k in 1:length(sites2)){
      s3 = sites2[k]
      sites3 = sites1[-k]
      
      tempsitesDF = sitesDF %>%
        filter(site %in% sites3)
      
      date_detection = c(date_detection, (tempsitesDF %>%
                                            filter(aegypti > 0) %>%
                                            filter(progYear == min(progYear)) %>%
                                            pull(progYear))[1])
    }
  }
}
toc() # quite long

m_delays = 12*(date_detection-correct_date_detection)

summary(m_delays)

# DOUBTS
# why are abundance data limited to 50? At least they look so
# is "collection_date" the starting date of sampling and "add_date" the ending date? 
# is it possible to have daily weather series (min temp, max temp, avg temp, cum rain, relativ humidity) for the 4 sites?


