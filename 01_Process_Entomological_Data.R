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

# These are test for the disease: therefore rows are repeated for each sample. Le

totDFmod <- totDF %>%
  group_by(pool_id, trap_type, site_code, collection_date, species, sex, sex_condition) %>%
  summarise(num_count = mean(num_count)) %>% 
  ungroup()

# if they belong to the same pool id, it means that they are the same mosquito that are tested over only multiple diseases.
# however, if they belong to different pool id, it means that they are different mosquitos - there fore they should be summed.

totDFmod <- totDF %>%
  group_by(trap_type, site_code, collection_date, species, sex, sex_condition) %>%
  summarise(num_count = sum(num_count)) %>% 
  ungroup()
# This need to be checked

summary(totDFmod$num_count) # There are never 0: is this normal?

# check #1: there should be only one trap/date/species/sex. Is it so?
totDFmod$colchk = paste0(totDFmod$trap_type, "_", totDFmod$site_code, "_", totDFmod$collection_date, "_", totDFmod$species, "_", totDFmod$sex_condition)
-sort(-table(totDFmod$colchk))[1]

x = totDFmod %>% filter(site_code == "2278",
                        collection_date == "2008-07-17",
                        species == "Culex quinquefasciatus",
                        sex_condition == "Females - Mixed")

x0 = totDF %>% filter(site_code == "2278",
                      collection_date == "2008-07-17",
                      species == "Culex quinquefasciatus",
                      sex_condition == "Females - Mixed")

# ok ther are no more duplicated lines
saveRDS(totDFmod, file = paste0(folderDataLocal, "/totDFmod_ElDorado_Sepulveda.rds"))
totDFmod <- readRDS(file = paste0(folderDataLocal, "/totDFmod_ElDorado_Sepulveda.rds"))

# Variables definition ----

species = unique(totDF$species)
sites = unique(totDF$site_code)
traps = unique(totDF$trap_type)

# let's plot the actovity period for each trap (in weeks)
year_start = year(min(totDF$collection_date))
year_end = year(max(totDF$collection_date))

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
for(i in 1:nrow(totDF)){
  mi = week(totDF$collection_date[i])
  yi = year(totDF$collection_date[i])
  si = totDF$site_code[i]
  ai = totDF$area[i]
  quinquefasciatusi = (totDF$species[i] == "Culex quinquefasciatus")
  tarsalisi = (totDF$species[i] == "Culex tarsalis")
  stigmatosomai = (totDF$species[i] == "Culex stigmatosoma")
  aegyptii = (totDF$species[i] == "Aedes aegypti")
  
  r = which(sitesDF$year == yi & sitesDF$week == mi & sitesDF$site == si)
  
  sitesDF$active[r] = 1
  sitesDF$area[r] = ai
  sitesDF$quinquefasciatus[r] = sitesDF$quinquefasciatus[r] + quinquefasciatusi
  sitesDF$tarsalis[r] = sitesDF$tarsalis[r] + tarsalisi
  sitesDF$stigmatosoma[r] =  sitesDF$stigmatosoma[r] + stigmatosomai
  sitesDF$aegypti[r]= sitesDF$aegypti[r] + aegyptii
}
toc() #26-31 sec

# Image preprocessing

sitesDF <- sitesDF %>%
  filter(!is.na(area))

sitesDF$progYear = round(sitesDF$week/52 + sitesDF$year,2)

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
  scale_fill_viridis_c(option = "D", direction = -1)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  facet_wrap(.~area, scales = "free_y", space = "free_y")+
  ggtitle("Detection of Ae. aegypti")

## C. quinquefasciatus----
 
ggplot(data = sitesDF, aes(x = progYear, y = site, fill = quinquefasciatus))+
  geom_tile()+
  scale_fill_viridis_c(option = "D", direction = -1)+
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
         filter(site_code == 2655,
                species == "Culex quinquefasciatus"),
       aes(x = collection_date,
           y = num_count)) +
  geom_point()

# what does it mean "50"? 

hist(totDFmod$num_count)

save(sitesDF, file = paste0(folderDataLocal, "/sitesDF_ElDorado_Sepulveda.rds"))

# can I compute diversity indices by week (how can I do that? Max abundance is limited to 50...)

# DOUBTS
# why are abundance data limited to 50? At least they look so
# is "collection_date" the starting date of sampling and "add_date" the ending date? 
# is it possible to have daily weather series (min temp, max temp, avg temp, cum rain, relativ humidity) for the 4 sites?
# What is NA trap? 
# why not zeros? Are these line sdeleted or just by channe? 

#With Mathilde
# is week a good scale? Waht abouth month?
# should indicators be diffeerent depending on traps?
# how to compute biodiversity indices depending on different traps?
