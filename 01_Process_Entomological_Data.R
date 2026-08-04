# preprocess entomological observation files

library(xlsx) #4.5
library(pracma)
library(lubridate)
library(tidyverse)
library(ISOweek)

# datafodler

folderData = "C:/Users/2024ar003/Desktop/Alcuni file permanenti/Post_doc_biodivecity/Dati/US_Los_Angeles"
folderDataLocal = "Data"

dir.create(folderDataLocal)
# Preprocessing ----

#  save data as RDS (once for all)
tic()
ElDoradoDF = read.xlsx(file = paste0(folderData, "/El Dorado Park data.xlsx"),
                       sheetName = "ElDoradoTrapData")
toc()
SepulvedaDF = read.xlsx(file = paste0(folderData, "/Sepulveda Basin data.xlsx"),
                        sheetName = "SepulvedaTrapData")
toc()

# read.xlsx is super slow:

# the two tables have the same items (rbind is ok)

ElDoradoDF$Area = "El Dorado"
ElDoradoDF$Landscape = "green area"
SepulvedaDF$Area = "Sepulveda"
SepulvedaDF$Landscape = "green area"

#rbind

totDF = rbind(ElDoradoDF, SepulvedaDF)

saveRDS(totDF, file = paste0(folderDataLocal, "/totDF_ElDorado_Sepulveda.rds"))
totDF = readRDS(file = paste0(folderDataLocal, "/totDF_ElDorado_Sepulveda.rds"))

# Genus adding; estimating daily total presence per trap per night per site

totDFmod <- totDF %>%
  rename(siteCode = Site.Code) %>%
  mutate(tot = pmax(Males, 0, na.rm = T) + pmax(Females, 0, na.rm = T)) %>%
  mutate(collectionWeek = ISOweek(CollectionDate)) %>%
  group_by(siteCode, collectionWeek, Species, Area, TrapType, Landscape) %>%
  summarize(totWeekly = sum(tot),
            totNightTraps = sum(pmax(X.Nights*X.Traps, 1, na.rm = T)), # there are some 0 nights
            avgDayTrap = totWeekly/totNightTraps, # just do the aveage per day per trap 
            AvgAbundance = mean(avgDayTrap)) %>% 
  ungroup() %>%
  mutate(Genus = case_when(Species %in% c("erythrothorax", "quinquefasciatus", "tarsalis", "incidens", "pipiens", "thriambus","restuans") ~ "Culex",
                           Species %in% c("stigmatosoma", "inornata", "particeps") ~ "Culiseta",
                           Species %in% c("aegypti", "sierrensis", "washinoi", "notoscriptus", "increpitus") ~ "Aedes",
                           Species %in% c("franciscanus", "freeborni", "hermsi") ~ "Anopheles",
                           Species %in% c("signifera") ~ "Orthopodomyia"))

# # Just to tell, here's the species:
# 
# Culex erythrothorax
# Culex quinquefasciatus
# Culex tarsalis
# Culex incidens
# Culex pipiens
# Culex thriambus
# Culex restuans
# Culiseta stigmatosoma
# Culiseta inornata
# Culiseta particeps
# Aedes aegypti
# Aedes sierrensis
# Aedes washinoi
# Aedes notoscriptus
# Aedes increpitus
# Anopheles franciscanus
# Anopheles freeborni
# Anopheles hermsi
# Orthopodomyia signifera

saveRDS(totDFmod, file = paste0(folderDataLocal, "/totDFmod_ElDorado_Sepulveda.rds"))
