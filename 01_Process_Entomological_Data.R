# preprocess entomological observation files

library(xlsx) #4.5
library(pracma)
library(tidyverse)

# datafodler

folderData = "C:/Users/2024ar003/Desktop/Alcuni file permanenti/Post_doc_biodivecity/Dati/US_Los_Angeles"
folderDataLocal = "Data"

dir.create(folderDataLocal)
# Preprocessing ----

#  save data as RDS (done once for all)
tic()
ElDoradoDF = read.xlsx(file = paste0(folderData, "/El Dorado Park abundance and arbovirus data.xlsx"),
                       sheetName = "Sheet2")
SepulvedaDF = read.xlsx(file = paste0(folderData, "/Sepulveda Basin abundance and arbovirus data.xlsx"),
                        sheetName = "Sheet2")
toc()

# read.xlsx is super slow:

# the two tables have the same items (rbind is ok)

ElDoradoDF$area = "El Dorado"
ElDoradoDF$type = "green area"
SepulvedaDF$area = "Sepulveda"
SepulvedaDF$type = "green area"

#rbind

totDF = rbind(ElDoradoDF, SepulvedaDF)

saveRDS(totDF, file = paste0(folderDataLocal, "/Abundance_ElDorado_Sepulveda.rds"))

totDF <- readRDS(file = paste0(folderDataLocal, "/Abundance_ElDorado_Sepulveda.rds"))


# Merging pools ----

# These are test for the disease: therefore rows are repeated for each sample. Le

totDFmod <- totDF %>%
  group_by(area, pool_id, trap_type, site_code, collection_date, species, sex, sex_condition) %>%
  summarise(num_count = mean(num_count)) %>%
  ungroup()

# if they belong to the same pool id, it means that they are the same mosquito that are tested over only multiple diseases.
# however, if they belong to different pool id, it means that they are different mosquitos - there fore they should be summed.

totDFmod <- totDF %>%
  group_by(area, trap_type, site_code, collection_date, species, sex, sex_condition) %>%
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

totDFmod$colchk = NULL

# ok ther are no more duplicated lines
saveRDS(totDFmod, file = paste0(folderDataLocal, "/totDFmod_ElDorado_Sepulveda.rds"))
