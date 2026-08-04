# the most important trap

# the one which, if removed, changes the result the most

# let's take from code MC

library(pracma)
library(tidyverse)
library(lubridate)
library(ggplot2)

folderDataLocal = "Data"
folderOutput = "Outputs"

# load data
totDFmod <- readRDS(file = paste0(folderDataLocal, "/totDFsel_ElDorado_Sepulveda.rds"))
siteWeeksDF <- readRDS(file = paste0(folderDataLocal, "/siteWeeksDFsel_ElDorado_Sepulveda.rds"))

# Recompute variables
species = unique(totDFmod$Species)
species = species[-which(is.na(species))]
sites = unique(totDFmod$siteCode)

# Benchmark definition ----
# Detection delay
correctDateDetection = (siteWeeksDF %>%
                          filter(aegypti > 0) %>%
                          filter(datesLabels == min(datesLabels)) %>%
                          pull(datesLabels))[1]

# Trend quiquefaciatus
correctTrend = siteWeeksDF %>%
  group_by(datesLabels) %>%
  summarise(mqf = mean(quinquefasciatus)) %>%
  ungroup()

# Alpha diversity: Shannon index (is already normalized. Here we consider tyhe whole time series

nV = c(sum(siteWeeksDF$quinquefasciatus), sum(siteWeeksDF$tarsalis), sum(siteWeeksDF$stigmatosoma), sum(siteWeeksDF$aegypti))
pV = nV/sum(nV)
correctShannon = - sum(pV*log2(pV))

# create a datafrem for each

# DA RIVEDERE!
trapDF = data.frame(trap = totDFmod %>%
  mutate(trap = paste0(siteCode, "_", TrapType)) %>%
  select(trap))
  
nSites = length(sites)

sitesDF = data.frame(site = as.factor(sites),
                     delay = NA,
                     mr = NA,
                     deltaShannon = NA)

tic()
for(i in 1:nSites){
    
    si = sites[i]
  
    tempsiteWeeksDF = siteWeeksDF %>%
      filter(!(site %in% si))
    
    # trend ----
    tempsitesQfDF = tempsiteWeeksDF %>%
      filter(quinquefasciatus > 0)
    
    if(nrow(tempsitesQfDF)>0){
      tempsitesQfDF <- tempsitesQfDF%>%
        group_by(datesLabels) %>%
        summarise(mqf = mean(quinquefasciatus)) %>%
        ungroup() 
      
      commonProYears = correctTrend$datesLabels %in% tempsitesQfDF$datesLabels
      
      tempR = cor(correctTrend %>% filter(datesLabels %in% tempsitesQfDF$datesLabels) %>% pull(mqf), tempsitesQfDF$mqf,  method = "spearman")
    } else {
      tempR  = NA
    }
    
    # delay ----
    tempsitesAegyptiDF = tempsiteWeeksDF %>%
      filter(aegypti > 0)
    
    if(nrow(tempsitesAegyptiDF) == 0){
      dateDetection = NA
    } else {
      dateDetection = (tempsitesAegyptiDF  %>%
                         filter(datesLabels == min(datesLabels)) %>%
                         pull(datesLabels)) [1]
    }
    
    weeksDelay = as.numeric((dateDetection-correctDateDetection)/7)
    
    # Alpha biodiversity ---- 
    nVi = c(sum(tempsiteWeeksDF$quinquefasciatus), sum(tempsiteWeeksDF$tarsalis), sum(tempsiteWeeksDF$stigmatosoma), sum(tempsiteWeeksDF$aegypti))
    N = sum(nVi) 
    pVi = nVi/N
    tempShannon = - sum(pVi*log2(pVi))
    

    # fill data frame
    sitesDF$delay[i] = weeksDelay
    sitesDF$mr[i] = 1 - tempR 
    sitesDF$deltaShannon[i] = abs(correctShannon - tempShannon)
  }
toc()

ggDelay = ggplot(sitesDF, aes(x = reorder(site, -delay), y = delay))+
  geom_col()+
  labs(title = "Delay due to trap removal", y = "weeks", x = "removed site")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggSpearman = ggplot(sitesDF, aes(x = reorder(site, -mr), y = mr))+
  geom_col()+
  labs(title = "Distance to overall C. quinquefasciatus trend", y = "1 - Spearman's r", x = "removed site")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggShannon = ggplot(sitesDF, aes(x = reorder(site, -deltaShannon), y = deltaShannon))+
  geom_col()+
  labs(title = "Distance to correct biodiversity indicator", y = "delta Shannon", x = "removed site")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

gtot = ggDelay / ggSpearman / ggShannon

ggsave(plot = gtot, filename = paste0(folderOutput, "/F - Site importance.png"),
       device = "png", width = 7, height = 5)

TRAP = SITE_TRAPTYPE