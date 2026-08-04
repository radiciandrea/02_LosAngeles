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
trapWeeksDF <- readRDS(file = paste0(folderDataLocal, "/trapWeeksDFsel_ElDorado_Sepulveda.rds"))

# Recompute variables
species = unique(totDFmod$Species)
species = species[-which(is.na(species))]
traps = unique(totDFmod$trap)

# Benchmark definition ----
# Detection delay
correctDateDetection = (trapWeeksDF %>%
                          filter(aegypti > 0) %>%
                          filter(datesLabels == min(datesLabels)) %>%
                          pull(datesLabels))[1]

worstDateDetection = (trapWeeksDF %>%
              filter(aegypti > 0) %>%
              filter(datesLabels == max(datesLabels)) %>%
              pull(datesLabels))[1]

# Trend quiquefaciatus
correctTrend = trapWeeksDF %>%
  group_by(datesLabels) %>%
  summarise(mqf = mean(quinquefasciatus)) %>%
  ungroup()

# Alpha diversity: Shannon index (is already normalized. Here we consider tyhe whole time series
VM = data.matrix(trapWeeksDF[,6:24])
nV = colSums(VM)
pV = nV/sum(nV)
correctShannon = - sum(pV*log2(pV))

ntraps = length(traps)

perfTrapsDF = data.frame(trap = as.factor(traps),
                     delay = NA,
                     mr = NA,
                     deltaShannon = NA,
                     shannonJensen = NA)

tic()
for(i in 1:ntraps){
    
    ti = traps[i]
    
    tempTrapWeeksDF = trapWeeksDF %>%
      filter((trap %in% ti))
  
    removedTrapWeeksDF = trapWeeksDF %>%
      filter(!(trap %in% ti))
    
    # trend ----
    temptrapsQfDF = removedTrapWeeksDF %>%
      filter(quinquefasciatus > 0)
    
    if(nrow(temptrapsQfDF)>0){
      temptrapsQfDF <- temptrapsQfDF%>%
        group_by(datesLabels) %>%
        summarise(mqf = mean(quinquefasciatus)) %>%
        ungroup() 
      
      commonProYears = correctTrend$datesLabels %in% temptrapsQfDF$datesLabels
      
      tempR = cor(correctTrend %>% filter(datesLabels %in% temptrapsQfDF$datesLabels) %>% pull(mqf), temptrapsQfDF$mqf,  method = "spearman")
    } else {
      tempR  = NA
    }
    
    # delay ----
    temptrapsAegyptiDF = tempTrapWeeksDF %>%
      filter(aegypti > 0)
    
    if(nrow(temptrapsAegyptiDF) == 0){
      dateDetection = NA
      weeksDelay = NA
    } else {
      dateDetection = (temptrapsAegyptiDF  %>%
                         filter(datesLabels == min(datesLabels)) %>%
                         pull(datesLabels)) [1]
      weeksDelay = as.numeric((dateDetection-worstDateDetection)/7)
    }
    
    
    
    # Alpha biodiversity ---- 
    # of the set without removed trap
    VMi = data.matrix(removedTrapWeeksDF[,6:24])
    nVi = colSums(VMi)
    Ni = sum(nVi) 
    pVi = nVi/Ni
    tempShannon = - sum(pVi*log2(pVi))
    
    # Beta biodiversity ---- 
    # of the set without removed trap
    M = colMeans(rbind(pV, pVi)) 
    KLD_PQ = sum(ifelse(pV >0, pV*log(pV/M), 0)) # Kullback–Leibler divergence, one side
    KLD_QP = sum(ifelse(pVi >0, pVi*log(pVi/M), 0)) # Kullback–Leibler divergence, other side
    tempJensenShannon = 0.5*KLD_PQ + 0.5*KLD_QP

    # fill data frame
    perfTrapsDF$delay[i] = weeksDelay
    perfTrapsDF$mr[i] = 1 - tempR 
    perfTrapsDF$deltaShannon[i] = abs(correctShannon - tempShannon)
    perfTrapsDF$shannonJensen[i] = tempJensenShannon
  }
toc()

ggDelay = ggplot(perfTrapsDF, aes(x = reorder(trap, delay), y = delay))+
  geom_col()+
  labs(title = "Delay looking only at trap (with respect with worst trap)", y = "weeks", x = "trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggSpearman = ggplot(perfTrapsDF, aes(x = reorder(trap, -mr), y = mr))+
  geom_col()+
  labs(title = "Distance to overall C. quinquefasciatus trend", y = "1 - Spearman's r", x = "removed trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

# ggShannon = ggplot(perfTrapsDF, aes(x = reorder(trap, -deltaShannon), y = deltaShannon))+
#   geom_col()+
#   labs(title = "Distance to correct biodiversity indicator", y = "delta Shannon", x = "removed trap")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         panel.background = element_rect(fill = "white"),
#         panel.grid = element_line(color = "gray90"))

ggShannonJensen = ggplot(perfTrapsDF, aes(x = reorder(trap, -shannonJensen), y = shannonJensen))+
  geom_col()+
  labs(title = "Distance to correct biodiversity indicator (beta)", y = "Shannon-Jensen", x = "of this trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

gtot = ggDelay / ggSpearman / ggShannonJensen

ggsave(plot = gtot, filename = paste0(folderOutput, "/G - trap importance.png"),
       device = "png", width = 7, height = 10)

