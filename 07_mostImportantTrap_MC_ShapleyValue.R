# Shapley value attribution MC algorithm
# https://en.wikipedia.org/wiki/Shapley_value

library(sensitivity)
library(pracma)
library(tidyverse)
library(Metrics)

folderDataLocal = "Data"
folderOutput = "Outputs"

# load data
totDFmod <- readRDS(file = paste0(folderDataLocal, "/totDFsel_ElDorado_Sepulveda.rds"))
trapWeeksDF <- readRDS(file = paste0(folderDataLocal, "/trapWeeksDFsel_ElDorado_Sepulveda.rds"))

# Recompute variables
species = unique(totDFmod$Species)
species = species[-which(is.na(species))]
traps = unique(totDFmod$trap)
nTraps = length(traps)

# Benchmark definition ----
# Detection delay
correctDateDetection = (trapWeeksDF %>%
                          filter(aegypti > 0) %>%
                          filter(datesLabels == min(datesLabels)) %>%
                          pull(datesLabels))[1]

# Trend quinquefasciatus
correctTrendDF = trapWeeksDF %>%
  group_by(datesLabels) %>%
  summarise(mqf = mean(quinquefasciatus)) %>%
  ungroup()

# Peak quinquefasciatus
correctPeakDF = correctTrendDF %>%
  mutate(year = year(datesLabels)) %>%
  group_by(year) %>%
  filter(mqf == max(mqf))%>%
  ungroup() #weird trend

# Season start quinquefasciatus
correctStartDF = correctTrendDF %>%
  mutate(year = year(datesLabels)) %>%
  filter(mqf > 0) %>%
  group_by(year) %>%
  filter(datesLabels == min(datesLabels))%>%
  ungroup() 

# Season end quinquefasciatus
correctEndDF = correctTrendDF %>%
  mutate(year = year(datesLabels)) %>%
  filter(mqf > 0) %>%
  group_by(year) %>%
  filter(datesLabels == max(datesLabels))%>%
  ungroup() 

# Alpha diversity: Shannon index (is already normalized. Here we consider tyhe whole time series
VM = data.matrix(trapWeeksDF[,6:24])
nV = colSums(VM)
pV = nV/sum(nV)
correctShannon = - sum(pV*log2(pV))

indicatorFun <- function(ti = traps){
  
  tempTrapWeeksDF = trapWeeksDF %>%
    dplyr::filter(trap %in% ti)

  # trend and seasonality  ----
  temptrapsQfDF = tempTrapWeeksDF %>%
    dplyr::filter(quinquefasciatus > 0)
  
  if(nrow(temptrapsQfDF)>0){
    # tend
    temptrapsQfDF <- temptrapsQfDF%>%
      group_by(datesLabels) %>%
      dplyr::summarise(mqf = mean(quinquefasciatus)) %>%
      ungroup() 
    
    tempR = cor(correctTrendDF %>%
                  dplyr::filter(datesLabels %in% temptrapsQfDF$datesLabels)
                %>% pull(mqf),
                temptrapsQfDF$mqf,  method = "spearman")
    
    # seasonality (rmse)
    
    tempPeakDF = temptrapsQfDF %>%
      mutate(year = year(datesLabels)) %>%
      group_by(year) %>%
      dplyr::filter(mqf == max(mqf))%>%
      ungroup() #weird trend
    
    tempPeakDateError =  rmse(as.numeric(correctPeakDF %>%
                                           dplyr::filter(year %in% tempPeakDF$year)
                                         %>% dplyr::pull(datesLabels))/7,
                              as.numeric(tempPeakDF$datesLabels)/7)
    
    tempStartDF = temptrapsQfDF %>%
      mutate(year = year(datesLabels)) %>%
      group_by(year) %>%
      filter(datesLabels == min(datesLabels))%>%
      ungroup() #weird trend
    
    tempSeasonStartError = rmse(as.numeric(correctStartDF %>%
                                             dplyr::filter(year %in% tempStartDF$year) %>%
                                             dplyr::pull(datesLabels))/7,
                                as.numeric(tempStartDF$datesLabels)/7)
    
    tempEndDF = temptrapsQfDF %>%
      mutate(year = year(datesLabels)) %>%
      group_by(year) %>%
      filter(datesLabels == max(datesLabels))%>%
      ungroup() #weird trend
    
    tempSeasonEndError = rmse(as.numeric(correctEndDF %>%
                                           dplyr::filter(year %in% tempEndDF$year) %>%
                                           dplyr::pull(datesLabels))/7,
                              as.numeric(tempEndDF$datesLabels)/7)
    
  } else {
    tempR  = 0
    tempPeakDateError = NA
    tempSeasonStartError = NA 
    tempSeasonEndError = NA
  }
  
  
  # delay ----
  tempTrapsAegyptiDF = tempTrapWeeksDF %>%
    dplyr::filter(aegypti > 0)
  
  if(nrow(tempTrapsAegyptiDF) == 0){
    weeksDelay = NA
  } else {
    dateDetection = (tempTrapsAegyptiDF  %>%
                       filter(datesLabels == min(datesLabels)) %>%
                       pull(datesLabels)) [1]
    weeksDelay = as.numeric((dateDetection-correctDateDetection)/7)
  }
  
  
  
  # Alpha biodiversity ---- 
  VM = data.matrix(tempTrapWeeksDF[,6:24])
  nVi = colSums(VM)
  # remove 0s
  Ni = sum(nVi) 
  pVi = nVi/Ni
  tempShannon = - sum(pVi[which(nVi>0)]*log2(pVi[which(nVi>0)]))
  
  
  return(c(weeksDelay, tempR, tempPeakDateError, tempSeasonStartError, tempSeasonEndError, tempShannon))
}

# # Shapley value attribution setup
nPerm = 100
contribDelay <- matrix(0, nPerm, nTraps)
contribR <- matrix(0, nPerm, nTraps)
contribPeakError <- matrix(0, nPerm, nTraps)
contribSeasonStart <- matrix(0, nPerm, nTraps)
contribSeasonEnd <- matrix(0, nPerm, nTraps)
contribShannon <- matrix(0, nPerm, nTraps)

tic()
# Shapley value attribution MC
for (p in seq_len(nPerm)) {
  includedTraps <- rep("", times = nTraps)
  prev_val <- indicatorFun(ti = includedTraps)   # indicator with no traps (baseline)
  perm <- sample(traps) # it is important t take them randomly, otherwise some will be always taken after others
  for (i in seq_len(nTraps)) {
    includedTraps = c(includedTraps, perm[i])
    new_val <- indicatorFun(ti = includedTraps)
    j = which(traps == perm[i])
    contribDelay[p, j] <- new_val[1] - prev_val[1]
    contribR[p, j] <- new_val[2] - prev_val[2]
    contribPeakError[p, j] <- new_val[3] - prev_val[3]
    contribSeasonStart[p, j] <- new_val[4] - prev_val[4]
    contribSeasonEnd[p, j] <- new_val[5] - prev_val[5]
    contribShannon[p, j] <- new_val[6] - prev_val[6]
    
    prev_val <- new_val
  }
  cat(p, "\n")
}
toc()


# elab
shapleyV <- colMeans(contribDelay, na.rm = T)

# Early detection
shapleyV <- colMeans(contribShannon, na.rm = T)

trap_delayShapley = data.frame(trap = traps,
                              trapType = sapply(X = traps, FUN = function(w){substr(w, 6, nchar(w))}),
                              shapley = shapleyV)

ggplot(trap_delayShapley, aes(x = reorder(trap, shapley), y = shapley, fill = trapType))+
  geom_col()+
  labs(title = "Importance of trap type in Ae. aegypti detection", y = "Per-trap Shapley value (delay in weeks)", x = "trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

# per Trap type 

trapType_delayShapley = trap_delayShapley %>%
  group_by(trapType) %>%
  summarise(averageShapley = mean(shapley))%>%
  ungroup()

ggplot(trapType_delayShapley, aes(x = reorder(trapType, averageShapley), y = averageShapley, fill = trapType))+
  geom_col()+
  labs(title = "Importance of trap type in Ae. aegypti detection", y = "Per-type average Shapley value (delay in weeks)", x = "trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))
