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

#Benchmark definition ----

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

# Alpha diversity: Shannon index 
VM = data.matrix(trapWeeksDF[,6:24])
nV = colSums(VM)
pV = nV/sum(nV)
correctShannon = - sum(pV*log2(pV))

indicatorFun <- function(ti = traps){
 
 # subset only ti traps
 tempTrapWeeksDF = trapWeeksDF %>%
 dplyr::filter(trap %in% ti)

 #Trend and seasonality----
 temptrapsQfDF = tempTrapWeeksDF %>%
 dplyr::filter(quinquefasciatus > 0)
 
 if(nrow(temptrapsQfDF)>0){
 
 ##Trend----
 
 temptrapsQfDF <- temptrapsQfDF%>%
  group_by(datesLabels) %>%
  dplyr::summarise(mqf = mean(quinquefasciatus)) %>%
  ungroup() 
 
 tempR = cor(correctTrendDF %>%
     dplyr::filter(datesLabels %in% temptrapsQfDF$datesLabels)
    %>% pull(mqf),
    temptrapsQfDF$mqf, method = "spearman")
 
 #Peak(rmse)---- 
 
 tempPeakDF = temptrapsQfDF %>%
  mutate(year = year(datesLabels)) %>%
  group_by(year) %>%
  dplyr::filter(mqf == max(mqf))%>%
  ungroup() #weird trend
 
 tempPeakDateError = rmse(as.numeric(correctPeakDF %>%
           dplyr::filter(year %in% tempPeakDF$year)
           %>% dplyr::pull(datesLabels))/7,
        as.numeric(tempPeakDF$datesLabels)/7)
 
 #Start(rmse)---- 
 tempStartDF = temptrapsQfDF %>%
  mutate(year = year(datesLabels)) %>%
  group_by(year) %>%
  filter(datesLabels == min(datesLabels))%>%
  ungroup() #weird trend
 
 tempSeasonStartError = rmse(as.numeric(correctStartDF %>%
            dplyr::filter(year %in% tempStartDF$year) %>%
            dplyr::pull(datesLabels))/7,
        as.numeric(tempStartDF$datesLabels)/7)
 
 #End(rmse)---- 
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
 tempR = 0
 tempPeakDateError = NA
 tempSeasonStartError = NA 
 tempSeasonEndError = NA
 }
 
 
 #Delay ---- 
 tempTrapsAegyptiDF = tempTrapWeeksDF %>%
 dplyr::filter(aegypti > 0)
 
 if(nrow(tempTrapsAegyptiDF) == 0){
 weeksDelay = NA
 } else {
 dateDetection = (tempTrapsAegyptiDF %>%
      filter(datesLabels == min(datesLabels)) %>%
      pull(datesLabels)) [1]
 weeksDelay = as.numeric((dateDetection-correctDateDetection)/7)
 }
 
 # Alpha biodiversity---- 
 VM = data.matrix(tempTrapWeeksDF[,6:24])
 nVi = colSums(VM)
 # remove 0s
 Ni = sum(nVi) 
 pVi = nVi/Ni
 tempShannon = - sum(pVi[which(nVi>0)]*log2(pVi[which(nVi>0)]))
 
 
 return(c(weeksDelay, tempR, tempPeakDateError, tempSeasonStartError, tempSeasonEndError, tempShannon))
}

# Shapley value attribution setup---
nPerm = 1000
contribDelay <- matrix(0, nPerm, nTraps)
contribR <- matrix(0, nPerm, nTraps)
contribPeakError <- matrix(0, nPerm, nTraps)
contribSeasonStartError <- matrix(0, nPerm, nTraps)
contribSeasonEndError <- matrix(0, nPerm, nTraps)
contribShannon <- matrix(0, nPerm, nTraps)

tic()
#Shapley value attribution loop----
for (p in seq_len(nPerm)) {
 includedTraps <- rep("", times = nTraps)
 prev_val <- indicatorFun(ti = includedTraps) # indicator with no traps (baseline)
 perm <- sample(traps) # it is important t take them randomly, otherwise some will be always taken after others
 for (i in seq_len(nTraps)) {
 includedTraps = c(includedTraps, perm[i])
 new_val <- indicatorFun(ti = includedTraps)
 j = which(traps == perm[i])
 contribDelay[p, j] <- new_val[1] - prev_val[1]
 contribR[p, j] <- new_val[2] - prev_val[2]
 contribPeakError[p, j] <- new_val[3] - prev_val[3]
 contribSeasonStartError[p, j] <- new_val[4] - prev_val[4]
 contribSeasonEndError[p, j] <- new_val[5] - prev_val[5]
 contribShannon[p, j] <- new_val[6] - prev_val[6]
 
 prev_val <- new_val
 }
 cat(p, "\n")
}
toc()
# with 1000 permutations: 2100 seconds

#sum up, save, and load----

ShapleyDF = data.frame(trap = traps,
      trapType = sapply(traps, function(w){substr(w, 6, nchar(w))}),
      shapleyDelay = colMeans(contribDelay, na.rm = T),
      shapleyR = colMeans(contribR, na.rm = T),
      shapleyPeakError = colMeans(contribPeakError, na.rm = T),
      shapleySeasonStartError = colMeans(contribSeasonStartError, na.rm = T),
      shapleySeasonEndError = colMeans(contribSeasonEndError, na.rm = T),
      shapleyShannon = colMeans(contribShannon, na.rm = T))

saveRDS(ShapleyDF, file = paste0(folderDataLocal, "/ShapleyDF_", nPerm, ".rds"))

ShapleyDF = readRDS(file = paste0(folderDataLocal, "/ShapleyDF_", nPerm, ".rds"))

# elaborate per trap type

trapType_ShapleyDF = ShapleyDF %>%
 group_by(trapType) %>%
 summarise(shapleyMeanDelay = mean(shapleyDelay),
   shapleyMeanR = mean(shapleyR),
   shapleyMeanPeakError = mean(shapleyPeakError),
   shapleyMeanSeasonStartError = mean(shapleySeasonStartError),
   shapleyMeanSeasonEndError = mean(shapleySeasonEndError),
   shapleyMeanShannon = mean(shapleyShannon),
   shapley05Delay = quantile(shapleyDelay, 0.05),
   shapley05R = quantile(shapleyR, 0.05),
   shapley05PeakError = quantile(shapleyPeakError, 0.05),
   shapley05SeasonStartError = quantile(shapleySeasonStartError, 0.05),
   shapley05SeasonEndError = quantile(shapleySeasonEndError, 0.05),
   shapley05Shannon = quantile(shapleyShannon, 0.05),
   shapley95Delay = quantile(shapleyDelay, 0.95),
   shapley95R = quantile(shapleyR, 0.95),
   shapley95PeakError = quantile(shapleyPeakError, 0.95),
   shapley95SeasonStartError = quantile(shapleySeasonStartError, 0.95),
   shapley95SeasonEndError = quantile(shapleySeasonEndError, 0.95),
   shapley95Shannon = quantile(shapleyShannon, 0.95))%>%
 ungroup()

# Plots----

## Early Ae. aegypti detection----
ggplot(ShapleyDF, aes(x = reorder(trap, shapleyDelay), y = shapleyDelay, fill = trapType))+
 geom_col()+
 labs(title = "Importance of traps in Ae. aegypti detection", y = "Per-trap Shapley value (delay in weeks)", x = "trap")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1),
  panel.background = element_rect(fill = "white"),
  panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, Ae. aegypti detection delay, ", nPerm, " reps, per trap.png"), device = "png", width = 10, height = 5)

#per Trap type + error
ggplot(trapType_ShapleyDF, aes(x = reorder(trapType, shapleyMeanDelay), y = shapleyMeanDelay, fill = trapType))+
 geom_col()+
 geom_errorbar(aes(ymin = shapley05Delay, ymax = shapley95Delay)) +
 labs(title = "Importance of trap type in Ae. aegypti detection", y = "Per-type average Shapley value (delay in weeks)", x = "trap")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1),
  panel.background = element_rect(fill = "white"),
  panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, Ae. aegypti detection delay, ", nPerm, " reps, per trap type.png"), device = "png", width = 5, height = 5)

## Biodiversity (Shannon)----
ggplot(ShapleyDF, aes(x = reorder(trap, -shapleyShannon), y = shapleyShannon, fill = trapType))+
  geom_col()+
  labs(title = "Importance of traps in determining Shannon index", y = "Per-trap Shapley value (rmse)", x = "trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, Shannon index, ", nPerm, " reps, per trap.png"), device = "png", width = 10, height = 5)

#per Trap type
ggplot(trapType_ShapleyDF, aes(x = reorder(trapType, -shapleyMeanShannon), y = shapleyMeanShannon, fill = trapType))+
  geom_col()+
  geom_errorbar(aes(ymin = shapley05Shannon, ymax = shapley95Shannon)) +
  labs(title = "Importance of trap type in determining Shannon index", y = "Per-type average Shapley value (rmse)", x = "trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, Shannon index, ", nPerm, " reps, per trap type.png"), device = "png", width = 5, height = 5)

## C. quinquefasciatus correct trend----
ggplot(ShapleyDF, aes(x = reorder(trap, -shapleyR), y = shapleyR, fill = trapType))+
 geom_col()+
 labs(title = "Importance of traps in C. quinquefasciatus trend", y = "Per-trap Shapley value (Spearman's r)", x = "trap")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1),
  panel.background = element_rect(fill = "white"),
  panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, C. quinquefasciatus trend, ", nPerm, " reps, per trap.png"), device = "png", width = 10, height = 5)

#per Trap type
ggplot(trapType_ShapleyDF, aes(x = reorder(trapType, -shapleyMeanR), y = shapleyMeanR, fill = trapType))+
 geom_col()+
  geom_errorbar(aes(ymin = shapley05R, ymax = shapley95R)) +
 labs(title = "Importance of trap type in C. quinquefasciatus trend", y = "Per-type average Shapley value (Spearman's r)", x = "trap")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1),
  panel.background = element_rect(fill = "white"),
  panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, C. quinquefasciatus trend, ", nPerm, " reps, per trap type.png"), device = "png", width = 5, height = 5)

## C. quinquefasciatus correct seasonal peak----
ggplot(ShapleyDF, aes(x = reorder(trap, shapleyPeakError), y = shapleyPeakError, fill = trapType))+
 geom_col()+
 labs(title = "Importance of traps in C. quinquefasciatus seasonal peak", y = "Per-trap Shapley value (rmse)", x = "trap")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1),
  panel.background = element_rect(fill = "white"),
  panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, C. quinquefasciatus seasonal peak, ", nPerm, " reps, per trap.png"), device = "png", width = 10, height = 5)

#per Trap type
ggplot(trapType_ShapleyDF, aes(x = reorder(trapType, shapleyMeanPeakError), y = shapleyMeanPeakError, fill = trapType))+
 geom_col()+
 geom_errorbar(aes(ymin = shapley05PeakError, ymax = shapley95PeakError)) +
 labs(title = "Importance of trap type in C. quinquefasciatus seasonal peak", y = "Per-type average Shapley value (rmse)", x = "trap")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1),
  panel.background = element_rect(fill = "white"),
  panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, C. quinquefasciatus seasonal peak, ", nPerm, " reps, per trap type.png"), device = "png", width = 5, height = 5)

## C. quinquefasciatus correct seasonal start----
ggplot(ShapleyDF, aes(x = reorder(trap, shapleySeasonStartError), y = shapleySeasonStartError, fill = trapType))+
  geom_col()+
  labs(title = "Importance of traps in C. quinquefasciatus seasonal start", y = "Per-trap Shapley value (rmse)", x = "trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, C. quinquefasciatus seasonal start, ", nPerm, " reps, per trap.png"), device = "png", width = 10, height = 5)

#per Trap type
ggplot(trapType_ShapleyDF, aes(x = reorder(trapType, shapleyMeanSeasonStartError), y = shapleyMeanSeasonStartError, fill = trapType))+
  geom_col()+
  geom_errorbar(aes(ymin = shapley05SeasonStartError, ymax = shapley95SeasonStartError)) +
  labs(title = "Importance of trap type in C. quinquefasciatus seasonal peak", y = "Per-type average Shapley value (rmse)", x = "trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, C. quinquefasciatus seasonal start, ", nPerm, " reps, per trap type.png"), device = "png", width = 5, height = 5)

## C. quinquefasciatus correct seasonal end----
ggplot(ShapleyDF, aes(x = reorder(trap, shapleySeasonEndError), y = shapleySeasonEndError, fill = trapType))+
  geom_col()+
  labs(title = "Importance of traps in C. quinquefasciatus seasonal end", y = "Per-trap Shapley value (rmse)", x = "trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, C. quinquefasciatus seasonal end, ", nPerm, " reps, per trap.png"), device = "png", width = 10, height = 5)

#per Trap type
ggplot(trapType_ShapleyDF, aes(x = reorder(trapType, shapleyMeanSeasonEndError), y = shapleyMeanSeasonEndError, fill = trapType))+
  geom_col()+
  geom_errorbar(aes(ymin = shapley05SeasonEndError, ymax = shapley95SeasonEndError)) +
  labs(title = "Importance of trap type in C. quinquefasciatus seasonal peak", y = "Per-type average Shapley value (rmse)", x = "trap")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))

ggsave(filename = paste0(folderOutput, "/G - Shapley value, C. quinquefasciatus seasonal end, ", nPerm, " reps, per trap type.png"), device = "png", width = 5, height = 5)
