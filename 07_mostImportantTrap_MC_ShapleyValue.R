library(sensitivity)
library(pracma)
library(tidyverse)

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

# Shapley value attribution MC algorithm
# https://en.wikipedia.org/wiki/Shapley_value

correctDateDetection = (trapWeeksDF %>%
                          filter(aegypti > 0) %>%
                          filter(datesLabels == min(datesLabels)) %>%
                          pull(datesLabels))[1]

indicatorFun <- function(ti = traps){
  
  tempTrapWeeksDF = trapWeeksDF %>%
    dplyr::filter(trap %in% ti)

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
  
  return(weeksDelay)
}


# # Shapley value attribution setup
contrib <- matrix(0, nPerm, nTraps)

# Shapley value attribution MC
for (p in seq_len(nPerm)) {
  perm <- sample(traps)
  includedTraps <- c()
  prev_val <- indicatorFun(ti = includedTraps)   # indicator with no traps (baseline)
  for (i in 1:nTraps) {
    includedTraps = c(includedTraps, perm[i])
    new_val <- indicatorFun(ti = includedTraps)
    contrib[p, which(traps == perm[i])] <- new_val - prev_val
    prev_val <- new_val
  }
  cat(p, "\n")
}



# run shapley_mc nPerm = 100

contrib = shapley_mc(indicatorFun = indicatorFun, traps = traps, nPerm = 10)

# elab
shapleyV <- colMeans(contrib, na.rm = T)

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
