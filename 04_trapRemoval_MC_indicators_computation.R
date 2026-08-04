# Surveillance: compute surveillance degradation (reduced number of traps)
# with MC approach

library(pracma)
library(tidyverse)
library(lubridate)

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

# to add: also mosquito detection ratio 
# (number of iterations when mosquito delay is below a given ration, such as 1 year)

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

# beta diversity: like if there were two different traps.
# Options : 
# Morisita–Horn index "two-community generalization of Simpson's index"
# Bray–Curtis dissimilarity "commonly reported beta-diversity metric in ecology"
# Jensen–Shannon divergence "information-theoretic analogue, tied to Shannon entropy"
# Hellinger distance ....


# Activity season for specific species (e.g. aedes aegypti, quinquefaciatus)
# still to do

#let's thake Jensen-Shannon
# https://medium.com/@vibhorkashyap/understanding-jensen-shannon-distance-a-friendly-guide-for-data-scientists-4cac664c3381
# article for biodiversity in gut microbiome
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10628149/

# Test on categorical distribution (to rethink)

# Exhaustively (too long) ---
# # delay with 3 trap less mechansitc one: very long.
# dateDetection = c()
# 
# tic()
# for(i in 1:length(traps)){
#   s1 = traps[i]
#   traps1 = traps[-i]
#   for(j in 1:length(traps1)){
#     s2 = traps1[j]
#     traps2 = traps1[-j]
#     for(k in 1:length(traps2)){
#       s3 = traps2[k]
#       traps3 = traps1[-k]
# 
#       temptrapWeeksDF = trapWeeksDF %>%
#         filter(trap %in% traps3)
# 
#       dateDetection = c(dateDetection, (temptrapWeeksDF %>%
#                                           filter(aegypti > 0) %>%
#                                           filter(datesLabels == min(datesLabels)) %>%
#                                           pull(datesLabels))[1])
#     }
#   }
# }
# toc() # 96 seconds; quite long
# 
# weeksDelay = 52*(dateDetection-correctDateDetection)
# 
# summary(weeksDelay)

# do it with monte carlo appraoches

# e.g., time we remove 1 < n < n_traps = 36
ntraps = length(traps)-1
nRep = 1000

indicatorDF = data.frame(idrep = rep(1:nRep, ntraps),
                      nRemtraps = rep(1:ntraps, each = nRep),
                      delay = NA,
                      spearmanR = NA,
                      shannon = NA,
                      jensenShannon = NA)

tic()
for(n in 1:ntraps){
  for(r in 1:nRep){
    idr = sample(1:ntraps, n)
    trapsr = traps[-idr]
    
    temptrapWeeksDF = trapWeeksDF %>%
      filter(trap %in% trapsr)
    
    # trend ----
    temptrapsQfDF = temptrapWeeksDF %>%
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
    temptrapsAegyptiDF = temptrapWeeksDF %>%
      filter(aegypti > 0)
    
    if(nrow(temptrapsAegyptiDF) == 0){
      dateDetection = NA
    } else {
      dateDetection = (temptrapsAegyptiDF  %>%
                         filter(datesLabels == min(datesLabels)) %>%
                         pull(datesLabels)) [1]
    }
    
    weeksDelay = as.numeric((dateDetection-correctDateDetection)/7)
    
    # Alpha biodiversity ---- 
    VM = data.matrix(temptrapWeeksDF[,6:24])
    nVi = colSums(VM)
    Ni = sum(nVi) 
    pVi = nVi/Ni
    tempShannon = - sum(pVi*log2(pVi))
    
    # Beta biodiversity ---- 
    M = colMeans(rbind(pV, pVi)) 
    KLD_PQ = sum(ifelse(pV >0, pV*log(pV/M), 0)) # Kullback–Leibler divergence, one side
    KLD_QP = sum(ifelse(pVi >0, pVi*log(pVi/M), 0)) # Kullback–Leibler divergence, other side
    tempJensenShannon = 0.5*KLD_PQ + 0.5*KLD_QP
    
    # fill data frame
    indicatorDF$delay[(n-1)*nRep + r] = weeksDelay
    indicatorDF$spearmanR[(n-1)*nRep + r] = tempR 
    indicatorDF$shannon[(n-1)*nRep + r] = tempShannon 
    indicatorDF$jensenShannon[(n-1)*nRep + r] = tempJensenShannon
    
  }
  cat(n, " over ", ntraps, "\n")
}
toc() # 12 per 10 sec, 1930 per 1000

saveRDS(indicatorDF, file = paste0(folderDataLocal, "/indicatorDF_", nRep, ".rds"))
indicatorDF = readRDS(file = paste0(folderDataLocal, "/indicatorDF_", nRep, ".rds"))

# Plots ----

plot(indicatorDF$nRemtraps, indicatorDF$delay)
plot(indicatorDF$nRemtraps, indicatorDF$MDR)
plot(indicatorDF$nRemtraps, indicatorDF$spearmanR)
plot(indicatorDF$nRemtraps, indicatorDF$shannon)
plot(indicatorDF$nRemtraps, indicatorDF$jensenShannon)
# plot(indicatorDF$nRemtraps, indicatorDF$pvalCategorical)
# to rethink delay (most)

## Plot delay ----

delayDFmod <- indicatorDF %>%
  group_by(nRemtraps) %>%
  summarise(d05 = quantile(delay, 0.05, na.rm = T),
            d25 = quantile(delay, 0.25, na.rm = T),
            d50 = quantile(delay, 0.55, na.rm = T),
            dAv = mean(delay, na.rm = T),
            d75 = quantile(delay, 0.75, na.rm = T),
            d95 = quantile(delay, 0.95, na.rm = T)) %>%
  ungroup()

# delay as a function of removed traps
ggplot(data = delayDFmod, aes(x = nRemtraps, y = dAv)) +
  geom_ribbon(aes(ymin = d05, ymax = d95), alpha = 0.2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  geom_point()+
  labs(
    title = "Delay in the detection of Ae. aegypti",
    x = "Number of removed traps", y = "Delay (weeks)"
  ) 

ggsave(filename = paste0(folderOutput, "/F - Ae. aegypti detection delay, ", nRep, " reps.png"), device = "png", width = 7, height = 5)

# ggplot(data = indicatorDF, aes(x = nRemtraps, y = delay, group = nRemtraps)) +
#   geom_boxplot(fill = "gray70")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         panel.background = element_rect(fill = "white"),
#         panel.grid = element_line(color = "gray90"))+
#   labs(
#     title = "Delay in the detection of Ae. aegypti",
#     x = "Number of removed traps", y = "Delay (weeks)") 

## Plot MDR----
# MDR for aegypti: detectionw within...

MDRDFmod <- indicatorDF %>%
  group_by(nRemtraps) %>%
  summarise(no_delay = 100*sum(delay == 0, na.rm = T)/nRep, # a month
            within_month_delay = 100*sum(delay <= 5, na.rm = T)/nRep, # a season
            within_year_delay = 100*sum(delay <= 52, na.rm = T)/nRep)%>% # a year
  ungroup()

MDRDFmodPV = pivot_longer(MDRDFmod, c("no_delay", "within_month_delay", "within_year_delay"), names_to = "delay", values_to = "mdr")

ggplot(data = MDRDFmodPV, aes(x = nRemtraps, y = mdr, color = delay)) +
  geom_point()+
  geom_line()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  labs(
    title = "Ae. aegypti detection ratio",
    x = "Number of removed traps", y = "% succesful surveillance")

ggsave(filename = paste0(folderOutput, "/F - Ae. aegypti detection ratio, ", nRep, " reps.png"), device = "png", width = 7, height = 5)

## Plot quinquefaciatus trend ----
# Spearman is ok but perhaps not informative on the overall dynamics
 
quiquefaciatusDFmod <-  indicatorDF %>%
  group_by(nRemtraps) %>%
  summarise(sr05 = quantile(spearmanR, 0.05, na.rm = T),
            sr25 = quantile(spearmanR, 0.25, na.rm = T),
            sr50 = quantile(spearmanR, 0.55, na.rm = T),
            srAv = mean(spearmanR, na.rm = T),
            sr75 = quantile(spearmanR, 0.75, na.rm = T),
            sr95 = quantile(spearmanR, 0.95, na.rm = T)) %>%
  ungroup()
  
# delay as a function of removed traps
ggplot(data = quiquefaciatusDFmod, aes(x = nRemtraps, y = srAv)) +
  geom_ribbon(aes(ymin = sr05, ymax = sr95), alpha = 0.2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  geom_point()+
  labs(
    title = "Correlation (Spearman) with complete C. quiquefaciatus observations",
    x = "Number of removed traps", y = "Spearman's rank r"
  ) 

ggsave(filename = paste0(folderOutput, "/F - correlation with C. quiquefaciatus series, ", nRep, " reps.png"), device = "png", width = 7, height = 5)

# ggplot(data = indicatorDF, aes(x = nRemtraps, y = spearmanR, group = nRemtraps)) +
#   geom_boxplot(fill = "gray70")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         panel.background = element_rect(fill = "white"),
#         panel.grid = element_line(color = "gray90"))+
#   labs(
#     title = "Correlation (Spearman) with complete C. quiquefaciatus observations",
#     x = "Number of removed traps", y = "Spearman's rank r"
#   ) 

## Plot shannon ----
# Alpha Biodiversity: what to take exactly?  shannon

alphaBiodiversityDFmod <-  indicatorDF %>%
  group_by(nRemtraps) %>%
  summarise(s05 = quantile(shannon, 0.05, na.rm = T),
            s25 = quantile(shannon, 0.25, na.rm = T),
            s50 = quantile(shannon, 0.55, na.rm = T),
            sAv = mean(shannon, na.rm = T),
            s75 = quantile(shannon, 0.75, na.rm = T),
            s95 = quantile(shannon, 0.95, na.rm = T)) %>%
  ungroup()

# shannon as a function of removed traps
ggplot(data = alphaBiodiversityDFmod, aes(x = nRemtraps, y = sAv)) +
  geom_ribbon(aes(ymin = s05, ymax = s95), alpha = 0.2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  geom_point()+
  labs(
    title = "Apparent alpha-biodiversity",
    x = "Number of removed traps", y = "Shannon index"
  ) 

ggsave(filename = paste0(folderOutput, "/F - Alpha-biodiversity (Shannon), ", nRep, " reps.png"), device = "png", width = 7, height = 5)

# ggplot(data = indicatorDF, aes(x = nRemtraps, y = shannon, group = nRemtraps)) +
#   geom_boxplot(fill = "gray70")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         panel.background = element_rect(fill = "white"),
#         panel.grid = element_line(color = "gray90"))+
#   labs(
#     title = "Apparent alpha-biodiversity",
#     x = "Number of removed traps", y = "Shannon index"
#   ) 

# actually work pretty well

## Plot jensen-Shannon ----
# beta Biodiversity:  jensen

betaBiodiversityDFmod <-  indicatorDF %>%
  group_by(nRemtraps) %>%
  summarise(js05 = quantile(jensenShannon, 0.05, na.rm = T),
            js25 = quantile(jensenShannon, 0.25, na.rm = T),
            js50 = quantile(jensenShannon, 0.55, na.rm = T),
            jsAv = mean(jensenShannon, na.rm = T),
            js75 = quantile(jensenShannon, 0.75, na.rm = T),
            js95 = quantile(jensenShannon, 0.95, na.rm = T)) %>%
  ungroup()

# js as a function of removed traps
ggplot(data = betaBiodiversityDFmod, aes(x = nRemtraps, y = jsAv)) +
  geom_ribbon(aes(ymin = js05, ymax = js95), alpha = 0.2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "gray90"))+
  geom_point()+
  labs(
    title = "Apparent beta-biodiversity",
    x = "Number of removed traps", y = "Jensen-Shannon divergence"
  ) 

ggsave(filename = paste0(folderOutput, "/F - Beta-biodiversity (Jensen-Shannon), ", nRep, " reps.png"), device = "png", width = 7, height = 5)

# ggplot(data = indicatorDF, aes(x = nRemtraps, y = jensenShannon, group = nRemtraps)) +
#   geom_boxplot(fill = "gray70")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         panel.background = element_rect(fill = "white"),
#         panel.grid = element_line(color = "gray90"))+
#   labs(
#     title = "Apparent beta-biodiversity",
#     x = "Number of removed traps", y = "Jensen-Shannon divergence"
#   ) 

