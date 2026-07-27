# Survaillance: compute survaillance degradation

library(pracma)
library(tidyverse)

folderDataLocal = "Data"

# load data
totDFmod <- readRDS(file = paste0(folderDataLocal, "/totDFmod_ElDorado_Sepulveda.rds"))
sitesDF <- readRDS(file = paste0(folderDataLocal, "/sitesDF_ElDorado_Sepulveda.rds"))

# Recompute variables
species = unique(totDF$species)
sites = unique(totDF$site_code)
traps = unique(totDF$trap_type)

# let's plot the actovity period for each trap (in weeks)
year_start = year(min(totDF$collection_date))
year_end = year(max(totDF$collection_date))

years = rep(year_start:year_end, each = 52)
weeks = rep(1:52, times = length(year_start:year_end))

# Detection delay

correctDateDetection = (sitesDF %>%
                          filter(aegypti > 0) %>%
                          filter(progYear == min(progYear)) %>%
                          pull(progYear))[1]

# delay with 3 site less mechansitc one: very long.
dateDetection = c()

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
      
      dateDetection = c(dateDetection, (tempsitesDF %>%
                                          filter(aegypti > 0) %>%
                                          filter(progYear == min(progYear)) %>%
                                          pull(progYear))[1])
    }
  }
}
toc() # 96 seconds; quite long

weeksDelay = 52*(dateDetection-correctDateDetection)

summary(weeksDelay)

# do it with monte carlo appraoches

# e.g., time we remove 1 < n < n_sites = 36
nSites = length(sites)-1
nRep = 100

delaysDF = data.frame(idrep = rep(1:nRep, nSites),
                      nRemSites = rep(1:nSites, each = nRep),
                      delay = NA)

tic()
for(n in 1:nSites){
  for(r in 1:nRep){
    idr = sample(1:nSites, n)
    sitesr = sites[-idr]
    
    tempsitesDF = sitesDF %>%
      filter(site %in% sitesr)%>%
      filter(aegypti > 0)
    
    if(nrow(tempsitesDF) == 0){
      dateDetection = NA
    } else {
      dateDetection = (tempsitesDF  %>%
                         filter(progYear == min(progYear)) %>%
                         pull(progYear)) [1]
    }
    
    weeksDelay = 52*(dateDetection-correctDateDetection)
    
    delaysDF$delay[(n-1)*nRep + r] = weeksDelay
    
  }
}
toc() # 7.6 sec per 100

plot(delaysDF$nRemSites, delaysDF$delay)

# to rethink

delaysDFmod <- delaysDF %>%
  group_by(nRemSites) %>%
  summarise(d25 = quantile(delay, 0.25, na.rm = T),
            d50 = quantile(delay, 0.55, na.rm = T),
            dAv = mean(delay, na.rm = T),
            d75 = quantile(delay, 0.75, na.rm = T)) %>%
  ungroup()

ggplot(data = delaysDFmod, aes(x = nRemSites, y = dAv)) +
  geom_point()

ggplot(data = delaysDF, aes(x = nRemSites, y = delay, group = nRemSites)) +
  geom_boxplot()