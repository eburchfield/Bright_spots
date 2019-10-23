library(tidyverse)
library(sp)
select <- dplyr::select
library(INLA)
library(spdplyr)
library(sf)
library(viridis)
library(rgeos)
library(RColorBrewer)
library(knitr)
library(randomForest)
library(MLmetrics)

county <- readRDS("./data/county.RDS")
null <- readRDS("./data/null.RDS") 
nulli <- readRDS("./out/null_indicators_p.RDS")
hadj <- readRDS("./data/neighborhood.RDS") # rook

# Data prep for models

# bin to nearest 10 (-2 does 100), resolution at which estimation occurs
null$TP<-round(null$TP,-1)  
null$GDD<-round(null$GDD,-1)
null$SDD<-round(null$SDD,-1)
null$EFP<-round(null$EFP,-1)
null$EXP<-round(null$EXP,-1)

nulli$TP<-round(nulli$TP,-1)  
nulli$GDD<-round(nulli$GDD,-1)
nulli$SDD<-round(nulli$SDD,-1)
nulli$EFP<-round(nulli$EFP,-1)
nulli$EXP<-round(nulli$EXP,-1)

# remove counties based on data availability

# (1) Drop counties with fewer than three years of data available over the five year period
yr5 <- null %>% group_by(GEOID) %>% summarize(n=n()) %>% filter(n>2)
yr5 <- unique(yr5$GEOID) # from 2203 counties to 1896 counties 
null <- null %>% filter(GEOID %in% yr5)

# (2) Drop regions with fewer than 5 counties (Cinner et al)
r <- null %>% group_by(LRR) %>% summarize(n = length(unique(GEOID)))
# no counties with thie limited number of counties though this provides an idea of sample size by LRR

# add ID
county_sub <- county %>% filter(GEOID %in% unique(null$GEOID)) %>%
  arrange(GEOID) %>% 
  mutate(STATE = as.factor(STATEFP)) %>%
  dplyr::select(GEOID, STATE) 
county_sub$ID <- 1:nrow(county_sub@data)

null <- merge(null, county_sub@data, by = "GEOID", all = T)

# regional shapefiles
lrr_shp <- readRDS("./out/lrr_shp.RDS") # build in BS_data_construction.html
lrr_shp <- sf::st_as_sf(lrr_shp)

rm(yr5)


