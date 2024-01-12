## ----setup, include=FALSE--------------------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo=T, warning = FALSE, message = FALSE) 


## ----loaddata--------------------------------------------------------------------------
alldata <- readRDS("nhanes_complete_data.rds")
dattable <- alldata$data # main dataset
colmeta <- alldata$coldata # metadata of columns


## --------------------------------------------------------------------------------------
dattable$svy_year <- factor(dattable$svy_year,
                                   levels=c("2011-2012", "2013-2014", "2015-2016", "2017-2020"))

summary(dattable$svy_year)

summary(dattable$demo_age_years)

summary(dattable$demo_gender)

summary(dattable$demo_race)

summary(dattable$INDFMPIR)


## ----message=FALSE---------------------------------------------------------------------
library(dplyr)
dattable_subset <- dattable %>% filter(bp_med_use == "No" & DRQSDIET != 1 &
                                         !is.na(DR1TKCAL) & !is.na(INDFMPIR))


## --------------------------------------------------------------------------------------

summary(dattable_subset$svy_year)

summary(dattable_subset$demo_age_years)

summary(dattable_subset$demo_gender)

summary(dattable$demo_race)

summary(dattable_subset$INDFMPIR)


## --------------------------------------------------------------------------------------
library(WeightIt)
library(cobalt)
library(CBPS)
init.bal <- bal.tab(svy_year ~ demo_age_years + demo_race + demo_gender + INDFMPIR,
                    data = dattable_subset, 
                    estimand = "ATE", 
                    thresholds = c(m = .05))
#init.bal

love.plot(init.bal)


## ----message=F-------------------------------------------------------------------------
W.out <- weightit(svy_year ~ demo_age_years + demo_race + demo_gender + INDFMPIR,
                  data = dattable_subset, estimand = "ATE", 
                  method = "cbps", over=F)

love.plot(W.out)

summary(W.out)


## --------------------------------------------------------------------------------------
bal.tab(W.out, stats = c("m", "v", "ks"), thresholds = c(m = .05))


## ----message=F-------------------------------------------------------------------------
dattable_subset$ps_weights <- W.out$weights
write.csv(dattable_subset, "selected_sample_IPW.csv",
          row.names=F)
saveRDS(dattable_subset, file="selected_sample_IPW.rds")

