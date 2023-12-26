library(dplyr)
# subsample data based on propensity score matching

alldata <- readRDS("nhanes_complete_data.rds")

dattable <- alldata$data
colmeta <- alldata$coldata

## remove individuals who are not on diet or anti-hypertension medication
dattable_subset <- dattable %>% filter(bp_med_use == "No" & DRQSDIET != 1 & 
                                         !is.na(INDFMPIR))
  
dattable_subset$svy_year <- factor(dattable_subset$svy_year,
                                   levels=c("2011-2012", "2013-2014", "2015-2016", "2017-2020"))

library(WeightIt)
library(cobalt)
## initial balance of covariates
init.bal <- bal.tab(svy_year ~ demo_age_years + demo_race + demo_gender + INDFMPIR,
                    data = dattable_subset, estimand = "ATE", thresholds = c(m = .05))
init.bal

## calculate the weights based on CBPS
library(WeightIt)
library(CBPS)
W.out <- weightit(svy_year ~ demo_age_years + demo_race + demo_gender + INDFMPIR,
                  data = dattable_subset, estimand = "ATE", 
                  method = "cbps", over=F)
summary(W.out)
## check if weighting is done
bal.tab(W.out, stats = c("m", "v"), thresholds = c(m = .05))
dattable_subset$psm_weights <- W.out$weights

saveRDS(dattable_subset, file="selected_sample_psm.rds")
