library(dplyr)

data <- readRDS("selected_sample_IPW.rds")
col_keys <- read.csv("nhanes_complete_key.csv")

## remove individuals who don't have any intake of protein/carbohydrate/fat
subset_data <- data %>% filter(DR1TPROT > 0 & DR1TCARB > 0 & DR1TTFAT > 0)

## impute zero mass with 0.001
subset_data$DR1TSUGR[subset_data$DR1TSUGR == 0] <- 0.001
subset_data$DR1TFIBE[subset_data$DR1TFIBE == 0] <- 0.001
subset_data$starch <- subset_data$DR1TCARB - subset_data$DR1TSUGR - subset_data$DR1TFIBE
subset_data$starch[subset_data$starch == 0] <- 0.001
subset_data$DR1TSFAT[subset_data$DR1TSFAT == 0] <- 0.001
subset_data$unsat_fatacid <- subset_data$DR1TMFAT + subset_data$DR1TPFAT
subset_data$unsat_fatacid[subset_data$unsat_fatacid == 0] <- 0.001
# subset_data$DR1TMFAT[subset_data$DR1TMFAT == 0] <- 0.001
# subset_data$DR1TPFAT[subset_data$DR1TPFAT == 0] <- 0.001
subset_data$DR1TALCO[subset_data$DR1TALCO == 0] <- 0.001

## calculate corresponding energy
subset_data <- subset_data %>% mutate(protein=DR1TPROT*4,
                                  sugar=DR1TSUGR*4,
                                  fiber=DR1TFIBE*2,
                                  starch=starch*4,
                                  sat_fatacid=DR1TSFAT*9,
                                  unsat_fatacid=unsat_fatacid*9,
                                  alcohol=DR1TALCO*7) %>%
  mutate(sum_cal = protein+sugar+fiber+starch+sat_fatacid+unsat_fatacid+alcohol) %>%
  mutate(protein=protein/sum_cal, sugar=sugar/sum_cal,
         fiber=fiber/sum_cal, starch=starch/sum_cal,
         sat_fatacid=sat_fatacid/sum_cal, unsat_fatacid=unsat_fatacid/sum_cal,
         alcohol=alcohol/sum_cal)

# orthornormal basis
alcohol_others <- c(rep(-sqrt(1/42), 6), sqrt(6/7))
protein_others <- c(sqrt(5/6), rep(-sqrt(1/30), 5), 0)
fat_carbo <- c(0, rep(-sqrt(2/15), 3), rep(sqrt(3/10), 2), 0)
fiber_othercarbo <- c(0, -sqrt(1/6), -sqrt(1/6), sqrt(2/3), rep(0, 3))
sugar_starch <- c(0, sqrt(1/2), -sqrt(1/2), rep(0, 4))
sat_unsat <- c(rep(0, 4), sqrt(1/2), -sqrt(1/2), 0)
effsize_mat <- cbind(alcohol_others, protein_others, fat_carbo, fiber_othercarbo, sugar_starch, sat_unsat)

log_macronutrient <- subset_data %>% select(protein, sugar, starch, fiber, sat_fatacid, unsat_fatacid,
                                            alcohol) %>% as.matrix() %>% log()

# calculate the mediators
ilr_macronutrient <- log_macronutrient %*% effsize_mat


# select the covariates
## propensity score weights
indv_weights <- subset_data$ps_weights

## ordinal survey periods (helmert coding)
subset_data$svy_year <- factor(subset_data$svy_year, order=T,
                               levels=c("2011-2012", "2013-2014", "2015-2016", "2017-2020"))
svy_period <- model.matrix(~svy_year, data=subset_data,
                           contrasts.arg = list(svy_year="contr.helmert"))
svy_period <- svy_period[, -1] %>% as.data.frame()

## centered age
centered_age <- subset_data$demo_age_years - mean(subset_data$demo_age_years)


## categorical race (white as baseline)
race <- model.matrix(~demo_race, data=subset_data)
race <- race[, -1] %>% as.data.frame()
colnames(race) <- c("Black", "Asian", "Hispanic", "Others")


## poverty
income_ratio <- subset_data$INDFMPIR
income_ratio[income_ratio == 0] <- 0.01 # replace zeros with 0.01
poverty <- log10(5/income_ratio)


## erengy intake
energy <- log(subset_data$sum_cal)
predictors <- cbind(indv_weights, svy_period, centered_age, race, poverty, energy)


## hypertension variable

hypertension <- subset_data$htn_accaha == "Yes"

# combine all the predictors, mediators and outcome
cleaned_data <- cbind(predictors, ilr_macronutrient, hypertension)
rownames(cleaned_data) <- subset_data$SEQN

write.csv(cleaned_data, "cleaned_data.csv")
saveRDS(cleaned_data, file="cleaned_data.rds")
