library(dplyr)

data <- readRDS("selected_sample_IPW.rds")
col_keys <- read.csv("nhanes_complete_key.csv")

## remove individuals who don't have any intake of protein/carbohydrate/fat
subset_data <- data %>% filter(DR1TPROT > 0 & DR1TCARB > 0 & DR1TTFAT > 0)

## impute zero mass with 0.001
subset_data$DR1TSUGR[subset_data$DR1TSUGR == 0] <- 0.001
subset_data$DR1TFIBE[subset_data$DR1TFIBE == 0] <- 0.001
subset_data$sugar <- subset_data$DR1TCARB - subset_data$DR1TFIBE
subset_data$DR1TSFAT[subset_data$DR1TSFAT == 0] <- 0.001
subset_data$unsat_fatacid <- subset_data$DR1TMFAT + subset_data$DR1TPFAT
subset_data$unsat_fatacid[subset_data$unsat_fatacid == 0] <- 0.001
# subset_data$DR1TMFAT[subset_data$DR1TMFAT == 0] <- 0.001
# subset_data$DR1TPFAT[subset_data$DR1TPFAT == 0] <- 0.001
subset_data$DR1TALCO[subset_data$DR1TALCO == 0] <- 0.001

## calculate corresponding energy
subset_data <- subset_data %>% mutate(protein=DR1TPROT*4,
                                  sugar=sugar*4,
                                  fiber=DR1TFIBE*2,
                                  sat_fatacid=DR1TSFAT*9,
                                  unsat_fatacid=unsat_fatacid*9,
                                  alcohol=DR1TALCO*7) %>%
  mutate(sum_cal = protein+sugar+fiber+sat_fatacid+unsat_fatacid+alcohol) %>%
  mutate(protein=protein/sum_cal, sugar=sugar/sum_cal,
         fiber=fiber/sum_cal,
         sat_fatacid=sat_fatacid/sum_cal, unsat_fatacid=unsat_fatacid/sum_cal,
         alcohol=alcohol/sum_cal)

# orthornormal basis
alcohol_others <- c(rep(-sqrt(1/30), 5), sqrt(5/6))
protein_others <- c(sqrt(4/5), rep(-sqrt(1/20), 4), 0)
fat_carbo <- c(0, rep(-1/2, 2), rep(1/2, 2), 0)
sugar_fiber <- c(0, sqrt(1/2), -sqrt(1/2), rep(0, 3))
sat_unsat <- c(rep(0, 3), sqrt(1/2), -sqrt(1/2), 0)
effsize_mat <- cbind(alcohol_others, protein_others, fat_carbo, sugar_fiber, sat_unsat)

log_macronutrient <- subset_data %>% select(protein, sugar, fiber, sat_fatacid, unsat_fatacid,
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
svy_period$svy_year1 <- svy_period$svy_year1 / 2
svy_period$svy_year2 <- svy_period$svy_year2 / 3
svy_period$svy_year3 <- svy_period$svy_year3 / 4

## centered age
centered_age <- subset_data$demo_age_years - mean(subset_data$demo_age_years)


## categorical race (white as baseline)
race <- model.matrix(~demo_race, data=subset_data)
race <- race[, -1] %>% as.data.frame()
colnames(race) <- c("Black", "Asian", "Hispanic", "Others")

## gender
gender <- subset_data$demo_gender == "Women"

## poverty
income_ratio <- subset_data$INDFMPIR
income_ratio[income_ratio == 0] <- 0.01 # replace zeros with 0.01
poverty <- log10(5/income_ratio)


## erengy intake
energy <- log(subset_data$sum_cal)
predictors <- cbind(indv_weights, svy_period, centered_age, race, gender, poverty, energy)


## hypertension variable
hypertension <- subset_data$htn_jnc7 == "Yes"

# combine all the predictors, mediators and outcome
cleaned_data <- cbind(predictors, ilr_macronutrient, hypertension)
rownames(cleaned_data) <- subset_data$SEQN

write.csv(cleaned_data, "cleaned_data.csv")
saveRDS(cleaned_data, file="cleaned_data.rds")



# mediation analysis
## regressions on the mediators
med1 <- function(mydata){
  # subset_mydata <- mydata[indices, ]
  # regression on alcohol
  lm_alcohol <- lm(alcohol_others ~ svy_year1 + svy_year2 + svy_year3 + centered_age +
                     Black + Asian + Hispanic + Others + gender+ poverty + energy, data=mydata,
                   weights=mydata$indv_weights)
  a_alcohol <- lm_alcohol$coefficients[-1]
  
  # regression on protein
  lm_protein <- lm(protein_others ~ svy_year1 + svy_year2 + svy_year3 + centered_age +
                     Black + Asian + Hispanic + Others + gender+ poverty + energy, data=mydata,
                   weights=mydata$indv_weights)
  a_protein <- lm_protein$coefficients[-1]
  
  # regression on fat vs carbohydrate
  lm_fat_carbo <- lm(fat_carbo ~ svy_year1 + svy_year2 + svy_year3 + centered_age +
                       Black + Asian + Hispanic + Others + gender+ poverty + energy, data=mydata,
                     weights=mydata$indv_weights)
  a_fat_carbo <- lm_fat_carbo$coefficients[-1]
  
  # regression on sugar vs fiber
  lm_sugar_fiber <- lm(sugar_fiber ~ svy_year1 + svy_year2 + svy_year3 + centered_age +
                         Black + Asian + Hispanic + Others + gender+ poverty + energy, data=mydata,
                       weights=mydata$indv_weights)
  a_sugar_fiber <-  lm_sugar_fiber$coefficients[-1]
  
  # regression on saturated vs unsaturated fat
  lm_sat_unsat <- lm(sat_unsat ~ svy_year1 + svy_year2 + svy_year3 + centered_age +
                       Black + Asian + Hispanic + Others + gender+ poverty + energy, data=mydata,
                     weights=mydata$indv_weights)
  a_sat_unsat <-  lm_sat_unsat$coefficients[-1]
  
  a_mat <- cbind(a_alcohol, a_protein, a_fat_carbo, a_sugar_fiber, a_sat_unsat)
  
  return(a_mat)
}

# utility function for bootstrap
boot_med1 <- function(mydata, i){
  subset_mydata <- mydata[i, ]
  return (med1(subset_mydata))
}


## regressions on the hypertension outcome
med2 <- function(mydata, direct=F){
  # subset_mydata <- mydata[indices, ]
  glm_hypertension <- glm(hypertension ~ svy_year1 + svy_year2 + svy_year3 + centered_age +
                            Black + Asian + Hispanic + Others + gender+ poverty + energy +
                            alcohol_others + protein_others + fat_carbo + sugar_fiber + sat_unsat, 
                          data=mydata,
                          weights=mydata$indv_weights, family="quasibinomial")
  c_hypertension <- glm_hypertension$coefficients[seq(2, 12)]
  b_hypertension <- glm_hypertension$coefficients[seq(13, 17)]
  if (direct){
    return(c_hypertension)
  } else{
    return(b_hypertension)
  }
}

# utility function for bootstrap
boot_direct <- function(mydata, i){
  subset_mydata <- mydata[i, ]
  return (med2(subset_mydata, direct=T))
}

# utility function for bootstrap
boot_med2 <- function(mydata, i){
  subset_mydata <- mydata[i, ]
  return (med2(subset_mydata, direct=F))
}


boot_indirect <- function(mydata, i){
  
  subset_mydata <- mydata[i, ]
  a_coefmat <- med1(subset_mydata)
  b_coef <- med2(subset_mydata, direct=F)
  indirect_eff_mat <- a_coefmat * t(replicate(nrow(a_coefmat), b_coef))
  return(indirect_eff_mat)
  
}


library(boot)
## bootstrap linear regressions on the nutrients
med1_bootstrap <- boot(cleaned_data, boot_med1, R=5000,
                       parallel='multicore')
med1_estimate <- med1_bootstrap$t0
med1_lower <- apply(med1_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.025)) %>%
  matrix(nrow=nrow(med1_estimate))
med1_upper <- apply(med1_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.975)) %>%
  matrix(nrow=nrow(med1_estimate))
med1_result <- (med1_lower > 0) * (med1_upper > 0) - (med1_lower < 0) * (med1_upper < 0)
rownames(med1_result) <- rownames(med1_estimate)
colnames(med1_result) <- colnames(med1_estimate)

## bootstrap on the direct effects between demographics and hypertension
direct_bootstrap <- boot(cleaned_data, boot_direct, R=5000,
                       parallel='multicore', ncpus=4)
direct_estimate <- direct_bootstrap$t0
direct_lower <- apply(direct_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.025))
direct_upper <- apply(direct_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.975))
direct_result <- (direct_lower > 0) * (direct_upper > 0) - 
  (direct_lower < 0) * (direct_upper < 0)
names(direct_result) <- names(direct_estimate)


## bootstrap logistic regression between nutrients and hypertension
med2_bootstrap <- boot(cleaned_data, boot_med2, R=5000,
                       parallel='multicore', ncpus=4)
med2_estimate <- med2_bootstrap$t0
med2_lower <- apply(med2_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.025))
med2_upper <- apply(med2_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.975))
med2_result <- (med2_lower > 0) * (med2_upper > 0) - (med2_lower < 0) * (med2_upper < 0)
names(med2_result) <- names(med2_estimate)



## indirect effect from demographics to nutrient to hypertension
indirect_bootstrap <- boot(cleaned_data, boot_indirect, R=5000,
                           parallel='multicore', ncpus=4)
indirect_estimate <- indirect_bootstrap$t0
indirect_lower <- apply(indirect_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.025)) %>%
  matrix(nrow=nrow(indirect_estimate))
indirect_upper <- apply(indirect_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.975)) %>%
  matrix(nrow=nrow(indirect_estimate))
indirect_result <- (indirect_lower > 0) * (indirect_upper > 0) - 
  (indirect_lower < 0) * (indirect_upper < 0)
rownames(indirect_result) <- rownames(indirect_estimate)
colnames(indirect_result) <- colnames(indirect_estimate)


