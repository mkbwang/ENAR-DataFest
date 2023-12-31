## ----setup, include=FALSE---------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----load and filter data, message=F----------------------------------------------
library(dplyr)
data <- readRDS("selected_sample_IPW.rds")
col_keys <- read.csv("nhanes_complete_key.csv")
subset_data <- data %>% filter(DR1TPROT > 0 & DR1TCARB > 0 & DR1TTFAT > 0)


## ---------------------------------------------------------------------------------
# fiber
subset_data$DR1TFIBE[subset_data$DR1TFIBE == 0] <- 0.001 
# sugar is carbohydrate - fiber
subset_data$sugar <- subset_data$DR1TCARB - subset_data$DR1TFIBE
# saturated fat
subset_data$DR1TSFAT[subset_data$DR1TSFAT == 0] <- 0.001
# unsaturated fat
subset_data$unsat_fatacid <- subset_data$DR1TMFAT + subset_data$DR1TPFAT
subset_data$unsat_fatacid[subset_data$unsat_fatacid == 0] <- 0.001
# alcohol
subset_data$DR1TALCO[subset_data$DR1TALCO == 0] <- 0.001


## ----calculate composition--------------------------------------------------------
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


## ----orthonormal bases------------------------------------------------------------
alcohol_others <- c(rep(-sqrt(1/30), 5), sqrt(5/6))
protein_others <- c(sqrt(4/5), rep(-sqrt(1/20), 4), 0)
fat_carbo <- c(0, rep(-1/2, 2), rep(1/2, 2), 0)
sugar_fiber <- c(0, sqrt(1/2), -sqrt(1/2), rep(0, 3))
sat_unsat <- c(rep(0, 3), sqrt(1/2), -sqrt(1/2), 0)
effsize_mat <- cbind(alcohol_others, protein_others, fat_carbo, sugar_fiber, sat_unsat)
t_effsize_mat <- t(effsize_mat)
rownames(t_effsize_mat) <- c("Alcohol VS Others", "Protein VS Fat and Carbohydrates", 
                             "Fat VS Carbohydrates", "Sugar VS Fiber", "Sat Fat VS Unsat Fat")
colnames(t_effsize_mat) <- c("Protein", "Sugar", "Fiber", 
                             "Sat Fat", "Unsat Fat", "Alcohol")
knitr::kable(as.data.frame(t_effsize_mat), digits=4)


## ----ilr transformation-----------------------------------------------------------

log_macronutrient <- subset_data %>% select(protein, sugar, fiber, sat_fatacid, unsat_fatacid,
                                            alcohol) %>% as.matrix() %>% log()

ilr_macronutrient <- log_macronutrient %*% effsize_mat



## ----survey period----------------------------------------------------------------
subset_data$svy_year <- factor(subset_data$svy_year, order=T,
                               levels=c("2011-2012", "2013-2014", "2015-2016", "2017-2020"))
svy_period <- model.matrix(~svy_year, data=subset_data,
                           contrasts.arg = list(svy_year="contr.helmert"))
svy_period <- svy_period[, -1] %>% as.data.frame()
svy_period$svy_year1 <- svy_period$svy_year1 / 2
svy_period$svy_year2 <- svy_period$svy_year2 / 3
svy_period$svy_year3 <- svy_period$svy_year3 / 4


## ----illustrate helmert coding, echo=F--------------------------------------------
svy1 <- c(-1/2, 1/2, 0,0)
svy2 <- c(-1/3, -1/3, 2/3, 0)
svy3 <- c(-1/4, -1/4, -1/4, 3/4)
svymat <- cbind(svy1, svy2, svy3)
rownames(svymat) <- c("2011-2012", "2013-2014", "2015-2016", "2017-2020")
colnames(svymat) <- c("svy_year1", "svy_year2", "svy_year3")
knitr::kable(as.data.frame(svymat), digits=4)


## ----center age-------------------------------------------------------------------
centered_age <- subset_data$demo_age_years - mean(subset_data$demo_age_years)


## ----race-------------------------------------------------------------------------
race <- model.matrix(~demo_race, data=subset_data)
race <- race[, -1] %>% as.data.frame()
colnames(race) <- c("Black", "Asian", "Hispanic", "Others")


## ----gender-----------------------------------------------------------------------
gender <- subset_data$demo_gender == "Women"


## ----poverty----------------------------------------------------------------------
income_ratio <- subset_data$INDFMPIR
income_ratio[income_ratio == 0] <- 0.01 # replace zeros with 0.01
poverty <- log10(5/income_ratio)


## ----energy-----------------------------------------------------------------------
energy <- log(subset_data$sum_cal)


## ---------------------------------------------------------------------------------
indv_weights <- subset_data$ps_weights
predictors <- cbind(indv_weights, svy_period, centered_age, 
                    race, gender, poverty, energy)


## ----hypertension output----------------------------------------------------------
hypertension <- subset_data$htn_jnc7 == "Yes"

# combine all the predictors, mediators and outcome
cleaned_data <- cbind(predictors, ilr_macronutrient, hypertension)
rownames(cleaned_data) <- subset_data$SEQN

# save the data
write.csv(cleaned_data, "cleaned_data.csv")
saveRDS(cleaned_data, file="cleaned_data.rds")


## ---- out.width = "400px", echo=F-------------------------------------------------
knitr::include_graphics("mediation/Mediation_Analysis.jpg")


## ----med1-------------------------------------------------------------------------
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
                       Black + Asian + Hispanic + Others + gender+ poverty + energy,
                     data=mydata,
                     weights=mydata$indv_weights)
  a_fat_carbo <- lm_fat_carbo$coefficients[-1]
  
  # regression on sugar vs fiber
  lm_sugar_fiber <- lm(sugar_fiber ~ svy_year1 + svy_year2 + svy_year3 + centered_age +
                         Black + Asian + Hispanic + Others + gender+ poverty + energy,
                       data=mydata,
                       weights=mydata$indv_weights)
  a_sugar_fiber <-  lm_sugar_fiber$coefficients[-1]
  
  # regression on saturated vs unsaturated fat
  lm_sat_unsat <- lm(sat_unsat ~ svy_year1 + svy_year2 + svy_year3 + centered_age +
                       Black + Asian + Hispanic + Others + gender+ poverty + energy,
                     data=mydata,
                     weights=mydata$indv_weights)
  a_sat_unsat <-  lm_sat_unsat$coefficients[-1]
  
  a_mat <- cbind(a_alcohol, a_protein, a_fat_carbo, a_sugar_fiber, a_sat_unsat)
  
  return(a_mat)
}


## ----utility med1-----------------------------------------------------------------
# utility function for bootstrap
boot_med1 <- function(mydata, i){
  subset_mydata <- mydata[i, ]
  return (med1(subset_mydata))
}


## ----direct effect----------------------------------------------------------------
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


## ---------------------------------------------------------------------------------
# utility function for bootstrap
boot_direct <- function(mydata, i){
  subset_mydata <- mydata[i, ]
  return (med2(subset_mydata, direct=T))
}


## ---------------------------------------------------------------------------------
boot_med2 <- function(mydata, i){
  subset_mydata <- mydata[i, ]
  return (med2(subset_mydata, direct=F))
}


## ---------------------------------------------------------------------------------
boot_indirect <- function(mydata, i){
  subset_mydata <- mydata[i, ]
  a_coefmat <- med1(subset_mydata)
  b_coef <- med2(subset_mydata, direct=F)
  indirect_eff_mat <- a_coefmat * t(replicate(nrow(a_coefmat), b_coef))
  return(indirect_eff_mat)
}


## ----fit linear regression--------------------------------------------------------
library(boot)
## bootstrap linear regressions on the nutrients
if (!file.exists("mediation/demo_nutrient.rds")){
  set.seed(2023)
  begin <- proc.time()
  med1_bootstrap <- boot(cleaned_data, boot_med1, R=5000,
                         parallel='multicore')
  duration <- proc.time() - begin
  print(duration)
  saveRDS(med1_bootstrap, "mediation/demo_nutrient.rds")
} else{
  med1_bootstrap <- readRDS("mediation/demo_nutrient.rds")
}


## ----linear regression inference, message=F---------------------------------------
med1_estimate <- med1_bootstrap$t0
med1_lower <- apply(med1_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.025)) %>%
  matrix(nrow=nrow(med1_estimate))
med1_upper <- apply(med1_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.975)) %>%
  matrix(nrow=nrow(med1_estimate))
med1_result <- (med1_lower > 0) * (med1_upper > 0) - (med1_lower < 0) * (med1_upper < 0)
rownames(med1_result) <- rownames(med1_estimate)
rownames(med1_result)[9] <- "genderFemale"
colnames(med1_result) <- c("Alcohol_Others", "Protein_Others", "Fat_Carbo", 
                           "Sugar_Fiber", "Sat_Unsat")


## ----visualize linear regression, message=F, fig.height=3-------------------------
library(tidyr)
med1_result <- as.data.frame(med1_result)
med1_result$X <- rownames(med1_result)
med1_result_long <- pivot_longer(med1_result, -c(X), values_to="Inference", names_to = "Mediator")
med1_result_long$X <- factor(med1_result_long$X, 
                             levels=c("svy_year1", "svy_year2", "svy_year3", 
                                      "Black", "Asian", "Hispanic", "Others",
                                      "genderFemale", "centered_age", "poverty", "energy"))
med1_result_long$Inference <- as.factor(med1_result_long$Inference)
library(ggplot2)
ggplot(med1_result_long, aes(Mediator, X)) + geom_tile(aes(fill=Inference), color="#222222") +
  scale_fill_manual(values=c("#3366ff", "#FFFFFF", "#cc0000")) + theme_bw() + 
  theme(text = element_text(size=12)) + theme(legend.position = "none")


## ----direct effect inference------------------------------------------------------
if (!file.exists("mediation/demo_hypertension_direct.rds")){
  set.seed(2023)
  begin <- proc.time()
  direct_bootstrap <- boot(cleaned_data, boot_direct, R=5000,
                           parallel='multicore', ncpus=4)
  duration <- proc.time() - begin
  print(duration)
  saveRDS(direct_bootstrap, "mediation/demo_hypertension_direct.rds")
} else{
  direct_bootstrap <- readRDS("mediation/demo_hypertension_direct.rds")
}
direct_estimate <- direct_bootstrap$t0
direct_lower <- apply(direct_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.025))
direct_upper <- apply(direct_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.975))
direct_result <- (direct_lower > 0) * (direct_upper > 0) - 
  (direct_lower < 0) * (direct_upper < 0)
names(direct_result) <- names(direct_estimate)
direct_inference <- cbind(direct_estimate, direct_lower, direct_upper, direct_result)
colnames(direct_inference) <- c("Estimate", "95% Lower", "95% Upper", "Significance")
rownames(direct_inference)[9] <- "Female"
knitr::kable(as.data.frame(direct_inference), digits=4)


## ----nutrient on hypertension-----------------------------------------------------
if (!file.exists("mediation/nutrient_hypertension.rds")){
  set.seed(2023)
  begin <- proc.time()
  med2_bootstrap <- boot(cleaned_data, boot_med2, R=5000,
                         parallel='multicore', ncpus=4)
  duration <- proc.time() - begin
  print(duration)
  saveRDS(med2_bootstrap, "mediation/nutrient_hypertension.rds")
} else{
  med2_bootstrap <- readRDS("mediation/nutrient_hypertension.rds")
}
med2_estimate <- med2_bootstrap$t0
med2_lower <- apply(med2_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.025))
med2_upper <- apply(med2_bootstrap$t, 2, 
                    function(values) quantile(values, probs=0.975))
med2_result <- (med2_lower > 0) * (med2_upper > 0) - (med2_lower < 0) * (med2_upper < 0)
names(med2_result) <- names(med2_estimate)

med2_inference <- cbind(med2_estimate, med2_lower, med2_upper, med2_result)
colnames(med2_inference) <- c("Estimate", "95% Lower", "95% Upper", "Significance")
knitr::kable(as.data.frame(med2_inference), digits=4)


## ---------------------------------------------------------------------------------
if (!file.exists("mediation/demo_nutrient_indirect.rds")){
  set.seed(2023)
  begin <- proc.time()
  indirect_bootstrap <- boot(cleaned_data, boot_indirect, R=5000,
                             parallel='multicore', ncpus=4)
  duration <- proc.time() - begin
  saveRDS(indirect_bootstrap, "mediation/demo_nutrient_indirect.rds")
} else{
  indirect_bootstrap <- readRDS("mediation/demo_nutrient_indirect.rds")
}
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
rownames(indirect_result)[9] <- "genderFemale"
colnames(indirect_result) <- c("Alcohol_Others", "Protein_Others", "Fat_Carbo", 
                           "Sugar_Fiber", "Sat_Unsat")


## ----visualize indirect effect, message=F, fig.height=3---------------------------
library(tidyr)
indirect_result <- as.data.frame(indirect_result)
indirect_result$X <- rownames(indirect_result)
indirect_result_long <- pivot_longer(as.data.frame(indirect_result),
                                     -c(X), values_to="Inference", names_to = "Mediator")
indirect_result_long$X <- factor(indirect_result_long$X, 
                             levels=c("svy_year1", "svy_year2", "svy_year3", 
                                      "Black", "Asian", "Hispanic", "Others",
                                      "genderFemale", "centered_age", "poverty", "energy"))
indirect_result_long$Inference <- as.factor(indirect_result_long$Inference)
library(ggplot2)
ggplot(indirect_result_long, aes(Mediator, X)) + 
  geom_tile(aes(fill=Inference), color="#222222") +
  scale_fill_manual(values=c("#3366ff", "#FFFFFF", "#cc0000")) + theme_bw() + 
  theme(text = element_text(size=12)) + theme(legend.position = "none")

