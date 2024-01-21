# ENAR DataFest 2024
The codes in this repository are owned by Mukai Wang (University of Michigan), Bulun Te (University of Michigan) and Suene Paloma Lima (Ohio University). They contain the data cleaning, exploratory data analysis and mediation analysis of [NHANES data](https://www.cdc.gov/nchs/nhanes/index.htm) for the [ENAR 2024 DataFest](https://www.enar.org/meetings/spring2024/program/datafest_submission.cfm). The project is titled "Mediation Analysis of Macronutrient Composition and Hypertension  Among US Adults".


* `nutrient_matching.R:` Retrieve the diet data from survey responses between 2011 and 2020 in [NHANES](https://wwwn.cdc.gov/nchs/nhanes/) and merge them with the hypertension diagnoses in [cardioStatsUSA](https://github.com/jhs-hwg/cardioStatsUSA) based on SEQN number.
* `subsample.R`: Filter out survey respondents who were on special medication or diet. Calculate weights for individual responses based on [covariate balancing propensity scores](https://imai.fas.harvard.edu/research/files/CBPS.pdf)  to balance the covariates between different survey periods (2011-2012, 2013-2014, 2015-2016, 2017-2020).
* `macronutrient_mediation.R`: Mediation analysis with macronutrients as mediators. The mediators are isometric log transformation of six macronutrients (sugar, fiber, saturated fat, unsaturated fat, protein, alcohol).
* `Exploratory Data Analysis/`: Codes for exploratory data analysis
* `Report.pdf`: Submitted report PDF.
