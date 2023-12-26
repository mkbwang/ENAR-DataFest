load("nhanes_data.rda")
load("nhanes_key.rda")

# change ID column name for merging with other datasets
colnames(nhanes_data)[1] <- "SEQN" 
svy_ids <- nhanes_data$SEQN

# read 2011-2012 data
library(haven)
library(dplyr)
demo_11 <- read_xpt(file.path("other_data", "nhanes11", "DEMO_G.XPT"))
## include extra demographics variable:
## country of birth, time staying in US, education level, marital status, ratio of family income to poverty 
demo_11_subset <- demo_11 %>% select(SEQN, DMDBORN4, DMDYRSUS, DMDEDUC2, 
                                     DMDMARTL, INDFMPIR)
diet_11 <- read_xpt(file.path("other_data", "nhanes11", "DR1TOT_G.XPT"))



# read 2013-2014 data
demo_13 <- read_xpt(file.path("other_data", "nhanes13", "DEMO_H.XPT"))
demo_13_subset <- demo_13 %>% select(SEQN, DMDBORN4, DMDYRSUS, DMDEDUC2, 
                                     DMDMARTL, INDFMPIR)
diet_13 <- read_xpt(file.path("other_data", "nhanes13", "DR1TOT_H.XPT"))



# read 2015-2016 data
demo_15 <- read_xpt(file.path("other_data", "nhanes15", "DEMO_I.XPT"))
demo_15_subset <- demo_15 %>% select(SEQN, DMDBORN4, DMDYRSUS, DMDEDUC2, 
                                     DMDMARTL, INDFMPIR)
diet_15 <- read_xpt(file.path("other_data", "nhanes15", "DR1TOT_I.XPT"))



# read 2017-2020 data
demo_17 <- read_xpt(file.path("other_data", "nhanes17", "P_DEMO.XPT"))
demo_17_subset <- demo_17 %>% select(SEQN, DMDBORN4, DMDYRUSZ, DMDEDUC2, 
                                     DMDMARTZ, INDFMPIR)
## some column names in 2017 data are slightly different from previous years
demo_17_subset <- demo_17_subset %>% 
  rename(DMDYRSUS=DMDYRUSZ,
         DMDMARTL=DMDMARTZ)
diet_17 <- read_xpt(file.path("other_data", "nhanes17", "P_DR1TOT.XPT"))
diet_17 <- diet_17 %>% rename(WTDRD1=WTDRD1PP, 
                              WTDR2D=WTDR2DPP,
                              DR1TWS=DR1TWSZ)

# intersection of diet columns
diet_columns <- Reduce(intersect, list(colnames(diet_11), colnames(diet_13),
                                       colnames(diet_15), colnames(diet_17)))

diet_filter <- grepl("DRQ", diet_columns) | grepl("DR1", diet_columns)
diet_columns <- c("SEQN", diet_columns[diet_filter])

diet_11 <- diet_11 %>% select(all_of(diet_columns))
diet_13 <- diet_13 %>% select(all_of(diet_columns))
diet_15 <- diet_15 %>% select(all_of(diet_columns))
diet_17 <- diet_17 %>% select(all_of(diet_columns))

extradata_11 <- demo_11_subset %>% inner_join(diet_11, by="SEQN")
extradata_13 <- demo_13_subset %>% inner_join(diet_13, by="SEQN")
extradata_15 <- demo_15_subset %>% inner_join(diet_15, by="SEQN")
extradata_17 <- demo_17_subset %>% inner_join(diet_17, by="SEQN")


extradata <- rbind(extradata_11, extradata_13, extradata_15, extradata_17)

# merge the available clinical data and the retrieved diet data
complete_data <- nhanes_data %>% inner_join(extradata, by="SEQN")
all_columns <- colnames(complete_data)
all_attributes <- lapply(complete_data, function(cdata) attr(cdata, "label"))
all_attributes$demo_race_black <- "Whether race is black"
coldata <- data.frame(Column=all_columns,
                      Label=unlist(all_attributes))


saveRDS(list(data=complete_data, coldata=coldata),
        file="nhanes_complete_data.rds")
write.csv(complete_data, 'nhanes_complete_data.csv',
          row.names=F)
write.csv(coldata, 'nhanes_complete_key.csv',
          row.names=F)

