# R script to analyze powdery mildew resistance in male hop lines
# published in Journal of Plant Registrations (DOI: https://doi.org/10.1002/plr2.70035)

# load packages
#install.packages("tidyverse")
library("tidyverse")
#install.packages("dplyr")
library("dplyr")
#install.packages("readxl")
library("readxl")
#install.packages("Skillings.Mack")
library(Skillings.Mack)
#install.packages("dunn.test")
library("dunn.test")
#install.packages("FSA")
library("FSA")
#install.packages("rcompanion")
library("rcompanion")
#install.packages("coin")
library("coin")
#install.packages("ggh4x")
library("ggh4x")

# note on phenotype data: edited to change typo in one male genotype identity
# 076-051M -> 0716_051M

#### Whole Plant #### 
# read in data
Mildew_WholePlant<- read_excel("/Users/kaylaaltendorf/Documents/PM Resistant Males/Data/2021 male germplasm PM screening Data.xlsx",
                       sheet = "Whole Plant Data",
                       na=".") %>%
  dplyr::select(Identifier, Mildew, Run, Bench, Node, Leaf,
         `Ordinal Rating Bine 1/colony counts`, `Ordinal Rating Bine 2/colony counts`, `Reassessment Bine 1`, `Reassessment Bine 2`) %>%
  rename(FirstRating_Bine1 = `Ordinal Rating Bine 1/colony counts`,
         FirstRating_Bine2 = `Ordinal Rating Bine 2/colony counts`,
         SecondRating_Bine1 = `Reassessment Bine 1`,
         SecondRating_Bine2 = `Reassessment Bine 2`) %>%
  rowwise() %>% 
  mutate(MedianBineRating = median(c(FirstRating_Bine1, FirstRating_Bine2, SecondRating_Bine1), na.rm = T),
         Block = paste(Bench, sep = "_"),
         Node = as.factor(Node),
         Leaf = as.factor(Leaf),
         Identifier = sub('\\-', '_', Identifier)) # find the median value across bines for each run/bench/node/leaf

unique(Mildew_WholePlant$Identifier)

#### Determining Resistant or Susceptible for Whole Plants ####
disease_summary <- Mildew_WholePlant %>% 
  group_by(Identifier, Mildew) %>% 
  summarise(median_disease_score = median(FirstRating_Bine1, FirstRating_Bine2, SecondRating_Bine1, na.rm = T), 
            max_disease_score = max(FirstRating_Bine1, FirstRating_Bine2, SecondRating_Bine1, na.rm = T)) 

disease_summary$susceptible <- NA

# separate out by isolate
disease_summary_pre2012 <- disease_summary %>% filter(Mildew == "Pre-2012")
disease_summary_v6 <- disease_summary %>% filter(Mildew == "V6")
disease_summary_cascadeadapted <- disease_summary %>% filter(Mildew == "Cascade-adapted")

# pre-2012
# for the pre-2012 mildew, if it had a median score >2, susceptible, or any leaf over >=3 
for (i in 1:nrow(disease_summary_pre2012)) {
  if (disease_summary_pre2012$median_disease_score[i] >= 2) {disease_summary_pre2012$susceptible[i] = "yes"}
  if (disease_summary_pre2012$max_disease_score[i] >=3) {disease_summary_pre2012$susceptible[i] = "yes"}
}

# v6 - same criteria as above
for (i in 1:nrow(disease_summary_v6)) {
  if (disease_summary_v6$median_disease_score[i] >= 2) {disease_summary_v6$susceptible[i] = "yes"}
  if (disease_summary_v6$max_disease_score[i] >=3) {disease_summary_v6$susceptible[i] = "yes"}
}

# cascade adapted - anything over 0 is susceptible
for (i in 1:nrow(disease_summary_cascadeadapted)) {
  if (disease_summary_cascadeadapted$max_disease_score[i] != 0) {disease_summary_cascadeadapted$susceptible[i] = "yes"}
}

# rbind them together
disease_summary <- rbind(disease_summary_pre2012, disease_summary_v6, disease_summary_cascadeadapted)

# change IDs
# if it starts with 02, 03, 07, 11, 86, 03, then add a W beforehand
for (i in 1:nrow(disease_summary)) {
  if(as.character(substr(disease_summary$Identifier[i], 1, 2)) %in% sprintf("%02d", c(02, 03, 07, 11, 86, 03))) {
    disease_summary$Identifier[i] <- paste("W", disease_summary$Identifier[i], sep = "")
  }
}


#### Pre-2012 Isolate ####

Mildew_WholePlant_Pre2012 <- Mildew_WholePlant %>%
  dplyr::select(Identifier, Mildew, Block, Node, Leaf,
         MedianBineRating) %>%
  filter(Mildew == "Pre-2012") %>%
  group_by(Identifier, Mildew, Block) %>%
  summarise(BlockMedian = median(MedianBineRating, na.rm = T))

xtabs( ~ Identifier + Block, data=Mildew_WholePlant_Pre2012)

Ski.Mack(Mildew_WholePlant_Pre2012$BlockMedian,
         groups = Mildew_WholePlant_Pre2012$Identifier,
         blocks = Mildew_WholePlant_Pre2012$Block) # P < 0.000000

Mildew_WholePlant_Pre2012_CI = groupwiseMedian(BlockMedian ~ Identifier,
                                          data       = Mildew_WholePlant_Pre2012,
                                          conf       = 0.95,
                                          R          = 100,
                                          boot = TRUE,
                                          wilcox = FALSE,
                                          percentile = TRUE,
                                          bca        = FALSE,
                                          exact = FALSE,
                                          digits     = 2) %>%
  mutate(Plot_Median = if_else(is.na(Boot.median), # Using bootstrapped medians unless the is NO variance
                               Median, 
                               Boot.median))

ggplot(Mildew_WholePlant_Pre2012_CI, aes(x = Identifier,
                                         y = Boot.median)) +
  geom_point()+
  geom_errorbar(aes(x=Identifier, ymin=Percentile.lower, ymax=Percentile.upper))

mildew1 <- Mildew_WholePlant_Pre2012_CI %>% mutate(mildew = "HPM 663")


#### V6 Isolate #### 
Mildew_WholePlant_V6<- Mildew_WholePlant %>%
  dplyr::select(Identifier, Mildew, Block, Node, Leaf,
         MedianBineRating) %>%
  filter(Mildew == "V6") %>%
  group_by(Identifier, Mildew, Block) %>%
  summarise(BlockMedian = median(MedianBineRating, na.rm = T))

xtabs( ~ Identifier + Block, data=Mildew_WholePlant_V6)

Ski.Mack(Mildew_WholePlant_V6$BlockMedian, # test for significance
         groups = Mildew_WholePlant_V6$Identifier,
         blocks = Mildew_WholePlant_V6$Block) # P = 0.012 

Mildew_WholePlant_V6_Varience<- Mildew_WholePlant_V6 %>% # filtering entry NAMES that have some variance in measures
  group_by(Identifier) %>%
  summarise(IdentifierMedian = var(BlockMedian, na.rm = T)) %>%
  filter(IdentifierMedian != 0)

Mildew_WholePlant_V6_NoVarience<- Mildew_WholePlant_V6 %>% # filtering those NAMES with NO variance
  group_by(Identifier) %>%
  summarise(IdentifierMedian = var(BlockMedian, na.rm = T)) %>%
  filter(IdentifierMedian == 0)

Mildew_WholePlant_V6_VarianceData<- Mildew_WholePlant_V6_Varience %>% # adding the data of those with variance
  dplyr::select(Identifier) %>%
  left_join(Mildew_WholePlant_V6, by = "Identifier")

Mildew_WholePlant_V6_VarianceData_CI = groupwiseMedian(BlockMedian ~ Identifier,    # calculating CI
                                                      data = Mildew_WholePlant_V6_VarianceData,
                                          conf       = 0.95,
                                          R          = 100,
                                          boot = TRUE,
                                          wilcox = FALSE,
                                          percentile = TRUE,
                                          bca        = FALSE,
                                          exact = FALSE,
                                          digits     = 2)

Mildew_WholePlant_V6_NoVarienceData<- Mildew_WholePlant_V6_NoVarience %>% # adding the data of those with NO variance
  dplyr::select(Identifier) %>%
  left_join(Mildew_WholePlant_V6, by = "Identifier") %>%
  group_by(Identifier) %>%
  summarise(Median = median(BlockMedian, na.rm = T)) %>%
  mutate(n = rep(NA),
         Median = Median,
         Boot.median = rep(NA),
         Conf.level = rep(0.95),
         Percentile.lower = rep(NA),
         Percentile.upper = rep(NA)) 
  
Mildew_WholePlant_V6_CombinedData<- bind_rows(Mildew_WholePlant_V6_VarianceData_CI,
                                              Mildew_WholePlant_V6_NoVarienceData) %>%
  mutate(Plot_Median = if_else(is.na(Boot.median), # Using bootstrapped medians unless the is NO varience
                               Median, 
                               Boot.median))
  
ggplot(Mildew_WholePlant_V6_CombinedData, aes(x = Identifier,
                                         y = Plot_Median)) +
  geom_point()+
  geom_errorbar(aes(x=Identifier, ymin=Percentile.lower, ymax=Percentile.upper))

mildew2 <- Mildew_WholePlant_V6_CombinedData %>% mutate(mildew = "HPM 609")


#### Cascade Adapted Isolate #### 
Mildew_WholePlant_CascadeAdapted<- Mildew_WholePlant %>%
  dplyr::select(Identifier, Mildew, Block, Node, Leaf,
         MedianBineRating) %>%
  filter(Mildew == "Cascade-adapted") %>%
  group_by(Identifier, Mildew, Block) %>%
  summarise(BlockMedian = median(MedianBineRating, na.rm = T))

xtabs( ~ Identifier + Block, data=Mildew_WholePlant_CascadeAdapted)

Ski.Mack(Mildew_WholePlant_CascadeAdapted$BlockMedian,
         groups = Mildew_WholePlant_CascadeAdapted$Identifier,
         blocks = Mildew_WholePlant_CascadeAdapted$Block) # P < 0.001 

Mildew_WholePlant_CascadeAdapted_CI = groupwiseMedian(BlockMedian ~ Identifier,
                                            data       = Mildew_WholePlant_CascadeAdapted,
                                            conf       = 0.95,
                                            R          = 1000,
                                            boot = TRUE,
                                            wilcox = FALSE,
                                            percentile = TRUE,
                                            bca        = FALSE,
                                            digits     = 2) %>%
  mutate(Plot_Median = if_else(is.na(Boot.median), # Using bootstrapped medians unless the is NO varience
                               Median, 
                               Boot.median))

ggplot(Mildew_WholePlant_CascadeAdapted_CI, aes(x = Identifier,
                                              y = Plot_Median)) +
  geom_point()+
  geom_errorbar(aes(x=Identifier, ymin=Percentile.lower, ymax=Percentile.upper))

mildew3 <- Mildew_WholePlant_CascadeAdapted_CI %>% mutate(mildew = "HPM 1084")


#### Leaf #### 
Mildew_LeafData<- read_excel("/Users/kaylaaltendorf/Documents/PM Resistant Males/Data/2021 male germplasm PM screening Data.xlsx",
                               sheet = "Detached Leaf Data",
                               na=".") %>%
  dplyr::select(Identifier, Mildew, Run, `Petri Dish`, Leaf,
         OrdinalFirstRating, CountFirstRating) %>%
  rowwise() %>%
  mutate(`Petri Dish` = as.factor(`Petri Dish`),
         Leaf = as.factor(Leaf),
         Identifier = sub('\\-', '_', Identifier)) #finding the median value across bines for each run/bench/node/leaf

#### Determining Resistant or Susceptible for Leaves ####
leaf_summary <- Mildew_LeafData %>% 
  group_by(Identifier, Mildew) %>% 
  summarise(max_ordinal_count = max(OrdinalFirstRating, CountFirstRating, na.rm = T))

leaf_summary$susceptible <- NA

# for the leaves, anything above zero is susceptible
leaf_summary_HPM200 <- leaf_summary %>% filter(Mildew == "HPM200")
leaf_summary_HPM204 <- leaf_summary %>% filter(Mildew == "HPM204")
leaf_summary_HPM1285 <- leaf_summary %>% filter(Mildew == "HPM1285")

# HPM200
for (i in 1:nrow(leaf_summary_HPM200)) {
  if (leaf_summary_HPM200$max_ordinal_count[i] != 0) {leaf_summary_HPM200$susceptible[i] = "yes"}
}

# HPM204
for (i in 1:nrow(leaf_summary_HPM204)) {
  if (leaf_summary_HPM204$max_ordinal_count[i] != 0) {leaf_summary_HPM204$susceptible[i] = "yes"}
}

# HPM1285
for (i in 1:nrow(leaf_summary_HPM1285)) {
  if (leaf_summary_HPM1285$max_ordinal_count[i] != 0) {leaf_summary_HPM1285$susceptible[i] = "yes"}
}

# rbind them together
leaf_summary <- rbind(leaf_summary_HPM200, leaf_summary_HPM204, leaf_summary_HPM1285)

# change IDs
# if it starts with 02, 03, 07, 11, 86, 03, then add a W beforehand
for (i in 1:nrow(leaf_summary)) {
  if(as.character(substr(leaf_summary$Identifier[i], 1, 2)) %in% sprintf("%02d", c(02, 03, 07, 11, 86, 03))) {
    leaf_summary$Identifier[i] <- paste("W", leaf_summary$Identifier[i], sep = "")
  }
}

#### HPM200 Isolate #### 
Mildew_LeafData_HPM200 <- Mildew_LeafData %>%
  dplyr::select(Identifier, Mildew, Run, `Petri Dish`, Leaf,
         CountFirstRating) %>%
  filter(Mildew == "HPM200") %>%
  group_by(Identifier, Mildew, Run, `Petri Dish`) %>%
  summarise(BlockMean = mean(CountFirstRating, na.rm = T),  # taking the MEAN of leaves per petri dish, assuming petri dishes are experimental units
            Identifier = as.factor(Identifier))

xtabs( ~ Identifier + Run + `Petri Dish`, data=Mildew_LeafData_HPM200)
summary(Mildew_LeafData_HPM200)

kruskal_test(BlockMean ~ Identifier, data = Mildew_LeafData_HPM200,
             distribution = approximate(nresample = 1000))

Mildew_LeafData_HPM200_CI = groupwiseMedian(BlockMean ~ Identifier,
                                            data       = Mildew_LeafData_HPM200,
                                            conf       = 0.95,
                                            R          = 1000,
                                            boot = TRUE,
                                            wilcox = FALSE,
                                            percentile = TRUE,
                                            bca        = FALSE,
                                            digits     = 2) %>%
  mutate(Plot_Median = if_else(is.na(Boot.median), # Using bootstrapped medians unless the is NO variance
                               Median, 
                               Boot.median))

ggplot(Mildew_LeafData_HPM200_CI, aes(x = Identifier,
                                                y = Plot_Median)) +
  geom_point()+
  geom_errorbar(aes(x=Identifier, ymin=Percentile.lower, ymax=Percentile.upper))

mildew4 <- Mildew_LeafData_HPM200_CI %>% mutate(mildew = "HPM 200")

#### HPM204 Isolate #### 
Mildew_LeafData_HPM204<- Mildew_LeafData %>%
  dplyr::select(Identifier, Mildew, Run, `Petri Dish`, Leaf,
         CountFirstRating) %>%
  filter(Mildew == "HPM204") %>%
  group_by(Identifier, Mildew, Run, `Petri Dish`) %>%
  summarise(BlockMean = mean(CountFirstRating, na.rm = T), # taking the MEAN of leaves per petri dish, assuming petri dishes are experimental units
            Identifier = as.factor(Identifier))

xtabs( ~ Identifier + Run + `Petri Dish`, data=Mildew_LeafData_HPM204)
summary(Mildew_LeafData_HPM204)

kruskal_test(BlockMean ~ Identifier, data = Mildew_LeafData_HPM204,
             distribution = approximate(nresample = 1000))

Mildew_LeafData_HPM204_CI = groupwiseMedian(BlockMean ~ Identifier,
                                             data       = Mildew_LeafData_HPM204,
                                             conf       = 0.95,
                                             R          = 1000,
                                             boot = TRUE,
                                             wilcox = FALSE,
                                             percentile = TRUE,
                                             bca        = FALSE,
                                             digits     = 2) %>%
  mutate(Plot_Median = if_else(is.na(Boot.median), # Using bootstrapped medians unless the is NO varience
                               Median, 
                               Boot.median))


ggplot(Mildew_LeafData_HPM204_CI, aes(x = Identifier,
                                      y = Plot_Median)) +
  geom_point()+
  geom_errorbar(aes(x=Identifier, ymin=Percentile.lower, ymax=Percentile.upper))


mildew5 <- Mildew_LeafData_HPM204_CI %>% mutate(mildew = "HPM 204")

#### HPM1285 Isolate #### 
Mildew_LeafData_HPM1285<- Mildew_LeafData %>%
  dplyr::select(Identifier, Mildew, `Petri Dish`, Leaf, # removed run, there seems to be some odd things going on with run
         OrdinalFirstRating) %>%
  filter(Mildew == "HPM1285") %>%
  group_by(Identifier, Mildew, `Petri Dish`) %>%
  summarise(BlockMean = mean(OrdinalFirstRating, na.rm = T), # taking the MEAN of leaves per petri dish, assuming petri dishes are experimental units
            Identifier = as.factor(Identifier))

xtabs( ~ Identifier + `Petri Dish`, data=Mildew_LeafData_HPM1285)
summary(Mildew_LeafData_HPM1285)

kruskal_test(BlockMean ~ Identifier, data = Mildew_LeafData_HPM1285,
             distribution = approximate(nresample = 1000)) # P < 0.001




Mildew_LeafData_HPM1285_CI = groupwiseMedian(BlockMean ~ Identifier,
                      data       = Mildew_LeafData_HPM1285,
                      conf       = 0.95,
                      R          = 1000,
                      boot = TRUE,
                      wilcox = FALSE,
                      percentile = TRUE,
                      bca        = FALSE,
                      digits     = 2) %>%
  mutate(Plot_Median = if_else(is.na(Boot.median), # Using bootstrapped medians unless the is NO varience
                               Median, 
                               Boot.median))

ggplot(Mildew_LeafData_HPM1285_CI, aes(x = Identifier,
                                      y = Plot_Median)) +
  geom_point()+
  geom_errorbar(aes(x=Identifier, ymin=Percentile.lower, ymax=Percentile.upper))

mildew6 <- Mildew_LeafData_HPM1285_CI %>% mutate(mildew = "HPM 1285")

### full figure
install.packages("stringr")
require('stringr')
install.packages("cowplot")
require('cowplot')
install.packages("car")
require('car')
install.packages("ggpubr")
require("ggpubr")
install.packages("glue")
require("glue")
install.packages("ggtext")
require("ggtext")


for (i in 1:nrow(mildew1)) {
  if (mildew1$Identifier[i] %in% mildew2$Identifier) {
    mildew1$Resistant[i] <- "yes"
  }
  else (mildew1$Resistant[i] <- "no")
}


for (i in 1:nrow(mildew2)) {
  if (mildew2$Identifier[i] %in% mildew3$Identifier) {
    mildew2$Resistant[i] <- "yes"
  }
  else (mildew2$Resistant[i] <- "no")
}


for (i in 1:nrow(mildew3)) {
  if (mildew3$Identifier[i] %in% mildew4$Identifier) {
    mildew3$Resistant[i] <- "yes"
  }
  else (mildew3$Resistant[i] <- "no")
}

for (i in 1:nrow(mildew4)) {
  if (mildew4$Identifier[i] %in% mildew5$Identifier) {
    mildew4$Resistant[i] <- "yes"
  }
  else (mildew4$Resistant[i] <- "no")
}

for (i in 1:nrow(mildew5)) {
  if (mildew5$Identifier[i] %in% mildew6$Identifier) {
    mildew5$Resistant[i] <- "yes"
  }
  else (mildew5$Resistant[i] <- "no")
}

for (i in 1:nrow(mildew6)) {
  if (mildew6$Median[i] == 0) {
    mildew6$Resistant[i] <- "yes"
  }
  else (mildew6$Resistant[i] <- "no")
}


# first make a backbone using the first dataset
# make sure it has all the checks
checks <- c("Cascade", "Nugget", "Zenith", "Symphony", "Target")

mildews <- data.frame(Identifier = mildew1$Identifier)
checks <- data.frame(Identifier = checks)

backbone <- rbind(mildews, checks)

mildew1 <- backbone %>% 
  left_join(mildew1 %>% 
              dplyr::select(Identifier, mildew, Plot_Median, Percentile.lower,  Percentile.upper, Resistant), by = "Identifier")

mildew2 <- backbone %>% 
  left_join(mildew2 %>% 
              dplyr::select(Identifier, mildew, Plot_Median, Percentile.lower,  Percentile.upper, Resistant), by = "Identifier")

mildew3 <- backbone %>% 
  left_join(mildew3 %>% 
              dplyr::select(Identifier, mildew, Plot_Median, Percentile.lower,  Percentile.upper, Resistant), by = "Identifier")

mildew4 <- backbone %>% 
  left_join(mildew4 %>% 
              dplyr::select(Identifier, mildew, Plot_Median, Percentile.lower,  Percentile.upper, Resistant), by = "Identifier")

mildew5 <- backbone %>% 
  left_join(mildew5 %>% 
              dplyr::select(Identifier, mildew, Plot_Median, Percentile.lower,  Percentile.upper, Resistant), by = "Identifier")

mildew6 <- backbone %>% 
  left_join(mildew6 %>% 
              dplyr::select(Identifier, mildew, Plot_Median, Percentile.lower,  Percentile.upper, Resistant), by = "Identifier")




### 
mildews <- rbind(mildew1, mildew2, mildew3, mildew4, mildew5, mildew6)
head(mildews)

mildews <- mildews %>% filter(! is.na(mildew)) 
unique(mildews$mildew)

# reorder for the figure
mildews$mildew = factor(mildews$mildew, levels=c('HPM 663','HPM 609','HPM 1084','HPM 200', 'HPM 204', 'HPM 1285'))

# if it starts with 02, 03, 07, 11, 86, 03, then add a W beforehand
for (i in 1:nrow(mildews)) {
  if(as.character(substr(mildews$Identifier[i], 1, 2)) %in% sprintf("%02d", c(02, 03, 07, 11, 86, 03))) {
    mildews$Identifier[i] <- paste("W", mildews$Identifier[i], sep = "")
  }
}

# now update resistant status
final_summary <- rbind(leaf_summary, disease_summary)
final_summary <- final_summary %>% dplyr::select(Identifier, Mildew, susceptible)

final_summary$resistant <- NA

for (i in 1:nrow(final_summary)) {
  if(is.na(final_summary$susceptible[i])) {final_summary$resistant[i] <- "yes"}
  if(! is.na(final_summary$susceptible[i])) {
    if(final_summary$susceptible[i] == "yes") {final_summary$resistant[i] <- "no"}
  }
}

# make the final summary mildew column match
mildew_merge <- data.frame(Mildew = unique(final_summary$Mildew), mildew = c("HPM 200", "HPM 204", "HPM 1285", "HPM 663", "HPM 609", "HPM 1084"))
final_summary <- final_summary %>% left_join(mildew_merge, by = "Mildew") 


#### ggplot for Genotypes ####
# write.csv(mildews, "/users/kayla.altendorf/OneDrive - USDA/Documents/2025/Male PM Germplasm Release/Median Disease Scores.csv", row.names = F)
genotypes <- mildews %>% filter(! Identifier %in% checks$Identifier)

genotypes %>% filter()

# add in virulence
virulence <- data.frame(mildew = unique(genotypes$mildew), virulence = c("pre-2012; Vb, V3, V5", "R6-virulent; Vb, V3, V4, V5, V6", "Cascade-adapted; V6, V3, V5", 
                                                            "Vb, V1, V3", "Vb, V1, V2, V3, V5", "Vb, V3, V4, V6, VWH18"))
genotypes <- genotypes %>% left_join(virulence, by = "mildew")

genotypes <- genotypes %>% mutate(merge_col = paste(Identifier, mildew, sep = "_")) %>% dplyr::select(-Resistant)
final_summary <- final_summary %>% mutate(merge_col = paste(Identifier, mildew, sep = "_")) %>% ungroup()

final_summary <- final_summary %>% dplyr::select(merge_col, resistant)

genotypes <- left_join(genotypes, final_summary, by = "merge_col")

### if resistant = "no" in previous mildew, eliminate upcoming mildew ratings, as these should not have been phenotyped but were by mistake
# separate out all mildews
genotypes_HPM663 <- genotypes %>% filter(mildew == "HPM 663")
genotypes_HPM609 <- genotypes %>% filter(mildew == "HPM 609")
genotypes_HPM1084 <- genotypes %>% filter(mildew == "HPM 1084")
genotypes_HPM200 <- genotypes %>% filter(mildew == "HPM 200")
genotypes_HPM204 <- genotypes %>% filter(mildew == "HPM 204")
genotypes_HPM1285 <- genotypes %>% filter(mildew == "HPM 1285")

sus <- genotypes_HPM663 %>% filter(resistant == "no")

# remove these from subsequent mildews
genotypes_HPM609 <- genotypes_HPM609 %>% filter(! Identifier %in% sus$Identifier)
genotypes_HPM1084 <- genotypes_HPM1084 %>% filter(! Identifier %in% sus$Identifier)

# identify the sus in HPM 200 that was sus and remove it from 204
sus <- genotypes_HPM200 %>% filter(resistant == "no") 
genotypes_HPM204 <- genotypes_HPM204 %>% filter(! Identifier %in% sus$Identifier)

genotypes <- rbind(genotypes_HPM663, genotypes_HPM609, genotypes_HPM1084, genotypes_HPM200, genotypes_HPM204, genotypes_HPM1285)

p1 <- ggplot(genotypes, aes(Plot_Median, Identifier)) + 
  geom_point()  + 
  facet_wrap(~mildew + virulence, nrow = 1, scales = "free_x") + 
  geom_pointrange(aes(xmin = Percentile.lower, xmax = Percentile.upper, color = resistant)) + 
  scale_color_manual(values = c("darkgray", "black")) + 
  #scale_color_manual(values = c("#CC79A7", "#009E73")) + 
  theme_bw() + 
  labs(x = "Bootstrapped Median", y = "Genotype", color = "Resistant*") + 
  theme(axis.title.x=element_blank()) + 
  scale_y_discrete(limits=rev) +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    strip.text.x = element_text(size = 10), 
    strip.background = element_rect(fill = "white"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 15), 
    axis.text.y = element_text(color = "black", size = 10), 
    axis.text.x = element_text(color = "black", size = 15), 
    legend.text=element_text(size=15), 
    legend.title=element_text(size = 15))

p1 <- p1 + facetted_pos_scales(x = list(
  mildew == "HPM 663" ~ scale_x_continuous(limits = c(0, 5)), 
  mildew == "HPM 609" ~ scale_x_continuous(limits = c(0, 5)), 
  mildew == "HPM 1084" ~ scale_x_continuous(limits = c(0, 20)), 
  mildew == "HPM 200" ~ scale_x_continuous(limits = c(0, 5)), 
  mildew == "HPM 204" ~ scale_x_continuous(limits = c(0, 15)),   
  mildew == "HPM 1285" ~ scale_x_continuous(limits = c(0, 5))
))
  
#### ggplot for Checks ####
# the for loop didn't work for these because they're in every mildew eval
checks <- mildews %>% filter(Identifier %in% checks$Identifier)
nug <- checks %>% filter(Identifier == "Nugget") %>% mutate(Resistant = "no")
zen <- checks %>% filter(Identifier == "Zenith") %>% mutate(Resistant = "yes")
symph <- checks %>% filter(Identifier == "Symphony") %>% mutate(Resistant = "no")
targ <- checks %>% filter(Identifier == "Target") %>% mutate(Resistant = "no")
cas <- checks %>% filter(Identifier == "Cascade") %>% mutate(Resistant = "no") %>% mutate(Identifier = "    Cascade") # add a few blank spaces to cultivar names so that the margins are the same in the figure

checks <- rbind(nug, zen, symph, cas, targ)

# edit
checks[2, 6] <- "yes" # nugget resistant to 200
checks[3, 6] <- "yes" # nugget resistant to 204
checks[4, 6] <- "no" # zenith sus to 200
checks[5, 6] <- "no" # zenith sus to 204

checks <- checks %>% left_join(virulence, by = "mildew")

p2 <- ggplot(checks, aes(Plot_Median, Identifier)) + 
  geom_point()  + 
  facet_wrap(~mildew, nrow = 1, scales = "free_x") + 
  geom_pointrange(aes(xmin = Percentile.lower, xmax = Percentile.upper, color = Resistant)) + 
  scale_color_manual(values = c("darkgray", "black")) + 
  #scale_color_manual(values = c("#CC79A7", "#009E73")) + 
  theme_bw() + 
  labs(x = "Bootstrapped Median", y = "Check") + 
  theme(axis.title.x=element_blank()) + 
  scale_y_discrete(limits=rev) +
  theme(
    legend.position = "none",
    strip.text.x = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.title.x = element_text(face = "bold", size = 15), 
    axis.title.y = element_text(face = "bold", size = 15), 
    axis.text.y = element_text(color = "black", size = 10), 
    axis.text.x = element_text(color = "black", size = 15), 
    legend.text=element_text(size=15), 
    legend.title=element_text(size = 15))

p2 <- p2 + facetted_pos_scales(x = list(
  mildew == "HPM 663" ~ scale_x_continuous(limits = c(0, 5)), 
  mildew == "HPM 609" ~ scale_x_continuous(limits = c(0, 5)), 
  mildew == "HPM 1084" ~ scale_x_continuous(limits = c(0, 55)), 
  mildew == "HPM 200" ~ scale_x_continuous(limits = c(0, 5)), 
  mildew == "HPM 204" ~ scale_x_continuous(limits = c(0, 25)),   
  mildew == "HPM 1285" ~ scale_x_continuous(limits = c(0, 5))
))

plot_grid(p1, p2, nrow = 2, ncol = 1, rel_widths = c(1,1), rel_heights = c(1,0.1))

# export as pdf 15 x 15

# summarize findings
genotypes %>% filter(mildew == "HPM 663") %>% filter(!Identifier %in% c("Cascade", "Zenith", "Symphony", "Nugget", "Target")) %>% tally()
# 102 in total
genotypes %>% filter(mildew == "HPM 663") %>% filter(!Identifier %in% c("Cascade", "Zenith", "Symphony", "Nugget", "Target")) %>% group_by(resistant) %>% tally()
50+52 
# 50 were susceptible, removed
genotypes %>% filter(mildew == "HPM 609") %>% filter(!Identifier %in% c("Cascade", "Zenith", "Symphony", "Nugget", "Target")) %>% group_by(resistant) %>% tally()
17+30 # 47 were subsequently tested, 30 were resistant 
47-17
genotypes %>% filter(mildew == "HPM 1084") %>% filter(!Identifier %in% c("Cascade", "Zenith", "Symphony", "Nugget", "Target")) %>% group_by(resistant) %>% tally()
15+14
genotypes %>% filter(mildew == "HPM 200") %>% filter(!Identifier %in% c("Cascade", "Zenith", "Symphony", "Nugget", "Target")) %>% group_by(resistant) %>% tally()
1+12
genotypes %>% filter(mildew == "HPM 204") %>% filter(!Identifier %in% c("Cascade", "Zenith", "Symphony", "Nugget", "Target")) %>% group_by(resistant) %>% tally()
6+7
genotypes %>% filter(mildew == "HPM 204") %>% filter(!Identifier %in% c("Cascade", "Zenith", "Symphony", "Nugget", "Target")) %>% filter(resistant == "yes")
genotypes %>% filter(mildew == "HPM 1285") %>% filter(!Identifier %in% c("Cascade", "Zenith", "Symphony", "Nugget", "Target")) %>% group_by(resistant) %>% tally()
