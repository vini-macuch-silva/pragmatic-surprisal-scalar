library(tidyverse)
library(forcats)
library(fields)
library(factoextra)
library(brms)
library(bayesplot)
options(mc.cores = parallel::detectCores ())

# For bootstrapping 95% confidence intervals
library(bootstrap)
theta <- function(x,xdata,na.rm=T) {mean(xdata[x],na.rm=na.rm)}
ci.low <- function(x,na.rm=T) {
  mean(x,na.rm=na.rm) - quantile(bootstrap(1:length(x),1000,theta,x,na.rm=na.rm)$thetastar,.025,na.rm=na.rm)}
ci.high <- function(x,na.rm=T) {
  quantile(bootstrap(1:length(x),1000,theta,x,na.rm=na.rm)$thetastar,.975,na.rm=na.rm) - mean(x,na.rm=na.rm)}

# Import data
d <- read_csv("results-25-11-sentence.csv", locale = locale(encoding = 'ISO-8859-1'))

# Check number of data points per participant in order to remove incomplete submissions
table(d$submission_id, d$trial_number)

# Remove incomplete submissions and practice trials
d <- d %>% 
  filter(submission_id != 7950, submission_id != 7954) %>% 
  filter(trial_name != "practice_task_three")

# Convert variable 'trial_type' to factor
d$trial_type <- as.factor(d$trial_type)

# Clean data set
# Remove unnecessary variables
d2 <- d %>%
  filter(is.na(trial_type)) %>%
  select(-trial_type, -property_1, -property_2, -shape_1, -shape_2, -response_1, -response_2, -response_3, -display, -experiment_id) %>% 
  gather(Region, RT, -submission_id, -quantifier, -picture_type, -condition, -sentence, -response, -listNumber,
         -trial_number, -picture, -trial_name)

# Convert variable 'Region' to factor
d2$Region = factor(x = d2$Region,
                   levels = c("QUANT", "der", "SHAPE", "auf",
                              "dem", "Bild", "sind", "CRIT_3", "der_1",
                              "Box"))

# Convert to factor variables to be included in the models as categorical predictors
d2$submission_id <- as.factor(d2$submission_id)
d2$quantifier <- as.factor(d2$quantifier)
d2$condition <- as.factor(d2$condition)

# Check for outliers
# Check for outlier trials (= RTs larger/ smaller than 2.5 SD of the mean)
# Check for outlier participants (= total number of outlier trials larger than a third of all trials)
d2 <- d2 %>% group_by(Region, quantifier, condition) %>%
  mutate(outlier_trial = RT < (mean(RT) - 2.5*sd(RT)) | RT > (mean(RT) + 2.5*sd(RT)) ) %>%
  ungroup() %>% 
  group_by(submission_id) %>% 
  mutate(outlier_participant = mean(outlier_trial, na.rm = TRUE) > 0.3) %>% 
  ungroup() 

# Show outliers
show(paste0("Excluded trials: ", sum(d2$outlier_trial, na.rm = TRUE) ))
show(paste0("Excluded participants: ", sum(d2$outlier_participant) ))

# Save data set with removed outliers
d2 <- d2 %>% filter(outlier_trial == FALSE) %>% 
  select(-outlier_trial, -outlier_participant)

# Create variable 'logCon' (= logical conditions)
# Variable consists of the relevant combinations of 'quantifier' and 'condition', as per the experimental conditions
d2Logic <- d2 %>%
  filter(trial_name != "practice_task_three") %>% 
  mutate(logCon = case_when(quantifier == "some" & condition == "unbiased" ~ "Einige (Unbiased)",
                            quantifier == "some" & condition == "biased" ~ "Einige (Biased)",
                            quantifier == "some" & condition == "false" ~ "Einige (False)",
                            quantifier == "some" & condition == "underspecified" ~ "Einige (Infelict)",
                            quantifier == "all" & condition == "unbiased" ~ "Alle (Unbiased)",
                            quantifier == "all" & condition == "biased" ~ "Alle (Biased)",
                            quantifier == "all" & condition == "false" ~ "Alle (False)",
                            TRUE ~ "other")) %>% 
  select(-quantifier, -condition) %>% 
  filter(picture_type != "A" | logCon != "Einige (False)")

# Convert variable 'logCon' to factor
d2Logic$logCon <- factor(x = d2Logic$logCon, levels = c("Einige (Biased)", "Einige (Unbiased)", "Einige (Infelict)", "Alle (Biased)", "Alle (Unbiased)", "Einige (False)", "Alle (False)"))

# Extract picture labels to be used in the random structure of the model
d2Logic$picture <- str_sub(d2Logic$picture, 12)
d2Logic$picture <- gsub("-", "_", d2Logic$picture)
d2Logic$picture <- gsub("/", "_", d2Logic$picture)
d2Logic$picture <- str_sub(d2Logic$picture, end = -5)

# Filter unique item labels
picture_levels <- sort(unique(d2Logic$picture))

# Save unique item labels as levels of the variable 'picture'
levels(d2Logic$picture) <- picture_levels

# Reorder levels of variable 'logCon'
# Set "Einige (Biased)" as reference level
d2Logic$logCon = factor(d2Logic$logCon, ordered = F, levels = c("Einige (Biased)", "Einige (Unbiased)", "Einige (Infelict)", "Einige (False)", "Alle (Biased)", "Alle (Unbiased)", "Alle (False)"))

# Convert variable 'picture' to factor
d2Logic$picture = factor(x = d2Logic$picture, levels = picture_levels)

# Extract strings with shape terms from variable 'sentence' and store in a new variable called 'critical'
d2Logic$critical <- str_extract(d2Logic$sentence, "Quadrate|Dreiecke|Dreicke|Kreise")
d2Logic$critical <- gsub("Dreicke", "Dreiecke", d2Logic$critical)

# Convert variables 'critical' and 'trial_name' to factor
d2Logic$critical <- as.factor(d2Logic$critical)
d2Logic$trial_name <- as.factor(d2Logic$trial_name)

# table(d2Logic$logCon, d2Logic$critical)

# Run statistical model predicting RTs at the SHAPE region as a function of the logical conditions 
# RT_SHAPE = brm(log(RT) ~ logCon +
#                           (1 + logCon | submission_id) +  
#                           (1 | picture),
#                         iter = 4000,
#                         filter(d2Logic, Region == "SHAPE"))
  
# Save output (model object)
# saveRDS(RT_SHAPE, "RT_SHAPE.RDS")

# Run statistical model predicting RTs at the SHAPE region as a function of the logical conditions 
# and shape term
# RT_SHAPE_crit = brm(log(RT) ~ logCon * critical +
#                        (1 + logCon * critical | submission_id) +
#                        (1 | picture),
#                      iter = 4000,
#                      filter(d2Logic, Region == "SHAPE"))

# Save output (model object)
# saveRDS(RT_SHAPE_crit, "RT_SHAPE_crit.RDS")

# Run statistical model predicting RTs at the SHAPE region as a function of the logical conditions
# and experimental block
RT_SHAPE_trial = brm(log(RT) ~ logCon * trial_name +
                          (1 + logCon + trial_name | submission_id) +
                          (1 + trial_name | picture),
                        iter = 4000,
                        filter(d2Logic, Region == "SHAPE"))

# Save output (model object)
saveRDS(RT_SHAPE_trial, "RT_SHAPE_trial.RDS")


#############################################################
########                                               ######
########  FROM HERE ON NOT RELEVANT FOR YOUR CHECK-UP  ######
########                                               ######
#############################################################


# Run statistical model predicting RTs at the SHAPE region as a function of the logical conditions, 
# shape term, and experimental block
# RT_SHAPE_crit_trial = brm(log(RT) ~ logCon * critical * trial_name +
#                        (1 + logCon + critical + trial_name | submission_id) +
#                        (1 + trial_name | picture),
#                      iter = 4000,
#                      filter(d2Logic, Region == "SHAPE"))

# Save output (model object)
# saveRDS(RT_SHAPE_crit_trial, "RT_SHAPE_crit_trial.RDS")

# Load model object
RT_SHAPE <- readRDS("RT_SHAPE.RDS")

# Extract summary of samples
# Show probabilities of relevant comparisons
samples_summary_SHAPE <-  posterior_samples(RT_SHAPE) %>%
summarize(ProbabilityEinigeBiasedFasterEinigeUnbiased = sum((b_logConEinigeUnbiased > 0)) / n(),
          ProbabilityEinigeBiasedFasterEinigeInfelict = sum((b_logConEinigeInfelict > 0)) / n(),
          ProbabilityEinigeBiasedFasterEinigeFalse = sum((b_logConEinigeFalse > 0)) / n(),
          ProbabilityAlleBiasedFasterAlleUnbiased = sum((b_logConAlleUnbiased > b_logConAlleBiased)) / n(),
          ProbabilityAlleBiasedFasterAlleFalse = sum((b_logConAlleFalse > b_logConAlleBiased)) / n() ) %>%
  gather(key = "Hypothesis", value = "Probability") %>%
  mutate(Hypothesis = recode(Hypothesis, "ProbabilityEinigeBiasedFasterEinigeUnbiased" = "[einige] Unbiased > Biased",
                             "ProbabilityEinigeBiasedFasterEinigeInfelict" = "[einige] Infelict > Biased",
                             "ProbabilityEinigeBiasedFasterEinigeFalse" = "[einige] False > Biased",
                             "ProbabilityAlleBiasedFasterAlleUnbiased" = "[alle] Unbiased > Biased",
                             "ProbabilityAlleBiasedFasterAlleFalse" = "[alle] False > Biased"))

# Produce clean model summary (model estimates, credibility interval)
RT_SHAPE_tidy <- broom.mixed::tidy(RT_SHAPE)

# Filter relevant information (fixed effects)
RT_SHAPE_tidy2 <- filter(RT_SHAPE_tidy, effect == "fixed") %>% 
  select(-effect, -component, -group, -std.error) %>% 
  mutate_if(is.numeric, ~ round(., 2))


###############################
##### Exploratory analysis ####
###############################

# Create respondent type
d2LogicResp <- d2Logic %>%
  group_by(submission_id) %>%
  mutate(response_type = case_when(logCon == "Einige (Infelict)" & (response == "7" | response == "6" | response == "5") ~ "Semantic",
                                   logCon == "Einige (Infelict)" & (response == "1" | response == "2" | response == "3") ~ "Pragmatic",
                                   logCon == "Einige (Infelict)" & response == "4" ~ "Neutral",
                                   TRUE ~ "Other")) %>%
  mutate(respond_type = case_when(sum(response_type == "Pragmatic") == 0 & sum(response_type == "Semantic") > 0 ~ "Semantic",
                                  sum(response_type == "Semantic") == 0 & sum(response_type == "Pragmatic") > 0 ~ "Pragmatic",
                                  TRUE ~ "Mixed")) %>%
  ungroup() 


# Convert variable 'respond_type' to factor
d2LogicResp$respond_type <- as.factor(d2LogicResp$respond_type)

# Reorder levels of variable 'respond_type'
# Set 'Semantic' as reference level
d2LogicResp$respond_type = factor(d2LogicResp$respond_type, ordered = F, levels = c("Semantic", "Pragmatic", "Mixed"))

# Run statistical model predicting RTs at the SHAPE region as a function of the logical conditions
# and respondent type
# RT_SHAPE_respond = brm(log(RT) ~ logCon * respond_type +
#                           (1 + logCon | submission_id) +  
#                           (1 + respond_type | picture),
#                         iter = 4000,
#                         filter(d2LogicResp, Region == "SHAPE"))

# Save output (model object)
# saveRDS(RT_SHAPE_respond, "RT_SHAPE_respond.RDS")

# Run statistical model predicting RTs at the SHAPE region as a function of the logical conditions,
# respondent type, and shape term
# RT_SHAPE_respond_crit = brm(log(RT) ~ logCon * respond_type * critical +
#                               (1 + logCon + critical | submission_id) +
#                               (1 + respond_type | picture),
#                             iter = 4000,
#                             filter(d2LogicResp, Region == "SHAPE"))

# Save output (model object)
# saveRDS(RT_SHAPE_respond_crit, "RT_SHAPE_respond_crit.RDS")

# Run statistical model predicting RTs at the SHAPE region as a function of the logical conditions,
# respondent type, and experimental block
# RT_SHAPE_respond_trial = brm(log(RT) ~ logCon * respond_type * trial_name +
#                           (1 + logCon + trial_name | submission_id) +
#                           (1 + respond_type | picture),
#                         iter = 4000,
#                         filter(d2LogicResp, Region == "SHAPE"))

# Save output (model object)
# saveRDS(RT_SHAPE_respond_trial, "RT_SHAPE_respond_trial.RDS")

# Run statistical model predicting RTs at the SHAPE region as a function of the logical conditions,
# respondent type, shape term, and experimental block
# RT_SHAPE_respond_crit_trial = brm(log(RT) ~ logCon * respond_type * critical * trial_name +
#                               (1 + logCon + critical + trial_name | submission_id) +
#                               (1 + respond_type + trial_name | picture),
#                             iter = 4000,
#                             filter(d2LogicResp, Region == "SHAPE"))

# Save output (model object)
# saveRDS(RT_SHAPE_respond_crit_trial, "RT_SHAPE_respond_crit_trial.RDS")

# Load model object
RT_SHAPE_respond <- readRDS("RT_SHAPE_respond.RDS")

# Extract summary of samples
# Show probabilities of relevant comparisons
samples_summary_SHAPE_respond <-  posterior_samples(RT_SHAPE_respond) %>%
  summarize(ProbabilityEinigeBiasedFasterEinigeUnbiasedSemantic = sum((b_logConEinigeUnbiased > 0)) / n(),
            ProbabilityEinigeBiasedFasterEinigeInfelictSemantic = sum((b_logConEinigeInfelict > 0)) / n(),
            ProbabilityEinigeBiasedFasterEinigeFalseSemantic = sum((b_logConEinigeFalse > 0)) / n(),
            ProbabilityAlleBiasedFasterAlleUnbiasedSemantic = sum((b_logConAlleUnbiased > b_logConAlleBiased)) / n(),
            ProbabilityAlleBiasedFasterAlleFalseSemantic = sum((b_logConAlleFalse > b_logConAlleBiased)) / n(),
            ProbabilityEinigeBiasedFasterEinigeUnbiasedPragmatic = sum((b_logConEinigeUnbiased + `b_logConEinigeUnbiased:respond_typePragmatic` > 0)) / n(),
            ProbabilityEinigeBiasedFasterEinigeInfelictPragmatic = sum((b_logConEinigeInfelict + `b_logConEinigeInfelict:respond_typePragmatic` > 0)) / n(),
            ProbabilityEinigeBiasedFasterEinigeFalsePragmatic = sum((b_logConEinigeFalse + `b_logConEinigeFalse:respond_typePragmatic` > 0)) / n(),
            ProbabilityAlleBiasedFasterAlleUnbiasedPragmatic = sum((b_logConAlleUnbiased + `b_logConAlleUnbiased:respond_typePragmatic` > b_logConAlleBiased + `b_logConAlleBiased:respond_typePragmatic`)) / n(),
            ProbabilityAlleBiasedFasterAlleFalsePragmatic = sum((b_logConAlleFalse + `b_logConAlleFalse:respond_typePragmatic` > b_logConAlleBiased + `b_logConAlleBiased:respond_typePragmatic`)) / n(),
            ProbabilityEinigeBiasedFasterEinigeUnbiasedMixed = sum((b_logConEinigeUnbiased + `b_logConEinigeUnbiased:respond_typeMixed` > 0)) / n(),
            ProbabilityEinigeBiasedFasterEinigeInfelictMixed = sum((b_logConEinigeInfelict + `b_logConEinigeInfelict:respond_typeMixed` > 0)) / n(),
            ProbabilityEinigeBiasedFasterEinigeFalseMixed = sum((b_logConEinigeFalse + `b_logConEinigeFalse:respond_typeMixed` > 0)) / n(),
            ProbabilityAlleBiasedFasterAlleUnbiasedMixed = sum((b_logConAlleUnbiased + `b_logConAlleUnbiased:respond_typeMixed` > b_logConAlleBiased + `b_logConAlleBiased:respond_typeMixed`)) / n(),
            ProbabilityAlleBiasedFasterAlleFalseMixed = sum((b_logConAlleFalse + `b_logConAlleFalse:respond_typeMixed` > b_logConAlleBiased + `b_logConAlleBiased:respond_typeMixed`)) / n()) %>%
  gather(key = "Hypothesis", value = "Probability") %>%
  mutate(Hypothesis = recode(Hypothesis, "ProbabilityEinigeBiasedFasterEinigeUnbiasedSemantic" = "[einige][semantic] Unbiased > Biased",
                             "ProbabilityEinigeBiasedFasterEinigeInfelictSemantic" = "[einige][semantic] Infelict > Biased",
                             "ProbabilityEinigeBiasedFasterEinigeFalseSemantic" = "[einige][semantic] False > Biased",
                             "ProbabilityAlleBiasedFasterAlleUnbiasedSemantic" = "[alle][semantic] Unbiased > Biased",
                             "ProbabilityAlleBiasedFasterAlleFalseSemantic" = "[alle][semantic] False > Biased",
                             "ProbabilityEinigeBiasedFasterEinigeUnbiasedPragmatic" = "[einige][pragmatic] Unbiased > Biased",
                             "ProbabilityEinigeBiasedFasterEinigeInfelictPragmatic" = "[einige][pragmatic] Infelict > Biased",
                             "ProbabilityEinigeBiasedFasterEinigeFalsePragmatic" = "[einige][pragmatic] False > Biased",
                             "ProbabilityAlleBiasedFasterAlleUnbiasedPragmatic" = "[alle][pragmatic] Unbiased > Biased",
                             "ProbabilityAlleBiasedFasterAlleFalsePragmatic" = "[alle][pragmatic] False > Biased",
                             "ProbabilityEinigeBiasedFasterEinigeUnbiasedMixed" = "[einige][mixed] Unbiased > Biased",
                             "ProbabilityEinigeBiasedFasterEinigeInfelictMixed" = "[einige][mixed] Infelict > Biased",
                             "ProbabilityEinigeBiasedFasterEinigeFalseMixed" = "[einige][mixed] False > Biased",
                             "ProbabilityAlleBiasedFasterAlleUnbiasedMixed" = "[alle][mixed] Unbiased > Biased",
                             "ProbabilityAlleBiasedFasterAlleFalseMixed" = "[alle][mixed] False > Biased"))


# Produce clean model summary (model estimates, credibility interval)
RT_SHAPE_respond_tidy <- broom.mixed::tidy(RT_SHAPE_respond)

# Filter relevant information (fixed effects)
RT_SHAPE_respond_tidy2 <- filter(RT_SHAPE_respond_tidy, effect == "fixed") %>% 
  select(-effect, -component, -group, -std.error) %>% 
  mutate_if(is.numeric, ~ round(., 2))
