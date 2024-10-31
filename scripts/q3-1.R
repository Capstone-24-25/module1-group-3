library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('../data/biomarker-clean.RData')

## training & testing partition ##
##################################

set.seed(1234)

data_split <- initial_split(biomarker_clean, prop = 0.8, strata = group)
data_train <- training(data_split)
data_test <- testing(data_split)

## training split multiple testing ##
#####################################

ttests_out_train <- data_train %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1_train <- ttests_out_train %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## training split random forest ##
##################################
predictors_train <- data_train %>%
  select(-c(group, ados))

response_train <- data_train %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_train <- randomForest(x = predictors_train, 
                       y = response_train, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_train$confusion

# compute importance scores
proteins_s2_train <- rf_train$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## training split logistic regression ##
########################################

# select subset of interest
proteins_sstar_train <- intersect(proteins_s1_train, proteins_s2_train)

biomarker_sstar_train <- data_train %>%
  select(group, any_of(proteins_sstar_train)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_sstar_test <- data_test %>%
  select(group, any_of(proteins_sstar_train)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# fit logistic regression model to training set
fit_train <- glm(class ~ ., 
           data = biomarker_sstar_train, 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

results_train <- biomarker_sstar_test %>%
  add_predictions(fit_train, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second') %>% 
  write.csv("../results/figures/q3-1_results.csv")

proteins_s1_train %>%  write.csv("../results/figures/q3-1_proteins_s1.csv")
proteins_s2_train %>%  write.csv("../results/figures/q3-1_proteins_s2.csv")
proteins_sstar_train %>% write.csv("../results/figures/q3-1_proteins_sstar.csv")

## Notes: ##
# different proteins selected
# same sensitivity
# lower specificity, accuracy, and roc_auc