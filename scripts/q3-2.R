library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('../data/biomarker-clean.RData')

## multiple testing - top 20 ##
###############################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out_top20 <- biomarker_clean %>%
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
proteins_s1_top20 <- ttests_out_top20 %>%
  slice_min(p.adj, n = 20) %>%
  pull(protein)

## random forest - top 20##
###########################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2_top20 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 20) %>%
  pull(protein)

## logistic regression - top 20##
#################################

# select subset of interest
proteins_sstar_top20 <- intersect(proteins_s1_top20, proteins_s2_top20)

biomarker_sstar_top20 <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_top20)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split_top20 <- biomarker_sstar_top20 %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit_top20 <- glm(class ~ ., 
           data = training(biomarker_split_top20), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split_top20) %>%
  add_predictions(fit_top20, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second') %>% 
  write.csv("../results/figures/q3-2_results.csv")

# storing protein selections
proteins_s1_top20 %>%  write.csv("../results/figures/q3-2_proteins_s1.csv")
proteins_s2_top20 %>%  write.csv("../results/figures/q3-2_proteins_s2.csv")
proteins_sstar_top20 %>% write.csv("../results/figures/q3-2_proteins_sstar.csv")

## Notes: ##
# lower sensitivity, accuracy and roc_auc
# slightly higher specificity