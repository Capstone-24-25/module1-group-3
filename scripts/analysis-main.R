#### Task 1

load("../data/biomarker-clean.RData")

# clean data without log-transformation
var_names <- read_csv('../data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

biomarker_raw <- read_csv('../data/biomarker-raw.csv', 
                          skip = 2,
                          col_select = -2L,
                          col_names = c('group', 
                                        'empty',
                                        pull(var_names, abbreviation),
                                        'ados'),
                          na = c('-', '')) %>%
  filter(!is.na(group))

#investigate the potential reason for log-transformation
# check normality for some protein

par(mfrow = c(1, 2))

qqnorm(main = "Raw Data Q-Q plot for `Mcl-1`", biomarker_raw$`Mcl-1`)
qqline(biomarker_raw$`Mcl-1`, col = "red")
qqnorm(main = "Clean Data Q-Q plot for `Mcl-1`",biomarker_clean$`Mcl-1`)
qqline(biomarker_clean$`Mcl-1`, col = "red")

par(mfrow = c(1, 2))
qqnorm(main = "Clean Data Q-Q plot for `DERM`",biomarker_raw$DERM)
qqline(biomarker_raw$DERM, col = "red")
qqnorm(main = "Clean Data Q-Q plot for `DERM`",biomarker_clean$DERM)
qqline(biomarker_clean$DERM, col = "red")

#### Task 2

z_score <- biomarker_clean %>%
  group_by(group) %>%
  mutate(across(where(is.numeric), ~ scale(.))) %>%
  ungroup()

outlier_counts <- z_score %>%
  mutate(across(where(is.numeric), ~ abs(.) > 3)) %>%
  group_by(group) %>%
  summarise(across(where(is.logical), sum, na.rm = T))

outliers_ge_4 <- outlier_counts %>%
  pivot_longer(cols = -group, names_to = "Column", values_to = "Outlier_Count") %>%
  filter(Outlier_Count >= 4) %>%
  arrange(group, desc(Outlier_Count))
outliers_ge_4

group_counts <- outliers_ge_4 %>%
  count(group)
group_counts

total_outliers_by_group <- outlier_counts %>%
  rowwise() %>%
  mutate(Total_Outliers = sum(c_across(where(is.numeric)), na.rm = T)) %>%
  select(group, Total_Outliers)
total_outliers_by_group


### Methodological variations

# **Original Methodology (from *inclass-analysis.R*):**
  

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
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
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

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
# rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

# **Task 3** **(Data split before variable selection)**
  
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

biomarker_sstar_test %>%
  add_predictions(fit_train, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')


# **Task 3 (larger panel of proteins)**

# select top 20 proteins from multiple testing
proteins_s1_top20 <- ttests_out %>%
  slice_min(p.adj, n = 20) %>%
  pull(protein)

# select top 20 proteins from random forest
proteins_s2_top20 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 20) %>%
  pull(protein)

# logistic regression
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

testing(biomarker_split_top20) %>%
  add_predictions(fit_top20, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')


# **Task 3 (Fuzzy Intersection)**

# MULTIPLE TESTING SELECTION
proteins_s1_fuzzy <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST SELECTION
proteins_s2_fuzzy <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## FUZZY INTERSECTION
#####################
#Combine scores from both methods and assign fuzzy scores
all_proteins_fuzzy <- unique(c(proteins_s1_fuzzy, proteins_s2_fuzzy))

fuzzy_scores <- tibble(protein = all_proteins_fuzzy) %>%
  mutate(score = (protein %in% proteins_s1_fuzzy) + (protein %in% proteins_s2_fuzzy))

fuzzy_threshold <- 1

## LOGISTIC REGRESSION

# select subset of interest
proteins_sstar_fuzzy <- fuzzy_scores %>%
  filter(score >= fuzzy_threshold) %>%
  pull(protein)

biomarker_sstar_fuzzy <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_fuzzy)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split_fuzzy <- biomarker_sstar_fuzzy %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit_fuzzy <- glm(class ~ ., 
           data = training(biomarker_split_fuzzy), 
           family = 'binomial')

testing(biomarker_split_fuzzy) %>%
  add_predictions(fit_fuzzy, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
              truth = tr_c, pred,
              event_level = 'second')

### Improved classifier

# **Task 4 (Improved panel of proteins)**
  
# selected top 15 from multiple testing to prevent over-fitting
proteins_s1_improved <- ttests_out %>%
  slice_min(p.adj, n = 15) %>% 
  pull(protein)

# selected top 15 from multiple testing
proteins_s2_improved <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 15) %>% 
  pull(protein)

# logistic regression
proteins_sstar_improved <- intersect(proteins_s1_improved, proteins_s2_improved)

biomarker_sstar_improved <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_improved)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split_improved <- biomarker_sstar_improved %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit_improved <- glm(class ~ ., 
                    data = training(biomarker_split_improved), 
                    family = 'binomial')

testing(biomarker_split_improved) %>%
  add_predictions(fit_improved, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')