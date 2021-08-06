# Load packages
library(tidyverse)
library(tidymodels)
# library(broom.mixed) # for converting bayesian models to tidy tibbles
# library(dotwhisker)  # for visualizing regression results
library(themis) # step_smote


############################################################################### I ### Exploratory Data Analysis ----
path <- fs::path("", "Volumes", "Peres_Research", "Ovarian - Radiomics")

# No NA in features

# Clinical
clinical <- read_rds("clinical.rds")

# recurrence by age or year
library(GGally)

clinical %>% 
  select(has_the_patient_recurred, year_of_diagnosis, treatment_type, debulking_status,
         age_at_diagnosis, primary_site) %>% 
  ggpairs(columns = 2:ncol(.), aes(color = has_the_patient_recurred, alpha = 0.5))

clinical %>% 
  select(has_the_patient_recurred, treatment_type, debulking_status,
         primary_site, year_of_diagnosis, age_at_diagnosis) %>% 
  mutate(year_of_diagnosis = as.character(year_of_diagnosis)) %>%
  mutate(age_at_diagnosis = as.character(age_at_diagnosis)) %>%
  pivot_longer(2:ncol(.)) %>% 
  ggplot(aes(y = value, fill = has_the_patient_recurred)) +
  geom_bar(position = "fill")+
  facet_wrap(vars(name), scales = "free")+
  labs(x= NULL, y = NULL)

clinical %>% 
  group_by(year_of_diagnosis) %>% 
  summarise(recurred = sum(rec_event)) %>%
  ggplot(aes(year_of_diagnosis, recurred))+
  geom_line()

library(naniar) # upsetplot

clinical %>% 
  select(has_the_patient_recurred, ecog_pretrt, ecog_posttrt, ecog_recurr, 
         Ethnicity, grade_differentiation) %>%
  gg_miss_upset()
# Will do unknown, remove ecog

clinical %>% 
  select(has_the_patient_recurred, germline_brca1_mutation, germline_brca2_mutation, somatic_brca1_mutation,
         somatic_brca2_mutation, any_unclassified_brca_mutation) %>%
  gg_miss_upset()
# Will do unknown, remove any_unclassified_brca_mutation

clinical %>% 
  select(hypertension, diabetes_mellitus, hypercholesterolemia,
         chronic_kidney_disease, cardiac_conditions_including_bu) %>%
  gg_miss_upset()

############################################################################### II ### Data prepping ----
# remove variables non meaningful or redundant to recurrence
clinical_ml <- read_rds("clinical.rds") %>% 
  select(-c(vital_new, 
            gender, race, ethnicity, 
            baseline_ctscan_outside_moffitt, tumor_sequence_number, 
            "histology", tnm_edition_number_must_use, grade_differentiation,
         date_of_birth, any_unclassified_brca_mutation,
         comment_for_cardiac_comorbidity,
         vital_date_new, date_of_last_followup, fwdate_most_recent,
         germline_brca1_mutation, germline_brca2_mutation, somatic_brca1_mutation, ##### To add
         somatic_brca2_mutation, any_unclassified_brca_mutation,
         "hypertension", "diabetes_mellitus", "hypercholesterolemia",
         "chronic_kidney_disease", "cardiac_conditions_including_bu",
         comment_for_cardiac_comorbidity,
         "complete",
         ecog_posttrt, ecog_recurr, "ecog_pretrt_date", "ecog_posttrt_date", "ecog_recurr_date",
         contrast_enhancement, # What is the difference with w_w...
         age_at_Dx,
         
         months_at_first_neoadjuvant_chem, months_at_first_adjuvant_chem, months_at_first_chemo,
         months_at_first_surgery, age_at_surgery, age_at_first_recurrence, month_at_first_recurrence_Dx,
         months_at_surg_followup, months_at_neo_followup, months_at_chem_followup, months_of_surg_rec_free,
         months_of_neo_rec_free, months_of_chem_rec_free,
         months_of_dx_rec_free, recurrence_time, months_of_treat_rec_free,
         os_time, rec_event, os_event,
         has_the_patient_recurredafter_surg
  )) %>% 
  
  select(-c( # Include as year, month?
    date_of_diagnosis, baseline_ct_scan_date,
            date_of_first_neoadjuvant_chemot,
            date_of_first_surgery, date_of_first_adjuvant_chemother,
            "date_of_first_recurrence", "date_of_surgery_abstracted",
    first_treatment_date, first_chemo_date, "months_at_first_treatment", "rec_event_date",
    "months_at_treat_followup", "recurrence_date_after_surgery"
            ))

# Select column name of stable features from the radiomics data
concordance <- read_rds("concordance.rds") %>% 
  filter(value_CCC >= 0.95) %>% 
  select(name) %>% 
  mutate(stable_features = str_match(name, "([a-z][:digit:]*)_*")[,2])

stable_features <- paste0(paste("nor_", concordance$stable_features, "[a-z]", sep = ""), collapse = "|")
# radiomics <- read_csv(
#   paste0(path,
#          "/data/Merge clinical and radiomics data/Ovarian_Normalized_Radiomics_Features_05242021.csv")) #%>% 
#   select(matches(stable_features))

mldata <- read_rds("radiomics.rds") %>% 
  select(mrn, matches(stable_features)) %>% 
  right_join(., clinical_ml, by = "mrn")

# Explore what will need to be changed
skimr::skim(mldata)

# Cleaning
# rm(concordance, stable_features)


############################################################################### II ### Build model ----
# Build the dataset
# warning for date variable
meaningful_dates <- 
  c("baseline_ct_scan_date", "date_of_first_neoadjuvant_chemot", 
    "date_of_first_surgery", "date_of_first_adjuvant_chemother")
mldata <- mldata %>% 
  drop_na(starts_with("nor")) %>% 
  select(-c(contains("date")#,
            # meaningful_dates
  )) %>% 
  mutate(tnm_cs_mixed_group_stage = as.character(tnm_cs_mixed_group_stage)) %>% 
  # it is important that the outcome variable for training a (logistic) regression model is a factor.
  mutate_if(is.character, as.factor)

skimr::skim(mldata)


set.seed(1234)

# 1. Splitting the data
# 3/4 of the data into the training set but split evenly winthin race
data_split <- initial_split(mldata, prop = 3/4, strata = w_wo_contrast)
# Create training and testing datasets:
train_data <- training(data_split)
test_data  <- testing(data_split)

# 2. Data pre processing and features engineering + imputation
# Recipe
mldata_recipe <-
  # 1.model formula
  recipe(has_the_patient_recurred ~ ., data = train_data)  %>% 
  # 2.keep these variables but not use them as either outcomes or predictors
  # Keep id but not include in the model as predictor
  update_role(mrn, new_role = "ID") %>%
  # update_role(w_wo_contrast, new_role = "Other") %>% 
  # remove variables that contain only a single value.
  step_zv(all_predictors()) %>% # or step_nzv
  # 3.If factor with too much levels, collapse lower levels
  # step_other(Histology, threshold = 0.05) %>% 
  # data_recipe %>% prep() %>% juice() %>% count(Histology)
  
  # 4.Imputation
  step_unknown(all_nominal_predictors()) %>% 
  step_impute_median(all_numeric_predictors()) %>%
  # 5.change all factor to indicator/dummy variables for model that cannot handle factor variables
  step_dummy(all_nominal(), -all_outcomes()) #%>%
  
  # Feature engineering on dates
  # step_date(all_of(meaningful_dates), features = c("year", "month")) %>% 
  # step_rm(meaningful_dates)
  
  # LAST.For imbalance, model memorize the few example and
  # step_smote(has_the_patient_recurred) # Use nearest neighbor to create new synthetic observation almost similar

mldata_recipe

# estimate the required parameters from a training set that can be later applied to other data sets.
# learn what the model should be with the training data
mldata_prep <- prep(mldata_recipe, verbose = TRUE, log_changes = TRUE)
# Extract Finalized Training Set
juiced <- juice(mldata_prep)

############################################################################### II ### Data Tuning ----
# train hyperparameter
set.seed(123)
# 10 fold cross validation
mldata_folds <- vfold_cv(train_data, strata = w_wo_contrast)

############################################################################################### RAMDOM FOREST ----
# Model specification
ranger_spec <- rand_forest(
  # tune right value for the number of predictors that will be randomly sampled at each split when creating the tree models
  mtry = tune(), 
  # tune right value for the minimum number of data points in a node that are required for the node to be split further.
  min_n = tune(),
  trees = 1000) %>% 
  set_mode("classification") %>% 
  set_engine("ranger")

# Set up a work flow
ranger_workflow <- 
  workflow() %>% 
  # Add preprocessor
  add_recipe(mldata_recipe) %>% # here is an unfit mode
  add_model(ranger_spec) # add our model specification

# Tuning
set.seed(54691)
doParallel::registerDoParallel()
# will tune mtry and min_m on a grid
ranger_tune <-
  tune_grid(ranger_workflow, # Will take our workflow and apply it on
            resamples = mldata_folds, # each fold of the data for 
            grid = 20) # How many candidate point do I want to try
# Try the 20 point with each min_n and mtry
# Then train on analysis cv and assessed on the test of cv
# will compute performance matrix on each of these of the 20 candidate model

# May need to add specificity and sensitivity to metrics if we have rare events. these 2 will tell us how the model did for our positive and negative cases
# If sens is low it means the model had a really hard time finding the rare case (could mean step_smote is bad idea for this model)


############################################################################### II ### Explore Tuning Results ----

show_best(ranger_tune, metric = "roc_auc")
show_best(ranger_tune, metric = "accuracy")

# Visualize tuned parameters
autoplot(ranger_tune)
# mtry is the number of predictor randomly selected -> needs to be 
# mim_n is the minimal node size -> needs to be small

ranger_tune %>% 
  collect_metrics() %>% 
  filter(.metric == "roc_auc") %>% 
  select(mean, min_n, mtry) %>% 
  pivot_longer(min_n:mtry,
               values_to = "value",
               names_to = "param") %>% 
  ggplot(aes(value, mean, color = param)) +
  geom_point(show.legend = FALSE) +
  facet_wrap( . ~ param, scales = "free_x")

ranger_tune %>% select_best("accuracy")

# # Tune more?
# # Can make a regular grid
# new_grid <- grid_regular(
#   mtry(range = c(25, 75)),
#   min_n(range = c(0, 25)),
#   levels = 10
# )
# 
# set.seed(123)
# sec_tune_results <- tune_grid( # will tune mtry and min_m on a grid
#   tune_wf, # tune worflow
#   resamples = data_folds, # on this data
#   grid = new_grid
# )
# 
# sec_tune_results %>% collect_metrics() %>% 
#   filter(.metric == "roc_auc") %>% 
#   mutate(min_n = factor(min_n)) %>% 
#   ggplot(aes(mtry, mean, color = min_n)) +
#   geom_line(alpha = 0.5, size = 1.5) +
#   geom_point()

# Step Finalize the model with our tuned parameters
# best_auc <- select_best(ranger_tune, # or sec_tune_results if tuned more,
#                         "roc_auc")
# 
# final_rf <- finalize_model(ranger_spec,
#                            best_auc)
# # Aka same
# final_rf <- ranger_spec %>% 
#   finalize_model(select_best(ranger_tune, "accuracy"))

final_rf <- ranger_workflow %>% 
  finalize_workflow(select_best(ranger_tune, "accuracy"))
final_rf


############################################################################################### LOGISTIC REGRESSION ----

# Tune lambda

# logistic reg, choose regularized log reg with mixture 1 for Lasso model
# Because lots of var for few rows (samples) but don't know which vars are important

# Find the optimal value of lambda that minimizes the cross-validation error
# Need to prep data differently
train_data_lambda <- training(data_split) %>% 
  drop_na() %>% select(-mrn)

# Dummy code categorical predictor variables
x <- model.matrix(has_the_patient_recurred ~ ., data = train_data_lambda)[,-1]
# Convert the outcome (class) to a numerical variable
y <- ifelse(train_data_lambda$has_the_patient_recurred == "Recurrence", 1, 0)

library(glmnet)
set.seed(123)
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cv.lasso)
cv.lasso$lambda.min
# The left dashed vertical line indicates that the log of the optimal value of lambda is approximately -4, 
# which is the one that minimizes the prediction error. This lambda value will give the most accurate model.
cv.lasso$lambda.1se
# Fit the final model on the training data
model <- glmnet(x, y, alpha = 1, family = "binomial",
                lambda = cv.lasso$lambda.min)
# Display regression coefficients
coef(model)

# library(usemodels) # Gives a scaffolding of the modeling code
# use_glmnet(has_the_patient_recurred ~ ., data = train_data)

# glmnet_recipe <- 
#   recipe(formula = has_the_patient_recurred ~ ., data = train_data) %>% 
#   step_novel(all_nominal(), -all_outcomes()) %>% 
#   step_dummy(all_nominal(), -all_outcomes()) %>% 
#   step_zv(all_predictors()) %>% 
#   step_normalize(all_predictors(), -all_nominal()) 

glmnet_spec <- 
  logistic_reg(penalty = tune(), mixture = tune()) %>% 
  set_mode("classification") %>% 
  set_engine("glmnet") 

glmnet_workflow <- 
  workflow() %>% 
  add_recipe(mldata_recipe) %>% 
  add_model(glmnet_spec) 

glmnet_grid <- tidyr::crossing(penalty = 10^seq(-6, -1, length.out = 20), mixture = c(0.05, 
                                                                                      0.2, 0.4, 0.6, 0.8, 1)) 

glmnet_tune <- 
  tune_grid(glmnet_workflow, resamples = mldata_folds, grid = glmnet_grid) 

############################################################################### II ### Explore Tuning Results ----

show_best(glmnet_tune, metric = "roc_auc")
show_best(glmnet_tune, metric = "accuracy")

# Visualize tuned parameters
autoplot(glmnet_tune)
# mtry is the number of predictor randomly selected -> needs to be 
# mim_n is the minimal node size -> needs to be small

glmnet_tune %>% 
  collect_metrics() %>% 
  filter(.metric == "accuracy") %>% 
  select(mean, penalty, mixture) %>% 
  pivot_longer(penalty:mixture,
               values_to = "value",
               names_to = "param") %>% 
  ggplot(aes(value, mean, color = param)) +
  geom_point(show.legend = FALSE) +
  facet_wrap( . ~ param, scales = "free_x")

glmnet_tune %>% 
  collect_metrics() %>% 
  filter(.metric == "roc_auc") %>% 
  select(mean, penalty, mixture) %>% 
  pivot_longer(penalty:mixture,
               values_to = "value",
               names_to = "param") %>% 
  ggplot(aes(value, mean, color = param)) +
  geom_point(show.legend = FALSE) +
  facet_wrap( . ~ param, scales = "free_x")

glmnet_tune %>% select_best("roc_auc")

# Step Finalize the model with our tuned parameters
# best_acc <- select_best(glmnet_tune, # or sec_tune_results if tuned more,
#                         "accuracy")
# 
# final_glmnet <- finalize_model(glmnet_spec,
#                                best_acc)
# # Aka same
# final_glmnet <- glmnet_spec %>% 
#   finalize_model(select_best(glmnet_tune, "roc_auc"))

final_glmnet <- glmnet_workflow %>% 
  finalize_workflow(select_best(glmnet_tune, "roc_auc"))
final_glmnet


# Calculate Performance Metrics again with our last tuned model
#Plot the ROC curve
glmnet_results <- final_glmnet %>% 
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
    control = control_resamples(save_pred = TRUE)
  )

rf_results <- final_rf %>% 
  fit_resamples( # is not doing any tuning, measuring performance on cross validation data
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )

###############################################################################  Support vector machine # not necessary
###############################################################################  neural network

############################################################################################### EVALUATE MODELS ----
# Explore Performance Metrics
collect_metrics(rf_results)
collect_metrics(glmnet_results)
# Accuracy is llow
# Sensitivity has 50% chance finding the minority class pop

rf_results %>% 
  conf_mat_resampled()
glmnet_results %>% 
  conf_mat_resampled()

# Visualization
rf_results %>% # is the saved predictions
  collect_predictions() %>% 
  conf_mat(truth = has_the_patient_recurred, estimate = .pred_class) %>% 
  autoplot()
glmnet_results %>% 
  collect_predictions() %>% 
  conf_mat(truth = has_the_patient_recurred, estimate = .pred_class) %>% 
  autoplot()

rf_results %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(has_the_patient_recurred, .pred_Recurrence) %>% 
  autoplot()
glmnet_results %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(has_the_patient_recurred, .pred_Recurrence) %>% 
  autoplot()

rf_results %>% # Compare both models
  collect_predictions() %>% 
  mutate(model = "rf") %>% 
  bind_rows(glmnet_results %>% 
              collect_predictions() %>% 
              mutate(model = "glmet")) %>% 
  group_by(model) %>% 
  roc_curve(has_the_patient_recurred, .pred_Recurrence) %>% 
  autoplot()



# Calculate prediction after the fact
rf_results %>% collect_predictions() %>% group_by(id) %>% 
  ppv(truth = has_the_patient_recurred, estimate = .pred_class) # do histog
glmnet_results %>% collect_predictions() %>% group_by(id) %>% 
  ppv(truth = has_the_patient_recurred, estimate = .pred_class)



############################################################################################### Step Explore of features importance ----
library(vip) 

# Need to train the model one more time but without tuning to go faster
importance_spec <- glmnet_spec %>% 
  finalize_model(select_best(glmnet_tune, "roc_auc")) %>% 
  set_engine("glmnet", importance = "permutation") # permutation based importance

# represents the mean decrease in node impurity (and not the mean decrease in accuracy)
workflow() %>% 
  add_recipe(mldata_recipe) %>% 
  add_model(importance_spec) %>% 
  fit(train_data) %>% 
  extract_fit_parsnip() %>% 
  vip(aesthetics = list(alpha = 0.5, fill = "midnightblue"),
      # geom = "point",
      num_features = 20
  )

workflow() %>% 
  add_recipe(mldata_recipe) %>% 
  add_model(importance_spec) %>% 
  fit(train_data) %>% 
  extract_fit_parsnip() %>% 
  vi() %>% 
  mutate(Importance = ifelse(Sign == "NEG", -Importance, Importance)) %>% 
  ggplot(aes(x= desc(fct_reorder(Variable, abs(Importance)))  , y = Importance, fill = Sign)) +
  geom_bar(stat = "identity")

workflow() %>% 
  add_recipe(mldata_recipe) %>% 
  add_model(importance_spec) %>% 
  fit(train_data) %>% 
  extract_fit_parsnip() %>% 
  vi(method = "permute", target = , metrics = , reference = , pred_wrapper = ,
     train = juice(prep_data), nsim = 10) #%>% 
  # mutate(Importance = ifelse(Sign == "NEG", -Importance, Importance)) %>% 
  # ggplot(aes(x= desc(fct_reorder(Variable, abs(Importance)))  , y = Importance, fill = Sign)) +
  # geom_bar(stat = "identity")

workflow() %>% 
  add_recipe(mldata_recipe) %>% 
  add_model(importance_spec) %>% 
  fit(train_data) %>% 
  extract_fit_parsnip() %>% 
  vi() %>% 
  mutate(Importance = ifelse(Sign == "NEG", -Importance, Importance))


###################################### Ranger

# Need to train the model one more time but without tuning to go faster
importance_spec <- ranger_spec %>% 
  finalize_model(select_best(ranger_tune, "roc_auc")) %>% 
  set_engine("ranger", importance = "permutation") # permutation based importance

# represents the mean decrease in node impurity (and not the mean decrease in accuracy)
workflow() %>% 
  add_recipe(mldata_recipe) %>% 
  add_model(importance_spec) %>% 
  fit(train_data) %>% 
  extract_fit_parsnip() %>% 
  vip(aesthetics = list(alpha = 0.5, fill = "midnightblue"),
      # geom = "point",
      num_features = 20
  )

workflow() %>% 
  add_recipe(mldata_recipe) %>% 
  add_model(importance_spec) %>% 
  fit(train_data) %>% 
  extract_fit_parsnip() %>% 
  vi() %>% 
  mutate(Importance = ifelse(Sign == "NEG", -Importance, Importance)) %>% 
  ggplot(aes(x= desc(fct_reorder(Variable, abs(Importance)))  , y = Importance, fill = Sign)) +
  geom_bar(stat = "identity")

workflow() %>% 
  add_recipe(mldata_recipe) %>% 
  add_model(importance_spec) %>% 
  fit(train_data) %>% 
  extract_fit_parsnip() %>% 
  vi(method = "permute", target = has_the_patient_recurred, metrics = , reference = , pred_wrapper = ,
     train = juice(mldata_prep), nsim = 10) #%>% 
# mutate(Importance = ifelse(Sign == "NEG", -Importance, Importance)) %>% 
# ggplot(aes(x= desc(fct_reorder(Variable, abs(Importance)))  , y = Importance, fill = Sign)) +
# geom_bar(stat = "identity")






############################################################################################### AT LAST
# Step finalize Fit : fitting to the whole training and evaluating on the testing data
# With the model of our choice
final_fit <- final_rf %>% 
  last_fit(data_split)


# Step Explore the modelon testing set
# Collect metrics and compare number with the metrics from the training
# Can see if lower or higher...overfit our data, etc
collect_metrics(final_fit)
# Compare to the training previous number
collect_metrics(rf_results) # as a meminder of previous results
# Test data is a littke lower with samll SD

collect_metrics(final_fit) %>% 
  mutate(set = "testing") %>% 
  rename(mean = .estimate) %>% 
  bind_rows(collect_metrics(rf_results) %>% 
              mutate(set = "training")) %>% 
  filter(str_detect(.metric, "accuracy|roc_auc")) %>% 
  ggplot(aes(x= .metric, y=mean, fill = set))+
  geom_bar(stat = "identity",
           position = position_dodge()) + 
  geom_errorbar(aes(xmin = mean - std_err,
                    xmax = mean + std_err),
                width = 0.2, alpha = 0.5)
collect_metrics(final_fit) %>% 
  mutate(set = "testing") %>% 
  rename(mean = .estimate) %>% 
  bind_rows(collect_metrics(rf_results) %>% 
              mutate(set = "training")) %>% 
  filter(str_detect(.metric, "accuracy|roc_auc")) %>% 
  ggplot(aes(x= .metric, y=mean, color = set))+
  geom_point() + 
  geom_errorbar(aes(xmin = mean - std_err,
                    xmax = mean + std_err),
                width = 0.2, alpha = 0.5)


final_fit %>% 
  collect_predictions() %>% 
  conf_mat(truth = has_the_patient_recurred, estimate = .pred_class)
final_fit %>% 
  collect_predictions() %>% 
  conf_mat(truth = has_the_patient_recurred, estimate = .pred_class) %>% 
  autoplot()

last_fit_rf %>%
  collect_predictions() %>% 
  conf_mat(truth = has_the_patient_recurred, estimate = .pred_class) %>% 
  autoplot(type = "heatmap")



# Step Explore the prediction
final_fit %>% collect_predictions() %>% 
  mutate(is_predicton_correct = case_when(
    has_the_patient_recurred == .pred_class     ~ "Cool!",
    TRUE                                        ~ ":("
  )) %>% 
  ggplot(aes(is_predicton_correct))+
  geom_bar() #### scale percent

final_fit %>% collect_predictions() %>% 
  ggplot(aes(has_the_patient_recurred, .pred_class))+
  geom_point() +
  geom_abline()

# Predict on 1 element or new data
predict(final_fit$.workflow[[1]], test_data[15,]) 
# Here to test on the 15th element of the test data by itself






# If choose glmnet, can see estimate on the testing
final_fit_sec <- final_glmnet %>% 
  last_fit(data_split)

final_fit_sec %>%
  pull(.workflow) %>% 
  pluck(1) %>% 
  tidy(exponentiate = TRUE) %>% 
  arrange(desc(abs(estimate))) # %>%
  # kableExtra::kable(digits = 3)

final_fit_sec %>% # WARNING it uses the scaled data
  pull(.workflow) %>% 
  pluck(1) %>% 
  tidy() %>% 
  filter(term != "(Intercept)") %>% 
  arrange(desc(abs(estimate))) %>% 
  filter(abs(estimate) >0) %>% 
  ggplot(aes(estimate, fct_reorder(term, desc(estimate))))+
  geom_vline(xintercept = 0, color = "lightgrey", lty = 2, size = 1.2) +
  # geom_errorbar(aes(xmin = estimate - std.error,
  #                   xmax = estimate + std.error),
  #               width = 0.2, alpha = 0.5)+
  geom_point() + 
  theme_classic2()


collect_metrics(final_fit)
# Compare to the training previous number
collect_metrics(glmnet_results)

# Step Explore the prediction
final_fit_sec %>% collect_predictions() %>% 
  mutate(is_predicton_correct = case_when(
    has_the_patient_recurred == .pred_class     ~ "Cool!",
    TRUE                                        ~ ":("
  )) %>% 
  ggplot(aes(is_predicton_correct))+
  geom_bar()


# # Measure how they performed
# 
# results_train %>% 
#   group_by(model) %>% 
#   rmse(truth = truth, estimate = .prod) # .estimate = rootmean square error is lower so do better
# 
# results_test %>% 
#   group_by(model) %>% 
#   rmse(truth = truth, estimate = .prod) # .est same as train so didn't overfit
# # if rf becore bigger than the one for training


final_fit_sec %>% collect_predictions() %>% 
  mutate(set = "glmnet") %>% 
  bind_rows(final_fit %>% collect_predictions() %>% 
              mutate(set = "rf")) %>% 
  group_by(set) %>% 
  summarize(recurrence = mean(.pred_Recurrence),
            no_recurrence = mean(`.pred_No Recurrence`)) %>% 
  ungroup() %>% 
  pivot_longer(cols = -set, names_to = "status", values_to = "probability") %>% 
  ggplot(aes(x= set, y=probability, fill = status))+
  geom_bar(stat = "identity",
           position = position_dodge())

final_fit_sec %>% collect_predictions() %>% 
  mutate(set = "glmnet") %>% 
  bind_rows(final_fit %>% collect_predictions() %>% 
              mutate(set = "rf")) %>% 
  mutate(is_predicton_correct = case_when(
    has_the_patient_recurred == .pred_class     ~ "Cool!",
    TRUE                                        ~ ":("
  )) %>% 
  group_by(set, is_predicton_correct) %>% 
  summarize(recurrence = mean(.pred_Recurrence),
            no_recurrence = mean(`.pred_No Recurrence`)) %>% 
  ungroup() %>% 
  pivot_longer(cols = -c(set, is_predicton_correct), names_to = "status", values_to = "probability") %>% 
  ggplot(aes(x= is_predicton_correct, y=probability, fill = status))+
  geom_bar(stat = "identity",
           position = position_dodge()) + 
  facet_wrap(. ~ set)


############################################################################################### DECISION TREE
# xgboost_recipe <- 
#   recipe(formula = has_the_patient_recurred ~ ., data = train_data) %>% 
#   step_novel(all_nominal(), -all_outcomes()) %>% 
#   step_dummy(all_nominal(), -all_outcomes(), one_hot = TRUE) %>% 
#   step_zv(all_predictors()) 
library(xgboost)
xgboost_spec <- 
  boost_tree(trees = tune(), min_n = tune(), tree_depth = tune(), learn_rate = tune(), 
             loss_reduction = tune(), sample_size = tune()) %>% 
  set_mode("classification") %>% 
  set_engine("xgboost") 

xgboost_workflow <- 
  workflow() %>% 
  add_recipe(mldata_recipe) %>% 
  add_model(xgboost_spec) 

set.seed(57585)
xgboost_tune <-
  tune_grid(xgboost_workflow, resamples = mldata_folds, grid = 10)

######################################################################### Explore xgboost Tuning Results ----

show_best(xgboost_tune, metric = "roc_auc")
show_best(xgboost_tune, metric = "accuracy")

# Visualize tuned parameters
autoplot(xgboost_tune)
# mtry is the number of predictor randomly selected -> needs to be 
# mim_n is the minimal node size -> needs to be small

xgboost_tune %>% 
  collect_metrics() %>% 
  filter(.metric == "accuracy") %>% 
  select(mean, trees:sample_size) %>% 
  pivot_longer(trees:sample_size,
               values_to = "value",
               names_to = "param") %>% 
  ggplot(aes(value, mean, color = param)) +
  geom_point(show.legend = FALSE) +
  facet_wrap( . ~ param, scales = "free_x")

xgboost_tune %>% 
  collect_metrics() %>% 
  filter(.metric == "roc_auc") %>% 
  select(mean, trees:sample_size) %>% 
  pivot_longer(trees:sample_size,
               values_to = "value",
               names_to = "param") %>% 
  ggplot(aes(value, mean, color = param)) +
  geom_point(show.legend = FALSE) +
  facet_wrap( . ~ param, scales = "free_x")

xgboost_tune %>% select_best("roc_auc")

final_xgboost <- xgboost_workflow %>% 
  finalize_workflow(select_best(xgboost_tune, "roc_auc"))
final_xgboost


# Calculate Performance Metrics again with our last tuned model
#Plot the ROC curve
xgboost_results <- final_xgboost %>% 
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
    control = control_resamples(save_pred = TRUE)
  )



# OR SIMPLE DECISION TREE

tree_spec <- decision_tree(
  cost_complexity = tune(),
  tree_depth = tune(),
  min_n = tune()) %>% 
  set_engine("rpart") %>% 
  set_mode("classification")

tree_workflow <- 
  workflow() %>% 
  add_recipe(mldata_recipe) %>% 
  add_model(tree_spec) 

tree_grid <- grid_regular(cost_complexity(),
                          tree_depth(),
                          min_n(), levels = 4)

tree_tune <- tree_workflow %>% 
  tune_grid(
  resamples = mldata_folds,
  grid = tree_grid,
  metrics = metric_set(roc_auc, accuracy, sensitivity, specificity)
)

# Explore
head(tree_tune %>% 
  collect_metrics())

autoplot(tree_tune) + theme_light()

tree_tune %>% select_best("roc_auc")

final_tree <- tree_workflow %>% 
  finalize_workflow(select_best(tree_tune, "roc_auc"))
final_tree

tree_results <- final_tree %>% 
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity),
    control = control_resamples(save_pred = TRUE)
  )

############################################################################################### EVALUATE MODELS ----
# Explore Performance Metrics
collect_metrics(rf_results)
collect_metrics(glmnet_results)
collect_metrics(xgboost_results)
collect_metrics(tree_results)

# Accuracy is llow
# Sensitivity has 50% chance finding the minority class pop

rf_results %>% 
  conf_mat_resampled()
glmnet_results %>% 
  conf_mat_resampled()
xgboost_results %>% 
  conf_mat_resampled()
tree_results %>% 
  conf_mat_resampled()

############################################################################################### DRAW TREE

# plot the first tree
xgb.plot.tree(model = xgboost_results, trees = 1)

xgb.plot.multi.trees(
  xgboost_workflow,
  feature_names = NULL,
  features_keep = 5,
  plot_width = NULL,
  plot_height = NULL,
  render = TRUE
)


final_fit$

# Fitting
final_fit <- last_fit(final_xgboost, data_split)

write_rds(final_fit)

library(vip)

final_fit %>% 
  vip(geom = "col", aesthetics = list(fill = , alpha = ))

library(parttree)

ex_fit <- fit(final_tree, formula but with only 2 predictor, train_data)
train_data %>% 
  ggplot(aes(x= , y = )) +
  geom_partree(data = ex_fit, aes(fill(outcome), alpha = 0.3)) +
  geom_point(alpha = 0.7, aes(color = outcome)) +
  scale_color_viridis_c(aesthetics = c("color", "fill"))

# testing data
collect_metrics(final_rs)
collect_predictions(final_rs) %>% 
  ggplot(aes(outcome, .pred)) +
  geom_abline(slope = 1) +
  geom_point()



















