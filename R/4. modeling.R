# Load packages
library(tidyverse)
library(tidymodels)
# library(broom.mixed) # for converting bayesian models to tidy tibbles
# library(dotwhisker)  # for visualizing regression results


############################################################################### I ### Data prepping
path <- fs::path("", "Volumes", "Peres_Research", "Ovarian - Radiomics")

clinical <- read_rds("clinical.rds") %>% select(-Gender) %>% 
  select(-c(tnm_cs_mixed_group_stage, ecog_pretrt, ecog_posttrt, ecog_recurr, 
            months_at_first_neoadjuvant_chem, months_at_first_adjuvant_chem, months_at_first_chemo, 
            months_at_first_surgery, age_at_surgery, age_at_first_recurrence, month_at_first_recurrence_Dx, 
            months_at_surg_followup, months_at_neo_followup, months_at_chem_followup, months_of_surg_rec_free, 
            months_of_neo_rec_free, months_of_chem_rec_free, baseline_ctscan_outside_moffitt,
            Ethnicity, grade_differentiation, germline_brca1_mutation, germline_brca2_mutation, somatic_brca1_mutation,
            somatic_brca2_mutation, any_unclassified_brca_mutation,
            contrast_enhancement, w_wo_contrast, chronic_kidney_disease,
            
            # Include as year, month?
            date_of_first_neoadjuvant_chemot,
            date_of_first_surgery, date_of_first_adjuvant_chemother,
            
            months_of_dx_rec_free, recurrence_time, months_of_treat_rec_free, os_time
            ))




# Select column name of stable features in the data
concordance <- read_rds("concordance.rds") %>% 
  filter(value_CCC >= 0.95) %>% 
  select(name) %>% 
  mutate(stable_features = str_match(name, "([a-z][:digit:]*)_*")[,2])

stable_features <- paste0(paste("nor_", concordance$stable_features, "[a-z]", sep = ""), collapse = "|")
# radiomics <- read_csv(
#   paste0(path,
#          "/data/Merge clinical and radiomics data/Ovarian_Normalized_Radiomics_Features_05242021.csv")) #%>% 
#   select(matches(stable_features))

# radiomics <- read_rds("radiomics.rds")
mldata <- read_rds("radiomics.rds") %>% 
  select(mrn, matches(stable_features)) %>% 
  left_join(., clinical, by = "mrn") %>% 
  drop_na(starts_with("nor")) %>% 
  # it is important that the outcome variable for training a (logistic) regression model is a factor.
  mutate_if(is.character, as.factor)

# Explore what will need to be changed
skimr::skim(mldata)
# warning for date variable
meaningful_dates <- 
  c("baseline_ct_scan_date", "date_of_first_neoadjuvant_chemot", 
    "date_of_first_surgery", "date_of_first_adjuvant_chemother")


mldata <- mldata %>% 
  # Clean not meaningful var
  select(everything(), -c(contains("date"))#, 
         # c(baseline_ct_scan_date, date_of_first_neoadjuvant_chemot, date_of_first_surgery, date_of_first_adjuvant_chemother)
         ) %>% 
  select(-c(subject_number, comment_for_cardiac_comorbidity), # tumor_sequence_number, recurrence_time
         TNM = "_tnm_edition_number_must_use_") 

# Cleaning
# rm(concordance, stable_features)


############################################################################### II ### Build model
set.seed(1234)

# 1. Splitting the data
# 3/4 of the data into the training set but split evenly winthin race
data_split <- initial_split(mldata, prop = 3/4, strata = Race) # w_wo_contrast

# Create data frames for the two sets:
train_data <- training(data_split)
test_data  <- testing(data_split)

# 2. Recipe for data pre processing
# Recipe
data_recipe <-
  recipe(has_the_patient_recurred_ ~ ., data = train_data)  %>% 
  # keep these variables but not use them as either outcomes or predictors
  update_role(mrn, new_role = "ID") %>% 
  # change all factor to dummy variables
  # step_impute_mode(all_nominal(),, all_predictors()) %>% 
  step_dummy(all_nominal(),  -all_outcomes()) # %>%
  
  # Feature engineering on dates
  # step_date(all_of(meaningful_dates), features = c("year", "month")) %>% 
  # step_rm(meaningful_dates)

summary(data_recipe)

# estimate the required parameters from a training set that can be later applied to other data sets.
# learn what the model should be with the training data
mldata_prep <- prep(data_recipe, verbose = TRUE, log_changes = TRUE)
# Extract Finalized Training Set
juiced <- juice(mldata_prep)

# Build model specification
tune_spec <- rand_forest(
  # tune right value for the number of predictors that will be randomly sampled at each split when creating the tree models
  mtry = tune(), 
  trees = 1000,
  # tune right value for the minimum number of data points in a node that are required for the node to be split further.
  min_n = tune()
) %>% 
  set_mode("classification") %>% 
  set_engine("ranger")

tune_wf <- workflow() %>% 
  add_recipe(data_recipe) %>% # add unpreped recipe
  add_model(tune_spec)

# train hyperparameter
set.seed(123)
# 10 fold cross validation
data_folds <- vfold_cv(train_data)

doParallel::registerDoParallel() # Because not patient, parallel processing
set.seed(123)
# 
tune_results <- tune_grid( # will tune mtry and min_m on a grid
  tune_wf, # tune worflow
  resamples = data_folds, # on this data
  grid = 20 # do 20 point
)

tune_results %>% 
  collect_metrics() %>% 
  filter(.metric == "roc_auc") %>% 
  select(mean, min_n, mtry) %>% 
  pivot_longer(min_n:mtry,
               values_to = "value",
               names_to = "param") %>% 
  ggplot(aes(value, mean, color = param)) +
  geom_point(show.legend = FALSE) +
  facet_wrap( . ~ param, scales = "free_x")

tune_results %>% select_best("accuracy")

# Tune more?
# Can make a regular grid
new_grid <- grid_regular(
  mtry(range = c(25, 75)),
  min_n(range = c(0, 25)),
  levels = 10
)

set.seed(123)
sec_tune_results <- tune_grid( # will tune mtry and min_m on a grid
  tune_wf, # tune worflow
  resamples = data_folds, # on this data
  grid = new_grid
)

sec_tune_results %>% collect_metrics() %>% 
  filter(.metric == "roc_auc") %>% 
  mutate(min_n = factor(min_n)) %>% 
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point()

best_auc <- select_best(sec_tune_results, "roc_auc")

final_rf <- finalize_model(tune_spec,
               best_auc)

library(vip)

final_rf %>% 
  set_engine("ranger", importance = "permutation") %>% 
  fit(has_the_patient_recurred_ ~ .,
      data = juice(mldata_prep) %>% select(-mrn)) %>% 
  vip(geom = "point")

final_wf <- workflow() %>% 
  add_recipe(data_recipe) %>% 
  add_model(final_rf)

final_results <- final_wf %>% 
  last_fit(data_split)

final_results %>% 
  collect_metrics()

final_results %>% collect_predictions() %>% 
  mutate(is_predicton_correct = case_when(
    has_the_patient_recurred_ == .pred_class     ~ "Cool!",
    TRUE                                        ~ ":("
  )) %>% 
  ggplot(aes(is_predicton_correct))+
  geom_bar()



############################################################################### II ### Data Tuning 


############################################################################### II ### Machine Learning


  



# # 2. Splitting the data
# # 3/4 of the data into the training set but split evenly winthin race
# data_split <- initial_split(radiomics1, prop = 3/4, strata = c(race, w_wo_contrast))
# 
# # Create data frames for the two sets:
# train_data <- training(data_split)
# test_data  <- testing(data_split)




# Do we need dummy variable aka %>% 
# step_dummy(all_nominal(), -all_outcomes())


# logistic reg
mod <-
  logistic_reg(penalty = 0.01, mixture = 1/3) %>%
  # now say how you want to fit the model and another other options
  set_engine("glmnet", nlambda = 10)
translate(mod, engine = "glmnet")
# Support vector machine
# ramdom forest
rf_mod <- rand_forest(mode = "classification") %>% 
  set_engine("ranger")
rf_mod <- rf_mod %>% #train_data %>% 
  fit(has_the_patient_recurred_ ~ ., data = train_data)
# neural network




# Evaluate model

results_train <- mod %>% 
  predict(new_data = train_data) %>% 
  mutate(truth = train_data$recurrence, model = "glmet") %>% 
  bind_rows(rf_mod %>% 
              predict(new_data = train_data) %>% 
              mutate(truth = train_data$recurrence, model = "rf"))

results_test <- mod %>% 
  predict(new_data = test_data) %>% 
  mutate(truth = train_data$recurrence, model = "glmet") %>% 
  bind_rows(rf_mod %>% 
              predict(new_data = test_data) %>% 
              mutate(truth = train_data$recurrence, model = "rf"))

# Measure ho they performed

results_train %>% 
  group_by(model) %>% 
  rmse(truth = truth, estimate = .prod) # .estimate = rootmean square error is lower so do better

results_test %>% 
  group_by(model) %>% 
  rmse(truth = truth, estimate = .prod) # .est same as train so didn't overfit
# if rf becore bigger than the one for training


results_test %>% 
  mutate(train = "testing") %>% 
  bind_rows(results_train %>% mutate(train = "training")) %>% 
  ggplot(aes(truth, .prod, color = model)) +
  geom_abline(lty = 2, color = "grey80", size = 1.5)+
  geom_point(alpha = 0.5)+
  facet_wrap( ~ train)








