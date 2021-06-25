# Load packages
library(tidyverse)
library(tidymodels)
# library(broom.mixed) # for converting bayesian models to tidy tibbles
# library(dotwhisker)  # for visualizing regression results

############################################################################### I ### Exploratory Data Analysis
path <- fs::path("", "Volumes", "Peres_Research", "Ovarian - Radiomics")

# No NA in features

# Clinical
clinical <- read_rds("clinical.rds")

# recurrence by age or year
library(GGally)

clinical %>% 
  select(has_the_patient_recurred_, year_of_diagnosis, treatment_type, debulking_status,
         age_at_Dx, primary_site) %>% 
  ggpairs(columns = 2:ncol(.), aes(color = has_the_patient_recurred_, alpha = 0.5))

clinical %>% 
  select(has_the_patient_recurred_, treatment_type, debulking_status,
         primary_site, year_of_diagnosis, age_at_Dx) %>% 
  mutate(year_of_diagnosis = as.character(year_of_diagnosis)) %>%
  mutate(age_at_Dx = as.character(age_at_Dx)) %>%
  pivot_longer(2:ncol(.)) %>% 
  ggplot(aes(y = value, fill = has_the_patient_recurred_)) +
  geom_bar(position = "fill")+
  facet_wrap(vars(name), scales = "free")+
  labs(x= NULL, y = NULL)

clinical %>% 
  group_by(year_of_diagnosis) %>% 
  summarise(recurred = sum(rec_event)) %>%
  ggplot(aes(year_of_diagnosis, recurred))+
  geom_line()

library(naniar)

clinical %>% 
  select(has_the_patient_recurred_, ecog_pretrt, ecog_posttrt, ecog_recurr, 
         Ethnicity, grade_differentiation) %>%
  gg_miss_upset()
# Will do unknown, remove ecog

clinical %>% 
  select(has_the_patient_recurred_, germline_brca1_mutation, germline_brca2_mutation, somatic_brca1_mutation,
         somatic_brca2_mutation, any_unclassified_brca_mutation) %>%
  gg_miss_upset()
# Will do unknown, remove any_unclassified_brca_mutation

clinical %>% 
  select(hypertension, diabetes_mellitus, hypercholesterolemia,
         chronic_kidney_disease, cardiac_conditions_including_bu) %>%
  gg_miss_upset()

############################################################################### II ### Data prepping
# remove variables non meaningful or redundant to recurrence
clinical_ml <- read_rds("clinical.rds") %>% 
  select(-c(Gender, subject_number, '_tnm_edition_number_must_use_', baseline_ctscan_outside_moffitt, 
         date_of_birth, any_unclassified_brca_mutation,
         comment_for_cardiac_comorbidity,
         has_the_patient_recurred_after_surg,
         months_at_first_neoadjuvant_chem, months_at_first_adjuvant_chem, months_at_first_chemo, 
         months_at_first_surgery, age_at_surgery, age_at_first_recurrence, month_at_first_recurrence_Dx, 
         months_at_surg_followup, months_at_neo_followup, months_at_chem_followup, months_of_surg_rec_free, 
         months_of_neo_rec_free, months_of_chem_rec_free,
         months_of_dx_rec_free, recurrence_time, months_of_treat_rec_free, 
         os_time, rec_event, os_event, 
         has_the_patient_recurred_after_surg
  )) %>% 
  
  select(-c(ecog_pretrt, ecog_posttrt, ecog_recurr, 
            
            germline_brca1_mutation, germline_brca2_mutation, somatic_brca1_mutation,
            somatic_brca2_mutation,

            contrast_enhancement, # What is the difference with w_w...
            tumor_sequence_number, 
            # Include as year, month?
            date_of_first_neoadjuvant_chemot,
            date_of_first_surgery, date_of_first_adjuvant_chemother,
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
  left_join(., clinical_ml, by = "mrn")

# Explore what will need to be changed
skimr::skim(mldata)

# Cleaning
# rm(concordance, stable_features)


############################################################################### II ### Build model
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

library(themis)
# 2. Data pre processing and features engineering + imputation
# Recipe
data_recipe <-
  # 1.model formula
  recipe(has_the_patient_recurred_ ~ ., data = train_data)  %>% 
  # 2.keep these variables but not use them as either outcomes or predictors
  # Keep id but not include in the model as predictor
  update_role(mrn, new_role = "ID") %>% 
  # remove variables that contain only a single value.
  step_zv(all_predictors()) %>% 
  # 3.If factor with too much levels, collapse lower levels
  step_other(Histology, threshold = 0.05) %>% 
  # data_recipe %>% prep() %>% juice() %>% count(Histology)
  
  # 4.Imputation
  step_unknown(all_nominal_predictors()) %>% 
  # data_recipe %>% prep() %>% juice() %>% count(grade_differentiation)
  
  # 5.change all factor to dummy variables for model that cannot handle factor variables
  step_dummy(all_nominal(), -all_outcomes()) #%>%
  
  # Feature engineering on dates
  # step_date(all_of(meaningful_dates), features = c("year", "month")) %>% 
  # step_rm(meaningful_dates)
  
  # LAST.For imbalance, model memorize the few example and
  # step_smote(Race) # Use nearest neighbor to create new synthetic observation almost similar 

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
  add_recipe(data_recipe) %>% # add unpreped recipe, is an unfit model at the end of this line
  add_model(tune_spec)

############################################################################### II ### Data Tuning 
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
)# May need to add specificity and sensitivity to metrics because we have rare events. these 2 will tell us how the model did for our positive and negative cases
# If sens is low it means the model had a really hard time finding the rare case (could mean step_smor is bad idea for this model)



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





############################################################################### II ### Machine Learning




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



############################################################################### II ### Model Evaluation 


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

# Measure how they performed

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








