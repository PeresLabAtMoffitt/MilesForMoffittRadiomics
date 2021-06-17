# Load packages

library(tidymodels)
# library(broom.mixed) # for converting bayesian models to tidy tibbles
# library(dotwhisker)  # for visualizing regression results


############################################################################### I ### Data prep
path <- fs::path("", "Volumes", "Peres_Research", "Ovarian - Radiomics")

# 1.1. Clinical
clinical <- read_rds("clinical.rds")





# 1.2. Radiomics features
# Select column name of stable features in the data
concordance <- read_rds("concordance.rds") %>% 
  filter(value_CCC >= 0.95) %>% 
  select(name) %>% 
  mutate(stable_features = str_match(concordance$name, "([a-z][:digit:]*)_*")[,2])

stable_features <- paste0(paste("nor_", concordance$stable_features, "[a-z]", sep = ""), collapse = "|")
# radiomics <- read_csv(
#   paste0(path,
#          "/data/Merge clinical and radiomics data/Ovarian_Normalized_Radiomics_Features_05242021.csv")) #%>% 
#   select(matches(stable_features))

# radiomics <- read_rds("radiomics.rds")
radiomics1 <- read_rds("radiomics.rds") %>% 
  select(mrn, matches(stable_features)) %>% 
  left_join(., clinical, by = "mrn")

# Cleaning
rm(concordance, stable_features)


############################################################################### II ### Machine Learning
set.seed(1234)

# 1. Preprocessing the data
radiomics1 <- radiomics1 %>% 
  drop_na(starts_with("nor")) %>% 
  # it is important that the outcome variable for training a (logistic) regression model is a factor.
  mutate_if(is.character, as.factor)



# 2. Splitting the data
# 3/4 of the data into the training set but split evenly winthin race
data_split <- initial_split(radiomics1, prop = 3/4, strata = c(race, w_wo_contrast))

# Create data frames for the two sets:
train_data <- training(data_split)
test_data  <- testing(data_split)

# Recipe
data_recipe <- 
  recipe(has_the_patient_recurred_ ~ ., data = train_data)  %>% 
  # keep these variables but not use them as either outcomes or predictors
  update_role(mrn, date, lesionid, new_role = "ID") 

summary(data_recipe)


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








