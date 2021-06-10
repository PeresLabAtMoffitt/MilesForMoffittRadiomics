# install.packages("tidymodels")

library(tidymodels)
# library(broom.mixed) # for converting bayesian models to tidy tibbles
# library(dotwhisker)  # for visualizing regression results

radiomics <- read_rds("radiomics.rds")

set.seed(1234)

# 1. Preprocessing the data
data %>% 
  na.omit() %>% 
  # it is important that the outcome variable for training a (logistic) regression model is a factor.
  mutate_if(is.character, as.factor) # 



# 2. Splitting the data
# 3/4 of the data into the training set but split evenly winthin race
data_split <- initial_split(data, prop = 3/4, strata = race)

# Create data frames for the two sets:
train_data <- training(data_split)
test_data  <- testing(data_split)

# Recipe
data_recipe <- 
  recipe(recurrence ~ ., data = train_data)  %>% 
  # keep these variables but not use them as either outcomes or predictors
  update_role(mrn, date, lesionid, new_role = "ID") 

summary(data_recipe)


# Do we need dummy variable aka %>% 
# step_dummy(all_nominal(), -all_outcomes())







