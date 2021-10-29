library(tidyverse)
library(tidymodels)
#############
clinical <- read_rds("/Users/colinccm/Documents/GitHub/Peres/MilesForMoffittRadiomics/clinical.rds") %>% 
  select(mrn, 
         has_the_patient_recurred, vital_new, 
         rec_event, recurrence_time, os_event, os_time,
         age_at_diagnosis, tnm_cs_mixed_group_stage, 
         treatment_type, debulking_status,
         raceeth)

# concordance <- read_rds("/Users/colinccm/Documents/GitHub/Peres/MilesForMoffittRadiomics/concordance1.rds") %>% 
#   filter(value_CCC >= 0.95) %>% 
#   select(name) %>% 
#   mutate(stable_features = str_match(name, "([a-z][:digit:]*)_*")[,2])
concordance <- read_rds("/Users/colinccm/Documents/GitHub/Peres/MilesForMoffittRadiomics/concordance.rds") %>% 
  mutate(stable_features = str_match(namee, "([a-z][:digit:]*)_*")[,2])

stable_features <- paste0(paste(#"nor_", 
  concordance$stable_features, "[a-z]", sep = ""), collapse = "|")

mldata <- read_rds("/Users/colinccm/Documents/GitHub/Peres/MilesForMoffittRadiomics/radiomics.rds") %>% 
  select(mrn, contrastenhancementyn, matches(stable_features)) %>% 
  right_join(., clinical %>% 
               select(mrn, treatment_type, recurrence_time, rec_event),
             by = "mrn") %>% 
  filter(!is.na(.[2])) %>% 
  filter(contrastenhancementyn == "yes") %>% select(-contrastenhancementyn) %>% 
  select(-mrn)

set.seed(123)
adj <-  mldata %>% filter(treatment_type == "Upfront Surgery") %>% select(-treatment_type)
neo <- mldata %>% filter(treatment_type == "Upfront Neoadjuvant") %>% select(-treatment_type)
data_split <- initial_split(adj, prop = 3/4)
# Create training and testing datasets:
train_data <- training(data_split)
test_data  <- testing(data_split)
################
library(survival)
library(rpart)
library(partykit)
library(rattle)
library(ROCR)
library(survivalROC)

# Create 1 model for each treatment
adj_fit <- rpart(Surv(recurrence_time, rec_event) ~ . , data = adj,
                 control = rpart.control(cp = 0.03))
neo_fit <- rpart(Surv(recurrence_time, rec_event) ~ . , data = neo,
                 control = rpart.control(cp = 0.03))
adj_train_fit <- rpart(Surv(recurrence_time, rec_event) ~ . , data = train_data,
                 control = rpart.control(cp = 0.03))
printcp(adj_fit)
# printcp(adj_train_fit)
# printcp(neo_fit)
rpart.plot::rpart.plot(adj_fit)
fancyRpartPlot(adj_fit)
fancyRpartPlot(neo_fit)
fancyRpartPlot(adj_train_fit)

tneo_fit <- as.party(neo_fit) #Transfer "rpart" to "party"
tadj_fit <- as.party(adj_fit)

# predtree<-predict(tneo_fit,newdata=mldata,type="prob") #Prediction
# predtree<-predict(tadj_fit,newdata=mldata,type="prob") #Prediction
# predROCR <- prediction(predict(tneo_fit, newdata = mldata, type = "prob")[, 2],labels=Surv(recurrence_time, rec_event))
predict(tadj_fit, newdata = adj, type = "prob")[1]
predict(tneo_fit, newdata = neo, type = "prob")[[1]]

# tree2 = neo_fit
adj_fit[["frame"]]$var
adj_fit[["where"]]
# adj_fit$frame$yval = as.numeric(rownames(adj_fit$frame))
# neo_fit$frame$yval = as.numeric(rownames(neo_fit$frame))

#Get the survival function value of each sample
Surv_value = data.frame(predict(adj_fit, newdata=adj,type = "matrix"))[,1]
Out_adj=data.frame()
Out_adj=cbind(adj,Surv_value)

Surv_value = data.frame(predict(adj_fit, newdata=neo,type = "matrix"))[,1]
Out_neo=data.frame()
Out_neo=cbind(neo,Surv_value)

Surv_value = data.frame(predict(adj_train_fit, newdata=test_data,type = "matrix"))[,1]
Out_test=data.frame()
Out_test=cbind(test_data,Surv_value)
Surv_value = data.frame(predict(adj_train_fit, newdata=neo,type = "matrix"))[,1]
Out_neo_test=data.frame()
Out_neo_test=cbind(neo,Surv_value)
#ROC
roc_adj=survivalROC(Stime=Out_adj$recurrence_time, status=Out_adj$rec_event, 
                    marker = Out_adj$Surv_value, predict.time =365, method="KM")
roc_neo=survivalROC(Stime=Out_neo$recurrence_time, status=Out_neo$rec_event, 
                    marker = Out_neo$Surv_value, predict.time =365, method="KM")
roc_test=survivalROC(Stime=Out_test$recurrence_time, status=Out_test$rec_event, 
                    marker = Out_test$Surv_value, predict.time =365, method="KM")
roc_neo_test=survivalROC(Stime=Out_neo_test$recurrence_time, status=Out_neo_test$rec_event, 
                    marker = Out_neo_test$Surv_value, predict.time =365, method="KM")
roc_adj$AUC #Get the AUC of ROC plot
roc_neo$AUC
roc_test$AUC
roc_neo_test$AUC

wilcox.test(roc_neo$AUC, roc_adj$AUC)
wilcox.test(roc_test$AUC, roc_neo_test$AUC)

# plot(roc_neo$FP, roc_neo$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#f8766d", 
#      xlab="False positive rate", ylab="True positive rate",
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# plot(roc_adj$FP, roc_adj$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#f8766d", 
#      xlab="False positive rate", ylab="True positive rate",
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)

ADJ <- tibble(true_positive_rate = roc_adj$TP,
              false_positive_rate = roc_adj$FP)
NEO <- tibble(true_positive_rate = roc_neo$TP,
              false_positive_rate = roc_neo$FP)
NEO %>% 
  ggplot(aes(x= false_positive_rate, y= true_positive_rate))+
  geom_line(color = "blue")+
  geom_line(data = ADJ, color = "red")+
  geom_abline(intercept = 0, slope = 1, linetype = 3)+
  theme_minimal()

ks.test(NEO$false_positive_rate,ADJ$false_positive_rate)

TEST <- tibble(true_positive_rate = roc_test$TP,
              false_positive_rate = roc_test$FP)
NEO_TEST <- tibble(true_positive_rate = roc_neo_test$TP,
              false_positive_rate = roc_neo_test$FP)
NEO_TEST %>% 
  ggplot(aes(x= false_positive_rate, y= true_positive_rate))+
  geom_line(color = "blue")+
  geom_line(data = TEST, color = "red")+
  geom_abline(intercept = 0, slope = 1, linetype = 3)+
  theme_minimal()

ks.test(TEST$false_positive_rate,NEO_TEST$false_positive_rate)

# A prognostic biomarker is a clinical or biological characteristic that provides 
# information on the likely patient health outcome (e.g. disease recurrence) 
# irrespective of the treatment. 
# On the other hand, a predictive biomarker indicates the likely benefit to 
# the patient from the treatment, compared to their condition at baseline (Ruberg and Shen, 2015).


# Keep risk group from rpart to plot cox
adj <- bind_cols(adj, risk_group = adj_fit$where) %>% 
  mutate(risk_group = case_when(
    risk_group == 3          ~ "low",
    risk_group == 4          ~ "mid",
    risk_group == 5          ~ "high"
  ))

km <- survfit(Surv(recurrence_time, rec_event) ~ risk_group, adj) 
plot(km, lty = 1:4, mark.time = FALSE,
     xlab = "Years", ylab = "Progression")
legend(10, 0.3, paste('node', c(4,5,6)), lty = 1:3)


Surv(recurrence_time, rec_event) ~ .

##### COX

# data_split <- initial_split(mldata, prop = 3/4)
# # Create training and testing datasets:
# train_data <- training(data_split)
# test_data  <- testing(data_split)

cph_mod <- coxph(Surv(recurrence_time, rec_event) ~ ., data = adj %>% select(-risk_group), x = TRUE)

myplot <- survfit(Surv(time = adj$recurrence_time, event = adj$rec_event)~risk_group, data = adj)
ggsurvplot(myplot)
summary(cph_mod)

# predict(cph_mod, type = "expected", newdata = test_data)
# predict(cph_mod, type = "risk", newdata = test_data)

plot(survfit(cph_mod, newdata = adj))
plot(survfit(cph_mod, newdata = neo))
plot(predict(cph_mod, type = "expected", newdata = test_data))
survfit(cph_mod, newdata = test_data)


# Both restricted mean survival time (RMST) and survival function probabilities at specific times can be predicted with summary.survit():
cph_mod <- coxph(Surv(recurrence_time, rec_event) ~ ., data =train_data, x = TRUE)
cph_surv <- summary(survfit(cph_mod, newdata = test_data)#, times = .times
                    ) 
plot(survfit(cph_mod, newdata = test_data))

# Event time prediction (RMST)
cph_surv$table
# Survival probabilities
cph_surv 
# Survival probabilities can also be predicted using the pec package;
library(pec)
# predictSurvProb(cph_mod, newdata = test_data, times = .times)
predictSurvProb(cph_mod, newdata = test_data, times = c(1, 50, 100))


#validation
library(survcomp)
cindec_val <- concordance.index(pred_val, surv_time = validation_data$survival, surv_event = event, method = "noether")






############################################################ Surv KM risk group
library(survival)
library(rpart)
# z <- Surv(mldata$recurrence_time, mldata$rec_event)
# z
# fitsurv <- rpart::rpart(z ~ f65flatness + f80com_x_pxl + f85com_z_mm, data = mldata, method = "exp")
# fitsurv
# plot(fitsurv)
# text(fitsurv)
# 
# plot(fitsurv)
# text(fitsurv, use.n = TRUE, cex = 0.75, all = TRUE)
# 
# rpart.plot::rpart.plot(fitsurv)
# 
# path.rpart(fitsurv, node =1)??????
# 
# plotcp(fitsurv)
# # z is dependent variable
# fitsurv$where
# km <- survfit(z ~ fitsurv$where, data = mldata)
# km
# plot(km)


library(party)
# fitctree <- ctree(z ~ f65flatness + f80com_x_pxl + f85com_z_mm, data = mldata)
# fitctree
# plot(fitctree)
# print(fitctree)
# where(fitctree)
# stree <- treeresponse(fitctree)


# RPART
## cp complexity param use to give constraint of overfitting only give the split that are < of the value
# minsplit 
rpart::rpart(control = rpart.control( cp= , minsplit = , minbucket = , maxdepth = , xval = 10))
plotcp
# when under the dotted line less error so picked one of that
printcp
# look at x error because is calculated on the average cross val
# take the lowest xerror+ sd (first standard dev rule), take the highest xerror that are lower than that


# Want to get all patient dataset in one of the final node
# https://www.youtube.com/watch?v=1rNclbWruI0
# https://stackoverflow.com/questions/23924051/find-the-data-elements-in-a-data-frame-that-pass-the-rule-for-a-node-in-a-tree-m?rq=1

summary(fitsurv)
rpart.plot::rpart.plot(extra = 3)



# for loop to choose best cp
# https://gist.github.com/rudeboybert/6468b2e7b929b6241bfb8da7b4927a85
# library(tidyverse)
# library(rpart)
# library(Metrics)
# 
# # Reload house prices data
# train <- read_csv("https://rudeboybert.github.io/SDS293/static/train.csv")
# test <- read_csv("https://rudeboybert.github.io/SDS293/static/test.csv")
# 
# # Set number of folds
# k <- 5
# 
# # Randomly set k folds to training data
# train <- train %>% 
#   sample_frac(size = 1) %>% 
#   mutate(fold = rep(1:k, length = n())) %>% 
#   arrange(fold)
# 
# cp_values_grid <- seq(from = 0, to = 0.0015, len = 101)
# error_estimates <- rep(0, times = length(cp_values_grid))
# 
# error_estimate_per_fold <- rep(0, k)
# 
# for(j in 1:length(cp_values_grid )){
#   
#   current_cp_value <- cp_values_grid[j]
#   
#   for(i in 1:k){
#     train_cv <- train %>% 
#       filter(fold != i)
#     test_cv <- train %>% 
#       filter(fold == i)
#     
#     # Fit model:
#     trained_model <- rpart(SalePrice ~ GrLivArea + HalfBath + YearBuilt, 
#                            data = train_cv,
#                            control = rpart.control(cp = current_cp_value))
#     
#     # Get predictions
#     y_hat <- predict(trained_model, type="vector", newdata = test_cv)
#     
#     # Get error
#     error_estimate_per_fold[i] <- rmsle(actual = test_cv$SalePrice, predicted = y_hat)
#     
#   }
#   error_estimates[j] <- mean(error_estimate_per_fold)
# }
# 
# blah <- tibble(
#   cp_value = cp_values_grid,
#   error_estimate = error_estimates
# )
# ggplot(blah, aes(x = cp_value, y = error_estimate)) +
#   geom_point() +
#   labs(x = "Complexity parameter", y = "Estimate of RMSLE")
# 
# 
# 
# # Bonus: Use optimal complexity parameter value to make submissions on Kaggle
# # Since there are multiple cp values that yield the lowest estimated RMSLE, use
# # the smallest value since it yields the least complex tree.
# cp_star <- blah %>% 
#   arrange(error_estimate, cp_value) %>% 
#   slice(1) %>% 
#   pull(cp_value)
# 
# # Fit/train model on all training data
# trained_model_all <- rpart(SalePrice ~ GrLivArea + HalfBath + YearBuilt, 
#                            data = train,
#                            control = rpart.control(cp = cp_star))
# 
# # Visualize this tree:
# plot(trained_model_all, margin = 0.25)
# text(trained_model_all, use.n = TRUE)
# title("Classification & Regression Tree")
# box()
# 
# # Predict on test set
# test <- test %>% 
#   mutate(SalePriceHat = predict(trained_model_all, type="vector", newdata = test))
# 
# # Write predictions to csv following exact format required by Kaggle here
# # https://www.kaggle.com/c/house-prices-advanced-regression-techniques/submit
# test %>% 
#   select(Id, SalePrice = SalePriceHat) %>% 
#   write_csv("submission.csv")
# 
# # This yields a RMSLE of 0.22065! 

























