# Import packages
library(tidyverse)
library(scales)
library(lubridate)
library(data.table)
library(gtsummary)
library(gplots)
library(heatmap.plus)
library(RColorBrewer)
# library(corrplot)
library(ggcorrplot)
library(survival)
library(survminer)

################################################################################################# I ### Load data----
radiomics <- read.csv(paste0("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/Ovarian_Radiomics_Features_01062021.csv")) %>%
  `colnames<-`(str_remove(colnames(.), "F[0-9]*\\.") %>% tolower() %>% str_replace_all(., "\\.", "_"))
colnames(radiomics)
clinical <- readxl::read_xls("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/radiomic_CT_final.xls") %>% 
  `colnames<-`(str_remove(colnames(.), "_Cancer_Registry|_CS__Mixed_Group") %>%
                 str_replace_all( "__", "_") %>%
                 tolower()) %>%
  mutate(mrn = as.character(mrn))

################################################################################################# II ### Radiomics----
summary(radiomics)
radiomics <- radiomics %>% mutate(mrn = as.character(mrn))
# scale data from -1 to 1
for(i in 1:length(colnames(radiomics))) {
  if(class(radiomics[,i]) == "numeric" | class(radiomics[,i]) == "integer") {
    radiomics[,i] <- rescale(radiomics[,i], to=c(-1,1)) }
}
summary(radiomics)

################################################################################################# III ### Clinical----
clinical <- clinical %>% 
  mutate(across(.fns = ~ str_replace(., "Unknown|unknown", NA_character_))) %>% 
  # For summary stats
  mutate(age_at_Dx = round(interval(start = date_of_birth, end = date_of_diagnosis)/
                             duration(n=1, units = "years"), 2)) %>% 
  mutate(months_at_first_neoadjuvant_chem = round(interval(start = date_of_diagnosis, end = date_of_first_neoadjuvant_chemot)/
                                                 duration(n=1, units = "months"), 2)) %>% 
  mutate(months_at_first_adjuvant_chem = round(interval(start = date_of_diagnosis, end = date_of_first_adjuvant_chemother)/
                                              duration(n=1, units = "months"), 2)) %>% 
  mutate(months_at_first_chemo = coalesce(months_at_first_neoadjuvant_chem, months_at_first_adjuvant_chem)) %>% 
  mutate(first_chemo_date = coalesce(date_of_first_neoadjuvant_chemot, date_of_first_adjuvant_chemother)) %>% 
  
  mutate(months_at_first_surgery = round(interval(start = date_of_diagnosis, end = date_of_first_surgery)/
                                        duration(n=1, units = "months"), 2)) %>% 
  mutate(age_at_surgery = round(interval(start = date_of_birth, end = date_of_surgery)/
                                  duration(n=1, units = "years"), 2)) %>% 
  mutate(months_at_first_treatment = 
           coalesce(months_at_first_neoadjuvant_chem, months_at_first_surgery, months_at_first_adjuvant_chem)) %>% 
  mutate(first_treatment_date = 
           coalesce(date_of_first_neoadjuvant_chemot, date_of_first_surgery, date_of_first_adjuvant_chemother)) %>% 
  
  mutate(age_at_first_recurrence = round(interval(start = date_of_birth, end = date_of_first_recurrence)/
                                           duration(n=1, units = "years"), 2)) %>% 
  mutate(month_at_first_recurrence = round(interval(start = date_of_diagnosis, end = date_of_first_recurrence)/
                                           duration(n=1, units = "months"), 2)) %>% 
  
  # For survivals
  mutate(os_event = ifelse((vital_new == "Alive"), 0, 1)) %>% 
  mutate(rec_event = ifelse((has_the_patient_recurred_ == "no"), 0, 1)) %>%
  
  mutate(months_at_dx_followup = round(interval(start = date_of_diagnosis, end = fwdate_most_recent)/
                                      duration(n=1, units = "months"), 2)) %>% 
  mutate(months_at_surg_followup = round(interval(start = date_of_first_surgery, end = fwdate_most_recent)/
                                      duration(n=1, units = "months"), 2)) %>% 
  mutate(months_at_neo_followup = round(interval(start = date_of_first_neoadjuvant_chemot, end = fwdate_most_recent)/
                                      duration(n=1, units = "months"), 2)) %>% 
  mutate(months_at_chem_followup = round(interval(start = first_chemo_date, end = fwdate_most_recent)/
                                           duration(n=1, units = "months"), 2)) %>% 
  mutate(months_at_treat_followup = round(interval(start = first_treatment_date, end = fwdate_most_recent)/
                                      duration(n=1, units = "months"), 2)) %>% 
  
  mutate(recurrence_date = coalesce(date_of_first_recurrence, fwdate_most_recent)) %>% 
  mutate(months_of_dx_rec_free = round(interval(start = date_of_diagnosis, end = recurrence_date)/
                                           duration(n=1, units = "months"), 2)) %>% 
  mutate(months_of_surg_rec_free = round(interval(start = date_of_first_surgery, end = recurrence_date)/
                                             duration(n=1, units = "months"), 2)) %>% 
  mutate(months_of_neo_rec_free = round(interval(start = date_of_first_neoadjuvant_chemot, end = recurrence_date)/
                                           duration(n=1, units = "months"), 2)) %>% 
  mutate(months_of_chem_rec_free = round(interval(start = first_chemo_date, end = recurrence_date)/
                                             duration(n=1, units = "months"), 2)) %>% 
  mutate(months_of_treat_rec_free = round(interval(start = first_treatment_date, end = recurrence_date)/
                                             duration(n=1, units = "months"), 2))

 



################################################################################################# IV ### Bind df----

radiomics <- full_join(radiomics, clinical, by = "mrn") %>% 
  filter(!is.na(lesion_id))

