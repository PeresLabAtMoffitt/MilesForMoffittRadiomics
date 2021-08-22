# Import packages
library(tidyverse)
# library(scales)
library(lubridate)
library(data.table)
library(gtsummary)
library(gplots)
library(heatmap.plus)
library(RColorBrewer)
library(psych)
# library(corrplot)
library(ggcorrplot)
library(survival)
library(survminer)

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
################################################################################################# I ### Load data----
# features <- read_csv(paste0("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/Ovarian_Radiomics_Features_01062021.csv")) %>%
#   `colnames<-`(str_remove(colnames(.), "F[0-9]*\\:") %>% tolower() %>% gsub("[()^]", "",.) %>% str_replace_all(., " ", "_")) # %>% 
  # mutate(rad = "radiom") %>%
  # group_by(mrn) %>%
  # mutate(id = cur_group_id()) %>% # %>% str_replace_all(., "\\^", "")
  # ungroup() %>%
  # mutate(len = nchar(id)) %>%
  # mutate(zero = 5 - len) %>%
  # mutate(ii = stringi::stri_dup("0", .$zero)) %>%
  # select(c(rad, ii, id, mrn, everything())) %>%
  # unite(patient_id, rad:id, sep = "") %>%
  # select(c(patient_id, mrn, everything(), -c(len, zero)))

# ID_linkage <- features %>% select(c(mrn, patient_id)) %>% distinct(patient_id, .keep_all = TRUE)
# write_csv(ID_linkage, "ID_linkage.csv")
# ID_linkage <- read_csv("ID_linkage.csv")
# features <- features %>%
#   left_join(ID_linkage, ., by = "mrn") %>% select(-mrn)
# colnames(features)[54]

path <- fs::path("", "Volumes", "Peres_Research", "Ovarian - Radiomics")
features <-
  read_csv(paste0(path, "/data/analysis dataset/Clinical_normradiomics_08182021.csv")) %>% 
  drop_na(lesionid) %>% 
  filter(isbiggesttumor == 1) %>%
  select(mrn, contrastenhancementyn, lesionid, matches("^f[0-9]"))
contrast_info <- features %>% 
  select(mrn, contrastenhancementyn)
  

clinical <- readxl::read_xlsx(
  paste0(path,
         "/data/analysis dataset/radiomic_CT_final_v3.xlsx"))


################################################################################################# II ### Features----
summary(features[4])
class(features) <- "data.frame"
# scale data from -1 to 1
for(i in 2:length(colnames(features))) {
  if(class(features[,i]) == "numeric" | class(features[,i]) == "integer") {
    features[,i] <- scales::rescale(features[,i], to=c(-1,1)) }
}
summary(features[4])


################################################################################################# III ### Clinical----
clinical <- clinical %>% 
  janitor::clean_names() %>% 
  # mutate(mrn = as.character(mrn)) %>% 
  `colnames<-`(
    str_remove(colnames(.), "_cancer_registry") %>% 
      tolower() %>% 
      str_replace_all(., "__", "_")
  ) %>% 
  mutate(across(where(is.character), .fns = ~ tolower(.))) %>%
  mutate(across(where(is.character), .fns = ~ capwords(.))) %>%
  mutate(across(where(is.character), .fns = ~ str_replace(., "NANA", NA_character_))) %>%
  # Remove low grade - well differentiated cases
  filter(grade_differentiation != "Well Differentiated" | is.na(grade_differentiation)) %>% 
  filter(str_detect(treatment_type, "Upfront")) %>%
  # rename(race = "race", ethnicity = "ethnicity") %>%
  mutate(race = case_when(
    race %in% 
      c("Other", "Asian", "Pacif", "Unko", 
        "Filip", "Pakis")                                ~ "Other",
    race == "Unkno"                                      ~ "Unknown",
    TRUE                                                 ~ race
  )) %>% 
  mutate(ethnicity = case_when(
    ethnicity == "Non-spanish"                           ~ "Non-Hispanic",
    ethnicity %in% 
      c("Spanish Nos", "Puerto Rican", 
        "Mexican", "Cuban", "South/centra")              ~ "Hispanic",
    TRUE                                                 ~ ethnicity
  )) %>% 
  mutate(raceeth = case_when(
    race == "White" &
      ethnicity == "Non-Hispanic"        ~ "White Non-Hispanic",
    race == "Black" &
      ethnicity == "Non-Hispanic"        ~ "Black Non-Hispanic",
    ethnicity == "Hispanic"              ~ "Hispanic",
    race %in% 
      c("Other", "Unknown") |
      ethnicity == "Unknown"             ~ "Other/Unknown"
  )) %>% 
  mutate(raceeth = factor(raceeth, 
                          levels = c("White Non-Hispanic", "Black Non-Hispanic", 
                                     "Hispanic","Other/Unknown")) ) %>% 
  mutate(tnm_cs_mixed_group_stage = factor(tnm_cs_mixed_group_stage)) %>% 
  mutate(preDx_hypertension = case_when(
    str_detect(hypertension, "pre|Pre|both")                      ~ "Yes",
    TRUE                                                          ~ "No"
  )) %>%
  mutate(postDx_hypertension = case_when(
    str_detect(hypertension, "post|Post|both")                    ~ "Yes",
    TRUE                                                          ~ "No"
  )) %>%
  mutate(preDx_diabetes_mellitus = case_when(
    str_detect(diabetes_mellitus, "pre|Pre|both")                 ~ "Yes",
    TRUE                                                          ~ "No"
  )) %>%
  mutate(postDx_diabetes_mellitus = case_when(
    str_detect(diabetes_mellitus, "post|Post|both")               ~ "Yes",
    TRUE                                                          ~ "No"
  )) %>%
  mutate(preDx_hypercholesterolemia = case_when(
    str_detect(hypercholesterolemia, "pre|Pre|both")              ~ "Yes",
    TRUE                                                          ~ "No"
  )) %>%
  mutate(postDx_hypercholesterolemia = case_when(
    str_detect(hypercholesterolemia, "post|Post|both")            ~ "Yes",
    TRUE                                                          ~ "No"
  )) %>%
  mutate(preDx_chronic_kidney_disease = case_when(
    str_detect(chronic_kidney_disease, "pre|Pre|both")            ~ "Yes",
    TRUE                                                          ~ "No"
  )) %>%
  mutate(postDx_chronic_kidney_disease = case_when(
    str_detect(chronic_kidney_disease, "post|Post|both")          ~ "Yes",
    TRUE                                                          ~ "No"
  )) %>%
  mutate(preDx_cardiac_conditions = case_when(
    str_detect(cardiac_conditions_including_bu, "pre|Pre|both")   ~ "Yes",
    TRUE                                                          ~ "No"
  )) %>%
  mutate(postDx_cardiac_conditions = case_when(
    str_detect(cardiac_conditions_including_bu, "post|Post|both") ~ "Yes",
    TRUE                                                          ~ "No"
  )) %>%
  
  mutate(preDx_comorbidities = case_when(
    str_detect(hypertension, "pre|Pre") |
      str_detect(diabetes_mellitus, "pre|Pre") |
      str_detect(hypercholesterolemia, "pre|Pre") |
      str_detect(chronic_kidney_disease, "pre|Pre") |
      str_detect(cardiac_conditions_including_bu, "pre|Pre")       ~ "Comorbidities",
    TRUE                                                           ~ "No comorbidities"
  )) %>%
  mutate(debulking_status =
           str_remove(debulking_status,
                      " \\(.*")) %>%
  # mutate(any_germline_brca_mutation = case_when(
  #   germline_brca1_mutation == "Yes" |
  #     germline_brca2_mutation == "Yes"       ~ "Yes",
  #   germline_brca1_mutation == "No"          ~ "No",
  #   germline_brca2_mutation == "No"          ~ "No",
  #   TRUE                                     ~ "Unknown"
  # )) %>% 
  # mutate(any_somatic_brca_mutation = case_when(
  #   somatic_brca1_mutation == "Yes" |
  #     somatic_brca2_mutation == "Yes"        ~ "Yes",
  #   somatic_brca1_mutation == "No"           ~ "No",
  #   somatic_brca2_mutation == "No"           ~ "No",
  #   TRUE                                     ~ "Unknown"
  # ))

  mutate(date_of_first_adjuvant_chemother = as.POSIXct(date_of_first_adjuvant_chemother, format = "%m/%d/%y")) %>% 
  mutate(debulking_status = case_when(
    debulking_status == "Incomplete Records"    ~ NA_character_,
    TRUE                                        ~ debulking_status
  )) %>% 
  # left_join(ID_linkage, ., by = "mrn") %>% 
  select(#-mrn, 
    -c(moffitt_patient, summary_of_rx_1st_course, summary_of_rx_1st_course_at_t, subject_number))

clinical_var <- function(data) {
  data <- data %>% 
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
    mutate(age_at_surgery = round(interval(start = date_of_birth, end = date_of_first_surgery)/ # First or abstracted?
                                    duration(n=1, units = "years"), 2)) %>% 
    mutate(months_at_first_treatment = 
             coalesce(months_at_first_neoadjuvant_chem, months_at_first_surgery, months_at_first_adjuvant_chem)) %>% 
    
    mutate(age_at_first_recurrence = round(interval(start = date_of_birth, end = date_of_first_recurrence)/
                                             duration(n=1, units = "years"), 2)) %>% 
    mutate(month_at_first_recurrence_Dx = round(interval(start = date_of_diagnosis, end = date_of_first_recurrence)/
                                                  duration(n=1, units = "months"), 2)) %>% 
    # For recurrence
    mutate(first_treatment_date = 
             coalesce(date_of_first_neoadjuvant_chemot, date_of_first_surgery#, date_of_first_adjuvant_chemother
                      )) %>% 
    mutate(rec_event_date = coalesce(date_of_first_recurrence, fwdate_most_recent)) %>% 
    
    mutate(recurrence_time = round(interval(start = first_treatment_date, end = rec_event_date)/
                                                  duration(n=1, units = "months"), 2)) %>% 
    mutate(rec_event = ifelse((has_the_patient_recurred == "No"), 0, 1)) %>%
    mutate(has_the_patient_recurred = case_when(
      str_detect(has_the_patient_recurred, "Yes|yes")          ~ "Recurrence",
      str_detect(has_the_patient_recurred, "No|no")           ~ "No Recurrence"
    )) %>%
    mutate(has_the_patient_recurred = factor(has_the_patient_recurred, 
                                              levels = c("Recurrence", "No Recurrence")) ) %>% 
    
    # For survivals
    mutate(os_time = round(interval(start = first_treatment_date, end = fwdate_most_recent)/
                                           duration(n=1, units = "months"), 2)) %>% 
    mutate(os_event = ifelse((vital_new == "Alive"), 0, 1)) %>% 
    
    # Others
    mutate(months_at_surg_followup = round(interval(start = date_of_first_surgery, end = fwdate_most_recent)/
                                             duration(n=1, units = "months"), 2)) %>% 
    mutate(months_at_neo_followup = round(interval(start = date_of_first_neoadjuvant_chemot, end = fwdate_most_recent)/
                                            duration(n=1, units = "months"), 2)) %>% 
    mutate(months_at_chem_followup = round(interval(start = first_chemo_date, end = fwdate_most_recent)/
                                             duration(n=1, units = "months"), 2)) %>% 
    mutate(months_at_treat_followup = round(interval(start = first_treatment_date, end = fwdate_most_recent)/
                                              duration(n=1, units = "months"), 2)) %>% 
    
    mutate(months_of_dx_rec_free = round(interval(start = date_of_diagnosis, end = rec_event_date)/
                                           duration(n=1, units = "months"), 2)) %>% 
    mutate(months_of_surg_rec_free = round(interval(start = date_of_first_surgery, end = rec_event_date)/
                                             duration(n=1, units = "months"), 2)) %>% 
    mutate(months_of_neo_rec_free = round(interval(start = date_of_first_neoadjuvant_chemot, end = rec_event_date)/
                                            duration(n=1, units = "months"), 2)) %>% 
    mutate(months_of_chem_rec_free = round(interval(start = first_chemo_date, end = rec_event_date)/
                                             duration(n=1, units = "months"), 2)) %>% 
    mutate(months_of_treat_rec_free = round(interval(start = first_treatment_date, end = rec_event_date)/
                                              duration(n=1, units = "months"), 2)) %>% 
    mutate(recurrence_date_after_surgery = case_when(
      date_of_first_recurrence > date_of_first_surgery    ~ date_of_first_recurrence,
      TRUE                                                ~ NA_POSIXct_
    )) %>%
    mutate(has_the_patient_recurredafter_surg = ifelse(!is.na(recurrence_date_after_surgery), "Recurrence", "No Recurrence"))
}

clinical <- clinical_var(data = clinical)
write_rds(clinical, "clinical.rds")

################################################################################################# IV ### Bind df----
radiomics <- left_join(features %>% 
                         select(mrn, contrastenhancementyn, 
                                (order(colnames(.)))) ,
                       clinical,
                       by = "mrn") %>% 
  filter(!is.na(has_the_patient_recurred))
write_rds(radiomics, "radiomics.rds")


################################################################################################# V ### New file---- Obsolete 

# merge_data <- 
#   read_csv(paste0("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/merge_clinical_normalized_radiomics_M4MOC_03192021.csv")) %>% 
#   select("mrn", "date", "lesionid", "segmentedby", first(starts_with("nor_")):ncol(.)) %>% 
#   mutate(across(
#     .cols = contains("date"),
#     .fns = ~ as.POSIXct(., format = "%m/%d/%Y"))) 

# final_data <- clinical_cleaning(data = merge_data) %>% 
#   filter(!is.na(lesionid))
# write_rds(final_data, "radiomics.rds")


