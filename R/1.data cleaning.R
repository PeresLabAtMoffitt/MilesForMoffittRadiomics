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
features <- read_csv(paste0("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/Ovarian_Radiomics_Features_01062021.csv")) %>%
  `colnames<-`(str_remove(colnames(.), "F[0-9]*\\:") %>% tolower() %>% gsub("[()^]", "",.) %>% str_replace_all(., " ", "_")) # %>% 
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
ID_linkage <- read_csv("ID_linkage.csv")
features <- features %>%
  left_join(ID_linkage, ., by = "mrn") %>% select(-mrn)
colnames(features)[54]

path <- fs::path("", "Volumes", "Peres_Research", "Ovarian - Radiomics")

clinical <- readxl::read_xlsx(
  paste0(path,
         "/data/radiomic_CT_final_v3.xlsx"))


################################################################################################# II ### Features----
summary(features)
class(features) <- "data.frame"
# scale data from -1 to 1
for(i in 1:length(colnames(features))) {
  if(class(features[,i]) == "numeric" | class(features[,i]) == "integer") {
    features[,i] <- scales::rescale(features[,i], to=c(-1,1)) }
}
summary(features)


################################################################################################# III ### Clinical----


clinical <- clinical %>% 
  `colnames<-`(
    str_remove(colnames(.), "_cancer_registry") %>% 
      tolower() %>% 
      str_replace_all(., "__", "_")
  ) %>% 
  mutate(across(where(is.character), .fns = ~ tolower(.))) %>%
  mutate(across(where(is.character), .fns = ~ capwords(.))) %>%
  mutate(across(where(is.character), .fns = ~ str_replace(., "NANA", NA_character_))) %>%
  rename(Race = "race", Ethnicity = "ethnicity", Histology = "histology") %>%
  mutate(Race = case_when(
    Race %in% 
      c("Other", "Asian", "Pacif", "Unko", 
        "Filip", "Pakis")                                ~ "Other",
    Race == "Unkno"                                      ~ "Unknown",
    TRUE                                                 ~ Race
  )) %>% 
  mutate(Ethnicity = case_when(
    Ethnicity == "Non-spanish"                           ~ "Non-Hispanic",
    Ethnicity %in% 
      c("Spanish Nos", "Puerto Rican", 
        "Mexican", "Cuban", "South/centra")              ~ "Hispanic",
    TRUE                                                 ~ Ethnicity
  )) %>% 
  mutate(raceeth = case_when(
    Race == "White" &
      Ethnicity == "Non-Hispanic"        ~ "White Non-Hispanic",
    Race == "Black" &
      Ethnicity == "Non-Hispanic"        ~ "Black Non-Hispanic",
    Ethnicity == "Hispanic"              ~ "Hispanic",
    Race %in% 
      c("Other", "Unknown") |
      Ethnicity == "Unknown"             ~ "Other/Unknown"
  )) %>% 
  
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
  # Remove low grade - well differentiated cases
  filter(grade_differentiation != "Well Differentiated" | is.na(grade_differentiation)) %>% 
  # left_join(ID_linkage, ., by = "mrn") %>% 
  select(#-mrn, 
    -c(complete_, moffitt_patient, summary_of_rx_1st_course, summary_of_rx_1st_course_at_t, subject_number))

clinical_cleaning <- function(data) {
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
             coalesce(date_of_first_neoadjuvant_chemot, date_of_first_surgery, date_of_first_adjuvant_chemother)) %>% 
    mutate(rec_event_date = coalesce(date_of_first_recurrence, fwdate_most_recent)) %>% 
    
    mutate(recurrence_time = round(interval(start = first_treatment_date, end = rec_event_date)/
                                                  duration(n=1, units = "months"), 2)) %>% 
    mutate(rec_event = ifelse((has_the_patient_recurred_ == "no"), 0, 1)) %>%
    mutate(has_the_patient_recurred_ = case_when(
      str_detect(has_the_patient_recurred_, "Yes|yes")          ~ "Recurrence",
      str_detect(has_the_patient_recurred_, "No|no")           ~ "No Recurrence"
    )) %>%
    mutate(has_the_patient_recurred_ = factor(has_the_patient_recurred_, 
                                              levels = c("Recurrence", "No Recurrence")) ) %>% 
    
    # For survivals
    mutate(os_time = round(interval(start = date_of_diagnosis, end = fwdate_most_recent)/
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
    mutate(has_the_patient_recurred_after_surg = ifelse(!is.na(recurrence_date_after_surgery), "Recurrence", "No Recurrence"))
}

clinical <- clinical_cleaning(data = clinical)
write_rds(clinical, "clinical.rds")

################################################################################################# IV ### Bind df----
radiomics <- full_join(features, clinical, by = "patient_id") %>% 
  filter(!is.na(lesion_id))


################################################################################################# V ### New file----

merge_data <- 
  read_csv(paste0("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/merge_clinical_normalized_radiomics_M4MOC_03192021.csv")) %>% 
  select("mrn", "date", "lesionid", "segmentedby", first(starts_with("nor_")):ncol(.)) %>% 
  mutate(across(
    .cols = contains("date"),
    .fns = ~ as.POSIXct(., format = "%m/%d/%Y"))) 

final_data <- clinical_cleaning(data = merge_data) %>% 
  filter(!is.na(lesionid))
write_rds(final_data, "radiomics.rds")


