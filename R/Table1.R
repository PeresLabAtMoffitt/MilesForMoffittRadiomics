library(tidyverse)
library(gtsummary)
theme_gtsummary_compact()

tbl <- analysis_data %>% 
  select(age_at_diagnosis, year_of_diagnosis,
         race, ethnicity, raceeth, gender,
         tnm_cs_mixed_group_stage, 
         primary_site, grade_differentiation, 
         debulking_status,
         vital_new, has_the_patient_recurred,
         starts_with("preDx"),
         treatment_type,
         contrastenhancementyn
         ) %>% 
  tbl_summary(by= treatment_type,
              type = list(starts_with("PreDx") ~ "categorical")) %>% 
  bold_labels() %>% add_p() %>% bold_p(t= .05) %>% 
  add_overall() %>% 
  as_gt() 
tbl %>% 
  gt::gtsave(zoom=1.5,
    filename = "Radiomics Table1.png"
  )








