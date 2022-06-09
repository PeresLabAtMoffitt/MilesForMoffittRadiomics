# Import packages
library(tidyverse)


################################################################################ I ### Load data----
path <- fs::path("", "Volumes", "Peres_Research", "Ovarian - Radiomics")
contrast <- readxl::read_xlsx(
  paste0(path,
         "/data/analysis dataset/WorkingFile_contrast_05052022.xlsx"))


clinical <- read_rds(paste0(here::here(), "/clinical.rds"))


################################################################################ I ### Merge data----
radiogenomics_pilot_study <- contrast %>% 
  filter(treatment__type == "upfront surgery") %>% 
  select(mrn, treatment__type) %>% 
  left_join(., 
            clinical %>% 
              select(mrn, date_of_diagnosis, date_of_first_neoadjuvant_chemot, date_of_first_adjuvant_chemother,
                     first_chemo_date))

write_csv(radiogenomics_pilot_study, "ids for radiogenomics pilot study.csv")
