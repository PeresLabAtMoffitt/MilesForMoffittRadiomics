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
  mutate(age_at_Dx = round(interval(start = date_of_birth, end = date_of_diagnosis)/
                             duration(n=1, units = "years"), 2)) %>% 
  mutate(age_at_first_adjuvant_chem = round(interval(start = date_of_birth, end = date_of_first_adjuvant_chemother)/
                                              duration(n=1, units = "years"), 2)) %>% 
  mutate(age_at_first_neoadjuvant_chem = round(interval(start = date_of_birth, end = date_of_first_neoadjuvant_chemot)/
                                                 duration(n=1, units = "years"), 2)) %>% 
  mutate(age_at_first_surgery = round(interval(start = date_of_birth, end = date_of_first_surgery)/
                                        duration(n=1, units = "years"), 2)) %>% 
  mutate(age_at_first_chemo = coalesce(age_at_first_neoadjuvant_chem, age_at_first_adjuvant_chem)) %>% 
  mutate(age_at_surgery = round(interval(start = date_of_birth, end = date_of_surgery)/
                                  duration(n=1, units = "years"), 2)) %>% 
  mutate(age_at_first_recurrence = round(interval(start = date_of_birth, end = date_of_first_recurrence)/
                                           duration(n=1, units = "years"), 2)) %>% 
  mutate(month_at_first_recurrence = round(interval(start = date_of_diagnosis, end = date_of_first_recurrence)/
                                           duration(n=1, units = "months"), 2))






################################################################################################# IV ### Bind df----

radiomics1 <- full_join(radiomics, clinical, by = "mrn") %>% 
  filter(!is.na(lesion_id))

