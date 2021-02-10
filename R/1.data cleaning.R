# Import packages
library(tidyverse)
library(scales)


################################################################################################# I ### Load data----
radiomics <- read.csv(paste0("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/Ovarian_Radiomics_Features_01062021.csv")) %>%
  `colnames<-`(str_remove(colnames(.), "F[0-9]*\\.") %>% tolower() %>% str_replace_all(., "\\.", "_"))
colnames(radiomics)
clinical <- readxl::read_xls("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/radiomic_CT_final.xls") %>% 
  mutate(MRN = as.character(MRN))


################################################################################################# II ### Radiomics----
summary(radiomics)
radiomics <- radiomics %>% mutate(mrn = as.character(mrn))
# scale data from -1 to 1
for(i in 1:length(colnames(radiomics))) {
  if(class(radiomics[,i]) == "numeric" | class(radiomics[,i]) == "integer") {
    radiomics[,i] <- rescale(radiomics[,i], to=c(-1,1)) }
}
summary(radiomics)







################################################################################################# III ### Bind df----

radiomics1 <- full_join(radiomics, clinical, by = c("mrn"= "MRN")) %>% 
  filter(!is.na(lesion_id))

