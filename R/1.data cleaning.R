# Import packages
library(tidyverse)



# Load data
radiomics <- read.csv(paste0("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/Ovarian_Radiomics_Features_01062021.csv"))
clinical <- readxl::read_xls("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/radiomic_CT_final.xls")

summary(radiomics)

