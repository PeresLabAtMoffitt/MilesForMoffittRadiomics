# MilesForMoffittRadiomics

# Cleaning
The "01.data cleaning.R" script was used to clean and combined clinical (recoding/cleaning/create new variables) and 
radiomics feature data (scale from -1 to 1).  
The code generate radiomics.rds and clinical.rds that are used in the later analyses/modeling.  
Other R scripts are not useful anymore.

# ML
Use the Rmd.  
"Stable features.Rmd" calculate the ICC / CCC of each feature and select the best features with CCC ≥ 0.95.  
"survival modeling.Rmd" is the last modeling script using a Decision tree.  
"modeling-outcomes-balanced" build multiple models and choosing the best to predict recurrence Yes/No with a balanced data.  
