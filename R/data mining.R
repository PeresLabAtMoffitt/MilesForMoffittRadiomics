################################################################################################# I ### Clinical mining----
# Table clinical
radiomics %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(dataset, race, ethnicity, primary_site, histology, tnm_stage) %>% 
  tbl_summary(by = dataset, 
              sort = list(everything() ~ "frequency", tnm_stage ~ "alphanumeric")) %>% 
  bold_labels()

radiomics %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(treatment_type, race, ethnicity, primary_site, histology, tnm_stage) %>% 
  tbl_summary(by = treatment_type, 
              sort = list(everything() ~ "frequency", tnm_stage ~ "alphanumeric")) %>% 
  bold_labels() %>% add_p()

# Table age at recurremce, diagnosis, drugs, surgery
radiomics %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(dataset, "age_at_Dx", "months_at_first_neoadjuvant_chem", "months_at_first_surgery", "months_at_first_adjuvant_chem", 
         "months_at_first_chemo", "age_at_surgery", "age_at_first_recurrence", "month_at_first_recurrence") %>% 
  tbl_summary(by = dataset,
              digits = list(everything()~ 2),
              missing = "no") %>% 
  bold_labels()

# Plot treatment
radiomics %>% distinct(mrn, .keep_all = TRUE) %>% 
  ggplot(aes(x= summary_of_rx_1st_course)) + 
  geom_bar()+
  coord_flip()+
  theme_minimal()

radiomics %>% distinct(mrn, .keep_all = TRUE) %>% 
  ggplot(aes(x= treatment_type)) + 
  geom_bar()+
  coord_flip()+
  theme_minimal()

radiomics %>% distinct(mrn, .keep_all = TRUE) %>% 
  ggplot(aes(x= debulking_status)) + 
  geom_bar()+
  coord_flip()+
  theme_minimal()

# Table treatment by recurrence
radiomics %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(treatment_type, debulking_status, has_the_patient_recurred_, tnm_stage, summary_of_rx_1st_course) %>% 
  tbl_summary(by = has_the_patient_recurred_, 
              sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% add_p()

tbl1 <- radiomics %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(treatment_type, month_at_first_recurrence) %>% 
  tbl_summary(by = treatment_type, 
              digits = list(everything() ~ 2)) %>% 
  bold_labels() %>% add_p()
tbl2 <- radiomics %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(debulking_status, month_at_first_recurrence) %>% 
  tbl_summary(by = debulking_status, 
              digits = list(everything() ~ 2)) %>% 
  bold_labels() %>% add_p()
tbl3 <- radiomics %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(tnm_stage, month_at_first_recurrence) %>% 
  tbl_summary(by = tnm_stage, 
              digits = list(everything() ~ 2)) %>% 
  bold_labels() %>% add_p()
tbl4 <- radiomics %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(summary_of_rx_1st_course, month_at_first_recurrence) %>% 
  tbl_summary(by = summary_of_rx_1st_course, 
              digits = list(everything() ~ 2)) %>% 
  bold_labels() %>% add_p()

tbl_merge(list(tbl1, tbl2, tbl3, tbl4),
          tab_spanner = c("**treatment_type**", "**debulking_status**", "**tnm_stage**", "**summary_of_rx_1st_course**")) %>%
  bold_labels() %>%
  italicize_levels()

# Table brca by recurrence
radiomics %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(has_the_patient_recurred_, germline_brca1_mutation, germline_brca2_mutation, 
         somatic_brca1_mutation, somatic_brca2_mutation, any_unclassified_brca_mutation) %>%
  tbl_summary(by = has_the_patient_recurred_, 
              sort = list(everything() ~ "frequency"),
              missing = "no") %>% 
  bold_labels() %>% add_p()

# Plot brca germline somatic
radiomics %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, germline_brca1_mutation, somatic_brca1_mutation, germline_brca2_mutation, somatic_brca2_mutation) %>% 
  pivot_longer(cols = c(germline_brca1_mutation, somatic_brca1_mutation, germline_brca2_mutation, somatic_brca2_mutation), 
               names_to = "name",values_to = "value") %>%
  filter(value == "Yes") %>% 
  mutate(brca = case_when(
    str_detect(name, "brca1") ~ "brca1",
    str_detect(name, "brca2") ~ "brca2"
  )) %>% 
  mutate(name = str_remove_all(name, "_brca1|_brca2|_mutation")) %>%
  group_by(name, brca) %>% 
  summarise(count=n()) %>%
  ggplot(aes(x = brca, y= count, fill= name)) + 
  geom_bar(stat = "identity")+
  geom_text(aes(label = paste0("n=", count)), size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal()

# table comorbidities
radiomics %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  mutate(pre_dx_hypertension = ifelse(str_detect(hypertension, "pre"), "hypertension", "no hypertension")) %>% 
  mutate(pre_dx_diabetes = ifelse(str_detect(diabetes_mellitus, "pre"), "diabetes", "no diabetes")) %>% 
  mutate(pre_dx_hypercholesterolemia = ifelse(str_detect(hypercholesterolemia, "pre"), 
                                              "hypercholesterolemia", "no hypercholesterolemia")) %>% 
  mutate(pre_dx_CKd = ifelse(str_detect(chronic_kidney_disease, "pre"), "CKd", "no CKd")) %>% 
  mutate(pre_dx_cardiac_conditions = ifelse(str_detect(cardiac_conditions_including_bu, "pre"), 
                                            "cardiac conditions", "no cardiac conditions")) %>% 
  select("has_the_patient_recurred_", "pre_dx_hypertension", "pre_dx_diabetes",
         "pre_dx_hypercholesterolemia", "pre_dx_CKd", "pre_dx_cardiac_conditions") %>%
  tbl_summary(by = has_the_patient_recurred_, 
              sort = list(everything() ~ "frequency"),
              missing = "no") %>%
  # modify_spanning_header(starts_with("stat_") ~ "**Randomization Assignment**")
  bold_labels() %>% add_p()
  show_header_names(a)

################################################################################################# II ### Radiomics mining----
# ICC
library(irr) # Does not count patient with NA
ICC_rad <- radiomics %>% select(mrn, statistical_mean)
ICC_rad <- dcast(setDT(ICC_rad), mrn ~ rowid(mrn), 
                 value.var = "statistical_mean")  %>% 
  select(-mrn)
  
icc(
  ICC_rad, model = "twoway", # # between 0.75 and 0.90: good
  type = "consistency", unit = "average"
)

library(psych) # Koo and Li (2016)
ICC(ICC_rad)  # between 0.75 and 0.90: good
ICC(ICC_rad)$results[4,2]

# Correlation ----
# Test normality
qqnorm(radiomics$statistical_mean)
# Pearson
a <- radiomics[, 5:310] %>% select(where(~ any(. != 0)))

mat <- cor(a, use = "pairwise.complete.obs")
# corrplot.mixed(mat, tl.pos = "n")

ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           # outline.col = "darkblue", # the outline of the circle or square
           # hc.method = "complete",
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           # colors = viridis::inferno(n=3),
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           # p.mat = pmat, # Add correlation significance
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)

a <- radiomics[, 5:310] %>% select(contains("statistical")) %>% select(where(~ any(. != 0)))
mat <- cor(a, use = "pairwise.complete.obs")
ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)
a <- radiomics[, 5:310] %>% select(contains("histogram")) %>% select(where(~ any(. != 0)))
mat <- cor(a, use = "pairwise.complete.obs")
ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)
a <- radiomics[, 5:310] %>% select(contains("fraction")) %>% select(where(~ any(. != 0)))
mat <- cor(a, use = "pairwise.complete.obs")
ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)
a <- radiomics[, 5:310] %>% select(contains("avgcoocurrence")) %>% select(where(~ any(. != 0)))
mat <- cor(a, use = "pairwise.complete.obs")
ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)
a <- radiomics[, 5:310] %>% select(contains("avg_3d")) %>% select(where(~ any(. != 0)))
mat <- cor(a, use = "pairwise.complete.obs")
ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)
a <- radiomics[, 5:310] %>% select(contains("glszm")) %>% select(where(~ any(. != 0)))
mat <- cor(a, use = "pairwise.complete.obs")
ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)
a <- radiomics[, 5:310] %>% select(contains("ngtdm")) %>% select(where(~ any(. != 0)))
mat <- cor(a, use = "pairwise.complete.obs")
ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)
a <- radiomics[, 5:310] %>% select(contains("features")) %>% select(where(~ any(. != 0)))
mat <- cor(a, use = "pairwise.complete.obs")
ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)
a <- radiomics[, 5:310] %>% select(contains("wavelet")) %>% select(where(~ any(. != 0)))
mat <- cor(a, use = "pairwise.complete.obs")
ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)
a <- radiomics[, 5:310] %>% select(-matches("statistical|histogram|fraction|avgcoocurrence|avg_3d|glszm|ngtdm|features|wavelet")) %>% select(where(~ any(. != 0)))
mat <- cor(a, use = "pairwise.complete.obs")
ggcorrplot(mat, hc.order = TRUE, method = "circle", 
           type = "lower", # show the top half panel
           lab = TRUE, # add correlation nbr
           title = "Correlation between radiomics features",
           show.legend = TRUE, legend.title = "Correlation", show.diag = TRUE,
           lab_col = "darkblue", lab_size = 3, # col and size of the correlation nbr
           sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 10, 
           tl.cex = 10, tl.col = "red", tl.srt = 40,
           digits = 2
)

# Heatmap
df1 <- radiomics %>% select(-date, -segmented_by) %>% 
  unite(mrn, c(mrn, lesion_id), sep = "_", remove = TRUE) %>% 
  `row.names<-`(.$mrn) %>% 
  select(-mrn) 
df1 <- as.matrix(df1)

heatmap.2(df1, main = "Immune Marker Presentation",
          
          trace = "none", density="none", col=bluered(20), cexRow=1, cexCol = 1, 
          margins = c(10,5), # bottom, right
          ColSideColors = ,
          scale = "row")

df2 <- t(df1)

heatmap.2(df2, main = "Immune Marker Presentation",
          
          trace = "none", density="none", col=bluered(20), cexRow=1, cexCol = 1, 
          margins = c(10,5), # bottom, right
          ColSideColors = ,
          scale = "column")

condition_colors <- unlist(lapply(rownames(df2), function(x){
  if(grepl("statistical", x)) "#000004FF"
  if(grepl("histogram", x)) "#1B0C42FF"
  if(grepl("fraction", x)) "#4B0C6BFF"
  if(grepl("avgcoocurrence", x)) "#781C6DFF"
  if(grepl("avg_3d", x)) "#A52C60FF"
  if(grepl("glszm", x)) "#CF4446FF"
  if(grepl("ngtdm", x)) "#ED6925FF"
  if(grepl("features", x)) "#FB9A06FF"
  if(grepl("wavelet", x)) "#F7D03CFF"
  else if(grepl(".", x)) "grey"
}))

condition_colors <- unlist(lapply(rownames(df2), function(x){
  ifelse(grepl("statistical", x), "#000004FF", ifelse(grepl("histogram", x), "#1B0C42FF", "grey"))
  # if(grepl("histogram", x)) "#1B0C42FF"
  # if(grepl("fraction", x)) "#4B0C6BFF"
  # if(grepl("avgcoocurrence", x)) "#781C6DFF"
  # if(grepl("avg_3d", x)) "#A52C60FF"
  # if(grepl("glszm", x)) "#CF4446FF"
  # if(grepl("ngtdm", x)) "#ED6925FF"
  # if(grepl("features", x)) "#FB9A06FF"
  # if(grepl("wavelet", x)) "#F7D03CFF"
  # else if(grepl(".", x)) "grey"
}))

heatmap.2(df2, main = "Immune Marker Presentation",
          
          trace = "none", density="none", col=bluered(20), cexRow=1, cexCol = 1, 
          margins = c(10,5), # bottom, right
          RowSideColors = condition_colors,
          scale = "column")



