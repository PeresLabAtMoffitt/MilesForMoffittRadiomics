################################################################################################# I ### Clinical mining----
# Table clinical
radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(dataset, Race, Ethnicity, primary_site, Histology, tnm_stage) %>% 
  tbl_summary(by = dataset, 
              sort = list(everything() ~ "frequency", tnm_stage ~ "alphanumeric")) %>% 
  bold_labels()

radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(treatment_type, Race, Ethnicity, primary_site, Histology, tnm_stage) %>% 
  tbl_summary(by = treatment_type, 
              sort = list(everything() ~ "frequency", tnm_stage ~ "alphanumeric")) %>% 
  bold_labels() %>% add_p()

# Table age at recurremce, diagnosis, drugs, surgery
radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(dataset, "age_at_Dx", "months_at_first_neoadjuvant_chem", "months_at_first_surgery", "months_at_first_adjuvant_chem", 
         "months_at_first_chemo", "age_at_surgery", "age_at_first_recurrence", "month_at_first_recurrence") %>% 
  tbl_summary(by = dataset,
              digits = list(everything()~ 2),
              missing = "no") %>% 
  bold_labels() %>% add_n()

radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(treatment_type, "age_at_Dx", "months_at_first_neoadjuvant_chem", "months_at_first_surgery", "months_at_first_adjuvant_chem", 
         "months_at_first_chemo", "age_at_surgery", "age_at_first_recurrence", "month_at_first_recurrence") %>% 
  tbl_summary(by = treatment_type,
              digits = list(everything()~ 2),
              missing = "no") %>% 
  bold_labels() %>% add_n()

# Table treatment
radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(treatment_type, summary_of_rx_1st_course, debulking_status) %>% 
  tbl_summary(by = treatment_type, 
              sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% add_p()

# Plot treatment
radiomics %>% distinct(patient_id, .keep_all = TRUE) %>% 
  ggplot(aes(x= summary_of_rx_1st_course)) + 
  geom_bar()+
  coord_flip()+
  theme_minimal()

radiomics %>% distinct(patient_id, .keep_all = TRUE) %>% 
  ggplot(aes(x= treatment_type)) + 
  geom_bar()+
  coord_flip()+
  theme_minimal()

radiomics %>% distinct(patient_id, .keep_all = TRUE) %>% 
  ggplot(aes(x= debulking_status)) + 
  geom_bar()+
  coord_flip()+
  theme_minimal()

# Table treatment by recurrence
radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(treatment_type, debulking_status, has_the_patient_recurred_, tnm_stage, summary_of_rx_1st_course) %>% 
  tbl_summary(by = has_the_patient_recurred_, 
              sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% add_p()

radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(treatment_type, has_the_patient_recurred_) %>% 
  tbl_summary(by = treatment_type, 
              sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% add_p()

radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(treatment_type, has_the_patient_recurred_after_surg) %>% 
  tbl_summary(by = treatment_type, 
              sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% add_p()


tbl1 <- radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(treatment_type, month_at_first_recurrence) %>% 
  tbl_summary(by = treatment_type, 
              digits = list(everything() ~ 2)) %>% 
  bold_labels() %>% add_p()
tbl2 <- radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(debulking_status, month_at_first_recurrence) %>% 
  tbl_summary(by = debulking_status, 
              digits = list(everything() ~ 2)) %>% 
  bold_labels() %>% add_p()
tbl3 <- radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(tnm_stage, month_at_first_recurrence) %>% 
  tbl_summary(by = tnm_stage, 
              digits = list(everything() ~ 2)) %>% 
  bold_labels() %>% add_p()
tbl4 <- radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(summary_of_rx_1st_course, month_at_first_recurrence) %>% 
  tbl_summary(by = summary_of_rx_1st_course, 
              digits = list(everything() ~ 2)) %>% 
  bold_labels() %>% add_p()

tbl_merge(list(tbl1, tbl2, tbl3, tbl4),
          tab_spanner = c("**treatment_type**", "**debulking_status**", "**tnm_stage**", "**summary_of_rx_1st_course**")) %>%
  bold_labels() %>%
  italicize_levels()

# Table brca by recurrence
radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(has_the_patient_recurred_, germline_brca1_mutation, germline_brca2_mutation, 
         somatic_brca1_mutation, somatic_brca2_mutation, any_unclassified_brca_mutation) %>%
  tbl_summary(by = has_the_patient_recurred_, 
              sort = list(everything() ~ "frequency"),
              missing = "no") %>% 
  bold_labels() %>% add_p()

# Table brca by treatment
radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(treatment_type, germline_brca1_mutation, germline_brca2_mutation, 
         somatic_brca1_mutation, somatic_brca2_mutation, any_unclassified_brca_mutation) %>%
  tbl_summary(by = treatment_type, 
              sort = list(everything() ~ "frequency"),
              missing = "no") %>% 
  bold_labels() %>% add_p()

# Plot brca germline somatic
radiomics %>% distinct(patient_id, .keep_all = TRUE) %>% 
  select(patient_id, germline_brca1_mutation, somatic_brca1_mutation, germline_brca2_mutation, somatic_brca2_mutation) %>% 
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
radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
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
  # show_header_names(a)

radiomics %>% mutate(dataset = "Patients") %>% distinct(patient_id, .keep_all = TRUE) %>% 
  mutate(pre_dx_hypertension = ifelse(str_detect(hypertension, "pre"), "hypertension", "no hypertension")) %>% 
  mutate(pre_dx_diabetes = ifelse(str_detect(diabetes_mellitus, "pre"), "diabetes", "no diabetes")) %>% 
  mutate(pre_dx_hypercholesterolemia = ifelse(str_detect(hypercholesterolemia, "pre"), 
                                              "hypercholesterolemia", "no hypercholesterolemia")) %>% 
  mutate(pre_dx_CKd = ifelse(str_detect(chronic_kidney_disease, "pre"), "CKd", "no CKd")) %>% 
  mutate(pre_dx_cardiac_conditions = ifelse(str_detect(cardiac_conditions_including_bu, "pre"), 
                                            "cardiac conditions", "no cardiac conditions")) %>% 
  select(treatment_type, "pre_dx_hypertension", "pre_dx_diabetes",
         "pre_dx_hypercholesterolemia", "pre_dx_CKd", "pre_dx_cardiac_conditions") %>%
  tbl_summary(by = treatment_type, 
              sort = list(everything() ~ "frequency"),
              missing = "no") %>%
  # modify_spanning_header(starts_with("stat_") ~ "**Randomization Assignment**")
  bold_labels() %>% add_p() %>%  add_n()
  
################################################################################################# II ### Radiomics mining----
# ICC
library(irr) # Does not count patient with NA
ICC_rad <- radiomics %>% select(patient_id, statistical_mean)
ICC_rad <- dcast(setDT(ICC_rad), patient_id ~ rowid(patient_id), 
                 value.var = "statistical_mean")  %>% 
  select(-patient_id)
  
icc(
  ICC_rad, model = "twoway", # # between 0.75 and 0.90: good
  type = "consistency", unit = "average"
)


ICC(ICC_rad)  # between 0.75 and 0.90: good
ICC(ICC_rad)$results[4,2]

ICC_rad <- radiomics %>% select(patient_id, `3d_wavelet_p2_l2_c15`)
ICC_rad <- dcast(setDT(ICC_rad), patient_id ~ rowid(patient_id), 
                 value.var = "3d_wavelet_p2_l2_c15")  %>% 
  select(-patient_id)

ICC(ICC_rad)$results[4,2]
####################
ICC_data <- data.frame(matrix(nrow=1, ncol=0)) #data.frame(x = c("variable", "ICC"))
for(i in 1:length(colnames(features))) {
  
  rad <- features %>% select(patient_id, lesion_id)
  
  if(class(features[,i]) == "numeric" | class(features[,i]) == "integer") {
    
    ICC_df <- cbind(rad, value = features[,i])
    ICC_df <- ICC_df %>%
      mutate(Lesion = "Lesion0") %>%
      group_by(patient_id) %>%
      mutate(n = row_number(patient_id)) %>%
      ungroup() %>%
      unite(lesion_id, Lesion:n, sep = "", remove = TRUE, na.rm = TRUE) %>%
      pivot_wider(names_from = lesion_id, values_from = value) %>%
      select(c(starts_with("Lesion")))
    
    ICC <- ICC(ICC_df)$results[4,2]
    ICC_data <- cbind(ICC_data,ICC)
    
  }
  
}
colnames(ICC_data) <- colnames(features)[5:ncol(features)]
rm(ICC_df, ICC, rad)

ICC_data %>% 
  pivot_longer(cols = everything()) %>% mutate(ICC = "ICC") %>% 
  ggplot(aes(x = ICC, y = value))+
  geom_jitter()

ICC_data %>% 
  pivot_longer(cols = everything()) %>% mutate(ICC = "ICC") %>% 
  ggplot(aes(x = reorder(name, value)))+
  geom_segment(aes(xend = name, y= 0, yend = value))

library(plotly)
# library(viridis)
# library(hrbrthemes)

plot <- ICC_data %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(ICC = "ICC") %>% 
  mutate(reliability = case_when(
    value < 0.5 ~ "poor",
    value >= 0.5 &
      value < 0.75 ~ "moderate",
    value >= 0.75 &
      value < 0.9 ~ "good",
    value >= 0.9 ~ "excellent"
  )) %>% 
  mutate(text = paste("Feature: ", name, "\nICC: ", round(value, 2), "\nStrict Reliability: ", reliability, sep="")) %>% 
  ggplot(aes(x = ICC, y = value, color = reliability, text=text))+
  labs(title = "Reliability of the radiomic features based on strict ICC value", x= NULL, y= "ICC value")+
  scale_color_brewer(palette = "Set2") +
  scale_x_discrete(labels = "")+
  theme_minimal()+
  geom_point(size = 1, position=position_jitter(w=0.8))
plot <- ggplotly(plot, tooltip="text")
plot








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

# Heatmap----
df1 <- features %>% select(-date, -segmented_by) %>% 
  unite(patient_id, c(patient_id, lesion_id), sep = "_", remove = TRUE) %>% 
  `row.names<-`(.$patient_id) %>% 
  select(-patient_id) 
df1 <- as.matrix(df1)

heatmap.2(df1, main = "",
          
          trace = "none", density="none", col=bluered(20), cexRow=1, cexCol = 1, 
          margins = c(10,5), # bottom, right
          ColSideColors = ,
          scale = "row")

df2 <- t(df1)

heatmap.2(df2, main = "",
          
          trace = "none", density="none", col=bluered(20), cexRow=1, cexCol = 1, 
          margins = c(10,5), # bottom, right
          ColSideColors = ,
          scale = "column")

condition_colors <- unlist(lapply(rownames(df2), function(x){
  if(grepl("statistical", x)) "black"
  if(grepl("histogram", x)) "snow"
  if(grepl("fraction", x)) "snow4"
  if(grepl("avgcoocurrence", x)) "steelblue1"
  if(grepl("avg_3d", x)) "tan1"
  if(grepl("glszm", x)) "thistle1"
  if(grepl("ngtdm", x)) "turquoise1"
  if(grepl("features", x)) "violetred"
  if(grepl("wavelet", x)) "springgreen"
  else if(grepl(".", x)) "grey"
}))

condition_colors <- unlist(lapply(rownames(df2), function(x){
  ifelse(grepl("statistical", x), "black", 
         ifelse(grepl("histogram", x), "snow", 
                ifelse(grepl("fraction", x), "snow4",
                       ifelse(grepl("avgcoocurrence", x), "steelblue1",
                              ifelse(grepl("avg_3d", x), "tan1",
                                     ifelse(grepl("glszm", x), "thistle1",
                                            ifelse(grepl("ngtdm", x), "turquoise1",
                                                   ifelse(grepl("features", x), "violetred",
                                                          ifelse(grepl("wavelet", x), "springgreen", "grey")
         ))))))))
}))

heatmap.2(df2, main = "",
          
          trace = "none", density="none", col=bluered(20), cexRow=1, cexCol = 1, 
          margins = c(10,5), # bottom, right
          RowSideColors = condition_colors,
          scale = "column")



