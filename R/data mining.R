################################################################################################# I ### Clinical mining----
# radiomics1 %>% 
#   ggplot(aes(x= Race_Cancer_Registry)) + 
#   geom_bar()+
#   coord_flip()+
#   theme_minimal()
# 
# radiomics1 %>% 
#   ggplot(aes(x= Ethnicity_Cancer_Registry)) + 
#   geom_bar()+
#   coord_flip()+
#   theme_minimal()
# 
# radiomics1 %>% 
#   ggplot(aes(x= Primary_Site)) + 
#   geom_bar()+
#   coord_flip()+
#   theme_minimal()
# 
# radiomics1 %>% 
#   ggplot(aes(x= Histology)) + 
#   geom_bar()+
#   coord_flip()+
#   theme_minimal()
# 
# radiomics1 %>% 
#   ggplot(aes(x= TNM_CS__Mixed_Group_Stage)) + 
#   geom_bar()+
#   coord_flip()+
#   theme_minimal()
radiomics1 %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(dataset, Race_Cancer_Registry, Ethnicity_Cancer_Registry, Primary_Site, Histology, TNM_CS__Mixed_Group_Stage) %>% 
  tbl_summary(by = dataset, 
              sort = list(everything() ~ "frequency", TNM_CS__Mixed_Group_Stage ~ "alphanumeric")) %>% 
  bold_labels()

radiomics1 %>% mutate(dataset = "Patients") %>% distinct(mrn, .keep_all = TRUE) %>% 
  select(treatment__type, Race_Cancer_Registry, Ethnicity_Cancer_Registry, Primary_Site, Histology, TNM_CS__Mixed_Group_Stage) %>% 
  tbl_summary(by = treatment__type, 
              sort = list(everything() ~ "frequency", TNM_CS__Mixed_Group_Stage ~ "alphanumeric")) %>% 
  bold_labels() %>% add_p()

radiomics1 %>% distinct(mrn, .keep_all = TRUE) %>% 
  ggplot(aes(x= Summary_of_Rx__1st_course)) + 
  geom_bar()+
  coord_flip()+
  theme_minimal()

radiomics1 %>% distinct(mrn, .keep_all = TRUE) %>% 
  ggplot(aes(x= treatment__type)) + 
  geom_bar()+
  coord_flip()+
  theme_minimal()

radiomics1 %>% distinct(mrn, .keep_all = TRUE) %>% 
  ggplot(aes(x= Debulking_status)) + 
  geom_bar()+
  coord_flip()+
  theme_minimal()

# table age at recurremce, diagnosis, drugs, surgery

# pie chart brca germline somatic + table

# table comorbidities



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


# Correlation ----
# Test normality
qqnorm(radiomics1$statistical_mean)
# Pearson
a <- radiomics1[, 5:310] %>% select(where(~ any(. != 0)))

mat <- cor(a, use = "pairwise.complete.obs")
corrplot.mixed(mat, tl.pos = "n")

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

a <- radiomics1[, 5:310] %>% select(contains("statistical")) %>% select(where(~ any(. != 0)))
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
a <- radiomics1[, 5:310] %>% select(contains("histogram")) %>% select(where(~ any(. != 0)))
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
a <- radiomics1[, 5:310] %>% select(contains("fraction")) %>% select(where(~ any(. != 0)))
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
a <- radiomics1[, 5:310] %>% select(contains("avgcoocurrence")) %>% select(where(~ any(. != 0)))
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
a <- radiomics1[, 5:310] %>% select(contains("avg_3d")) %>% select(where(~ any(. != 0)))
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
a <- radiomics1[, 5:310] %>% select(contains("glszm")) %>% select(where(~ any(. != 0)))
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
a <- radiomics1[, 5:310] %>% select(contains("ngtdm")) %>% select(where(~ any(. != 0)))
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
a <- radiomics1[, 5:310] %>% select(contains("features")) %>% select(where(~ any(. != 0)))
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
a <- radiomics1[, 5:310] %>% select(contains("wavelet")) %>% select(where(~ any(. != 0)))
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
a <- radiomics1[, 5:310] %>% select(-matches("statistical|histogram|fraction|avgcoocurrence|avg_3d|glszm|ngtdm|features|wavelet")) %>% select(where(~ any(. != 0)))
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



