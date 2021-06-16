################################################################################# Select stable features

# ICC
first_seg <- read_csv("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/Ovarian_Radiomics_Features_TestRetest_first segmentation.csv") %>%
  `colnames<-`(str_remove(colnames(.), " _first_segmentation") %>% tolower() %>% gsub(":", "_",.)) %>% 
  mutate(date = as.Date(date, format = "%m/%d/%y"))

second_seg <- read_csv("/Users/colinccm/Documents/GitHub/Peres/data/Radiomics/Ovarian_Radiomics_Features_TestRetest_second segmentation.csv") %>%
  `colnames<-`(str_remove(colnames(.), " _second_seg") %>% tolower() %>% gsub(":", "_",.)) %>% 
  mutate(date = as.Date(as.character(date), format = "%Y%m%d"))


segmentation <- bind_rows(first_seg, second_seg) %>% 
  mutate(seg = "segmentation0") %>% 
  unite(segmentationcohort, c(seg, segmentationcohort), sep = "") %>% 
  mutate(mrn = as.character(mrn)) %>% 
  # select(-c(date, 'lesion id', 'segmented by')) %>% 
  as.data.frame(.)

ICC_data <- data.frame(matrix(nrow=1, ncol=0)) #data.frame(x = c("variable", "ICC"))
for(i in 1:length(colnames(segmentation))) {
  
  rad <- segmentation %>% select(mrn, segmentationcohort)
  
  if(class(segmentation[,i]) == "numeric" | class(segmentation[,i]) == "integer") {
    
    ICC_df <- cbind(rad, value = segmentation[,i])
    ICC_df <- ICC_df %>%
      pivot_wider(names_from = segmentationcohort, values_from = value) %>%
      select(c(starts_with("segmentation0")))
    
    # ICC <- ICC(ICC_df)$results[3,2] # 2 way mixed effect for single rater
    ICC <- icc(
      ICC_df, model = "twoway", 
      type = "agreement", unit = "single"
    )$value
    ICC_data <- cbind(ICC_data,ICC)
    
  }
  
}
colnames(ICC_data) <- colnames(segmentation)[6:ncol(segmentation)]
rm(ICC_df, ICC, i, rad)

library(plotly)

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

ICC_data <- ICC_data %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(reliability = case_when(
    value < 0.5 ~ "poor",
    value >= 0.5 &
      value < 0.75 ~ "moderate",
    value >= 0.75 &
      value < 0.9 ~ "good",
    value >= 0.9 ~ "excellent"
  )) %>% arrange(reliability)
ICC_data



# CCC
library(DescTools)

CCC_data <- data.frame(matrix(nrow=1, ncol=0)) #data.frame(x = c("variable", "CCC"))
for(i in 1:length(colnames(segmentation))) {
  
  rad <- segmentation %>% select(mrn, segmentationcohort)
  
  if(class(segmentation[,i]) == "numeric" | class(segmentation[,i]) == "integer") {
    
    CCC_df <- cbind(rad, value = segmentation[,i])
    CCC_df <- CCC_df %>%
      pivot_wider(names_from = segmentationcohort, values_from = value) %>%
      select(c(starts_with("segmentation0")))
    
    CCC <- CCC(x=CCC_df$segmentation01, y=CCC_df$segmentation02, ci = "z-transform", conf.level = 0.95, na.rm = FALSE)$rho.c$est
    CCC_data <- cbind(CCC_data,CCC)
    
  }
  
}
colnames(CCC_data) <- colnames(segmentation)[6:ncol(segmentation)]
rm(CCC_df, CCC, i, rad)

library(plotly)

plot <- CCC_data %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(CCC = "CCC") %>% 
  mutate(reliability = case_when(
    value < 0.5 ~ "poor",
    value >= 0.5 &
      value < 0.75 ~ "moderate",
    value >= 0.75 &
      value < 0.9 ~ "good",
    value >= 0.9 ~ "excellent"
  )) %>% 
  mutate(text = paste("Feature: ", name, "\nCCC: ", round(value, 2), "\nStrict Reliability: ", reliability, sep="")) %>% 
  ggplot(aes(x = CCC, y = value, color = reliability, text=text))+
  labs(title = "Reliability of the radiomic features based on strict CCC value", x= NULL, y= "CCC value")+
  scale_color_brewer(palette = "Set2") +
  scale_x_discrete(labels = "")+
  theme_minimal()+
  geom_point(size = 1, position=position_jitter(w=0.8))
plot <- ggplotly(plot, tooltip="text")
plot

CCC_data <- CCC_data %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(reliability = case_when(
    value >= 0.95      ~ "stable",
    TRUE               ~ "not considred as stable"
  )) %>% arrange(desc(reliability))
CCC_data

# tmp.lab <- data.frame(lab = paste("CCC: ", 
#                                   round(CCC_data$rho.c[,1], digits = 2), " (95% CI ", 
#                                   round(CCC_data$rho.c[,2], digits = 2), " - ",
#                                   round(CCC_data$rho.c[,3], digits = 2), ")", sep = ""))
# 
# z <- lm(method2 ~ method1)
# alpha <- summary(z)$coefficients[1,1]
# beta <-  summary(z)$coefficients[2,1]
# tmp.lm <- data.frame(alpha, beta)
path <- fs::path("", "Volumes", "Peres_Research", "Ovarian - Radiomics", "Christelle")
write_csv(CCC_data %>% filter(value > 0.90), paste0(path, "/stable features/selected stable features at 90p.csv"))

concordance <- full_join(CCC_data, ICC_data, by = "name", suffix = c("_CCC", "_ICC")) %>% arrange(desc(value_CCC))
write_csv(concordance, paste0(path, "/stable features/concordance CCC and ICC.csv"))




