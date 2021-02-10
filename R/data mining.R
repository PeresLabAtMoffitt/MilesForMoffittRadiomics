################################################################################################# I ### Clinical mining----




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


# Heatmap
library(gplots)
library(heatmap.plus)
library(RColorBrewer)

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



