library(sf)
library(dplyr)

sf<-st_read(file.choose())

##Parse canopy relative height metrics
# Remove brackets and extra spaces
sf$val_clean <- gsub("\\[|\\]", "", sf$canopy_h_metrics)
sf$val_clean <- trimws(sf$val_clean)

# Split into list of values
sf$val_list <- strsplit(sf$val_clean, "\\s+")

# Convert list to data frame with separate columns
df_parsed <- as.data.frame(do.call(rbind, sf$val_list))

col=c("rh10","rh15","rh20","rh25","rh30","rh35","rh40","rh45","rh50","rh55", "rh60","rh65", 
      "rh70", "rh75", "rh80", "rh85", "rh90", "rh95","rh98","rh99")

# Rename columns (optional)
colnames(df_parsed) <- col


##Merge and rename columns
sf.2=select(sf, -val_clean,-val_list)
# Extract centroid coordinates
centroids <- st_coordinates(st_centroid(sf))
sf.3 <- cbind(sf.2, centroids)
st_geometry(sf.3) <- NULL
sf.4=cbind(sf.3,df_parsed)

#rename column
sf.4 <- sf.4 %>%
  rename(class=SAMPLE_1 ,
         vv=SAMPLE_2,
         vh=SAMPLE_3,
         b2=SAMPLE_4,
         b3=SAMPLE_5,
         b4=SAMPLE_6,
         b8=SAMPLE_7,
         b11=SAMPLE_8,
         b12=SAMPLE_9)

write.table(sf.4,"H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\dataset_v10_fin.csv",sep=";")