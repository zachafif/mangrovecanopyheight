library(sf)
library(dplyr)
library(terra)
library(randomForest)
library(ModelMetrics)
library(ggplot2)
library(tidyr)
library(gridExtra)

###Data Preparation
##Parse canopy relative height metrics
# Remove brackets and extra spaces
inp<-st_read("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\dataset_v10.geojson")
inp$val_clean <- gsub("\\[|\\]", "", inp$canopy_h_metrics)
inp$val_clean <- trimws(inp$val_clean)

# Split into list of values
inp$val_list <- strsplit(inp$val_clean, "\\s+")

# Convert list to data frame with separate columns
df_parsed <- as.data.frame(do.call(rbind, inp$val_list))

col=c("rh10","rh15","rh20","rh25","rh30","rh35","rh40","rh45","rh50","rh55", "rh60","rh65", 
      "rh70", "rh75", "rh80", "rh85", "rh90", "rh95","rh98","rh99")

# Rename columns (optional)
colnames(df_parsed) <- col


#Merge and rename columns
inp.2=select(inp, -val_clean,-val_list)
# Extract centroid coordinates
centroids <- st_coordinates(st_centroid(inp))
inp.3 <- cbind(inp.2, centroids)
st_geometry(inp.3) <- NULL
inp.4=cbind(inp.3,df_parsed)

#rename column
inp.4 <- inp.4 %>%
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

###Prediction Modelling
##Load dataset and filter data
d<-read.csv2("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\dataset_v10_fin.csv")
d= type.convert(d, as.is = TRUE) #Convert string columns (number) to numeric
sf=st_read("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\Trat_BioHeight_30112025-20251209T123704Z-1-001\\Trat_BioHeight_30112025\\DataRevised_30112025\\GCP_selected.geojson")
sf$ID=1:nrow(sf)
vect=vect(sf)

d.filt <- subset(d, d$b2 > -1000 # Filter dataset based on sentinel-2 that has nan value
                 & d$b3 > -1000
                 & d$b4 > -1000
                 & d$b8 > -1000
                 & d$b11 > -1000
                 & d$b12 > -1000, )
d.filt2<-d.filt %>% filter(vh < 0 | vv < 0) # Filter dataset based on sentinel-1 that has nan value
d.filt3<-d.filt2%>% filter(h_canopy>0,h_canopy<30) # Filter dataset based on icesat-2 canopy height
d.filt3$ort_h=d.filt3$h_te_median+22.021 # Calculate orthometric height
d.filt4<-d.filt3%>% filter(ort_h<20,class !=6) # Filter dataset based on icesat-2 orthometric height and class
d.final<-d.filt4 #all features
d.final$X<-1:dim(d.final)[1] 

#Full coverage dataset
r<-terra::rast("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\stacked_w2w_v4.tif")
w2w <- as.data.frame(r, xy = TRUE)
colnames(w2w)[3]<-'class'
colnames(w2w)[4]<-'vv'
colnames(w2w)[5]<-'vh'
colnames(w2w)[6]<-'b2'
colnames(w2w)[7]<-'b3'
colnames(w2w)[8]<-'b4'
colnames(w2w)[9]<-'b8'
colnames(w2w)[10]<-'b11'
colnames(w2w)[11]<-'b12'

w2w.filt<-w2w%>%filter(!is.na(class),class!=0)


##Train-Test Split using Stratified Sampling based on species class##
a <- table(d.final$class) #Count frequency of each class
prop <- 0.7 #Set the proportion for training data
freqclass <- as.numeric(a) 
trfreq <- round(freqclass * prop, digits = 0) #Calculate how many samples per class
d.fin <- d.final[order(d.final$class),] 
sizeArg <- trfreq

# Create stratified sample data frame
rows <- dim(d.final)[1]+1
iter<-50
smpdf = data.frame(matrix(NA, nrow=round((rows * prop)),ncol=0))

#Run stratified sampling in a loop
for (i in 1:50){
 smp <- sampling:::strata(d.fin,'class', sizeArg, method=c('srswor'))
 smpind <- smp$ID_unit
 smpdf <- cbind(smpdf, smpind)
}
colnames(smpdf)  <- c(paste('smp_',c(1:iter),sep=''))

##Model Training and Testing##
predictor_cols <- c( "vv","vh","b2","b3","b3","b4","b8","b11","b12")

# Specify the independent columns (relative height metrics)
target_cols <- c("rh10","rh15","rh20","rh25","rh30","rh35","rh40","rh45","rh50","rh55", "rh60","rh65", 
                 "rh70", "rh75", "rh80", "rh85", "rh90", "rh95","rh98","rh99")

# Initialize a list to store models
models_list <- list()
eval_results <- data.frame(rh=character(),RMSE.md=numeric(),RMSE.vd=numeric(),Rsquared=numeric(),mape=numeric())

#Train-test dataset split
smp <- smpdf[,1]
train <- d.final[smp,]
test <- d.final[!d.final$X %in% train$X ,]

# Iterate over multiple relative height metrics to fit Random Forest models
for (target in target_cols) {
  print(paste0("Modelling attempt ", target," Started!"))
  print(start_time <- Sys.time())
  formula_str <- paste(target, "~", paste(predictor_cols, collapse = " + "))
  model_formula <- as.formula(formula_str)
  
  # Fit the Random Forest model
  set.seed(123)
  model <- randomForest(model_formula, data = train)
  
  # Store the model in the list
  models_list[[target]] <- model
  
  # Predict on test set
  pred <- predict(model, test)
  rmse <- sqrt(mean((test[[target]] - pred)^2))
  #rsq <- cor(test[[target]], pred)^2
  #eval_results[[target]] <- list(RMSE = rmse, Rsquared = rsq)
  
  cat("\nRandom Forest model for target:", target, "\n")
  print(models_list[[target]])
  print(paste0("Modelling attempt ", target," Finished!"))

  # Predict on full coverage dataset
  print(paste0("Modelling W2W attempt ", target," Started!"))
  print(start_time <- Sys.time())
  set.seed(123)
  prdw2w <- predict(model, newdata = w2w.filt)
  w2w.pred<-cbind(w2w.filt,prdw2w)
  w2w.pred.fin=w2w.pred%>%dplyr::select(x,y,prdw2w)
  w2w.pred.raster <- terra::rast(w2w.pred.fin, type = "xyz")
  sample.rst <- terra::extract(w2w.pred.raster, vect)
  sample.rst.join=sample.rst%>%left_join(sf,by="ID")
  sample.rst.join<- na.omit(sample.rst.join)
  rmse.tr <- sqrt(mean((sample.rst.join$mean_x - sample.rst.join$prdw2w)^2))
  mape=mean(abs((sample.rst.join$mean_x - sample.rst.join$prdw2w) / sample.rst.join$mean_x)) * 100
  rsq=summary(lm(sample.rst.join$mean_x ~ sample.rst.join$prdw2w))$r.squared
  #eval_results[[target]] <- list(RMSE.md = rmse,RMSE.vd = rmse.tr, Rsquared = rsq)
  cat("\nEvaluation for target:", target, "\n")
  
  # Record model validation result
  new<-data.frame(rh=target,RMSE.md=rmse,RMSE.vd=rmse.tr,Rsquared=rsq,mape=mape)
  eval_results <- rbind(eval_results,new)
  print(eval_results%>%filter(rh==target))
  
  # Export Canopy Height Model (raster)
  writeRaster(w2w.pred.raster,paste0("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\w2w_res_20251209_",target,".tif"),overwrite=TRUE)
  print(paste0("Modelling W2W attempt ", target," Finished!"))
  
  # Export Scatterplot of predicted vs reference
  text=paste0("RMSE: ",round(rmse.tr,digit=2)," m")
  plot <- ggplot(sample.rst.join,aes(x=mean_x, y=prdw2w)) +
    geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
    geom_abline(a=0,b=0,linetype=1,colour="red")+
    labs(title = target,x = "Actual (m)", y = "Predicted (m)") +
    coord_fixed()+  # Set aspect ratio t  o be equal
    xlim(0,30) +
    ylim(0,30)
  plot+theme_bw()+theme(plot.title = element_text(hjust = 0.5))+ annotate("text", x = 15, y = 25, label = text,col="red")
  ggsave(paste0("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\plot_20251209_",target,".png"), plot = plot)
  assign(paste0("plot.",target),plot)
              
}

## Export Model Validation Table
write.table(eval_results,"H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\20251209_rhmodel.csv",sep=";",row.names = FALSE,quote = FALSE)

## Optional - Create a line graph from the evaluation result
cur.eval <- eval_results[27:46, ]

# create rh.num
cur.eval$rh.num <- as.numeric(substr(cur.eval$rh, start = 3, stop = 5))
cur.eval$rh.num[19] <- 100
cur.eval$rh.num[20] <- 105

# pivot longer for plotting
cur.long <- tidyr::pivot_longer(
  cur.eval,
  cols = c(RMSE.vd, RMSE.md),
  names_to = "set",
  values_to = "RMSE"
)

cur.long$set <- factor(cur.long$set,
                       levels = c("RMSE.vd", "RMSE.md"),
                       labels = c("Field data", "Test set"))

ggplot(cur.long, aes(x = rh.num, y = RMSE, color = set)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  theme_bw(base_size = 12) +
  labs(
    title = "Model Performance across RH Metrics",
    x = "Relative Height Metrics",
    y = "RMSE (m)",
    color = "Dataset"
  ) +
  scale_x_continuous(
    breaks = sort(unique(cur.eval$rh.num)),
    labels = cur.eval$rh
  ) + 
  scale_y_continuous(
    breaks = seq(0, max(cur.long$RMSE) + 1, by = 1),
    limits = c(0, NA)
  )+
  theme(
    plot.title = element_text(hjust=0.5,face = "bold"),
    # strong rectangular border around panel
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    # remove grey background explicitly (keeps minimal style)
    panel.background = element_rect(fill = "white", colour = NA),
    # darker axis lines for stern look
    axis.line = element_line(colour = "black", linewidth = 0.4),
    # optional: tone down grid
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

###Comparison with Global Products

#Load CHM products
tolan=rast("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\tolan_chm.tif")
potapov=rast("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\potapov_chm_30m.tif")
simard=rast("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\simard_chm.tif")
lang=rast("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\lang_chm.tif")

#Extract raster based in field plot location
sample.tolan <- terra::extract(tolan, vect)
sample.potapov <- terra::extract(potapov, vect)
sample.simard <- terra::extract(simard, vect)
sample.lang <- terra::extract(lang, vect)

#Combine result into one dataframe
globchm=cbind(sample.tolan,sample.potapov,sample.simard,sample.lang)
colnames(globchm)[2]="tolan_chm"
colnames(globchm)[4]="potapov_chm"
colnames(globchm)[3]="ID_2"
colnames(globchm)[5]="ID_3"
colnames(globchm)[7]="ID_4"
colnames(globchm)[8]="lang_chm"
globchm=globchm %>% left_join(sample.rst.join,by="ID") %>%  na.omit()
colnames(globchm)[12]="mH"
globchm$zachary_chm=globchm$prdw2w

#Calculate Evaluation Metrics
rmse.tolan=sqrt(mean((globchm$mH - globchm$tolan_chm)^2))
rmse.potapov=sqrt(mean((globchm$mH - globchm$potapov_chm)^2))
rmse.simard=sqrt(mean((globchm$mH - globchm$simard_chm)^2))
rmse.zachary=sqrt(mean((globchm$mH - globchm$prdw2w)^2))
rmse.lang=sqrt(mean((globchm$mH - globchm$lang_chm)^2))

r2.tolan=cor(globchm$mH,globchm$tolan_chm)^2
r2.potapov=cor(globchm$mH,globchm$potapov_chm)^2
r2.simard=cor(globchm$mH,globchm$simard_chm)^2
r2.zachary=cor(globchm$mH,globchm$prdw2w)^2
r2.lang=cor(globchm$mH,globchm$lang_chm)^2

chm_cols <- c("zachary_chm","tolan_chm","potapov_chm","simard_chm","lang_chm")

#arrange results
globchm.res=data.frame(
  model=chm_cols,
  rmse=c(rmse.zachary,rmse.tolan,rmse.potapov,rmse.simard,rmse.lang),
  r2=c(r2.zachary,r2.tolan,r2.potapov,r2.simard,r2.lang)
)

#Save result
write.table(globchm.res,"H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\globalcmh_20251224.csv",sep=";",row.names = FALSE,quote = FALSE)

#Plot Graph
plot.lang_chm <- ggplot(globchm,aes(x=mH, y=lang_chm)) +
  geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
  geom_abline(a=0,b=0,linetype=1,colour="red")+
  annotate(geom="text", x=15, y=30, label=paste0("RMSE: ",round(rmse.lang,2)," m"),
           color="black")+
  labs(title ="Lang",x = "Actual (m)", y = "Predicted (m)") +
  coord_fixed()+  # Set aspect ratio t  o be equal
  xlim(0,30) +
  ylim(0,30)

plot.simard_chm <- ggplot(globchm,aes(x=mH, y=simard_chm)) +
  geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
  geom_abline(a=0,b=0,linetype=1,colour="red")+
  annotate(geom="text", x=15, y=30, label=paste0("RMSE: ",round(rmse.simard,2)," m"),
           color="black")+
  labs(title ="Simard",x = "Actual (m)", y = "Predicted (m)") +
  coord_fixed()+  # Set aspect ratio t  o be equal
  xlim(0,30) +
  ylim(0,30)

plot.potapov_chm <- ggplot(globchm,aes(x=mH, y=potapov_chm)) +
  geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
  geom_abline(a=0,b=0,linetype=1,colour="red")+
  annotate(geom="text", x=15, y=30, label=paste0("RMSE: ",round(rmse.potapov,2)," m"),
           color="black")+
  labs(title ="Potapov",x = "Actual (m)", y = "Predicted (m)") +
  coord_fixed()+  # Set aspect ratio t  o be equal
  xlim(0,30) +
  ylim(0,30)

plot.tolan_chm <- ggplot(globchm,aes(x=mH, y=tolan_chm)) +
  geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
  geom_abline(a=0,b=0,linetype=1,colour="red")+
  annotate(geom="text", x=15, y=30, label=paste0("RMSE: ",round(rmse.tolan,2)," m"),
           color="black")+
  labs(title ="Tolan",x = "Actual (m)", y = "Predicted (m)") +
  coord_fixed()+  # Set aspect ratio t  o be equal
  xlim(0,30) +
  ylim(0,30)

plot.zachary_chm <- ggplot(globchm,aes(x=mH, y=zachary_chm)) +
  geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
  geom_abline(a=0,b=0,linetype=1,colour="red")+
  annotate(geom="text", x=15, y=30, label=paste0("RMSE: ",round(rmse.zachary,2)," m"),
           color="black")+
  labs(title ="Afif",x = "Actual (m)", y = "Predicted (m)") +
  coord_fixed()+  # Set aspect ratio t  o be equal
  xlim(0,30) +
  ylim(0,30)

plots<- map(chm_cols, ~ {
  ggplot(globchm,aes(x=mH, y=.data[[.x]])) +
    geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
    geom_abline(a=0,b=0,linetype=1,colour="red")+
    labs(title = paste0(.x),x = "Actual (m)", y = "Predicted (m)") +
    coord_fixed()+  # Set aspect ratio t  o be equal
    xlim(0,30) +
    ylim(0,30)+
    theme_minimal()
})

#arrange graph into grid
grid.arrange(plots[[1]],
             plots[[2]],
             plots[[3]],
             plots[[4]],
             plots[[5]],ncol = 5)

grid.arrange(plot.zachary_chm,
             plot.lang_chm,
             plot.potapov_chm,
             plot.simard_chm,
             plot.tolan_chm,
              ncol=5)

#Extract raster based on prediction result pixel location
pred_sf.cl$ID = 1:nrow(pred_sf.cl)
tolan.w2w <- terra::extract(tolan, vect(pred_sf.cl))
potapov.w2w <- terra::extract(potapov, vect(pred_sf.cl))
simard.w2w <- terra::extract(simard, vect(pred_sf.cl))
lang.w2w <- terra::extract(lang, vect(pred_sf.cl))

#merge into one dataframe
globchm.w2w=cbind(tolan.w2w,potapov.w2w,simard.w2w,lang.w2w)
colnames(globchm.w2w)[2]="tolan_chm"
colnames(globchm.w2w)[4]="potapov_chm"
colnames(globchm.w2w)[3]="ID_2"
colnames(globchm.w2w)[5]="ID_3"
colnames(globchm.w2w)[7]="ID_4"
colnames(globchm.w2w)[8]="lang_chm"
globchm.w2w=globchm.w2w %>% left_join(pred_sf.cl,by="ID") 

#Pivot longer for plottinh
globalchm.w2w.longer=globchm.w2w %>%
  pivot_longer(
    cols = c(prdw2w,tolan_chm, potapov_chm,simard_chm,lang_chm),
    names_to = "chm",
    values_to = "height"
  )

#Column and row renaming
globalchm.w2w.longer=globalchm.w2w.longer[globalchm.w2w.longer$class!=6,]
globalchm.w2w.longer$chm[globalchm.w2w.longer$chm == "prdw2w"] <- "zachary_chm"
globalchm.w2w.longer$chm[globalchm.w2w.longer$chm == "zachary_chm"] <- "afif_chm"
globalchm.w2w.longer$class_name[globalchm.w2w.longer$class_name == "Class 7"] <-  "Class 6"
globalchm.w2w.longer$class_name=paste0("Class ",globalchm.w2w.longer$class)

#Plot barplot of height distribution based on species class
order.bar<- c("Class 1", "Class 2", "Class 3", "Class 4", "Class 5", "Class 6")
ggplot(globalchm.w2w.longer, aes(fill = chm, y = height, x = factor(class_name, levels = order.bar))) + 
  geom_bar(position = position_dodge(0.8), 
           stat = "summary", 
           fun = "median", 
           alpha = 0.8, 
           width = 0.7) +
  scale_fill_viridis_d(option = "plasma", name = "Model") +
  labs(
    title = "Height Distribution by Class",
    x = "Species Communities Class",
    y = "Height (m)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 15)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey50"),
    
    # Professional panel styling
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.3),
    
    # Axis styling
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    
    # Legend
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.margin = margin(b = 10),
    
    # Margins
    plot.margin = margin(15, 15, 15, 15)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
