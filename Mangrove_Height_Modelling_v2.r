library(sf)
library(dplyr)
library(terra)
library(randomForest)
library(ModelMetrics)

##Load dataset and filter data
d<-read.csv2("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\dataset_v10_fin.csv")
d= type.convert(d, as.is = TRUE)
sf=st_read("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\transect.shp")
sf$ID=1:nrow(sf)
vect=vect(sf)

d <- d %>% ##Calculate indices (subject to change)
  mutate(
    CIg = (b8 / b3) - 1,
    NDVI = (b8 - b4) / (b8 + b4),
    GNDVI = (b8 - b3) / (b8 + b3),
    SAVI = ((b8 - b4) / (b8 + b4 + 0.5)) * (1 + 0.5),
    NDII = (b8 - b11) / (b8 + b11),
    MDI1 = (b8 - b11) / b11,
    MDI2 = (b8 - b12) / b12,
    MNDWI = (b3 - b11) / (b3 + b11),
    NDI = (vv - vh) / (vv + vh),
    RVI = (4 * vh) / (vv + vh),
    RI = vh / vv
  )

d.filt <- subset(d, d$b2 > -1000 # Filter dataset based on sentinel-2 that has nan value
                 & d$b3 > -1000
                 & d$b4 > -1000
                 & d$b8 > -1000
                 & d$b11 > -1000
                 & d$b12 > -1000, )
d.filt2<-d.filt %>% filter(vh < 0 | vv < 0) # Filter dataset based on sentinel-1 that has nan value
d.filt3<-d.filt2%>% filter(h_canopy>0,h_canopy<30,NDVI != -9999) 
d.filt3$ort_h=d.filt3$h_te_median+22.021
d.filt4<-d.filt3%>% filter(ort_h<20,class !=6)
d.final<-d.filt4[c(2,22:30,53:78)] #all features
d.final$X<-1:dim(d.final)[1]


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

# Specify the names of dependent var columns
predictor_cols <- c( "vv","vh","b2","b3","b3","b4","b8","b11","b12",
"CIg","NDVI","SAVI","NDII","MDI1","MDI2","MNDWI","RVI","RI")

# Specify the independent columns (relative height metrics)
target_cols <- c("rh10","rh15","rh20","rh25","rh30","rh35","rh40","rh45","rh50","rh55", "rh60","rh65", 
                 "rh70", "rh75", "rh80", "rh85", "rh90", "rh95","rh98","rh99")

# Initialize a list to store models
models_list <- list()
eval_results <- data.frame(rh=character(),RMSE.md=numeric(),RMSE.vd=numeric(),Rsquared=numeric())

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
  rsq <- cor(test[[target]], pred)^2
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
  rmse.tr <- sqrt(mean((sample.rst.join$mH - sample.rst.join$prdw2w)^2))
  #eval_results[[target]] <- list(RMSE.md = rmse,RMSE.vd = rmse.tr, Rsquared = rsq)
  cat("\nEvaluation for target:", target, "\n")
  
  # Record model validation result
  new<-data.frame(rh=target,RMSE.md=rmse,RMSE.vd=rmse.tr,Rsquared=rsq)
  eval_results <- rbind(eval_results,new)
  print(eval_results%>%filter(rh==target))
  
  # Export Canopy Height Model (raster)
  writeRaster(w2w.pred.raster,paste0("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\w2w_res_20251127_",target,".tif"))
  print(paste0("Modelling W2W attempt ", target," Finished!"))
  
  # Export Scatterplot of predicted vs reference
  text=paste0("RMSE: ",round(rmse.tr,digit=2)," m")
  plot <- ggplot(sample.rst.join,aes(x=mH, y=prdw2w)) +
    geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
    geom_abline(a=0,b=0,linetype=1,colour="red")+
    labs(title = target,x = "Actual (m)", y = "Predicted (m)") +
    coord_fixed()+  # Set aspect ratio t  o be equal
    xlim(0,30) +
    ylim(0,30)
  plot+theme_bw()+theme(plot.title = element_text(hjust = 0.5))+ annotate("text", x = 15, y = 25, label = text,col="red")
  ggsave(paste0("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\plot_20251127_",target,".png"), plot = plot)
  assign(paste0("plot.",target),plot)
              
}

## Export Model Validation Table
write.table(eval_results,"H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\20251127_rhmodel.csv",sep=";",row.names = FALSE,quote = FALSE)