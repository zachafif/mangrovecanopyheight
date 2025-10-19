library(sf)
library(dplyr)
library(terra)
library(caret)
library(randomForest)
library(ModelMetrics)
library(parallel)
library(doParallel)

##Load dataset and filter data
d<-read.csv2("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\dataset_v10.csv")
d= type.convert(d, as.is = TRUE)
d$h_canopy_num<-as.numeric(d$h_canopy)
d.filt <- subset(d, d$b2 > -1000 # Filter dataset based on sentinel-2 that has nan value
                 & d$b3 > -1000
                 & d$b4 > -1000
                 & d$b8 > -1000
                 & d$b11 > -1000
                 & d$b12 > -1000, )
d.filt2<-d.filt %>% filter(lo_vh < 0, lo_vv < 0,psar_hh >0 ,psar_hv> 0) # Filter dataset based on sentinel-1 that has nan value
d.filt3<-d.filt2%>% filter(h_canopy_num>2,h_canopy_num<25,NDVI != -9999) 
d.filt4<-d.filt3%>% filter(chm_meta>2,abs(chm_meta-h_canopy_num)<5) %>% select(-chm)

##Train-Test Split using Stratified Sampling based on species class

a <- table(d.filt4$class) #Count frequency of each class
prop <- 0.7 #Set the proportion for training data
freqclass <- as.numeric(a) 
trfreq <- round(freqclass * prop, digits = 0) #Calculate how many samples per class
d.fin <- d.filt4[order(d.filt4$class),] 
sizeArg <- trfreq

# Create stratified sample data frame
rows <- dim(d.filt4)[1]
iter<-50
smpdf = data.frame(matrix(NA, nrow=round((rows * prop)),ncol=0))

#Run stratified sampling in a loop
for (i in 1:50){
 smp <- sampling:::strata(d.fin,'class', sizeArg, method=c('srswor'))
 smpind <- smp$ID_unit
 smpdf <- cbind(smpdf, smpind)
}
colnames(smpdf)  <- c(paste('smp_',c(1:iter),sep=''))

##Model Training and Testing
#Create empty dataframe to store result
mdlRes <- data.frame(mdlNme=character(),
                     iter=numeric(),
                     smp=numeric(),
                     rmse_cal=numeric(),
                     mae_cal=numeric(),
                     mape_cal=numeric(),
                     rmse_val=numeric(),
                     mae_val=numeric(),
                     mape_val=numeric())
#Loop model training 
for (i in 1:3) {
  smp <- smpdf[,i]
  calPnt <- d.filt4[smp,]
  
  valPnt <- d.filt4[!d.filt4$X %in% calPnt$X ,]
  
  caldat<-calPnt%>%select(-h_canopy)%>%select(-h_canopy_num)%>%select(-X)
  
  valdat<-valPnt%>%select(-h_canopy)%>%select(-h_canopy_num)%>%select(-X)
  
  cntr <- trainControl(method='repeatedcv',
                       number=3,
                       repeats=3)
  
  print(paste0("Modelling attempt ", i," Started!"))
  print(start_time <- Sys.time())
  RF <- caret::train(caldat,calPnt$h_canopy,method='rf', importance=T, trControl = cntr)
  print(end_time <- Sys.time())
  print(paste0("Modelling attempt ", i," Finished!"))
      
  cal.prdRF <- predict(RF)
  val.prdRF <- predict(RF,newdata=valPnt) #Predict Testing dataset
  
  #Calculate model metrics
  MAPE = (1/nrow(calPnt)) * sum(abs(calPnt$h_canopy_num - cal.prdRF)/ calPnt$h_canopy_num) * 100
  MAPE.val = (1/nrow(valPnt)) * sum(abs(valPnt$h_canopy_num - val.prdRF)/ valPnt$h_canopy_num) * 100
  rmse.cal=rmse(calPnt$h_canopy_num,cal.prdRF)
  rmse.val=rmse(valPnt$h_canopy_num,val.prdRF)
  mae.cal=mae(calPnt$h_canopy_num,cal.prdRF)
  mae.val=mae(valPnt$h_canopy_num,val.prdRF)
  
  
  #saveRDS(RF, file = "H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\rf_model_v7.rds")
  
  smpNum <- prop
  # write results to data frame
  new<-data.frame(mdlname='RF',iter=i,smp=smpNum,
    rmse_cal=round(rmse.cal,3),mae_cal=round(mae.cal,3),mape_cal=round(MAPE,3),
    rmse_val=round(rmse.val,3),mae_val=round(mae.val,3),mape_val=round(MAPE.val,3))
  mdlRes <- rbind(mdlRes,new)
  
  print(i)
  print(round(rmse.cal,3))
  print(round(rmse.val,3))
  print(round(MAPE,3))
  print(round(MAPE.val,3))
  
}

##predict wall to wall dataset

r<-terra::rast("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\stacked_w2w_v3.tif")
w2w <- as.data.frame(r, xy = TRUE)

colnames(w2w)[3]<-'chm'
colnames(w2w)[4]<-'psar_hh'
colnames(w2w)[5]<-'psar_hv'
colnames(w2w)[6]<-'lo_vv'
colnames(w2w)[7]<-'lo_vh'
colnames(w2w)[8]<-'b2'
colnames(w2w)[9]<-'b3'
colnames(w2w)[10]<-'b4'
colnames(w2w)[11]<-'b8'
colnames(w2w)[12]<-'b11'
colnames(w2w)[13]<-'b12'
colnames(w2w)[31]<-'chm_meta'
colnames(w2w)[32]<-'class'

w2w_c <- w2w %>%
  mutate(
    CIg = (b8 / b3) - 1,
    DVI = b8 - b4,
    EVI = 2.5 * (b8 - b4) / (b8 + 6 * b4 - 7.5 * b2 + 1),
    FDI = b8 - (b3 + b4),
    NDVI = (b8 - b4) / (b8 + b4),
    SR = b8 / b4,
    TNDVI = sqrt((b8 - b4) / (b8 + b4) + 0.5),
    GNDVI = (b8 - b3) / (b8 + b3),
    SAVI = ((b8 - b4) / (b8 + b4 + 0.5)) * (1 + 0.5),
    NDII = (b8 - b11) / (b8 + b11),
    MDI1 = (b8 - b11) / b11,
    MDI2 = (b8 - b12) / b12,
    MNDWI = (b3 - b11) / (b3 + b11),
    NDI = (lo_vv - lo_vh) / (lo_vv + lo_vh),
    RVI = (4 * lo_vh) / (lo_vv + lo_vh),
    RI = lo_vh / lo_vv
  )

w2w.filt<-w2w_c%>%filter(!is.na(class))

prdw2w <- predict(RF, newdata = w2w.filt)

w2w.pred<-cbind(w2w.filt,prdw2w)

w2w.pred.fin=w2w.pred%>%select(x,y,prdw2w)

w2w.pred.raster <- terra::rast(w2w.pred.fin, type = "xyz")

writeRaster(w2w.pred.raster,"H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\w2w_res_20251018.tif")

plot(w2w.pred.raster,main="Mangrove Height Model - Trat, Thailand")



## Optional - Create Scatterplot

text=paste0("RMSE: ",round(rmse.val,digit=2)," m")

res<-data.frame(ref=valPnt$h_canopy_num,pred=val.prdRF)

plot <- ggplot(res,aes(x=ref, y=pred)) +
  geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
  geom_abline(a=0,b=0,linetype=1,colour="red")+
  labs(title = "Predicted vs Actual",x = "Actual (m)", y = "Predicted (m)") +
  coord_fixed()+  # Set aspect ratio t  o be equal
  xlim(0,30) +
  ylim(0,30)

plot+theme_bw()+theme(plot.title = element_text(hjust = 0.5))+ annotate("text", x = 15, y = 25, label = text,col="red")