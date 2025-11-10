library(sf)
library(dplyr)
library(terra)
library(caret)
library(randomForest)
library(ModelMetrics)
library(parallel)
library(doParallel)

##Load dataset and filter data
d<-read.csv2("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\dataset_v10_fin.csv")
d= type.convert(d, as.is = TRUE)
d$h_canopy_num<-as.numeric(d$h_canopy)

d$NDVI=(d$b8-d$b4)/(d$b8+d$b4)
d$GNDVI=(d$b8-d$b3)/(d$b8+d$b3)
d$SAVI=(d$b8-d$b4)/(d$b8+d$b4+0.5)*(1+0.5)
d$NDII=(d$b8-d$b11)/(d$b8+d$b11)


d.filt <- subset(d, d$b2 > -1000 # Filter dataset based on sentinel-2 that has nan value
                 & d$b3 > -1000
                 & d$b4 > -1000
                 & d$b8 > -1000
                 & d$b11 > -1000
                 & d$b12 > -1000, )
d.filt2<-d.filt %>% filter(vh < 0 | vv < 0) # Filter dataset based on sentinel-1 that has nan value
d.filt3<-d.filt2%>% filter(h_canopy_num>2,h_canopy_num<25,NDVI != -9999) 

vars=c("class","vv","vh","b2","b3","b4","b8","b11","b12","NDVI","GNDVI","SAVI","NDII","h_canopy_num")

d.final<-d.filt3[vars]
d.final$X<-1:dim(d.final)[1]


##Train-Test Split using Stratified Sampling based on species class

a <- table(d.final$class) #Count frequency of each class
prop <- 0.7 #Set the proportion for training data
freqclass <- as.numeric(a) 
trfreq <- round(freqclass * prop, digits = 0) #Calculate how many samples per class
d.fin <- d.final[order(d.final$class),] 
sizeArg <- trfreq

# Create stratified sample data frame
rows <- dim(d.final)[1]
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
for (i in 1:5) {
  smp <- smpdf[,i]
  trainPnt <- d.final[smp,]
  
  testPnt <- d.final[!d.final$X %in% trainPnt$X ,]
  
  traindat<-trainPnt%>%select(-h_canopy_num,-X)
  
  testdat<-testPnt%>%select(-h_canopy_num,-X)
  
  cntr <- trainControl(method='repeatedcv',
                       number=3,
                       repeats=3)
  
  print(paste0("Modelling attempt ", i," Started!"))
  print(start_time <- Sys.time())
  RF <- caret::train(traindat,trainPnt$h_canopy_num,method='rf', importance=T, trControl = cntr)
  print(end_time <- Sys.time())
  print(paste0("Modelling attempt ", i," Finished!"))
      
  train.prdRF <- predict(RF)
  test.prdRF <- predict(RF,newdata=testPnt) #Predict Testing dataset
  
  #Calculate model metrics
  MAPE = (1/nrow(trainPnt)) * sum(abs(trainPnt$h_canopy_num - train.prdRF)/ trainPnt$h_canopy_num) * 100
  MAPE.test = (1/nrow(testPnt)) * sum(abs(testPnt$h_canopy_num - test.prdRF)/ testPnt$h_canopy_num) * 100
  rmse.train=rmse(trainPnt$h_canopy_num,train.prdRF)
  rmse.test=rmse(testPnt$h_canopy_num,test.prdRF)
  mae.train=mae(trainPnt$h_canopy_num,train.prdRF)
  mae.test=mae(testPnt$h_canopy_num,test.prdRF)
  
  
  saveRDS(RF, file = paste0("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\rf_model_tr_",i,".rds"))
  
  smpNum <- prop
  # write results to data frame
  new<-data.frame(mdlname='RF',iter=i,smp=smpNum,
    rmse.train=round(rmse.train,3),mae.train=round(mae.train,3),mape.train=round(MAPE,3),
    rmse.test=round(rmse.test,3),mae.test=round(mae.test,3),mape.test=round(MAPE.test,3))
  mdlRes <- rbind(mdlRes,new)
  
  print(paste0("Attempts:" ,i))
  print(paste0("RMSE Train: ",round(rmse.train,3)))
  print(paste0("RMSE Test: ",round(rmse.test,3)))
  print(paste0("MAPE Train: ",round(MAPE,3)))
  print(paste0("MAPE Test: ",round(MAPE.test,3)))
  
}

##predict wall to wall dataset

filename <- file.choose()
RF <- readRDS(filename)

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

w2w_c <- w2w %>%
  mutate(
    NDVI = (b8 - b4) / (b8 + b4),
    GNDVI = (b8 - b3) / (b8 + b3),
    SAVI = ((b8 - b4) / (b8 + b4 + 0.5)) * (1 + 0.5),
    NDII = (b8 - b11) / (b8 + b11),
  )

w2w.filt<-w2w_c%>%filter(!is.na(class),class!=0)


prdw2w <- predict(RF, newdata = w2w.filt)

w2w.pred<-cbind(w2w.filt,prdw2w)

w2w.pred.fin=w2w.pred%>%select(x,y,prdw2w)

w2w.pred.raster <- terra::rast(w2w.pred.fin, type = "xyz")

writeRaster(w2w.pred.raster,"H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\w2w_res_20251106.tif")

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