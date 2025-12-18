
tolan=rast("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\tolan_chm.tif")
potapov=rast("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\potapov_chm_30m.tif")
simard=rast("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\simard_chm.tif")

sample.tolan <- terra::extract(tolan, vect)
sample.potapov <- terra::extract(potapov, vect)
sample.simard <- terra::extract(simard, vect)

globchm=cbind(sample.tolan,sample.potapov,sample.simard)
colnames(globchm)[2]="tolan_chm"
colnames(globchm)[4]="potapov_chm"
colnames(globchm)[3]="ID_2"
colnames(globchm)[5]="ID_3"
colnames(globchm)[10]="mH"

globchm=globchm %>% left_join(sample.rst.join,by="ID") %>%  na.omit()


sqrt(mean((globchm$mH - globchm$tolan_chm)^2))
sqrt(mean((globchm$mH - globchm$potapov_chm)^2))
sqrt(mean((globchm$mH - globchm$simard_chm)^2))
sqrt(mean((globchm$mH - globchm$prdw2w)^2))

cor(globchm$mH,globchm$prdw2w)^2
cor(globchm$mH,globchm$tolan_chm)^2
cor(globchm$mH,globchm$potapov_chm)^2
cor(globchm$mH,globchm$simard_chm)^2

cor(globchm.cln$Elevation,globchm.cln$mH)

varImpPlot(model,main="Variable Importance")

cor(sample.rst.join$prdw2w,sample.rst.join$mH)

chm_cols <- c("zachary_chm","tolan_chm","potapov_chm","simard_chm")

globchm$zachary_chm=globchm$prdw2w

plot.zachary_chm <- ggplot(globchm,aes(x=mH, y=prdw2w)) +
  geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
  geom_abline(a=0,b=0,linetype=1,colour="red")+
  labs(title = "zachary_chm",x = "Actual (m)", y = "Predicted (m)") +
  coord_fixed()+  # Set aspect ratio t  o be equal
  xlim(0,30) +
  ylim(0,30)

library(gridExtra)
grid.arrange(plot.zachary_chm,
             plot.tolan_chm,
             plot.potapov_chm,
             plot.simard_chm,ncol = 4)

write.table(globchm,"H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\globalcmh_20251214.csv",sep=";",row.names = FALSE,quote = FALSE)