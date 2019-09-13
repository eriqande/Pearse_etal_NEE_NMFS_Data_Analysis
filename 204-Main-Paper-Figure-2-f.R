library("ggplot2")
data<-read.csv("data/SNP_Survey_Table1-1.csv", header=TRUE, sep=",")
data$Latitude<-as.numeric(as.character(data$Latitude))
data$Longitude<-as.numeric(as.character(data$Longitude))

d1<-subset(data, data$Migratory_Access=="Below barrier")
d2<-subset(d1, Population!="Kamchatka (anadromous)")
d3<-subset(d2, Population!="Kamchatka (nonanadromous)")
d4<-subset(d3, Population!="Columbia, Snake, Salmon, Lemhi")
dataSub<-d4
x<-dataSub$Long
y<-dataSub$Lat
xy<-cbind(x,y)

variance<-(dataSub$InversionFreq*(1-dataSub$InversionFreq))/(2*dataSub$N)
SEp<-sqrt((dataSub$InversionFreq*(1-dataSub$InversionFreq))/(2*dataSub$N) )

pdf("Figure2E.pdf", width=(11/2), height=(8.5/2))

fit<-lm(dataSub$InversionFreq ~ dataSub$Latitude, weights=dataSub$N)

summary<-summary(fit)
summary(fit)

coeff=coefficients(fit)
proportion=(dataSub$N/max(dataSub$N))

ggplot(data=dataSub, aes(dataSub$Latitude,dataSub$InversionFreq))+
  #geom_point(color="red", alpha=0.5, size=5*proportion)+
  geom_point(aes(size=N),alpha=0.5, col="red")+
  geom_abline(intercept=fit$coefficients[1],slope=fit$coefficients[2],linetype="solid",size=1,color="black", alpha=0.5)+
  theme_classic()+
  labs(title="Plot of Inversion Frequency as Function of Latitude" 
       # , subtitle=paste("y = ",round(coeff[2],digits=2),"x",round(coeff[1],digits=2),", Adj. R² = ",
       #                round(summary$adj.r.squared,digits=2),sep="")
  )+
  theme(plot.title = element_text(hjust = 0.5, size=rel(1.4)),
        #plot.subtitle = element_text(hjust = 0.5, size=rel(1)),
        axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2)),
        legend.title.align=0.5 )+
  xlab("Latitude")+
  ylab("Inversion Frequency")+
  geom_errorbar(aes(x=dataSub$Latitude, ymin=dataSub$InversionFreq-SEp, ymax=dataSub$InversionFreq+SEp), alpha=0.5)+
  scale_size( breaks=c(5,20,40,60))

dev.off()

