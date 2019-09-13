library("ggplot2")
library("ggrepel")
library("rgdal")
library("raster")
library("dplyr")

data<-read.csv("data/SNP_Survey_Table1-1.csv", header=TRUE, sep=",")
data$Latitude<-as.numeric(as.character(data$Latitude))
data$Longitude<-as.numeric(as.character(data$Longitude))

d1<-subset(data, data$Migratory_Access=="Below barrier")
d2<-subset(d1, Population!="Kamchatka (anadromous)")
d3<-subset(d2, Population!="Kamchatka (nonanadromous)")
d4<-subset(d3, Population!="Columbia, Snake, Salmon, Lemhi")
dataSub<-d4
df<-tbl_df(dataSub)
df<-select(df, Pop_ID, InversionFreq, N, Latitude, Longitude)


#### DOWNLOAD WORLD CLIM DATA TO CONTINUE 
# Download tmean_*.bil from https://www.worldclim.org/
months<-c("January","February","March","April","May","June","July","August","September","October","November","December")
merged=raster(paste("./tmean_","1",".bil",sep=""))
i<-1
df1<-df %>% mutate(Temp= extract(merged, cbind(df$Longitude,df$Latitude),method='bilinear')) %>% mutate(Month = paste(i)) 
dfFull<-df1

for (i in 2:12) {
  merged=raster(paste("./tmean_",i,".bil",sep=""))
  df2<-df %>% mutate(Temp=extract(merged, cbind(df$Longitude,df$Latitude), method='bilinear')) %>% mutate(Month = paste(i))
  dfFull<-full_join(dfFull,df2, by=c("Pop_ID","InversionFreq","N", "Latitude","Longitude","Temp","Month"))
}

dfFull<-arrange(dfFull, as.numeric(Month))

dat<-dfFull
dat$Group<-dfFull$Month
Group<-dat$Group
#dat <- within(dat, Group <- factor(Group, Group))
#Changed 11/5/2018
#dat <- within(dat, Group <- factor(Group))
dat$Group<-factor(dat$Group, levels=as.factor(sort(unique(Group))))

model<-lm(InversionFreq ~ Temp + factor(Month), data=dat, weights=N)
predicted<-stats::predict(model, dat)

#To get Rsquared
#modded from http://stackoverflow.com/questions/19699858/ggplot-adding-regression-line-equation-and-r2-with-facet
lm_eqn=function(df){
  m=lm(InversionFreq ~ Temp, weights=N,data=df)
  eq <- substitute("Adj. R²="~r2,
                   #italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(
                     #a = format(coef(m)[1], digits = 2), 
                     #b = format(coef(m)[2], digits = 2), 
                     r2 = format(summary(m)$adj.r.squared, digits = 2)))
  as.character(as.expression(eq));                 
}

#lm_eqn(subset(dat, dat$Group==1))
rsquares<-months;
for (i in 1:12) {
  rsquares[i]<-paste(rsquares[i],lm_eqn(subset(dat, dat$Group==i)))
}

sub1<-gsub("\"","",rsquares)
sub2<-gsub(" ~ ","",sub1)
#Convert temp
dat$Temp<-(dat$Temp/10)

pdf("E9-temp-plot.pdf")

ggplot(dat, aes(x=Temp, y=InversionFreq))+geom_point(aes(size=N),alpha=0.5)+
  #geom_line(aes(x=dfFull$Temp, y=predicted))+
  facet_wrap( ~ Group, ncol=3, labeller=as_labeller(setNames(sub2, unique(Group))))+
  geom_smooth(aes(Temp, InversionFreq, weight=N), method=lm,se=FALSE,fullrange=TRUE)+
  ylim(0,1)+
  labs(title="Inversion Frequency as a Function of Mean Monthly Temperature", x="Mean Monthly Temperature °C", y="Inversion Frequency")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=rel(1.3)),
        axis.title.x = element_text(size = rel(1.1)),
        axis.title.y = element_text(size = rel(1.1)),
        legend.title.align=0.5 )
dev.off()

#Mean Monthly Temperature °C", Inversion Frequency"