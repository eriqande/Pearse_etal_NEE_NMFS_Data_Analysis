library("adegenet")
library("ape")
library("dplyr")
library("plyr")
library("ggplot2")
library("poppr")
library("mmod")
library("ggrepel")
library("ggtree")

data<-read.genepop("data/surveyEdited.gen")
toKeep=c("R14383", "R38346", "R30551", "R19198", "R30355", "R21023", "R40252", "R08985", "R23044", "R04944", "R14589", "R30589", "R02783", "R34286", "R12611", "SH114448", "R02676", "R14511", "R23813", "R33010", "R19633", "R11238", "R04760", "R16543", "R27050", "R30015", "R04203", "R21561", "R24370", "R26255", "R39365", "R13730", "R24068", "R40342", "R24493", "R40580", "SH121006", "R03418", "R19731", "R33562", "OMS00144", "R10226", "R20353", "R33411", "R11790", "R02551", "R03042", "R40094", "R36639", "R40319", "SH127236", "R30678", "R12730", "R13191", "R37560", "R18251", "R02785", "OMS00169", "R04143", "R28680", "R32692", "R12528", "R00586", "R03924", "R21330", "R35324", "R09768", "R40152", "R23330", "R37205", "R14225", "R11353", "R30220", "R34385", "R19868")
data<-data[loc=toKeep]
#data<-data[data@pop!="nelsoniM"]

load("data/grp.rda")
load("data/dapc.rda")
#inversionTable
csv<-read.csv("data/omy-5-survey-83-pops.csv", header=TRUE)
sizes<-subset(csv, select=c(Pop_ID,N,InversionFreq,trial_map_order,Longitude,Latitude))

#grp3 is 469 presumably AR individuals

#grp 1 is presumbaly AA type, 297 individuals
d1<-data[dapc$assign==1]
d1<-d1[d1@pop!="nelsoni"]
d1<-d1[d1@pop!="KamchaA"]
d1<-d1[d1@pop!="FHVirg"]
d1<-d1[d1@pop!="FHWyom"]
d1<-d1[d1@pop!="ARHShas"]
d1<-d1[d1@pop!="ARHMocc"]

X1<-scaleGen(d1,NA.method="mean")
pca1<-dudi.pca(X1,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
s.label(pca1$li)


#grp2 is 851 inds and presumably RR
d2<-data[dapc$assign==2]
d2<-d2[d2@pop!="nelsoni"]
d2<-d2[d2@pop!="KamchaA"]
d2<-d2[d2@pop!="KamchaR"]
d2<-d2[d2@pop!="FHVirg"]
d2<-d2[d2@pop!="FHWyom"]
d2<-d2[d2@pop!="HCHKmlp"]
d2<-d2[d2@pop!="ARHEagl"]

X2<-scaleGen(d2,NA.method="mean")
pca2<-dudi.pca(X2,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
s.label(pca2$li)

#get grp3 to count AR
d3<-data[dapc$assign==3]
d3<-d3[d3@pop!="nelsoni"]

df3<-as.data.frame(table(d3@pop))  %>% rename_(AR="Freq") 

df1<-as.data.frame(table(d1@pop)) %>% rename_(AA="Freq")
df2<-as.data.frame(table(d2@pop)) %>% rename_(RR="Freq")

dataPop<-as.data.frame(table(data@pop)) %>% rename_(GenePopN="Freq")
genos<-full_join(df1,df3) %>% full_join(df2, by="Var1") %>% full_join(dataPop) %>% rename_(Pop_ID="Var1") %>% full_join(sizes)
genos[is.na(genos)]<-0

write.table(genos, file="popSizes.txt",sep="\t")
#grp2<-find.clusters(d2, max.n.clust=90)
#dapc2 <- dapc(d2, grp2$grp)
#dapc2<-dapc(d2, d2@pop)
#scatter(dapc2)

#Where do these odd individuals come from?
#d1[grp1$grp==3]@pop
#d2[grp2$grp==12]@pop

#> d1[dapc1$assign==5]@pop
#[1] MFEelM   MFEelM   MFEelM   MFEelM   MFEelM   EelVArsM EelVArsM EelVArsM BigCrM   BigCrM   BigCrM   BigCrM   BigCrM   BigCrM  
#[15] BigCrM   BigCrM   AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGrandeM AGrandeM
#[29] AGrandeM AGrandeM AGrandeM SisquocM UpMatilM VenSAntM EFSnGabM SnLuisRM
#Levels: MFEelM EelVArsM BigCrM AGLopezM AGrandeM SisquocM UpMatilM VenSAntM EFSnGabM SnLuisRM

#> d1[grp1$grp==5]@pop
#[1] MFEelM   MFEelM   MFEelM   MFEelM   MFEelM   EelVArsM EelVArsM EelVArsM BigCrM   BigCrM   BigCrM   BigCrM   BigCrM   BigCrM  
#[15] BigCrM   BigCrM   AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGLopezM AGrandeM AGrandeM
#[29] AGrandeM AGrandeM AGrandeM SisquocM UpMatilM VenSAntM EFSnGabM SnLuisRM
#Levels: MFEelM EelVArsM BigCrM AGLopezM AGrandeM SisquocM UpMatilM VenSAntM EFSnGabM SnLuisRM

#trees
#colfunc <- colorRampPalette(c("yellow", "red"))
colfunc <- colorRampPalette(c("cornflowerblue","darkblue"))

#pdf("./AARRtrees.pdf")


d1pops<-unique(d1@pop)
col1=colfunc(length(d1@pop))[rank(d1@pop)]
key1=colfunc(length(d1pops))[rank(d1pops)]

#col1=rainbow(length(d1@pop))[rank(d1@pop)]
#key1=rainbow(length(d1pops))[rank(d1pops)]

#D1 <- dist(tab(d1))
#tre1 <- nj(D1)
tre1<-aboot(d1, tree="nj",distance="edwards.dist",missing="mean", cutoff=50, sample=100) %>% ladderize()
plot(tre1, type="unrooted", show.tip.label=FALSE,main="AA Survey Data, n=297")
tiplabels(pch=21, col="black", bg=col1, cex=1)
legend("topleft", legend=d1pops, cex=0.4,fill=key1)

d2pops<-unique(d2@pop)
#key2=rainbow(length(d2pops))[rank(d2pops)]
#col2= rainbow(length(d2@pop))[rank(d2@pop)]
col2=colfunc(length(d2@pop))[rank(d2@pop)]
key2=colfunc(length(d2pops))[rank(d2pops)]

#D2 <- dist(tab(d2))
#tre2 <- nj(D2)
#tre2<-aboot(d2, tree="nj",distance="edwards.dist",missing="mean", cutoff=50, sample=100) %>% ladderize()
#plot(tre2, type="unrooted", show.tip.label=FALSE,main="RR Survey Data, n=850")
#tiplabels(pch=21, col="black", bg=col2, cex=1)
#legend("topleft", legend=d2pops, cex=0.4,fill=key2)

#dev.off()

#color only these pops
#FortRocM UKDeepM UKSprngM UKSpragM
#pdf("./interiorGroup.pdf")
#RRcols<-rep("grey",length(d2@pop))
#intCol<-as.character(d2@pop)
#RRcols[which(intCol=="FortRocM")]<-"green"
#RRcols[which(intCol=="UKDeepM")]<-"green"
#RRcols[which(intCol=="UKSpringM")]<-"green"
#RRcols[which(intCol=="UKSpragM")]<-"green"

#keysRR<-rep("grey",length(d2pops))
#keyint<-as.character(d2pops)
#keysRR[which(keyint=="FortRocM")]<-"green"
#keysRR[which(keyint=="UKDeepM")]<-"green"
#keysRR[which(keyint=="UKSprngM")]<-"green"
#keysRR[which(keyint=="UKSpragM")]<-"green"

#plot(tre2, type="unrooted", show.tip.label=FALSE,main="RR Survey Data, n=850")
#tiplabels(pch=21, col="black", bg=RRcols, cex=1)
#legend("topleft", legend=d2pops, cex=0.4,fill=keysRR)

#dev.off()

#PCA of types
matrix1<-pca1$li[,1:2]
matrix2<-pca2$li[,1:2]
(kip1<-100*pca1$eig/sum(pca1$eig))
(kip2<-100*pca2$eig/sum(pca2$eig))

#pdf("./AARRPCA.pdf")

ggplot(matrix1)+geom_point(aes(Axis1,Axis2), col=col1)+xlab(paste("Axis 1 ",round(kip1[1],2),"%",sep=""))+
  ylab(paste("Axis 2 ",round(kip1[2],2),"%",sep=""))+ggtitle("AA Survey SNP PCA")+theme_classic()

ggplot(matrix2)+geom_point(aes(Axis1,Axis2), col=col2)+xlab(paste("Axis 1 ",round(kip2[1],2),"%",sep=""))+
  ylab(paste("Axis 2 ",round(kip2[2],2),"%",sep=""))+ggtitle("RR Survey SNP PCA")+theme_classic()

#dev.off()

#distances

#we want to replace d1@pop with population numbers
#sizes[sizes$Pop_ID=="SnLuisR",]$trial_map_order
popNum <- function(x) {result<-sizes[sizes$Pop_ID==x,]$trial_map_order
return(result)}
popLat<-function(x) {result<-sizes[sizes$Pop_ID==x,]$Latitude
return(result)}

AAtree<-aboot(d1, strata=d1@pop, tree="nj",sample=1000,distance="edwards.dist",missing="mean",cutoff = 50) %>% ladderize()
RRtree<-aboot(d2, strata=d2@pop, tree="nj",sample=1000,distance="edwards.dist",missing="mean",cutoff = 50) %>% ladderize()

vector1<-as.character(d1@pop)
pop1<-lapply(AAtree$tip.label, popNum)

vector2<-as.character(d2@pop)
pop2<-lapply(RRtree$tip.label, popNum)

#colfunc <- colorRampPalette(c("yellow", "red"))
colfunc <- colorRampPalette(c("white","#005A9C"))
colfunc2<-colorRampPalette(c("white", "orange"))
popLat1 <-lapply(AAtree$tip.label, popLat)
popLat2 <-lapply(RRtree$tip.label, popLat)

#col1=colfunc(length(unlist(popLat1)))[rank(unlist(popLat1))]
#col2=colfunc2(length(unlist(popLat2)))[rank(unlist(popLat2))]
col1="royalblue4"
col2="#E31A1C"
AAtree$tip.label<-pop1
RRtree$tip.label<-pop2

AAtree$node.label<-signif(AAtree$node.label, digits=2)
RRtree$node.label<-signif(RRtree$node.label, digits=2)

pdf("8cd.pdf", width=8, height=8)
#cutoff = 50
plot(AAtree, type="unrooted",lab4ut="axial", show.node.label=TRUE, 
     show.tip.label=FALSE,cex=0.8,
     edge.width=2)
#tiplabels(pch=21, col="black", bg=col1, cex=1.5)
tiplabels(AAtree$tip.label, col="white", bg=col1, adj=c(0.5, 0.5), cex=.9, frame="circle")
#legend("topleft", legend=unlist(pop1), cex=0.45,fill=col1)
add.scale.bar(length = 0.01,lwd=2)
#dev.off()
#plot(AAtree, show.node.label=TRUE, main="1,000 bootstrap of AA Data")
#add.scale.bar(length = 0.01)
#pdf("RR.pdf")
plot(RRtree, type="unrooted",lab4ut="axial", show.node.label=TRUE,
     show.tip.label=FALSE, cex=0.7,
     edge.width=2)
#tiplabels(pch=21, col="black", bg=col2, cex=1.5)
tiplabels(RRtree$tip.label, col="white", bg=col2, adj=c(0.5, 0.5), 
          cex=.9, frame="circle")

#legend("topleft", legend=unlist(pop2), cex=0.45,fill=col2)

add.scale.bar(length = 0.01,lwd=2)

#plot(RRtree, show.node.label=TRUE, main="1,000 bootstrap of RR Data",cex=0.7)
#add.scale.bar(length = 0.01)
dev.off()




