#/usr/bin/Rscript

#We are going to read in omy05 inversion data and non inversion SNPs and make an inertia plot.

library("adegenet")
library("ape")
library("pegas")
library("vcfR")
library("ggplot2")
library("ggrepel")


# Multiple plot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

nameList<-c("Shasta","North CA","South CA","Idaho","Klamath","Klamath","Idaho","Klamath","Central CA","Shasta","Central CA","Shasta","Kamloops","Kamchatka","Columbia","Klamath","Columbia","North CA","Central CA","North CA","Klamath","North CA","North CA","Idaho","Central CA","Kamchatka","Central CA","WA Coast","WA Coast","WA Coast","WA Coast","Columbia","Columbia","Columbia","WA Coast","WA Coast","WA Coast","Upper Columbia","Upper Columbia","Upper Columbia","Upper Columbia","WA Coast","WA Coast")

load(file="data/genindI.rda")
load(file="data/genind.rda")

vector<-c("AA","RR","AA","RR","RR","RR","AA","AA","AA","AA","RR","AA","RR","AA","RR","RR","RR","RR","RR","RR","RR","RR","RR","AA","RR","RR","RR","RR","RR","RR","RR")
fvector=factor(vector)

genindI@pop<-fvector
genind@pop<-fvector

#Make a tree
D <- dist(tab(genind))
tre <- nj(D)
#pdf("basicTreeGenomeWide.pdf")
par(xpd=TRUE)
plot(tre, type="unrooted", edge.w=2)
#dev.off()

#replace missing data
XI <- scaleGen(genindI, NA.method="mean")
#pca
pca1<-dudi.pca(XI,cent=FALSE,scale=FALSE,scannf=FALSE,nf=4)

s.label(pca1$li)
title("PCA of AA amd RR Types\naxes 1-2") 
add.scatter.eig(pca1$eig[1:20], posi="bottomright", 3,1,2)


#replace missing data
X <- scaleGen(genind, NA.method="mean")
#pca
pca2<-dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=4)
s.label(pca2$li)
title("PCA of Genome-Wide Data\naxes 1-2") 
add.scatter.eig(pca2$eig[1:20], posi="bottomright", 3,1,2)

#co1 <- coinertia(pca2, pca1)
#co1
#s.match(co1$mX, co1$mY, clab=.75)
#plot(co1)

#plot as blue and red
#pdf("GenomePaperPCA.pdf")
s.class(pca2$li, pop(genind), xax=1, yax=3, col=transp(c("blue","red"), .7), pch=17, axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)
par(new=TRUE)
s.class(pca1$li, pop(genind), xax=1, yax=3, col=c("blue","red"), axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)
title("PCA of Genome-Wide Data\nand Inversion Region Data\naxes 1-2") 
#dev.off()

#pdf("GenomePaperPCAV2.pdf")
s.label(pca1$li)
par(new=TRUE)
s.label(pca2$li)
#dev.off()

s.label(pca1$li, label=paste(nameList,genindI@pop,sep="-"), cpoint=2)
s.label(pca2$li, label=paste(nameList,genindI@pop,sep="-"), cpoint=2)

#s.label(pca2$li, label=paste(indNames(genindI),genindI@pop,sep="-"), cpoint=2)

#What about as discrete pops
#plot(pca2$li[, 1:2])
points1<-vector
points1[points1=="AA"]<-"17"
points1[points1=="RR"]<-"17"
points2<-vector
points2[points2=="AA"]<-"royalblue4"
points2[points2=="RR"]<-"#E31A1C"
points3<-vector
points3[points3=="AA"]<-"royalblue4"
points3[points3=="RR"]<-"#E31A1C"
#plot(pca2$li[, 1:2],pch=17, col=points2)
#plot(pca1$li[, 1:2],pch=19, col=points2)
#ggplot
genotype<-vector

matrix1<-pca1$li[, 1:4]
#compute %explained
(kip1<-100*pca1$eig/sum(pca1$eig))
(kip2<-100*pca2$eig/sum(pca2$eig))
#ggplot(matrix1)+geom_point(aes(Axis1,Axis2), color=points2)+geom_text_repel(aes(Axis1,Axis2, label=nameList))+theme_classic(base_size=16)
p1<-ggplot(matrix1)+geom_point(aes(Axis1,Axis2, size=3, color=points2, shape=factor(vector), fill=points2),alpha=0.75)+
  # geom_label_repel(
  #  aes(Axis1,Axis2, fill=genotype, label=nameList),
  # fontface = 'bold', color = 'white',
  #  box.padding = unit(0.35, "lines"),
  # point.padding = unit(0.5, "lines"),
  #  segment.color = 'grey50'
  #)+
  xlab(paste("PC 1 - ",round(kip1[1],2),"%",sep=""))+
  ylab(paste("PC 2 - ",round(kip1[2],2),"%",sep=""))+
  theme_classic(base_size=30)+
  ggtitle("Inversion Region")+
  theme(legend.position='none', plot.title = element_text(hjust = 0.5, size=rel(1.5)))+
  scale_color_manual(values=c("#E31A1C", "royalblue4"))


matrix2<-pca2$li[, 1:4]
#ggplot(matrix2)+geom_point(aes(Axis1,Axis2), color=points2)+geom_text_repel(aes(Axis1,Axis2, label=nameList))+theme_classic(base_size=16)

p2<-ggplot(matrix2)+geom_point(aes(Axis1,Axis2, color=points2, size=3, shape=factor(vector), fill=points2),alpha=0.75)+
  #scale_shape(solid = TRUE)+
  #geom_label_repel(
  #aes(Axis1,Axis2, fill=genotype, label=nameList),
  #fontface = 'bold', color = 'white',
  #box.padding = unit(0.35, "lines"),
  #point.padding = unit(0.5, "lines"),
  #segment.color = 'grey50'
  #)+
  xlab(paste("PC 1 - ",round(kip2[1],2),"%",sep=""))+
  ylab(paste("PC 2 - ",round(kip2[2],2),"%",sep=""))+
  theme_classic(base_size=30)+
  ggtitle("Genome-Wide")+
  theme(legend.position='none', plot.title = element_text(hjust = 0.5, size=rel(1.5)))+
  scale_color_manual(values=c("#E31A1C", "royalblue4"))



#Now to add numbers
nameList<-c("63","83","39","42","28","67","67","39","62","63","42","62","28","67","67","26","26","26","26","33","33","33","24","24","24","29","29","29","29","20","20")

p1<-ggplot(matrix1)+geom_point(aes(Axis1,Axis2, size=3, color=points2, shape=factor(vector), fill=points2),alpha=0.75)+
  geom_label_repel(
    aes(Axis1,Axis2, fill=genotype, label=nameList),
    fontface = 'bold', color = 'white',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50',
    fill=points2
  )+
  xlab(paste("PC 1 - ",round(kip1[1],2),"%",sep=""))+
  ylab(paste("PC 2 - ",round(kip1[2],2),"%",sep=""))+
  theme_classic(base_size=30)+
  ggtitle("Inversion Region")+
  theme(legend.position='none', plot.title = element_text(hjust = 0.5, size=rel(1.5)))+
  scale_color_manual(values=c("#E31A1C", "royalblue4"))

p2<-ggplot(matrix2)+geom_point(aes(Axis1,Axis2, color=points2, size=3, shape=factor(vector), fill=points2),alpha=0.75)+
  geom_label_repel(
    aes(Axis1,Axis2, label=nameList),
    fontface = 'bold', color = 'white',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50',
    fill=points2
  )+
  xlab(paste("PC 1 - ",round(kip2[1],2),"%",sep=""))+
  ylab(paste("PC 2 - ",round(kip2[2],2),"%",sep=""))+
  theme_classic(base_size=30)+
  ggtitle("Genome-Wide")+
  theme(legend.position='none', plot.title = element_text(hjust = 0.5, size=rel(1.5)))+
  scale_color_manual(values=c("#E31A1C", "royalblue4"))


##Now we need to get PC3 and 4 out there.

p3<-ggplot(matrix1)+geom_point(aes(Axis3,Axis4, size=3, color=points2, 
                                   shape=factor(vector), fill=points2),alpha=0.75)+
  geom_label_repel(
    aes(Axis3,Axis4, fill=genotype, label=nameList),
    fontface = 'bold', color = 'white',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50',
    fill=points2
  )+
  xlab(paste("PC 3 - ",round(kip1[3],2),"%",sep=""))+
  ylab(paste("PC 4 - ",round(kip1[4],2),"%",sep=""))+
  theme_classic(base_size=30)+
  ggtitle("Inversion Region")+
  theme(legend.position='none', plot.title = element_text(hjust = 0.5,
                                                          size=rel(1.5)))+
  scale_color_manual(values=c("#E31A1C", "royalblue4"))

p4<-ggplot(matrix2)+geom_point(aes(Axis3,Axis4, color=points2, size=3, 
                                   shape=factor(vector), fill=points2),alpha=0.75)+
  geom_label_repel(
    aes(Axis3,Axis4, label=nameList),
    fontface = 'bold', color = 'white',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50',
    fill=points2
  )+
  xlab(paste("PC 3 - ",round(kip2[3],2),"%",sep=""))+
  ylab(paste("PC 4 - ",round(kip2[4],2),"%",sep=""))+
  theme_classic(base_size=30)+
  ggtitle("Genome-Wide")+
  theme(legend.position='none', plot.title = element_text(hjust = 0.5, size=rel(1.5)))+
  scale_color_manual(values=c("#E31A1C", "royalblue4"))


###Let's do a better job plotting this 07/02/2018
matrix1$Locations<-nameList
matrix2$Locations<-nameList
matrix1$Names<-rownames(pca1$l1)
matrix2$Names<-rownames(pca2$l1)
library(tidyr)

pI<-ggplot(matrix1)+geom_point(aes(Axis1,Axis2, size=3, color=points2, shape=factor(vector), fill=points2),alpha=0.75)+
  geom_label_repel(
    aes(Axis1,Axis2, fill=genotype, label=Locations),
    fontface = 'bold', color = 'white',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50',
    fill=points2
  )+
  xlab(paste("PC 1 - ",round(kip1[1],2),"%",sep=""))+
  ylab(paste("PC 2 - ",round(kip1[2],2),"%",sep=""))+
  theme_classic(base_size=30)+
  ggtitle("Inversion Region")+
  theme(legend.position='none', plot.title = element_text(hjust = 0.5, size=rel(1.5)))+
  scale_color_manual(values=c("#E31A1C", "royalblue4"))

pII<-ggplot(matrix2)+geom_point(aes(Axis1,Axis2, color=points2, size=3, shape=factor(vector), fill=points2),alpha=0.75)+
  geom_label_repel(
    aes(Axis1,Axis2, label=Locations),
    fontface = 'bold', color = 'white',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50',
    fill=points2
  )+
  xlab(paste("PC 1 - ",round(kip2[1],2),"%",sep=""))+
  ylab(paste("PC 2 - ",round(kip2[2],2),"%",sep=""))+
  theme_classic(base_size=30)+
  ggtitle("Genome-Wide")+
  theme(legend.position='none', plot.title = element_text(hjust = 0.5, size=rel(1.5)))+
  scale_color_manual(values=c("#E31A1C", "royalblue4"))

pdf( file = "E8ab.pdf", width=11, height=7.0)
multiplot(pII,pI, cols=2)
dev.off()
