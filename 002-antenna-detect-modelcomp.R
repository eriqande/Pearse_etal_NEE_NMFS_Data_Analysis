# This script fits a series of glms to beach antenna detection data to eluciate the effects of sex and 
# genotype on emmigration rate.

rm(list = ls())
library('RODBC')
library('mgcv')
library('formula.tools')

dataf <- readRDS(file="BCant.rds")

# genotype 1 = Resident, genotype 3 = Anadromous
# the coding scheme for the sex*genotype interaction is like this:
# sg1 -sex-dependnet dominance reversal: female AA= female RR; male AR=male RR. 
# sg2 -unspecified dominance; each sex-geno combo has a different smooth
# sg3 -resident dominance with sex differences
# sg4 -anadromous dominance with sex differences


# Fit the models. A priori, I expected a unimodal hump near 150mm for groups prone to migrating, 
# and a smooth delcine for groups not migrating but subject to self-thinning. Hence, GAMs-- emigration
# rate is a smooth function of length at release, and the question is whether that function 
# differs among different groups of fish. 


# the base model-- detections a smooth function of length at release
bc.gam.fl <- gam(detect ~ s(fl), family=binomial, data=dataf)

# sex effect-- males and females might be different
bc.gam.fl.sex <- gam(detect ~ s(fl, by=sex), family=binomial, data=dataf)

# genotype effect
bc.gam.fl.geno <- gam(detect ~ s(fl, by=geno) +geno, family=binomial, data=dataf)

# full interaction-- six different smooths
bc.gam.fl.genofull <- gam(detect ~ s(fl, by=sg2) +sg2, family=binomial, data=dataf)

# resident dominance- four smooths
bc.gam.fl.resdom <- gam(detect ~ s(fl, by=sg3) +sg3, family=binomial, data=dataf)

# anadromous dominance- four smooths
bc.gam.fl.anadom <- gam(detect ~ s(fl, by=sg4) +sg4, family=binomial, data=dataf)

# sex dependent dominance reversal- four differnt smooths
bc.gam.fl.genosdr <- gam(detect ~ s(fl, by=sg1) +sg1, family=binomial, data=dataf)



# pack up the results
model<-c(bc.gam.fl$formula, bc.gam.fl.sex$formula, bc.gam.fl.geno$formula, 
           bc.gam.fl.genofull$formula, bc.gam.fl.resdom$formula, bc.gam.fl.anadom$formula,
           bc.gam.fl.genosdr$formula)

aic<-c(bc.gam.fl$aic, bc.gam.fl.sex$aic, bc.gam.fl.geno$aic, 
         bc.gam.fl.genofull$aic, bc.gam.fl.resdom$aic, bc.gam.fl.anadom$aic,
         bc.gam.fl.genosdr$aic)

aic.min<-min(aic)
delta.aic <- aic-aic.min

model.comp.table<-data.frame(as.character(model), aic, delta.aic)
ordered.model.comp.table <- model.comp.table[order(delta.aic),]

write.csv(ordered.model.comp.table, 'modelcomptable.csv')


