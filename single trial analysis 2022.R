setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
library(asreml)
library(reshape)
library(doBy)

#read data
data<- read.csv('YT_2022_phenotype_download.csv', as.is=TRUE)

############################
## Data conversions
############################
#convert yld to bu/acre
convYld<- function(y){
  x<- y/c(60 * 0.453592 * 2.47105)
  return(x)
}

#convert tw to lbs/bu
convTwt<- function(y){
  x<- y/1000 *2.2046 *35.2391
  return(x)
}  

#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#change -9 to NA
data[which(data[,'Grain.yield...kg.ha.CO_321.0001218']<0),'Grain.yield...kg.ha.CO_321.0001218']<- NA
data[which(data[,'Grain.test.weight...g.l.CO_321.0001210']<0),'Grain.test.weight...g.l.CO_321.0001210']<- NA

#convert yield and test weight to common units
data[,'Grain.yield...kg.ha.CO_321.0001218']<- convYld(data[,'Grain.yield...kg.ha.CO_321.0001218'])
data[,'Grain.test.weight...g.l.CO_321.0001210']<- convTwt(data[,'Grain.test.weight...g.l.CO_321.0001210'])
data[,"Plant.height...cm.CO_321.0001301"]<- data[,"Plant.height...cm.CO_321.0001301"] *0.393701
colnames(data)<- gsub('kg.ha.CO_321.0001218', "bu.ac", colnames(data))
colnames(data)<- gsub('g.l.CO_321.0001210', "lbs.bu", colnames(data))
colnames(data)<- gsub('..cm.CO_321.0001301', "inches", colnames(data))

#order by row, column, and trial
data$rowNumber<- as.numeric(data$rowNumber)
data$colNumber<- as.numeric(data$colNumber)
data<- data[order(data$rowNumber),]
data<- data[order(data$colNumber),]
data<- data[order(data$studyName),]

#convert to factors
data$colNumber<- as.factor(data$colNumber)
data$rowNumber<- as.factor(data$rowNumber)
data$germplasmName<- as.factor(data$germplasmName)
data$replicate<- as.factor(data$replicate)
data$studyName<- as.factor(data$studyName)

ustud<- unique(as.character(data$studyName))
for(i in 1:length(ustud)){
  
  #select a single study
  sub<- droplevels.data.frame(data[which(data$studyName == ustud[i]),])
  
  #traits
  trt<- colnames(sub)[c(41:44)]
  
  for(j in 1:length(trt)){
    if(!is.na(var(sub[,trt[j]], na.rm=TRUE))){
      
      #fit random model
      amod1<- asreml(fixed=as.formula(paste(trt[j], "~1", sep="")),
                     random= ~germplasmName,
                     residual= ~ar1(colNumber):ar1(rowNumber),
                     data=sub, na.action = na.method(y='include', x='include'))
      amod1<- mkConv(amod1)
      blups<- predict(amod1, classify='germplasmName', ignore=c('(Intercept)'))$pvals
      pev<- blups[,'std.error']^2
      Vg<- summary(amod1)$varcomp['germplasmName','component']
      rel<- 1-(pev/Vg) 
      blups<- data.frame(blups, rel, trait=trt[j], study=ustud[i])
      
      #fit fixed model
      amod2<- asreml(fixed=as.formula(paste(trt[j], "~1+germplasmName", sep="")),
                     residual= ~ar1(colNumber):ar1(rowNumber),
                     data=sub, na.action = na.method(y='include', x='include'))
      amod2<- mkConv(amod2)
      blues<- predict(amod2, classify='germplasmName', pworkspace=64e7)
      blues<- data.frame(blues, trait=trt[j], study=ustud[i])
      
      if(i==1 & j==1){
        blupsAll<- blups
        bluesAll<- blues
      }else{
        blupsAll<- rbind(blupsAll, blups)
        bluesAll<- rbind(bluesAll, blues)
      }
    }
  }
}

summary_by(blupsAll, rel~study+trait)
summary_by(bluesAll, pvals.std.error~study+trait, na.rm=TRUE)
write.csv(bluesAll, file = '~/Documents/GitHub/Wheat-Selection-Decisions-2022/YT blues 2022.csv')
