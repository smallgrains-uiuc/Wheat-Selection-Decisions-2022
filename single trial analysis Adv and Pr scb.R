setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
library(asreml)
library(asremlPlus)
library(reshape)
data<- read.csv('All scb data 2021.csv', as.is=TRUE)

#function to check model convergence and update until converged (tolerate a 1.5% change in components)
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >2, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#change variables to factors
data$germplasmName<- as.factor(as.character(data$germplasmName))
data$blockNumber<- as.factor(as.character(data$blockNumber))
data$entryType<- as.factor(as.character(data$entryType))
data$replicate<- as.factor(as.character(data$replicate))
data$rowNumber<- as.factor(as.character(data$rowNumber))
data$colNumber<- as.factor(as.character(data$colNumber))

#variance components row names for mv analysis
rnm<- c("trait:germplasmName!trait_FHB.grain.incidence.....CO_321.0001155:FHB.grain.incidence.....CO_321.0001155",
        "trait:germplasmName!trait_FHB.DON.content...ppm.CO_321.0001154",
        "trait:germplasmName!trait_FHB.incidence.....CO_321.0001149:FHB.incidence.....CO_321.0001149",
        "trait:germplasmName!trait_FHB.severity.....CO_321.0001440:FHB.severity.....CO_321.0001440",
        "trait:germplasmName!trait_Heading.time...Julian.date..JD..CO_321.0001233:Heading.time...Julian.date..JD..CO_321.0001233")

############################
#         Adv_Scb_21       #
############################

Adv_Scb_21<- droplevels.data.frame(data[which(data$studyName=='Adv_Scb_21'),])

#fit model UV for SEV
cols<- 1:length(unique(Adv_Scb_21$colNumber))
rows<- 1:length(unique(Adv_Scb_21$rowNumber))
colNumber<- sort(rep(cols, length(unique(Adv_Scb_21$rowNumber))))
rowNumber<- rep(rows, length(unique(Adv_Scb_21$colNumber)))
crdf<- data.frame(colNumber, rowNumber)
Adv_Scb_21<- merge(crdf, Adv_Scb_21, all.x=TRUE)
Adv_Scb_21$colNumber<- as.factor(Adv_Scb_21$colNumber)
Adv_Scb_21$rowNumber<- as.factor(Adv_Scb_21$rowNumber)

amod1<- asreml(fixed=FHB.grain.incidence.....CO_321.0001155~1+blockNumber,
               random= ~germplasmName+rowNumber+colNumber,
               data=Adv_Scb_21, na.action = na.method(y='include', x='include'))
amod1<- update(amod1)
coefficients(amod1)$fixed
blups<- predict(amod1, classify='germplasmName', ignore=c('(Intercept)', 'blockNumber'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg)
mean(rel) #0.476
plot(varioGram(amod1))
summary(amod1)$varcomp

#Fit multivariate model for blups
amod1<- asreml(fixed=cbind(FHB.DON.content...ppm.CO_321.0001154, FHB.grain.incidence.....CO_321.0001155, FHB.incidence.....CO_321.0001149,FHB.severity.....CO_321.0001440, Heading.time...Julian.date..JD..CO_321.0001233)~1+trait+at(trait):blockNumber,
               random= ~us(trait):germplasmName+ at(trait):rowNumber+ at(trait):colNumber,
               residual = ~id(units):us(trait),
               data=Adv_Scb_21,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify='trait:germplasmName', ignore=c('trait',
                                                            'at(trait, FHB.grain.incidence.....CO_321.0001155):blockNumber',
                                                            'at(trait, Heading.time...Julian.date..JD..CO_321.0001233):blockNumber',
                                                            'at(trait, FHB.severity.....CO_321.0001440):blockNumber',
                                                            'at(trait, FHB.incidence.....CO_321.0001149):blockNumber',
                                                            'at(trait, FHB.DON.content...ppm.CO_321.0001154):blockNumber'), as.is=TRUE)
blups<- p$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp[rnm[c(2,1,3,4,5)],'component']
rel<- 1-(pev/rep(Vg, 176)) 
blups<- data.frame(blups, rel)
mean(blups[which(blups$trait=='FHB.DON.content...ppm.CO_321.0001154'),'rel']) #0.81
mean(blups[which(blups$trait=='FHB.grain.incidence.....CO_321.0001155'),'rel']) #0.68
mean(blups[which(blups$trait=='FHB.incidence.....CO_321.0001149'),'rel']) #0.76
mean(blups[which(blups$trait=='FHB.severity.....CO_321.0001440'),'rel']) #0.60
mean(blups[which(blups$trait=='Heading.time...Julian.date..JD..CO_321.0001233'),'rel']) #0.937


#Fit multivariate model for blues
amod1<- asreml(fixed=cbind(FHB.DON.content...ppm.CO_321.0001154, FHB.grain.incidence.....CO_321.0001155, FHB.incidence.....CO_321.0001149,
                           FHB.severity.....CO_321.0001440, Heading.time...Julian.date..JD..CO_321.0001233)~1+trait+us(trait):germplasmName,
               random= ~at(trait):rowNumber+ at(trait):colNumber+at(trait):blockNumber,
               residual = ~id(units):us(trait),
               data=Adv_Scb_21,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify='trait:germplasmName')$pvals
blues<- data.frame(study='Adv_Scb_21', p)
write.csv(blues, file='Adv_Scb_21blues.csv')

############################
#         Pr_Scb_21       #
############################
Pr_Scb_21<- droplevels.data.frame(data[which(data$studyName=='Pr_Scb_21'),])

#fit model UV for SEV
cols<- 1:length(unique(Pr_Scb_21$colNumber))
rows<- 1:length(unique(Pr_Scb_21$rowNumber))
colNumber<- sort(rep(cols, length(unique(Pr_Scb_21$rowNumber))))
rowNumber<- rep(rows, length(unique(Pr_Scb_21$colNumber)))
crdf<- data.frame(colNumber, rowNumber)
Pr_Scb_21<- merge(crdf, Pr_Scb_21, all.x=TRUE)
Pr_Scb_21$colNumber<- as.factor(Pr_Scb_21$colNumber)
Pr_Scb_21$rowNumber<- as.factor(Pr_Scb_21$rowNumber)

amod1<- asreml(fixed=FHB.grain.incidence.....CO_321.0001155~1+blockNumber,
               random= ~germplasmName+rowNumber+colNumber,
               data=Pr_Scb_21, na.action = na.method(y='include', x='include'))
amod1<- update(amod1)
coefficients(amod1)$fixed
blups<- predict(amod1, classify='germplasmName', ignore=c('(Intercept)', 'blockNumber'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg)
mean(rel) #0.70
plot(varioGram(amod1))
summary(amod1)$varcomp

#Fit multivariate model for blups
amod1<- asreml(fixed=cbind(FHB.DON.content...ppm.CO_321.0001154, FHB.grain.incidence.....CO_321.0001155, FHB.incidence.....CO_321.0001149,FHB.severity.....CO_321.0001440, Heading.time...Julian.date..JD..CO_321.0001233)~1+trait,
               random= ~us(trait):germplasmName+ at(trait):rowNumber+ at(trait):colNumber+at(trait):blockNumber,
               residual = ~id(units):us(trait),
               data=Pr_Scb_21,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify='trait:germplasmName')
blups<- p$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp[rnm[c(2,1,3,4,5)],'component']
rel<- 1-(pev/rep(Vg, nrow(blups)/5)) 
blups<- data.frame(blups, rel)
mean(blups[which(blups$trait=='FHB.DON.content...ppm.CO_321.0001154'),'rel']) #0.785
mean(blups[which(blups$trait=='FHB.grain.incidence.....CO_321.0001155'),'rel']) #0.71
mean(blups[which(blups$trait=='FHB.incidence.....CO_321.0001149'),'rel']) #0.65
mean(blups[which(blups$trait=='FHB.severity.....CO_321.0001440'),'rel']) #0.59
mean(blups[which(blups$trait=='Heading.time...Julian.date..JD..CO_321.0001233'),'rel']) #0.89

#Fit multivariate model for blues
amod1<- asreml(fixed=cbind(FHB.DON.content...ppm.CO_321.0001154, FHB.grain.incidence.....CO_321.0001155, 
                           FHB.incidence.....CO_321.0001149,FHB.severity.....CO_321.0001440, 
                           Heading.time...Julian.date..JD..CO_321.0001233)~1+trait+trait:germplasmName,
               random= ~at(trait):rowNumber+ at(trait):colNumber+at(trait):blockNumber,
               residual = ~id(units):us(trait),
               data=Pr_Scb_21,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify='trait:germplasmName', pworkspace=64e6)$pvals
blues<- data.frame(study='Pr_Scb_21', p)
write.csv(blues, file='Pr_Scb_21blues.csv')


############################
#         IL_Scb_22       #
############################
headers<- read.csv('IL_Scb_22.csv', skip = 3, header = F)[1,]
df<- read.csv('IL_Scb_22.csv', skip = 4, header = F)
colnames(df)<- headers

#make missing data NA
df[which(df$`FHB grain incidence - %|CO_321:0001155` <0),'FHB grain incidence - %|CO_321:0001155']<- NA

#convert variables to factors
df$germplasmName<- as.factor(as.character(df$germplasmName))
df$blockNumber<- as.factor(as.character(df$blockNumber))
df$entryType<- as.factor(as.character(df$entryType))
df$replicate<- as.factor(as.character(df$replicate))
df$rowNumber<- as.factor(as.character(df$rowNumber))
df$colNumber<- as.factor(as.numeric(df$colNumber)) 

colnames(df)[31]<- 'FHB.grain.incidence.....CO_321.0001155'
colnames(df)[34]<- 'Heading.time...Julian.date..JD..CO_321.0001233'

#Fit multivariate model for blues
amod1<- asreml(fixed=cbind(FHB.grain.incidence.....CO_321.0001155, Heading.time...Julian.date..JD..CO_321.0001233)~1+trait+trait:germplasmName+trait:replicate,
               random= ~at(trait):ar1(colNumber),
               residual = ~id(units):us(trait),
               data=df,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify='trait:germplasmName', pworkspace=64e6)$pvals
blues<- data.frame(study='IL_Scb_22', p)
write.csv(blues, file='IL_Scb_22blues.csv')

#Fit multivariate model for blups
amod1<- asreml(fixed=cbind(FHB.grain.incidence.....CO_321.0001155, Heading.time...Julian.date..JD..CO_321.0001233)~1+trait+trait:replicate,
               random= ~us(trait):germplasmName+at(trait):ar1(colNumber),
               residual = ~id(units):us(trait),
               data=df,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify='trait:germplasmName', ignore=c('trait', 'trait:replicate'), pworkspace=64e6)$pvals
Vg<- summary(amod1)$varcomp[paste('trait:germplasmName!trait_', p$trait, ":", p$trait, sep=""), 'component']
r2<- 1-(p$std.error^2)/Vg
p<- data.frame(p, r2)
