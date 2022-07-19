setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2021")
library(asreml)
library(reshape)
data<- read.csv('TrialsForAnalysis_July12.2021.csv', as.is=TRUE)

#suppress bad yield trial plots on the Maxwell field
pltex<- read.csv('maxwell_excludedplots.csv')
data[match(pltex[,1], data$observationUnitName),'Grain.yield...kg.ha.CO_321.0001218']<- NA
data[match(pltex[,1], data$observationUnitName),'Grain.test.weight...g.l.CO_321.0001210']<- NA
data[match(pltex[,1], data$observationUnitName),'Grain.moisture.content.....CO_321.0001198']<- NA
row.names(data)<- data$observationUnitName

#suppress suspicious values
data[match('Pr_Urb_21-plot902', data$observationUnitName),'Grain.test.weight...g.l.CO_321.0001210']<- NA
data[match('Adv_Urb_21-plot454', data$observationUnitName),'Grain.test.weight...g.l.CO_321.0001210']<- NA

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
  while(any(pctchg >2, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#convert yield and test weight to common units
data[,'Grain.yield...kg.ha.CO_321.0001218']<- convYld(data[,'Grain.yield...kg.ha.CO_321.0001218'])
data[,'Grain.test.weight...g.l.CO_321.0001210']<- convTwt(data[,'Grain.test.weight...g.l.CO_321.0001210'])
data[,"Plant.height...cm.CO_321.0001301"]<- data[,"Plant.height...cm.CO_321.0001301"] *0.393701
colnames(data)<- gsub('kg.ha.CO_321.0001218', "bu.ac", colnames(data))
colnames(data)<- gsub('g.l.CO_321.0001210', "lbs.bu", colnames(data))
colnames(data)<- gsub('..cm.CO_321.0001301', "inches", colnames(data))

#exclude the moisture data from analysis
data<- data[,-c(33)]

#order rows and columns
data<- data[order(data$rowNumber),]
data<- data[order(data$colNumber),]

#change variables to factors
data$germplasmName<- as.factor(as.character(data$germplasmName))
data$blockNumber<- as.factor(as.character(data$blockNumber))
data$entryType<- as.factor(as.character(data$entryType))
data$replicate<- as.factor(as.character(data$replicate))
data$rowNumber<- as.factor(as.character(data$rowNumber))
data$colNumber<- as.factor(as.character(data$colNumber))

#variance components row names for mv analysis
rnm<- c("trait:germplasmName!trait_Grain.yield...bu.ac:Grain.yield...bu.ac",
"trait:germplasmName!trait_Heading.time...Julian.date..JD..CO_321.0001233:Heading.time...Julian.date..JD..CO_321.0001233",
"trait:germplasmName!trait_Grain.test.weight...lbs.bu:Grain.test.weight...lbs.bu",
"trait:germplasmName!trait_Plant.height.inches:Plant.height.inches")

############################
#         Adv_Urb_21       #
############################
Adv_Urb_21<- droplevels.data.frame(data[which(data$studyName=="Adv_Urb_21"),])

#fit model UV for yield
amod1<- asreml(fixed=Grain.yield...bu.ac~1,
               random= ~germplasmName+blockNumber+rowNumber+colNumber,
               data=Adv_Urb_21, na.action = na.method(y='include', x='include'))
coefficients(amod1)$fixed
blups<- predict(amod1, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg)
mean(rel)
plot(varioGram(amod1))

#fit model MV (spatial analysis is not possible for this trial)
amod1<- asreml(fixed=cbind(Grain.yield...bu.ac, Heading.time...Julian.date..JD..CO_321.0001233,
                           Grain.test.weight...lbs.bu, Plant.height.inches)~1+trait,
               random= ~us(trait):germplasmName+ at(trait):blockNumber+at(trait):rowNumber+at(trait):colNumber,
               residual= ~id(units):diag(trait),
               data=Adv_Urb_21, na.action = na.method(y='include', x='include'))
coefficients(amod1)$fixed
blups<- predict(amod1, classify='trait:germplasmName', ignore=c('trait'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp[rnm,'component']
rel<- 1-(pev/Vg)
blups<- data.frame(blups, rel)


#fit model MV to get blues
amod1<- asreml(fixed=cbind(Grain.yield...bu.ac, Heading.time...Julian.date..JD..CO_321.0001233,
                           Grain.test.weight...lbs.bu, Plant.height.inches)~1+trait+trait:germplasmName,
               random= ~at(trait):blockNumber+at(trait):rowNumber+at(trait):colNumber,
               residual= ~id(units):diag(trait),
               data=Adv_Urb_21, na.action = na.method(y='include', x='include'))
p<- predict(amod1, classify = 'germplasmName:trait', pworkspace=64e7)$pvals
blues<- data.frame(study='Adv_Urb_21', p)
write.csv(blues, file='~/Documents/GitHub/Wheat-Selection-Decisions-2022/Adv_Urb_21blues.csv')


############################
#         Pr_Urb_21       #
############################
Pr_Urb_21<- droplevels.data.frame(data[which(data$studyName=="Pr_Urb_21"),])

#fit model UV for yield
amod1<- asreml(fixed=Grain.yield...bu.ac~1,
               random= ~germplasmName+blockNumber+rowNumber+colNumber,
               data=Pr_Urb_21, na.action = na.method(y='include', x='include'))
amod1<- update(amod1)
coefficients(amod1)$fixed
blups<- predict(amod1, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg)
mean(rel)
plot(varioGram(amod1))

#fit model MV (spatial analysis is not possible for this trial)
amod1<- asreml(fixed=cbind(Grain.yield...bu.ac, Heading.time...Julian.date..JD..CO_321.0001233,
                           Grain.test.weight...lbs.bu, Plant.height.inches)~1+trait,
               random= ~us(trait):germplasmName+ at(trait):blockNumber+at(trait):rowNumber+at(trait):colNumber,
               residual= ~id(units):diag(trait),
               data=Pr_Urb_21, na.action = na.method(y='include', x='include'))
amod1<- update(amod1)
coefficients(amod1)$fixed
blups<- predict(amod1, classify='trait:germplasmName', ignore=c('trait'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp[rnm,'component']
rel<- 1-(pev/Vg)
blups<- data.frame(blups, rel)

#fit model MV to get blues
amod1<- asreml(fixed=cbind(Grain.yield...bu.ac, Heading.time...Julian.date..JD..CO_321.0001233,
                           Grain.test.weight...lbs.bu, Plant.height.inches)~1+trait+us(trait):germplasmName,
               random= ~at(trait):blockNumber+at(trait):rowNumber+at(trait):colNumber,
               residual= ~id(units):diag(trait),
               data=Pr_Urb_21, na.action = na.method(y='include', x='include'))
amod1<- update(amod1)
p<- predict(amod1, classify = 'germplasmName:trait', pworkspace=64e7)$pvals
blues<- data.frame(study='Pr_Urb_21', p)
write.csv(blues, file='~/Documents/GitHub/Wheat-Selection-Decisions-2022/Pr_Urb_21blues.csv')


############################
#         Adv_Neo_21       #
############################
Adv_Neo_21<- droplevels.data.frame(data[which(data$studyName=="Adv_Neo_21"),])
Adv_Neo_21<- Adv_Neo_21[order(Adv_Neo_21$rowNumber),]
Adv_Neo_21<- Adv_Neo_21[order(Adv_Neo_21$colNumber),]


#fit model UV for yield
amod1<- asreml(fixed=Grain.yield...bu.ac~1,
               random= ~germplasmName+blockNumber+rowNumber+colNumber,
               residual= ~ar1(colNumber):ar1(rowNumber),
               data=Adv_Neo_21, na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
coefficients(amod1)$fixed
blups<- predict(amod1, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg) #0.7033814
mean(rel)
plot(varioGram(amod1))


#fit model MV with spatial analysis
amod1<- asreml(fixed=cbind(Grain.yield...bu.ac, Heading.time...Julian.date..JD..CO_321.0001233,
                           Grain.test.weight...lbs.bu, Plant.height.inches)~1+trait,
               random= ~us(trait):germplasmName+ at(trait):rowNumber+ at(trait):colNumber,
               residual = ~id(units):diag(trait),
               data=Adv_Neo_21,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify='trait:germplasmName', ignore=c('trait'))
blups<- p$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp[rnm,'component']
rel<- 1-(pev/Vg) #0.7283098
blups<- data.frame(blups, rel)

#fit model MV to get blues
amod1<- asreml(fixed=cbind(Grain.yield...bu.ac, Heading.time...Julian.date..JD..CO_321.0001233,
                           Grain.test.weight...lbs.bu, Plant.height.inches)~1+trait+trait:germplasmName,
               random= ~ at(trait):rowNumber+ at(trait):colNumber,
               residual = ~id(units):diag(trait),
               data=Adv_Neo_21,  na.action = na.method(y='include', x='include'))
p<- predict(amod1, classify = 'germplasmName:trait', pworkspace=64e7)$pvals
blues<- data.frame(study='Adv_Neo_21', p)
write.csv(blues, file='~/Documents/GitHub/Wheat-Selection-Decisions-2022/Adv_Neo_21blues.csv')


############################
#         Pr_Neo_21       #
############################
Pr_Neo_21<- droplevels.data.frame(data[which(data$studyName=="Pr_Neo_21"),])
Pr_Neo_21<- Pr_Neo_21[order(Pr_Neo_21$rowNumber),]
Pr_Neo_21<- Pr_Neo_21[order(Pr_Neo_21$colNumber),]

#fit model UV for yield
amod1<- asreml(fixed=Grain.yield...bu.ac~1,
               random= ~germplasmName+blockNumber+rowNumber+colNumber,
               residual= ~ar1(colNumber):ar1(rowNumber),
               data=Pr_Neo_21, na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
coefficients(amod1)$fixed
blups<- predict(amod1, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg) #0.6626489
mean(rel)
plot(varioGram(amod1))


#fit model MV with spatial analysis
amod1<- asreml(fixed=cbind(Grain.yield...bu.ac, Heading.time...Julian.date..JD..CO_321.0001233,
                           Grain.test.weight...lbs.bu, Plant.height.inches)~1+trait,
               random= ~us(trait):germplasmName+ at(trait):ar1(rowNumber)+ at(trait):ar1(colNumber),
               residual = ~id(units):diag(trait),
               data=Pr_Neo_21,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify='trait:germplasmName', ignore=c('trait'))
blups<- p$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp[rnm,'component']
rel<- 1-(pev/Vg) 
blups<- data.frame(blups, rel)
mean(blups[which(blups$trait=='Grain.yield...bu.ac'),'rel']) #0.6859805

#fit model MV to get blues
amod1<- asreml(fixed=cbind(Grain.yield...bu.ac, Heading.time...Julian.date..JD..CO_321.0001233,
                           Grain.test.weight...lbs.bu, Plant.height.inches)~1+trait+us(trait):germplasmName,
               random= ~ at(trait):ar1(rowNumber)+ at(trait):ar1(colNumber),
               residual = ~id(units):diag(trait),
               data=Pr_Neo_21,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify = 'germplasmName:trait', pworkspace=64e7)$pvals
blues<- data.frame(study='Pr_Neo_21', p)
write.csv(blues, file='~/Documents/GitHub/Wheat-Selection-Decisions-2022/Pr_Neo_21blues.csv')

############################
#         Adv_Stj_21       #
############################
rnm<- c("trait:germplasmName!trait_Grain.yield...bu.ac:Grain.yield...bu.ac",
        "trait:germplasmName!trait_Grain.test.weight...lbs.bu:Grain.test.weight...lbs.bu")

Adv_Stj_21<- droplevels.data.frame(data[which(data$studyName=="Adv_Stj_21"),])
Adv_Stj_21<- Adv_Stj_21[order(Adv_Stj_21$rowNumber),]
Adv_Stj_21<- Adv_Stj_21[order(Adv_Stj_21$colNumber),]

cols<- 1:length(unique(Adv_Stj_21$colNumber))
rows<- 1:length(unique(Adv_Stj_21$rowNumber))
colNumber<- rep(cols, length(unique(Adv_Stj_21$rowNumber)))
rowNumber<- rep(rows, length(unique(Adv_Stj_21$colNumber)))
crdf<- data.frame(colNumber, rowNumber)
Adv_Stj_21<- merge(crdf, Adv_Stj_21, all.x=TRUE)
Adv_Stj_21$colNumber<- as.factor(Adv_Stj_21$colNumber)
Adv_Stj_21$rowNumber<- as.factor(Adv_Stj_21$rowNumber)

#fit model UV for yield
amod1<- asreml(fixed=Grain.yield...bu.ac~1,
               random= ~germplasmName+blockNumber+rowNumber+colNumber,
               residual= ~ar1(colNumber):ar3(rowNumber),
               data=Adv_Stj_21, na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
coefficients(amod1)$fixed
blups<- predict(amod1, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg) 
mean(rel) #0.7326799
plot(varioGram(amod1))

#fit model MV with spatial analysis
amod1<- asreml(fixed=cbind(Grain.yield...bu.ac,Grain.test.weight...lbs.bu)~1+trait,
               random= ~us(trait):germplasmName+at(trait):ar1(blockNumber)+ at(trait):ar3(rowNumber)+ at(trait):ar1(colNumber),
               residual = ~id(units):us(trait),
               data=Adv_Stj_21,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify='trait:germplasmName', ignore=c('trait'))
blups<- p$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp[rnm,'component']
rel<- 1-(pev/Vg) 
blups<- data.frame(blups, rel)
mean(blups[which(blups$trait=='Grain.yield...bu.ac'),'rel']) #0.5546

#fit model UV models to get blues
amod1<- asreml(fixed=Grain.yield...bu.ac~1+germplasmName,
               random= ~blockNumber+rowNumber+colNumber,
               residual= ~ar1(colNumber):ar3(rowNumber),
               data=Adv_Stj_21, na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify = 'germplasmName', pworkspace=64e7)$pvals
bluesY<- data.frame(study='Adv_Stj_21', trait='Grain.yield...bu.ac', p)

amod1<- asreml(fixed=Grain.test.weight...lbs.bu~1+germplasmName,
               random= ~blockNumber+rowNumber+colNumber,
               residual= ~ar1(colNumber):ar3(rowNumber),
               data=Adv_Stj_21, na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify = 'germplasmName', pworkspace=64e7)$pvals
bluesT<- data.frame(study='Adv_Stj_21', trait='Grain.test.weight...lbs.bu', p)
blues<- rbind(bluesY, bluesT)
blues<- blues[,c(1,3,2,4:6)]
write.csv(blues, file='~/Documents/GitHub/Wheat-Selection-Decisions-2022/Adv_Stj_21blues.csv')


###########################
#         Pr_Stj_21       #
###########################
rnm<- c("trait:germplasmName!trait_Grain.yield...bu.ac:Grain.yield...bu.ac",
        "trait:germplasmName!trait_Grain.test.weight...lbs.bu:Grain.test.weight...lbs.bu")

Pr_Stj_21<- droplevels.data.frame(data[which(data$studyName=="Pr_Stj_21"),])
Pr_Stj_21<- Pr_Stj_21[order(Pr_Stj_21$rowNumber),]
Pr_Stj_21<- Pr_Stj_21[order(Pr_Stj_21$colNumber),]

cols<- 1:length(unique(Pr_Stj_21$colNumber))
rows<- 1:length(unique(Pr_Stj_21$rowNumber))
colNumber<- sort(rep(cols, length(unique(Pr_Stj_21$rowNumber))))
rowNumber<- rep(rows, length(unique(Pr_Stj_21$colNumber)))
crdf<- data.frame(colNumber, rowNumber)

Pr_Stj_21<- merge(crdf, Pr_Stj_21, all.x=TRUE)
Pr_Stj_21$colNumber<- as.factor(Pr_Stj_21$colNumber)
Pr_Stj_21$rowNumber<- as.factor(Pr_Stj_21$rowNumber)

#fit model UV for yield
amod1<- asreml(fixed=Grain.yield...bu.ac~1,
               random= ~germplasmName+blockNumber+rowNumber+colNumber,
               residual= ~ar1(colNumber):ar1(rowNumber),
               data=Pr_Stj_21, na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
coefficients(amod1)$fixed
blups<- predict(amod1, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg) 
mean(rel) #0.5952482
plot(varioGram(amod1))

#fit model MV with spatial analysis
amod1<- asreml(fixed=cbind(Grain.yield...bu.ac,Grain.test.weight...lbs.bu)~1+trait,
               random= ~us(trait):germplasmName+at(trait):blockNumber+ at(trait):ar1(rowNumber)+ at(trait):ar1(colNumber),
               residual = ~id(units):us(trait),
               data=Pr_Stj_21,  na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify='trait:germplasmName', ignore=c('trait'))
blups<- p$pvals
pev<- blups[,'std.error']^2
Vg<- summary(amod1)$varcomp[rnm,'component']
rel<- 1-(pev/Vg) 
blups<- data.frame(blups, rel)
mean(blups[which(blups$trait=='Grain.yield...bu.ac'),'rel']) #0.4416881

#fit model UV models to get blues
amod1<- asreml(fixed=Grain.yield...bu.ac~1+germplasmName,
               random= ~blockNumber+rowNumber+colNumber,
               residual= ~ar1(colNumber):ar1(rowNumber),
               data=Pr_Stj_21, na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify = 'germplasmName', pworkspace=64e7)$pvals
bluesY<- data.frame(study='Pr_Stj_21', trait='Grain.yield...bu.ac', p)

amod1<- asreml(fixed=Grain.test.weight...lbs.bu~1+germplasmName,
               random= ~blockNumber+rowNumber+colNumber,
               residual= ~ar1(colNumber):ar1(rowNumber),
               data=Pr_Stj_21, na.action = na.method(y='include', x='include'))
amod1<- mkConv(amod1)
p<- predict(amod1, classify = 'germplasmName', pworkspace=64e7)$pvals
bluesT<- data.frame(study='Pr_Stj_21', trait='Grain.test.weight...lbs.bu', p)
blues<- rbind(bluesY, bluesT)
blues<- blues[,c(1,3,2,4:6)]
write.csv(blues, file='~/Documents/GitHub/Wheat-Selection-Decisions-2022/Pr_Stj_21blues.csv')








