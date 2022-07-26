setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
library(asreml)
library(reshape)
library(gaston)
library(rrBLUP)
all<- read.csv('ILTrainingSetPhenotypesJuly.21.2022.csv', row.names=1)

#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#get data from 2020 to 2022
allsub<- droplevels.data.frame(all[which(all$year>19),])

#change year to a factor
allsub$year<- as.factor(as.character(allsub$year))

#subset traits
traits_sub<- c('Grain.yield...bu.ac', 'Grain.test.weight...lbs.bu',
                 'Plant.height.inches','Heading.time...Julian.date..JD..')
allsub2<- droplevels.data.frame(allsub[which(allsub$trait %in% traits_sub),])

#Exclude trials without data on stage 4 lines and exclude YT_War_22 and YT_Msn_22
allsub2<- allsub2[-which(allsub2$study %in% c('Adv_Urb_20','Adv_StJ_20','Adv_Car_20',
                                              'Adv_Neo_20', 'Pr_Neo_21', 'Pr_Urb_21', 'Pr_Stj_21', 'AdvHY_Urb_20',
                                              'YT_War_22','YT_Msn_22')),]

#include only the lines in all locations
s4lines<- read.csv('s4lines22Notes.csv', row.names=1)[,1]
s4lines<- s4lines[-1]
allsub2<- allsub2[which(allsub2$germplasmName %in% s4lines),]

#subset the yield trials
if(length(grep('Scb', unique(allsub2$study)))>0){
  trials_sub<- unique(allsub2$study)[-grep('Scb', unique(allsub2$study))]
  allsub2<- droplevels.data.frame(allsub2[which(allsub2$study %in% trials_sub),])
}

#cast and melt data to include missing data
q<- melt(cast(allsub2, site+loc+year+study+germplasmName~trait, value='predicted.value'), 
         variable_name='trait', 
         id.vars=c('site','loc','year','study','germplasmName'), na.rm=FALSE)
q2<- melt(cast(allsub2, site+loc+year+study+germplasmName~trait, value='wt'), 
          variable_name='trait', 
          id.vars=c('site','loc','year','study','germplasmName'), na.rm=FALSE)
allsub3<- merge(q, q2, by=c('site','loc','year','study','germplasmName', 'trait'), 
                all.x=TRUE, all.y=TRUE)

#convert variables to factors
ix<- which(colnames(allsub3) %in%c('site','loc','year','study','germplasmName', 'trait'))
for(i in ix){
  allsub3[,i ]<- as.factor(as.character(allsub3[, i]))
}

#change NA weights to 0
if(length(which(is.na(allsub3$value.y)))>0){
  allsub3[which(is.na(allsub3$value.y)),'value.y']<- 0
}

#is a reliability estimate needed
relNeeded<- TRUE

#Multi-trait model without the relationship matrix
mod<- asreml(fixed=value.x~trait, 
             random=~trait:site+us(trait):germplasmName+trait:site:germplasmName,
             weights=value.y, data=allsub3, 
             na.action = na.method(y='include', x='include'), 
             family=asr_gaussian(dispersion = 1), workspace=64e6)
mod<- mkConv(mod)

if(relNeeded){
  blups<- predict(mod, classify='germplasmName:trait',
                  ignore=c('trait', '(Intercept)'))$pvals
  pev<- blups[,'std.error']^2
  Vg<- summary(mod)$varcomp[paste('trait:germplasmName!trait_', blups$trait, ":", blups$trait, sep=""),'component']
  rel<- 1-(pev/Vg) 
  blups_centered<- data.frame(blups, rel)
}

blups_Agronomic<- predict(mod, classify='germplasmName:trait')$pvals
write.csv(blups_Agronomic, file='blupsStage4_AgronomicJul21.csv')

################
#add net merit
################

blups_AgronomicW<- cast(blups_Agronomic, germplasmName~trait, value='predicted.value')
blups_Scab<- read.csv('BLUPscabJul18.csv')
blups_ScabW<- cast(blups_Scab, germplasmName~trait, value='predicted.value')

#combine scab with agronomic
blups<- merge(blups_AgronomicW, blups_ScabW, by='germplasmName', all.x=TRUE)

##For net merit calcualtion
wheat_price0<- mean(c(9.9128, 7.0402, 5.4621, 4.9414, 4.9757, 4.4014, 4.3945))
soybean_price<- mean(c(16.1773, 13.6890, 9.5344, 8.9298, 9.3456, 9.7820, 9.8753))

#wheat price fcn
wheatPrice<- function(fdk, don, twt, wheat_price0){
  if(don==0){
    donDiscount<- 0
  }else{
    donDiscount<- don*-0.2
  }
  if(fdk==0){
    fdkDiscount<- 0
  }else{
    fdkDiscount<- fdk*(-0.04/5)
  }
  twtDiscount<- c(58-twt)*-.2
  twtDiscount[which(twtDiscount>0)]<- 0
  wheat_price<- wheat_price0+donDiscount+fdkDiscount+twtDiscount
  return(wheat_price)
}

#net merit function
netMerit<- function(headings, yields, dons, fdks, twt, wheat_price0, soybean_price){
  wheat_price1<- wheatPrice(fdks, dons, twt, wheat_price0)
  soy_yld_gain<- 0.5* (137.7407-headings)
  soy_profit_gain<- soy_yld_gain*soybean_price
  wheat_profit<- yields*wheat_price1
  total_profit<- wheat_profit + soy_profit_gain
  return(total_profit)
}
blups$zero<- 0
blups<- data.frame(blups, net_merit=netMerit(blups$`Heading.time...Julian.date..JD..`,
                                              blups$`Grain.yield...bu.ac`,
                                              blups$`FHB.DON.content...ppm.`,
                                              blups$`FHB.grain.incidence.....`, 
                                              blups$`Grain.test.weight...lbs.bu`,
                                              wheat_price0, soybean_price))
blups<- data.frame(blups, net_merit_noScab=netMerit(blups$`Heading.time...Julian.date..JD..`,
                                             blups$`Grain.yield...bu.ac`,
                                             blups$zero,
                                             blups$zero, 
                                             blups$`Grain.test.weight...lbs.bu`,
                                             wheat_price0, soybean_price))
blups<- blups[order(blups$net_merit, decreasing=TRUE),]
blups<- blups[,-c(8)]
colnames(blups)[c(2:7)]<- c('Test weight', 'Yield', 'Heading date', 'Height', 'Vomitoxin', "Kernel damage")
for(i in 2:ncol(blups)){
  blups[,i]<- round(blups[,i], 1)
}
write.csv(blups, file='stage 4 multi-trait smry.csv')
