setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
library(asreml)
library(reshape)
library(gaston)
library(rrBLUP)
library(Matrix)

#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#read in means from stage 1 of analysis
all<- read.csv('ILTrainingSetPhenotypesJuly.21.2022.csv', row.names=1)

#Get the marker data to include in analysis
geno_file  <- "IL_2022_all_regions_samp_filt_fullnames_dedup_imp.vcf.gz"
geno <- read.vcf(geno_file)
geno@ped$id<- sub("^.*:", "", geno@ped$id)

#correct the line names
yrseries<- gsub("20", "", as.character(2000:2019))
for(i in 1:length(yrseries)){
  geno@ped$id<- gsub(paste("IL", yrseries[i], sep=""), yrseries[i], geno@ped$id)
}
geno@ped$id<- gsub('IL20', "2020", geno@ped$id)
geno@ped$id<- gsub('IL21', "IL2021", geno@ped$id)
geno@ped$id<- gsub('16LCSDH', "IL16LCSDH", geno@ped$id)
geno@ped$id<-gsub("PIO-25R74", "Pio25R74", geno@ped$id)
geno@ped$id<-gsub("KASKASKIA", "Kaskaskia", geno@ped$id)

###########use marker data to estimate an error variance in un replicated trials

###############Lexington KY
lex<- read.csv('YT_Lex_22.csv')
genoLex<- select.inds(geno, id %in% unique(lex$germplasmName))
inter_lines <- intersect(unique(lex$germplasmName), genoLex@ped$id)
lex<- droplevels.data.frame(lex[which(lex$germplasmName %in% inter_lines),])
## prep relationship matrix
genoLex<- select.snps(genoLex, maf > 0.01)
genoLex2<- as.matrix(genoLex)-1
K<- A.mat(genoLex2)
K2<- nearPD(K)
K2<- K2$mat

utrait<- unique(lex$trait)
for(i in 1:length(utrait)){
  sub<- droplevels.data.frame(lex[which(lex$trait==utrait[i]),])
  sub$germplasmName<- as.factor(sub$germplasmName)
  #fit models
  mod<- asreml(fixed= predicted.value~1, random=~ vm(germplasmName, K2), data=sub)
  mod<- mkConv(mod)
  sub$std.error<- sqrt(summary(mod)$varcomp['units!R','component'])
  sub$status<- 'Estimable'
  sub$wt<- 1/summary(mod)$varcomp['units!R','component']
  all<- rbind(all, sub)
}
#remove Truman
all<- all[-which(all$germplasmName =='TRUMAN'),]

##########Purdue 
Lay<- read.csv('YT_Lay_22.csv')
genoLay<- select.inds(geno, id %in% unique(Lay$germplasmName))
inter_lines <- intersect(unique(Lay$germplasmName), genoLay@ped$id)
Lay<- droplevels.data.frame(Lay[which(Lay$germplasmName %in% inter_lines),])
## prep relationship matrix
genoLay<- select.snps(genoLay, maf > 0.01)
genoLay2<- as.matrix(genoLay)-1
K<- A.mat(genoLay2)
K2<- nearPD(K)
K2<- K2$mat

utrait<- unique(Lay$trait)
for(i in 1:length(utrait)){
  sub<- droplevels.data.frame(Lay[which(Lay$trait==utrait[i]),])
  sub$germplasmName<- as.factor(sub$germplasmName)
  #fit models
  mod<- asreml(fixed= predicted.value~1, random=~ vm(germplasmName, K2), data=sub)
  mod<- mkConv(mod)
  sub$std.error<- sqrt(summary(mod)$varcomp['units!R','component'])
  sub$status<- 'Estimable'
  sub$wt<- 1/summary(mod)$varcomp['units!R','component']
  all<- rbind(all, sub)
}

##########Ohio 
OH<- read.csv('YT1_OH_22.csv')
genoOH<- select.inds(geno, id %in% unique(OH$germplasmName))
inter_lines <- intersect(unique(OH$germplasmName), genoOH@ped$id)
OH<- droplevels.data.frame(OH[which(OH$germplasmName %in% inter_lines),])
## prep relationship matrix
genoOH<- select.snps(genoOH, maf > 0.01)
genoOH2<- as.matrix(genoOH)-1
K<- A.mat(genoOH2)
K2<- nearPD(K)
K2<- K2$mat

utrait<- unique(OH$trait)
for(i in 1:length(utrait)){
  sub<- droplevels.data.frame(OH[which(OH$trait==utrait[i]),])
  sub$germplasmName<- as.factor(sub$germplasmName)
  #fit models
  mod<- asreml(fixed= predicted.value~1, random=~ vm(germplasmName, K2), data=sub)
  mod<- mkConv(mod)
  sub$std.error<- sqrt(summary(mod)$varcomp['units!R','component'])
  sub$status<- 'Estimable'
  sub$wt<- 1/summary(mod)$varcomp['units!R','component']
  all<- rbind(all, sub)
}


#get data from 2020 to 2022
allsub<- droplevels.data.frame(all[which(all$year>19),])

#change year to a factor
allsub$year<- as.factor(as.character(allsub$year))

#traits to use
traits_sub<- c('Grain.yield...bu.ac', 'Grain.test.weight...lbs.bu',
                 'Plant.height.inches','Heading.time...Julian.date..JD..')
allsub2<- droplevels.data.frame(allsub[which(allsub$trait %in% traits_sub),])

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

## Preliminary intersection of geno and pheno data
inter_lines <- intersect(unique(allsub3$germplasmName), geno@ped$id)
lines_diff<- setdiff(unique(allsub3$germplasmName), geno@ped$id)

#all lines tested in 2021
finalDesign<- read.csv('FinalDesignPlan_Aug24.csv')[,c('germplasmName')]
newLines<- finalDesign[finalDesign %in% geno@ped$id]

## subset the phenotypic data to include all those with genotypic data
allsub4<- droplevels.data.frame(allsub3[which(allsub3$germplasmName %in% inter_lines),])
#write.csv(allsub4, file='allsub4.csv')
#allsub4<- read.csv(file='allsub4.csv', row.names=1, stringsAsFactors = TRUE)

## subset the genotypic data to include all those with phenotypic data + new lines
#training and validation lines
geno<- select.inds(geno, id %in% unique(c(inter_lines, newLines)))

## subset the genotypic data to include only polymorphic snps
geno<- select.snps(geno, maf > 0.01)

#make relationship matrix
geno2<- as.matrix(geno)-1
K<- A.mat(geno2)

#make K positive semidefinite 
K2<- nearPD(K)
K2<- K2$mat

#remove extra objects
rm(K)
rm(geno2)

#Multi-trait model with the relationship matrix
#mod_test<- asreml(fixed=value.x~trait, 
#            random=~trait:study+us(trait):germplasmName+trait:study:germplasmName,
#            weights=value.y, data=allsub4, 
#            na.action = na.method(y='include', x='include'), 
#            family=asr_gaussian(dispersion = 1))
#mod_test<- mkConv(mod_test)
#blups_Agronomic_July28<- predict(mod_test, classify='germplasmName:trait')$pvals
#write.csv(blups_Agronomic_July28, file='blups_Agronomic_July28.csv')

#Multi-trait model with the relationship matrix
mod_gebv<- asreml(fixed=value.x~trait, 
             random=~trait:study+us(trait):vm(germplasmName, K2)
             +trait:study:germplasmName,
             weights=value.y, data=allsub4, 
             na.action = na.method(y='include', x='include'), 
             family=asr_gaussian(dispersion = 1), workspace='96gb')
mod_gebv<- mkConv(mod_gebv)
gebvs_Agronomic<- predict(mod_gebv, classify='germplasmName:trait', pworkspace='8gb')$pvals
write.csv(gebvs_Agronomic, 'gebvsIL22_AgronomicJuly30.csv')

#Fit each trait separately
yld<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Grain.yield...bu.ac'),]))
tw<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Grain.test.weight...lbs.bu'),]))
hd<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Heading.time...Julian.date..JD..'),]))
ht<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Plant.height.inches'),]))

#Yield
mod_yld<- asreml(fixed=value.x~1, 
                  random=~study+vm(germplasmName, K2),
                  weights=value.y, data=yld, 
                  na.action = na.method(y='include', x='include'), 
                  family=asr_gaussian(dispersion = 1), workspace='8gb')
mod_yld2<- update(mod_yld, random.= ~.+germplasmName:study)
yldGEBV<- predict(mod_yld2, classify='germplasmName')$pvals


#Test Weight
mod_tw<- asreml(fixed=value.x~1, 
                 random=~study+vm(germplasmName, K2),
                 weights=value.y, data=tw, 
                 na.action = na.method(y='include', x='include'), 
                 family=asr_gaussian(dispersion = 1), workspace='8gb')
mod_tw2<- update(mod_tw, random.= ~.+germplasmName:study)
twGEBV<- predict(mod_tw2, classify='germplasmName')$pvals

#Heading
mod_hd<- asreml(fixed=value.x~1, 
                random=~study+vm(germplasmName, K2),
                weights=value.y, data=hd, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='8gb')
mod_hd2<- update(mod_hd, random.= ~.+germplasmName:study)
hdGEBV<- predict(mod_hd2, classify='germplasmName')$pvals

#Height
mod_ht<- asreml(fixed=value.x~1, 
                random=~study+vm(germplasmName, K2),
                weights=value.y, data=ht, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='8gb')
mod_ht2<- update(mod_ht, random.= ~.+germplasmName:study)
htGEBV<- predict(mod_ht2, classify='germplasmName')$pvals



#########################################
### multi-trait BLUP for scab resistance
#########################################

#get scab resistance trait data
traits_sub<- c('FHB.DON.content...ppm.', 'FHB.grain.incidence.....')
allsub2<- droplevels.data.frame(allsub[which(allsub$trait %in% traits_sub),])

#subset the scab nursery
trials_sub<- unique(allsub2$study)[grep('Scb', unique(allsub2$study))]
allsub3<- droplevels.data.frame(allsub2[which(allsub2$study %in% trials_sub),])

#cast and melt data to include missing data
q<- melt(cast(allsub3, site+loc+year+study+germplasmName~trait, value='predicted.value'), 
         variable_name='trait', 
         id.vars=c('site','loc','year','study','germplasmName'), na.rm=FALSE)
q2<- melt(cast(allsub3, site+loc+year+study+germplasmName~trait, value='wt'), 
          variable_name='trait', 
          id.vars=c('site','loc','year','study','germplasmName'), na.rm=FALSE)
allsub4<- merge(q, q2, by=c('site','loc','year','study','germplasmName', 'trait'), 
                all.x=TRUE, all.y=TRUE)

#convert variables to factors
ix<- which(colnames(allsub4) %in%c('site','loc','year','study','germplasmName', 'trait'))
for(i in ix){
  allsub4[,i ]<- as.factor(as.character(allsub4[, i]))
}

#change NA weights to 0
allsub4[which(is.na(allsub4$value.y)),'value.y']<- 0

#is a reliability estimate needed
relNeeded<- FALSE

########### Model with the genomic relationship matrix

# subset the phenotypic data to include all those with genotypic data
allsub4<- droplevels.data.frame(allsub4[which(allsub4$germplasmName %in% unique(c(inter_lines, newLines))),])

#Multi-trait genomic blup model for scab resistance
modGEBVscab<- asreml(fixed=value.x~trait, 
                     random=~us(trait):vm(germplasmName, K2)+trait:site,
                     weights=value.y, data=allsub4, 
                     na.action = na.method(y='include', x='include'), 
                     family=asr_gaussian(dispersion = 1), workspace="8gb")

GEBVscab<- predict(modGEBVscab, classify='germplasmName:trait', pworkspace="4gb")$pvals
write.csv(GEBVscab, file='GEBVscabJul26.csv')
