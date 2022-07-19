setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
library(asreml)
library(reshape)
library(gaston)
library(rrBLUP)
all<- read.csv('ILTrainingSetPhenotypesJuly.18.2022.csv', row.names=1)

#consider two mega-environments?
twoME<-  TRUE

#get data from 2020 to 2022
allsub<- droplevels.data.frame(all[which(all$year>19),])

#change year to a factor
allsub$year<- as.factor(as.character(allsub$year))

if(twoME){
  #Consider yield in north and south as different traits
  south_locs<- unique(allsub$loc)[-grep('Urb|Scb', unique(allsub$loc))]
  north_locs<- c('Urb')
  allsub[which(allsub$loc %in% south_locs & allsub$trait=='Grain.yield...bu.ac'),'trait']<- paste("South.", 
                                                                                                  allsub[which(allsub$loc %in% south_locs & allsub$trait=='Grain.yield...bu.ac'),'trait'],sep="")
  allsub[which(allsub$loc %in% north_locs & allsub$trait=='Grain.yield...bu.ac'),'trait']<- paste("North.", 
                                                                                                  allsub[which(allsub$loc %in% north_locs & allsub$trait=='Grain.yield...bu.ac'),'trait'],sep="")
  #get data on agronomic traits
  traits_sub<- c('North.Grain.yield...bu.ac', 'South.Grain.yield...bu.ac', 'Grain.test.weight...lbs.bu',
                 'Plant.height.inches','Heading.time...Julian.date..JD..')
}else{
  traits_sub<- c('Grain.yield...bu.ac', 'Grain.test.weight...lbs.bu',
                 'Plant.height.inches','Heading.time...Julian.date..JD..')
}

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

#is a reliability estimate needed
relNeeded<- FALSE

#Multi-trait model without the relationship matrix
mod<- asreml(fixed=value.x~trait, 
             random=~trait:site+us(trait):germplasmName+trait:site:germplasmName,
             weights=value.y, data=allsub3, 
             na.action = na.method(y='include', x='include'), 
             family=asr_gaussian(dispersion = 1), workspace=64e6)

if(relNeeded){
  blups<- predict(mod, classify='germplasmName:trait',
                  ignore=c('trait', '(Intercept)'))$pvals
  pev<- blups[,'std.error']^2
  Vg<- summary(mod)$varcomp[paste('trait:germplasmName!trait_', blups$trait, ":", blups$trait, sep=""),'component']
  rel<- 1-(pev/Vg) 
  blups_centered<- data.frame(blups, rel)
}


blups_Agronomic<- predict(mod, classify='germplasmName:trait')$pvals
write.csv(blups_Agronomic, file='blups_AgronomicJul18.csv')

#get the genetic correlations
gencor<- matrix(nrow=length(traits_sub), ncol=length(traits_sub))
for(i in 1:length(traits_sub)){
  for(j in 1:length(traits_sub)){
    num1<-summary(mod)$varcomp[paste('trait:germplasmName!trait_', traits_sub[i], ":", traits_sub[j], sep=""),'component']
    num2<-summary(mod)$varcomp[paste('trait:germplasmName!trait_', traits_sub[j], ":", traits_sub[i], sep=""),'component']
    num<- unique(na.omit(c(num1,num2)))
    denom_a<- summary(mod)$varcomp[paste('trait:germplasmName!trait_', traits_sub[i], ":", traits_sub[i], sep=""),'component']
    denom_b<- summary(mod)$varcomp[paste('trait:germplasmName!trait_', traits_sub[j], ":", traits_sub[j], sep=""),'component']
    cr<- num/c(sqrt(denom_a)*sqrt(denom_b))
    gencor[i,j]<- cr
  }
}
colnames(gencor)<- traits_sub
row.names(gencor)<- traits_sub

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

## Preliminary intersection of geno and pheno data
inter_lines <- intersect(unique(allsub3$germplasmName), geno@ped$id)
lines_diff<- setdiff(unique(allsub3$germplasmName), geno@ped$id)

#all lines tested in 2021
finalDesign<- read.csv('FinalDesignPlan_Aug24.csv')[,c('germplasmName')]
newLines<- finalDesign[finalDesign %in% geno@ped$id]

## subset the phenotypic data to include all those with genotypic data
allsub4<- droplevels.data.frame(allsub3[which(allsub3$germplasmName %in% inter_lines),])

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
#mod_gebv<- asreml(fixed=value.x~trait, 
#             random=~trait:site+us(trait):vm(germplasmName, K2)
#             +trait:site:germplasmName,
#             weights=value.y, data=allsub4, 
#             na.action = na.method(y='include', x='include'), 
#             family=asr_gaussian(dispersion = 1), workspace='96gb')
#gebvs_Agronomic<- predict(mod_gebv, classify='germplasmName:trait', pworkspace='8gb')$pvals
gebvs_Agronomic<- read.csv('gebvsIL22_Agronomic.csv', row.names=1)

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

#Multi-trait model without the relationship matrix
mod<- asreml(fixed=value.x~trait, 
             random=~us(trait):germplasmName+trait:site,
             weights=value.y, data=allsub4, 
             na.action = na.method(y='include', x='include'), 
             family=asr_gaussian(dispersion = 1), workspace=64e6)

if(relNeeded){
  blups<- predict(mod, classify='germplasmName:trait',
                  ignore=c('trait', '(Intercept)'))$pvals
  pev<- blups[,'std.error']^2
  Vg<- summary(mod)$varcomp[paste('trait:germplasmName!trait_', blups$trait, ":", blups$trait, sep=""),'component']
  rel<- 1-(pev/Vg) 
  blups_centered<- data.frame(blups, rel)
}
blups_Scab<- predict(mod, classify='germplasmName:trait')$pval
head(blups_Scab)
write.csv(blups_Scab, file='BLUPscabJul18.csv')
###########Include the genomic relationship matrix

# subset the phenotypic data to include all those with genotypic data
allsub4<- droplevels.data.frame(allsub4[which(allsub4$germplasmName %in% unique(c(inter_lines, newLines))),])

#Multi-trait genomic blup model for scab resistance
modGEBVscab<- asreml(fixed=value.x~trait, 
                     random=~us(trait):vm(germplasmName, K2)+trait:site,
                     weights=value.y, data=allsub4, 
                     na.action = na.method(y='include', x='include'), 
                     family=asr_gaussian(dispersion = 1), workspace="8gb")

GEBVscab<- predict(modGEBVscab, classify='germplasmName:trait', pworkspace="4gb")$pvals
write.csv(GEBVscab, file='GEBVscabJul18.csv')
