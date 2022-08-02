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

#get data from 2022
allsub<- droplevels.data.frame(all[which(all$year>21),])

#change year to a factor
allsub$year<- as.factor(as.character(allsub$year))

######################################
##Factor analytic model for yield
######################################
allsub2<- droplevels.data.frame(allsub[which(allsub$trait %in% 'Grain.yield...bu.ac'),])

#Pick the trials to include
allsub2<- allsub2[which(allsub2$study %in% c("YT_Blk_22","YT_Neo_22","YT_Stj_22",
                                             "YT_Stp_22","YT_Urb_22" ,"YT_War_22")),]
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

#change Stj 22 to Addieville
allsub3$study<- as.character(allsub3$study)
allsub3$site<- as.character(allsub3$site)
allsub3$loc<- as.character(allsub3$loc)
allsub3[which(allsub3$study == 'YT_Stj_22'),'site']<- 'Add_22'
allsub3[which(allsub3$study == 'YT_Stj_22'),'loc']<- 'Add'
allsub3[which(allsub3$study == 'YT_Stj_22'),'study']<- 'YT_Add_22'
allsub3$study<- as.factor(allsub3$study)
allsub3$site<- as.factor(allsub3$site)
allsub3$loc<- as.factor(allsub3$loc)

#factor analytic model without the relationship matrix
mod<- asreml(fixed=value.x~1, 
             random=~site+fa(site, 2):germplasmName,
             weights=value.y, data=allsub3, 
             na.action = na.method(y='include', x='include'), 
             family=asr_gaussian(dispersion = 1), workspace=64e6)
mod<- mkConv(mod)

#site correlation matrox
met.corr <-function(object,site,faN=2){ # ,faRS=1
  
  n<-nlevels(site)
  
  varcomp<-summary(object)$varcomp['component']
  vcn<-row.names(varcomp)
  aimn<-vcn[grep('fa\\(.*,.*\\)',vcn)]
  varcomp1<-varcomp[aimn,]
  vect1<-varcomp1[1:n]
  w.var<-diag(vect1)
  vect2<-varcomp1[(n+1):((1+faN)*n)]
  L.var<-matrix(vect2,nrow=n)
  wL.var<-L.var%*%t(L.var)+w.var
  df<-wL.var
  for(i in 1:(n-1)){
    for(j in 2:n){
      if(i<j){df[i,j]<-df[j,i]/(sqrt(df[i,i]*df[j,j]))
      j<-j+1}
    }
    i<-i+1
  }
  rownames(df)<-levels(site)
  colnames(df)<-levels(site)
  
  df.2<-df
  
  for(i in 1:(n-1)){
    for(j in 2:n){
      if(i<j){df[j,i]<-df[i,j]
      j<-j+1}
    }
    i<-i+1
  }
  diag(df)<-1
  return(df)
}
df<- met.corr(mod, unique(allsub3$site, 2))
heatmap(df)
#df<- round(df, 2)
#stord<- c('Blv_22','War_22', 'Msn_22', 'Blk_22', 'Lex_22', 'Add_22','Neo_22','Stp_22','Urb_22','Stj_21','Neo_21','Urb_21','Urb_20','Stj_20','Car_20','Neo_20')
#write.csv(df[stord, stord], 'siteCorrelations.csv')

faBlups<- predict(mod, classify='germplasmName:site', pworkspace="4gb")$pvals
faBlups$predicted.value<- round(faBlups$predicted.value, 2)
faBlups_smry<- cast(faBlups, germplasmName~site, value='predicted.value')
faBlups_smry<- faBlups_smry[,c('germplasmName', colnames(faBlups_smry)[order(colMeans(faBlups_smry[,-1]))+1])]
AvgYld<- round(rowMeans(faBlups_smry[,-c(1)]), 2)
faBlups_smry<- cbind(faBlups_smry, AvgYld)

#subset the stage 2 lines
MasterEntryListILYT <- read.csv("~/Documents/Wheat/2022/MasterEntryListILYT.csv")
s2lines<- MasterEntryListILYT[which(MasterEntryListILYT$Stage ==2 | MasterEntryListILYT$Stage ==0),c('germplasmName', 'Source.observationUnitName')]

faBlups_smry_sub<- faBlups_smry[which(faBlups_smry$germplasmName %in% s2lines$germplasmName),]
faBlups_smry_sub<- merge(s2lines, faBlups_smry_sub, by='germplasmName')

#add small increase information
S600 <- read.csv("~/Documents/Wheat/2022/mirusFile_S600_seedlots.csv", row.names=1)
colnames(S600)[4]<- 'germplasmName'
faBlups_smry_sub<- merge(faBlups_smry_sub, S600, by='germplasmName', all.x=TRUE)
keepvdiscard<- read.csv('~/Documents/Wheat/2022/S600_KeepvDiscard_July5_table.csv', row.names=1)[,c('seedlot_name', 'action','Found')]

faBlups_smry_sub<- cbind(faBlups_smry_sub, keepvdiscard[match(faBlups_smry_sub$plot_name.,keepvdiscard$seedlot_name),])
#write.csv(faBlups_smry_sub, file='faBlups_smryJul25_stage2.csv')

#add gebvs to the summary
faBlups_smry_sub<- read.csv('faBlups_smryJul25_stage2.csv')
gebv<- read.csv('gebvsIL22_Agronomic.csv', row.names=1)
gebv<- gebv[which(gebv$germplasmName %in% faBlups_smry_sub$germplasmme),]
gebvwide<- cast(gebv, germplasmName~trait, value='predicted.value')
colnames(faBlups_smry_sub)[2]<- 'germplasmName'
faBlups_smry_sub2<- merge(faBlups_smry_sub, gebvwide, by='germplasmName', all.x=TRUE)

gebvscab<- read.csv('GEBVscabJul18.csv', row.names=1)
gebvScb_wide<- cast(gebvscab, germplasmName~trait, value='predicted.value')
faBlups_smry_sub2<- merge(faBlups_smry_sub2, gebvScb_wide, by='germplasmName', all.x=TRUE)


##################################
## add net merit to the summary  #
##################################

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

colnames(faBlups_smry_sub2)

faBlups_smry_sub3<- data.frame(faBlups_smry_sub2, merit=netMerit(faBlups_smry_sub2$`Heading.time...Julian.date..JD..`,
                                      faBlups_smry_sub2$`Grain.yield...bu.ac`,
                                      faBlups_smry_sub2$`FHB.DON.content...ppm.`,
                                      faBlups_smry_sub2$`FHB.grain.incidence.....`, 
                                      faBlups_smry_sub2$`Grain.test.weight...lbs.bu`,
                                         wheat_price0, soybean_price))
for(i in 21:27){
  faBlups_smry_sub3[,i]<- round(faBlups_smry_sub3[,i], 2)
}
#write.csv(faBlups_smry_sub3, file='stage2_selectionfile2022.csv')

#add stage 2 notes
notes<- read.csv('stage 2 notes.csv')
notes$notes<- trimws(notes$notes)
notes$notes2<- paste(notes$studyName, "[", notes$notes, "]", sep="")
ugid<- unique(notes$germplasmName)
notes$notes2<- gsub("_", "", notes$notes2)
notes3<- c()
for(i in 1:length(ugid)){
  notes3<- append(notes3, paste(notes[which(ugid[i] == notes$germplasmName),'notes2'], collapse=", "))
}
notesdf<- data.frame(ugid, notes3)
notesdf<- notesdf[match(faBlups_smry_sub3$germplasmName, notesdf$ugid),]
write.csv(notesdf, file='stage 2 notes formatted.csv')
