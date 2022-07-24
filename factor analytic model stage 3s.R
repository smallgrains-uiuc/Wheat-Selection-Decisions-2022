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

######################################
##Factor analytic model for yield
######################################
allsub2<- droplevels.data.frame(allsub[which(allsub$trait %in% 'Grain.yield...bu.ac'),])

#Pick the trials to include
allsub2<- allsub2[which(allsub2$study %in% c("YT_Blk_22","YT_Blv_22","YT_Neo_22","YT_Stj_22",
                                              "YT_Stp_22","YT_Urb_22" ,"YT_War_22","YT_Lex_22","YT_Msn_22",
                                              "Pr_Neo_21","Pr_Stj_21","Pr_Urb_21")),]
#include only the stage 3 and 4 lines
allsub2<- allsub2[which(allsub2$germplasmName %in% unique(allsub2[which(allsub2$study == "YT_Blv_22"),'germplasmName'])),]

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
AvgYld<- rowMeans(faBlups_smry[,-c(1)])
faBlups_smry<- cbind(faBlups_smry, AvgYld)
write.csv(faBlups_smry, file='faBlups_smryJul22_stage3and4.csv')

