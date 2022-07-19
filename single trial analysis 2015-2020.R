setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
library(asreml)
library(asremlPlus)
library(reshape)
data<- read.csv('trials 2015-2020.csv', as.is=TRUE)
row.names(data)<- data$observationUnitName

#Exclude data-points declared as outliers
data['A5S_Urb_20-plot216',"Heading.time...Julian.date..JD..CO_321.0001233"]<- NA
data['A6S_Urb_19-plot105',"Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot1048',"Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot1049', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot1096', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot2025', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot2114', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_Neo_20-plot518', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_Neo_20-plot559', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_Scb_20-plot3019', "FHB.grain.incidence.....CO_321.0001155"] <- NA
data['Adv_Stj_18-plot2009', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_StJ_20-plot317', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_StJ_20-plot322', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_StJ_20-plot520', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_StJ_20-plot539', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['AdvHY_Urb_20-plot30163', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['ORG_Urb_20-plot112', "Septoria.nodorum.leaf.blotch.severity...0.9.percentage.scale.CO_321.0501145"] <- NA
data['Pr_Sbmv_20-plot216', "Soil.borne.mosaic.plant.response...0.9.Response.Scale.CO_321.0501140"] <- NA
data['Pr_Urb_20-plot763', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Pr2_Neo_18-plot201', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Pr3_Urb_18-plot214', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Pr5_Urb_19-plot162', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Pr6_Urb_19-plot189', "Grain.test.weight...g.l.CO_321.0001210"] <- NA

data['A6S_Urb_19-plot103', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Brt_15-plot1013', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_17-plot1001', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot2004', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot2024', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot2028', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_20-plot246', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_Scb_20-plot2144', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Stj_16-plot1008', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_Stj_16-plot3035', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_Stj_17-plot3035', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_Urb_15-plot1007', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Pr2_Urb_16-plot233', "Grain.yield...kg.ha.CO_321.0001218"] <- NA

#sort plot no within trial
data<- data[order(data$plotNumber),]
data<- data[order(data$studyName),]

############################
##Data corrections
############################

pltsGy<- c('AdvHY_Urb_20-plot20162', 'AdvHY_Urb_20-plot30106', 'Adv_Neo_18-plot1095',
           'Adv_Neo_18-plot1097', 'Adv_Neo_18-plot1120', 'Adv_Neo_20-plot267',
           'Adv_Neo_20-plot340', 'Adv_StJ_20-plot118', 'Adv_StJ_20-plot139', 'Adv_StJ_20-plot511')
data[which(data$observationUnitName %in% pltsGy), "Grain.yield...kg.ha.CO_321.0001218"]<- NA

pltsTw<- c('Adv_Neo_18-plot1024', 'Adv_Neo_18-plot1073', 'Adv_Neo_18-plot1120', 'Adv_Neo_18-plot2001',
           'Adv_Neo_18-plot2021', 'Adv_Neo_18-plot2045', 'Adv_Neo_20-plot565', 'Pr_Car_20-plot135')
data[which(data$observationUnitName %in% pltsTw), "Grain.test.weight...g.l.CO_321.0001210"]<- NA

#test weight data, convert zero to missing
data[which(data[,"Grain.test.weight...g.l.CO_321.0001210"]==0),"Grain.test.weight...g.l.CO_321.0001210"]<- NA

############################
##Functions to be used later
############################

#mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#function to check model convergence and update until converged (tolerate a 1.5% change in components)
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >2, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#create vector of traits recorded in the trial
selectTraitcols<- function(trl, trtnms, thresh=0.2){
  ppres<- c()
  for(j in 1:length(trtnms)){
    vec<-trl[,trtnms[j]]
    ppres<- append(ppres, length(na.omit(vec))/length(vec))
  }
  return(trtnms[which(ppres>thresh)])
}

#function to add means, MSE, LSD, and CV
addRows<- function(df, varNm, label, vec){
  lenvec<- length(vec)
  df[nrow(df)+1,]<-df[nrow(df),]
  df[nrow(df), c(1:varNm)]<- rep("", varNm)
  df[nrow(df), varNm]<- label
  df[nrow(df), -c(1:varNm)]<- round(vec, 5)
  return(df)
}

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

############################
## Data curation
############################

#treat prelims from before 2020 as one trial with separate blocks in the un-replicated sites
studgrpExp<- c("Pr[0-9]_Car_19","Pr[0-9]_Stj_19", "Pr[0-9]_Urb_19", "Pr[0-9]_Scb_19",
 "Pr[0-9]_Urb_18","Pr[0-9]_Neo_18","Pr[0-9]_Stj_18", "Pr[0-9]_Rid_18",
               "Pr[0-9]_Neo_17", "Pr[0-9]_Stj_17", "Pr[0-9]_Rid_17", "Pr[0-9]_Urb_17",
               "Pr[0-9]_Brt_16", "Pr[0-9]_Car_16", "Pr[0-9]_Scb_16", "Pr[0-9]_Urb_16", "Pr[0-9]_Stj_16",
               "Pr[0-9]_Brt_15", "Pr[0-9]_Car_15", "Pr[0-9]_Scb_15", "Pr[0-9]_Stj_15", "Pr[0-9]_Urb_15")

for(x in 1:length(studgrpExp)){
  
  #get the study names and row positions of the study group
  studGrp<- unique(data[grep(c(studgrpExp[x]), data$studyName, fixed=FALSE),'studyName'])
  ixgrp<- which(data$studyName %in% studGrp)
  stdNms<- as.character(data[ixgrp,'studyName'])
  
  #make block the prelim number
  block<- matrix(unlist(strsplit(stdNms, split="_")), nrow=3)[1,]
  data[ixgrp,'blockNumber']<- block
  data[ixgrp,'studyName']<- studgrpExp[x]
  
  #determine checks and change entryType info
  mtck<- unique(data.frame(nm=data[ixgrp,'germplasmName'], block))
  cks<- names(which(table(mtck$nm)==length(unique(block))))
  ixgrpCk<- which(data[ixgrp, 'germplasmName'] %in% cks)
  data[ixgrp,'entryType'][ixgrpCk]<- 'check'
  
  #edit design info
  data[ixgrp,'studyDesign']<- 'Augmented RCBD'
}


#convert yield and test weight to common units
data[,'Grain.yield...kg.ha.CO_321.0001218']<- convYld(data[,'Grain.yield...kg.ha.CO_321.0001218'])
data[,'Grain.test.weight...g.l.CO_321.0001210']<- convTwt(data[,'Grain.test.weight...g.l.CO_321.0001210'])
data[,"Plant.height...cm.CO_321.0001301"]<- data[,"Plant.height...cm.CO_321.0001301"] *0.393701
colnames(data)<- gsub('kg.ha.CO_321.0001218', "bu.ac", colnames(data))
colnames(data)<- gsub('g.l.CO_321.0001210', "lbs.bu", colnames(data))
colnames(data)<- gsub('..cm.CO_321.0001301', "inches", colnames(data))

#shorten the trait names to avoid errors in model fitting
ixCOs<- grep('CO_', colnames(data))
colnames(data)[ixCOs]<- matrix(unlist(strsplit(colnames(data)[ixCOs], split="CO_")), nrow=2)[1,]
colnames(data)<- gsub("Soil.borne.mosaic.plant.response...0.9.Response.Scale." , "SBMV", colnames(data))

############################
##Subset each trial
############################

#get unique study names
stdnms<- unique(data$studyName)

#exclude study with singularities
#stdnms<- stdnms[-c(81)]

#empty vector trial names with no data will be added
nodata<-c()

#empty vector of trait and trial combinations with only one replication
trtnonrep<-c()

#loop begins
for(i in 1:length(stdnms)){
  cat(i, '\n')
  #subset single trial
  trl<- data[which(data$studyName==stdnms[i]),]
 
  ############################
  ##Extract design information
  ############################
  
  #get vector of all possible traits
  trtnms<- colnames(data)[c(40:ncol(data)-1)]
  
  #exclude plots with low plant stands
  ps<- trl[,'Plant.stand...0.9.density.scale.']
  if(any(ps<5, na.rm=T)){
    trl<- trl[-which(ps<5),]
  }
  
  #exclude milling and baking traits, notes,and plant stand as response variables
  trtnms<-setdiff(trtnms, c("Flour.protein.content.....","Flour.yield.score.....",
                            "Grain.hardness...skcs.index.","Lactic.Acid.SRC.score.....",
                            "Softness.equivalent.score.....","Sucrose.SRC.score.....",
                            'Plant.stand...0.9.density.scale.', 'notes')) 
  
  #vector of traits used in the selected trial
  ttrt<- selectTraitcols(trl, trtnms)
  
  if(length(ttrt)==0){
    nodata<- append(nodata, stdnms[i])
  }else{
  
  #convert missing lodging scores to zero
  if("Lodging.incidence...0.9.percentage.scale." %in% ttrt){
    #convert missing lodging scores to zero
    vrldg<- var(trl$Lodging.incidence...0.9.percentage.scale., na.rm=T)
    if(vrldg>0){
      trl[which(is.na(trl$Lodging.incidence...0.9.percentage.scale.)),
          'Lodging.incidence...0.9.percentage.scale.']<- 0
    }
  }
  
  #check if there is blocking for the traits measured
    if(length(unique(trl$replicate))>1){
      blkfac<- 'replicate'
    }else{
      blkfac<- 'blockNumber'
    }
    minBlknos<- c()
    for(a in 1:length(ttrt)){
      minBlkno<- length(unique(na.omit(trl[,c(blkfac, ttrt[a])])[,blkfac]))
      minBlknos<- append(minBlknos, minBlkno)
    }
    minBlkno<- max(minBlknos)
    
    #remove traits with no replication
    repTF<- minBlknos==1
    if(any(repTF)){
      trtnonrep<- append(trtnonrep, paste(stdnms[i], ttrt[which(repTF)], sep="-"))
      ttrt<- ttrt[-which(repTF)]
    }
    
  
  #single-trait or multitrait model?
  #if(length(ttrt)==1){
  #  uvvmv<- "UV"
  #  clasfy<- 'germplasmDbId'
  #}else{
  #  uvvmv<- "MV"
  #  clasfy<- 'germplasmDbId:trait'
  #}
  #if(trl$studyDesign[1]=='Augmented RCBD'){
  #  uvvmv<- "UV"
  #  clasfy<- 'germplasmDbId'
  #}
    #Use the UV model only
    uvvmv<- "UV"
    clasfy<- 'germplasmDbId'
    
  #############################
  ## Create the fixed formula
  #############################
  if(uvvmv== "UV"){
    fxform<- paste(ttrt, "~1+germplasmDbId", sep="")
  }
  if(uvvmv== "MV"){
    fxform<- paste('cbind(', paste(ttrt, collapse=", "), ")~1+trait+trait:germplasmDbId", sep="")
  }
  
  
  #############################
  ## Create the random formula
  #############################
  rform<-NA
  #add the blocking factor if any
  if(minBlkno>=1){
    if(uvvmv=='MV'){
      rform<- "~at(trait):blockNumber"
    }else{
      rform<- "~blockNumber"
    }
      
    #for augmented designs use the checks to estimate the block effect
    if(trl$studyDesign[1]=='Augmented RCBD'){ 
      rform<- "~at(entryType, 'check'):blockNumber"
      clasfy<- paste(clasfy, ":entryType", sep="")
    }
  }
  
  #############################
  ## Convert variables to factors
  #############################
  trl$germplasmDbId<- as.factor(as.character(trl$germplasmDbId))
  trl$blockNumber<- as.factor(as.character(trl$blockNumber))
  trl$entryType<- as.factor(as.character(trl$entryType))
  trl$replicate<- as.factor(as.character(trl$replicate))
  
  #################################
  ## Fit model and extract results
  #################################
  
  if(uvvmv=='MV'){
    if(class(rform)=='logical'){
      mod<- suppressWarnings(asreml(fixed=as.formula(fxform), residual=~id(units):us(trait), data=trl, trace=FALSE, aom=T, workspace=64e6,  na.action = na.method(y='include', x='include')))
    }else{
      mod<- suppressWarnings(asreml(fixed=as.formula(fxform), random=as.formula(rform), residual= ~id(units):us(trait), data=trl, trace=FALSE, aom=T, workspace=64e6,  na.action = na.method(y='include', x='include')))
    }
    mod<- mkConv(mod)
    
    #ignr<- row.names(coefficients(mod)$fixed)[grep('block', row.names(coefficients(mod)$fixed))]
    #ignr<- unique(matrix(unlist(strsplit(ignr, "_")), nrow=2)[1,])
    
    p<- suppressWarnings(predict(mod, classify = clasfy, pworkspace=64e7)) #try this
    if(trl$studyDesign[1]=='Augmented RCBD'){
      p$pvals<- p$pvals[which(p$pvals[,2]=='test'),]
      p$pvals<-p$pvals[,-c(2)]
    }
    blues<- p$pvals
  }

  
  if(uvvmv=='UV'){
    for(z in 1:length(fxform)){
      if(class(rform)=='logical'){
        mod<- suppressWarnings(asreml(fixed=as.formula(fxform[z]), data=trl, trace=FALSE, aom=T,workspace=64e6,na.action = na.method(y='include', x='include')))
      }else{
        mod<- suppressWarnings(asreml(fixed=as.formula(fxform[z]), random=as.formula(rform), data=trl, trace=FALSE, aom=T,workspace=64e6,na.action = na.method(y='include', x='include')))
      }
      mod<- mkConv(mod)
      p<- suppressWarnings(predict(mod, classify = clasfy, pworkspace=64e6)) ## try this
      if(trl$studyDesign[1]=='Augmented RCBD'){
        p$pvals<- p$pvals[which(p$pvals[,2]=='test'),]
        p$pvals<-p$pvals[,-c(2)]
      }
      if(z==1){
        blues<- data.frame(studyName=stdnms[i], trait=ttrt[z], p$pvals)
      }
      if(z>1){
        blues<- rbind(blues, data.frame(studyName=stdnms[i], trait=ttrt[z], p$pvals))
      }
    }
  
  }
  if('entryType' %in% colnames(blues)){
    blues<- blues[,-match('entryType', colnames(blues))]
    blues<- unique(blues)
  }
  
  #add study name to the blues table
  if(uvvmv=='UV'){
    #df<- data.frame(studyName=stdnms[i], trait=ttrt, blues)
    df<- blues
    df<- df[,c(1,3,2, 4:ncol(df))]
  }
  if(uvvmv=='MV'){
    df<- data.frame(studyName=stdnms[i], blues)
  }
  
  #get residuals
  jpeg(file=paste(stdnms[i], "-residuals.jpeg", sep=""))
  resids<- resid(mod, type="stdCond")
  plot(resid(mod, type="stdCond"), main=stdnms[i])
  dev.off()
  
  #make potential outlier table
  mltTrl<- melt(trl, id.vars=c('observationUnitName','germplasmDbId'), measure.vars=ttrt)
  mltTrl<-mltTrl[order(mltTrl$variable),]
  mltTrl<-mltTrl[order(mltTrl$observationUnitName),]
  residsTab<- cbind(mltTrl, resids)
  outTab<- residsTab[which(sqrt(resids^2)>3),]
  
  
  #################################
  ## combine all analysis results into one table
  #################################    
  if(exists('dfall')){
    dfall<- rbind(dfall, df)
    outTabs<- rbind(outTabs, outTab)
  }else{
    dfall<- df
    outTabs<- outTab
  }

  }
}

name<- data[match(dfall$germplasmDbId, data$germplasmDbId),'germplasmName']
dfall<- data.frame(name, dfall)
write.csv(dfall, file='univariate predicted values table.csv')
write.csv(outTabs, file='univariate possible outliers.csv')




