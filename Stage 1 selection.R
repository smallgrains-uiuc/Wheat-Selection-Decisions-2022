setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
library(reshape)

#subset the stage 1 lines and checks
MasterEntryListILYT <- read.csv("~/Documents/Wheat/2022/MasterEntryListILYT.csv")
s1lines<- MasterEntryListILYT[which(MasterEntryListILYT$Stage ==1 | MasterEntryListILYT$Stage ==0),c('germplasmName','purdy.pedigree')]

#add small increase information
S600 <- read.csv("~/Documents/Wheat/2022/mirusFile_S600_seedlots.csv", row.names=1)
colnames(S600)[4]<- 'germplasmName'
s1lines<- merge(s1lines, S600, by='germplasmName', all.x=TRUE)
keepvdiscard<- read.csv('~/Documents/Wheat/2022/S600_KeepvDiscard_July5_table.csv', row.names=1)[,c('seedlot_name', 'action','Found')]
keepvdiscard<- keepvdiscard[match(s1lines$plot_name., keepvdiscard$seedlot_name),]
s1lines<- cbind(s1lines, keepvdiscard)

#add notes
#add stage 2 notes
notes<- read.csv('stage 1 line notes.csv')
notes$notes<- trimws(notes$notes)
notes$notes2<- paste(notes$studyName, "[", notes$notes, "]", sep="")
ugid<- unique(notes$germplasmName)
notes$notes2<- gsub("_", "", notes$notes2)
notes3<- c()
for(i in 1:length(ugid)){
  notes3<- append(notes3, paste(notes[which(ugid[i] == notes$germplasmName),'notes2'], collapse=", "))
}
notesdf<- data.frame(ugid, notes3)
notesdf<- notesdf[match(s1lines$germplasmName, notesdf$ugid),]
s1lines<- cbind(s1lines, notesdf)

#read in the gebvs for agronomic traits
gebv<- read.csv('gebvsIL22_AgronomicJuly30.csv', row.names=1)
gebvW<- cast(gebv, germplasmName~trait, value='predicted.value')

#read gebvs for scab
gebvscb<- read.csv('GEBVscabJul26.csv', row.names=1)
gebvscbW<- cast(gebvscb, germplasmName~trait, value='predicted.value')

#all gebv
allgebv<- merge(gebvW, gebvscbW, by='germplasmName',all=TRUE)

#merge files
s1lines2<- merge(allgebv, s1lines, by='germplasmName', all.x=FALSE, all.y=TRUE)

#add merit
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
s1lines2<- data.frame(s1lines2, merit=netMerit(s1lines2$`Heading.time...Julian.date..JD..`,
                                        s1lines2$`Grain.yield...bu.ac`,
                                        s1lines2$`FHB.DON.content...ppm.`,
                                        s1lines2$`FHB.grain.incidence.....`, 
                                        s1lines2$`Grain.test.weight...lbs.bu`,
                                        wheat_price0, soybean_price))

s1lines2<- s1lines2[,c(1,22,3,2,4:7,21, 8:20)]
for(i in 2:8){
  s1lines2[,i]<- round(s1lines2[,i], 1)
}
#write.csv(s1lines2, file='stage1to2_selection.csv')

#make the list of seed needed from S600
ILYT2023 <- read.csv("~/Documents/Wheat/2023/ILYT_designplan2023.csv")
seed_from_S600<- merge(ILYT2023, S600, by='germplasmName')[,c(1,19:23, 25,26)]
#write.csv(seed_from_S600, file='selected_seedlotsS600.csv')

#make list for pure rows
fb<- read.csv('~/Documents/Wheat/2022/wSI_Urb_22.csv')
sikeep<- read.csv('~/Documents/Wheat/2023/SmallIncrease2022_keep.csv')
colnames(fb)[1]<- 'plotName'
sikeep<- merge(sikeep, fb, by='plotName')
View(sikeep)
write.csv(sikeep, file='~/Documents/Wheat/2023/SmallIncrease2022_keep.WithPRno.csv')
