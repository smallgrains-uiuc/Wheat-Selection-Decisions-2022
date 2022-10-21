library(reshape)
setwd("~/Documents/Wheat/2022")
parentLines0<- read.csv('CrossingParents2022.csv')[1:37,]
parentLines<- read.csv('CrossingParents2022.csv')[1:37,c(1,2)]
colnames(parentLines)[2]<- 'germplasmName'

#read in the gebvs for agronomic traits
setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
gebv<- read.csv('gebvsIL22_AgronomicJuly30.csv', row.names=1)
gebvW<- cast(gebv, germplasmName~trait, value='predicted.value')

#read gebvs for scab
gebvscb<- read.csv('GEBVscabJul26.csv', row.names=1)
gebvscbW<- cast(gebvscb, germplasmName~trait, value='predicted.value')

#all gebv
allgebv<- merge(gebvW, gebvscbW, by='germplasmName',all=TRUE)

#merge files
parentLines<- merge(allgebv, parentLines, by='germplasmName', all.x=FALSE, all.y=TRUE)

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

#make all combinations and their reciprocals
comb<- combn(c(1:37), 2)
for(i in 1:ncol(comb)){
    pair<- comb[,i]
    a<- parentLines[which(parentLines$ParentSel==pair[1]),]
    b<- parentLines[which(parentLines$ParentSel==pair[2]),]
    f1<- colMeans(rbind(a[,-1], b[,-1]))
    f1["FHB.grain.incidence....."]<- 0
    nm<- netMerit(f1["Heading.time...Julian.date..JD.."],
              f1["Grain.yield...bu.ac"],
              f1["FHB.DON.content...ppm."],
              f1["FHB.grain.incidence....."], 
              f1["Grain.test.weight...lbs.bu"],
              wheat_price0, soybean_price)[[1]]
    aname<- as.character(a$germplasmName)
    bname<- as.character(b$germplasmName)
    cross<- paste(aname, bname, sep=" x ")
    recip<- paste(bname, aname, sep=" x ")
    crossN<- paste(pair[1], pair[2], sep=" x ")
    recipN<- paste(pair[2], pair[1], sep=" x ")
    f1<- f1[1:6]
    f1<- t(data.frame(f1))
    df<- data.frame(p1=pair[1], p2=pair[2], crossN, reciprocalNm=recipN, netMerit=nm, cross, reciprocal=recip, f1)
    if(i==1){
      dfall<- df
    }else{
      dfall<- rbind(dfall, df)
    }
}
#make copy of original table before subsetting
dfall0<- dfall

#eliminate below average combinations
dfall<- dfall[which(dfall$netMerit> mean(dfall$netMerit)),]

#eliminate tall families
dfall<- dfall[which(dfall$`Plant.height.inches`< 38),]

#remove crosses already made
setwd("~/Documents/Wheat/2023")
prev<- read.csv('2021_ParentsCross_list.xlsx - Sheet1.csv')
prev_f1<- read.csv('2022-05-11-03-29-44_WWX1.22.F21seed.inventory_table.csv')[1:55,]
crosses<- paste(prev[match(prev_f1$P1, prev$Parent21no.),'Name.'], prev[match(prev_f1$P2, prev$Parent21no.),'Name.'], sep=" x ")
crossesRecip<- paste(prev[match(prev_f1$P2, prev$Parent21no.),'Name.'], prev[match(prev_f1$P1, prev$Parent21no.),'Name.'], sep=" x ")
crossesPrev<- c(crosses, crossesRecip)
dfall<- dfall[-which(as.character(dfall$cross) %in% crossesPrev),]

#trim by family
dfall$RankInfamily<- 0

dfall$p1uped<- parentLines0[match(dfall$p1, parentLines0$ParentSel),'uped']
dfall$p2uped<- parentLines0[match(dfall$p2, parentLines0$ParentSel),'uped']

for(i in 1:27){
  ix<- c(which(dfall$p1uped==i), which(dfall$p2uped==i))
  dfall[ix,'netMerit']
  dfall[ix,'RankInfamily']<- rank(-dfall[ix,'netMerit'])
  sub<- dfall[ix,][which(dfall[ix,'RankInfamily']<9),]
  if(i==1){
    suball<- sub
  }else{
    suball<- rbind(suball, sub)
  }
}
suball<- suball[,c(1:13)]
suball<- unique(suball)

#check how many times each one is used
dfall0$select<- 0
dfall0[match(suball$crossN, dfall0$crossN),'select']<- 1
dfall0$rank<- 667
dfall0[which(dfall0$select==1),'rank']<- rank(-dfall0[which(dfall0$select==1),'netMerit'])
dfall0[which(dfall0$select==0),'rank']<- rank(-dfall0[which(dfall0$select==0),'netMerit'])+178
setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
for(i in c(5,8:12)){
  dfall0[,i]<- round(dfall0[,i],1)
}
dfall0$relValue<- round((dfall0$netMerit/mean(dfall0$netMerit)) *100, 0)
write.csv(dfall0, file='cross ranknigs fall 2022.csv')

#parentLines<- data.frame(parentLines, merit=netMerit(parentLines$`Heading.time...Julian.date..JD..`,
#                                            parentLines$`Grain.yield...bu.ac`,
#                                            parentLines$`FHB.DON.content...ppm.`,
#                                            parentLines$`FHB.grain.incidence.....`, 
#                                            parentLines$`Grain.test.weight...lbs.bu`,
#                                            wheat_price0, soybean_price))
#parentLines$relativeMerit<- parentLines$merit- mean(parentLines$merit)

