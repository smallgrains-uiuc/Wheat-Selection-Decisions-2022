setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
a<- read.csv('Adv_Neo_21blues.csv', row.names=1)
c<- read.csv('Adv_Stj_21blues.csv', row.names=1)
d<- read.csv('Adv_Urb_21blues.csv', row.names=1)
a2<- read.csv('Pr_Neo_21blues.csv', row.names=1)
c2<- read.csv('Pr_Stj_21blues.csv', row.names=1)
d2<- read.csv('Pr_Urb_21blues.csv', row.names=1)
f<- read.csv('YT blues 2022 lexky.csv', row.names=1)
f<- f[,c('study', 'pvals.name', 'trait','pvals.predicted.value','pvals.std.error','pvals.status')]
colnames(f)<- colnames(d)
g<- read.csv('YT blues 2022 masonmi.csv', row.names=1)
g<- g[,c('study', 'pvals.Line.Name', 'trait','pvals.predicted.value','pvals.std.error','pvals.status')]
colnames(g)<- colnames(d)

setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2022")
b<- read.csv('Adv_Scb_21blues.csv', row.names=1)
b2<- read.csv('Pr_Scb_21blues.csv', row.names=1)
b3<- read.csv('IL_Scb_22blues.csv', row.names=1)
e<- read.csv('univariate predicted values table.csv', row.names=1)
e<- e[,c('studyName', 'name', 'trait','predicted.value','std.error','status')]
e<- e[-which(e$studyName %in% c("Adv_Neo_21","Adv_Scb_21","Adv_Stj_21","Adv_Urb_21")),]
colnames(e)<- colnames(d)
YT22<- read.csv('YT blues 2022.csv', row.names=1)[,c('study', 'pvals.germplasmName', 'trait','pvals.predicted.value','pvals.std.error','pvals.status')]
colnames(YT22)<- colnames(d)
all<- rbind(a, b, c, d, a2, b2, b3, c2, d2,  e, YT22, f, g)

#add weight to the file
var<- all$std.error^2
wt<- 1/var
all<- data.frame(all, wt)

#add more factors
loc<- matrix(unlist(strsplit(as.character(all$study), split="_")), nrow=3)[2,]
loc[which(loc=='StJ')]<- 'Stj'
year<- matrix(unlist(strsplit(as.character(all$study), split="_")), nrow=3)[3,]
site<- paste(loc, year, sep="_")
all<- data.frame(site, loc, year, all)
all$trait<- gsub("CO_321.0001154", "", all$trait)
all$trait<- gsub("CO_321.0001155", "", all$trait)
all$trait<- gsub("CO_321.0001149", "", all$trait)
all$trait<- gsub("CO_321.0001440", "", all$trait)
all$trait<- gsub("CO_321.0001233", "", all$trait)
traitsOfinterest<- unique(all$trait)[1:6]
all2<- all[which(all$trait %in% traitsOfinterest),]
write.csv(all2, 'ILTrainingSetPhenotypesJuly.21.2022.csv')

#remove an outlier 
all2<- read.csv('ILTrainingSetPhenotypesJuly.21.2022.csv', row.names=1)
all2<- all2[-which(all2$predicted.value>129 & all2$study=='YT_Stj_22'),]
write.csv(all2, 'ILTrainingSetPhenotypesJuly.21.2022.csv')


