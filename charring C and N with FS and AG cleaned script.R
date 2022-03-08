data<-read.csv("~/OneDrive - Nexus365/Feedsax/Charring results/Charring experiment/Experimental_charring/CharringdataOct29N1alldata.csv")
library(plyr)
library(nlme)
library(beeswarm)
library(lm.beta)

# 
# summary <- data.frame(ddply(data, c("Species", "temp", "time"), nrow))
# 
# summary2 <- ddply(data, c("Species", "temp", "time"), function(x) c(d15N=mean(x$normd15N), sd=sd(x$normd15N), pcN=mean(x$pcN), d13C=mean(x$normd13C), sd=sd(x$normd13C), pcN=mean(x$pcC), CN=mean(x$CN_Crun)))
#summary2 <- ddply(data, c("Species", "temp", "time"), function(x) c( d13C=mean(x$normd13C), sd=sd(x$normd13C), pcC=mean(x$pcC), CN=mean(x$CN_Crun)))
data2<-data

data2$Species<-droplevels(data2$Species)

TT2 <- paste(data$Species,data2$temp, data2$time, sep= "")
data2<-data.frame(data2,TT2)

#####work out which Rye to use from the repeats of 245 24h
RyeR<-data2[data2$TT2=="rye24524",]
RyeR$ID<-gsub('A','',RyeR$ID)
RyeR$ID<-gsub('B','',RyeR$ID)
RyeR$ID<-gsub('C','',RyeR$ID)

boxplot(RyeR$normd13C~RyeR$ID)

###Remove bad Rye samples 
RYE<-data2[data2$ID=="RCSPARERYE24524A"|data2$ID=="RCSPARERYE24524B"|data2$ID=="RCSPARERYE24524C",]
data2<-data2[data2$TT!="rye24524",]
data2<-rbind(data2,RYE)

data2<-data2[data2$Species!="RyeX",]
data2<-data2[data2$Species!="Wheat",]
data2<-data2[data2$ID!="RYE2304C",]


###work out if barley and wheat can be compared
bar_whe$normd15N<-as.numeric(bar_whe$normd15N)


bar_whe<-data2[data2$TT2=="Barley00"|data2$TT2=="Free-T wheat00"|data2$TT2=="BW00"|data2$TT2=="HBH00",]

boxplot(bar_whe$normd13C~bar_whe$TT2)
boxplot(as.numeric(bar_whe$normd15N)~bar_whe$TT2)
bar<-bar_whe[bar_whe$TT2=="Barley00"|bar_whe$TT2=="HBH00",]
whe<-bar_whe[bar_whe$TT2=="Free-T wheat00"|bar_whe$TT2=="BW00",]
t.test(bar$normd13C~bar$TT2)
t.test(as.numeric(bar$normd15N)~bar$TT2)
t.test(whe$normd13C~whe$TT2)
t.test(whe$normd15N~whe$TT2)

###work out if barley and wheat can be compared
bar$normd15N<-as.numeric(bar$normd15N)
bar3<-bar[bar$ID!="BAR0F",]
t.test(bar3$normd13C~bar3$TT2)
t.test(bar3$normd15N~bar3$TT2)

bar4<-bar3[bar3$ID!="BAR0D",]
t.test(bar4$normd13C~bar4$TT2)
t.test(bar4$normd15N~bar4$TT2)

# bar5<-bar4[bar4$ID!="BAR0C",]
#  t.test(bar5$normd13C~bar5$TT2)
#  t.test(bar5$normd15N~bar5$TT2)

 bar5<-bar4[bar4$ID!="BAR0E",]
t.test(bar5$normd13C~bar5$TT2)
t.test(bar5$normd15N~bar5$TT2)

### Only using rye, oat, BW and HB
data3<- data2[data2$Species=="BW"|data2$Species=="HBH"|data2$Species=="rye"|data2$Species=="oat"|data2$Species=="Free-T wheat"|data2$Species=="Barley",]

data3$Species<-gsub('Free-T wheat','BW',data3$Species)

data3$Species<-gsub('HBH','Barley',data3$Species)

data3$normd15N<-as.numeric(data3$normd15N)
data3$d15Nsd<-as.numeric(data3$d15Nsd)



summary <- data.frame(ddply(data3, c("Species", "temp", "time"), nrow))

########Carbon versions of all the figures ###############################

## Compare random slopes and intercetps models to just multiple regression
data3<-data3[data3$ID!="BAR0F",]
data3<-data3[data3$ID!="BAR0E",]
data3<-data3[data3$ID!="BAR0D",]
summary <- summarise(data3, mean(d13Csd, na.rm= TRUE))
summaryn <- summarise(data3, mean(d15Nsd, na.rm= TRUE))

mean(data3$d15Nsd, na.rm= TRUE)

summarycharred<- ddply(data3, c("Species","char"), function(x) c(d15N=mean(x$normd15N), sd=sd(x$normd15N), max=max(x$normd15N), min=min(x$normd15N), d13C=mean(x$normd13C), sd=sd(x$normd13C), max=max(x$normd13C), min=min(x$normd13C)))

summarycharred<- ddply(data3, c("Species","TT2"), function(x) c(d15N=mean(x$normd15N), sd=sd(x$normd15N), max=max(x$normd15N), min=min(x$normd15N), d13C=mean(x$normd13C), sd=sd(x$normd13C), max=max(x$normd13C), min=min(x$normd13C)))
# These are the values reported in table 5
###All temperatures 215-300
lm1 <- lm(normd13C ~ Species, data = data3)
summary(lm1)
lm2 <- lm(normd13C ~ char + Species, data=data3)#### and this one for CI - this was are saying predict d13c from charring and species - so species has a direct effect 
summary(lm2)
confint(lm2)

### without 300 data
no300data3<-data3[data3$temp!="300",]

lm1 <- lm(normd13C ~ Species, data = no300data3)
summary(lm1)
lm2 <- lm(normd13C ~ char + Species, data=no300data3)### data in paper
summary(lm2)
confint(lm2)#### and this one for CI

### without 215 data
no215data3<-data3[data3$temp!="215",]
lm1 <- lm(normd13C ~ Species, data = no215data3)
summary(lm1)
lm2 <- lm(normd13C ~ char + Species, data=no215data3)### data in paper
summary(lm2)
confint(lm2)#### and this one for CI

### without 215 and 300 data
no215300data3<-no215data3[no215data3$temp!="300",]
lm1 <- lm(normd13C ~ Species, data = no215300data3)
summary(lm1)
lm2 <- lm(normd13C ~ char + Species, data=no215300data3)### data in paper
summary(lm2)
confint(lm2)

##############
# 
# dataFC<- data3
# TT <- paste(dataFC$temp, dataFC$time, sep= "")
# dataFC <- data.frame(dataFC, TT)
# 
# taxon.1<-dataFC[dataFC$Species==unique(dataFC$Species)[1],]
# BAR.testFC<-aov(normd13C~TT, data=taxon.1) ##### no diff between any samples p =0.152
# posthoc<-TukeyHSD(x=BAR.testFC, 'TT')
# 
# lmBAR<-lm(normd13C~char, data=taxon.1)
# lmBAR<-lm(normd13C~time+temp, data=taxon.1)
# 
# 
# 
# taxon.1_215<-taxon.1[taxon.1$temp!="215",]
# BAR.testFC<-aov(normd13C~TT, data=taxon.1_215) ##### no diff between any samples p =0.152
# posthoc<-TukeyHSD(x=BAR.testFC, 'TT')
# 
# lmBAR<-lm(normd13C~charr, data=taxon.1_215) # not significant
# 
# taxon.2<-dataFC[dataFC$Species==unique(dataFC$Species)[2],]
# FT.testFC<-aov(normd13C~TT, data=taxon.2) ##### sig diff p = 0.00719
# posthoc<-TukeyHSD(x=FT.testFC, 'TT')#### 30024-00,
# lmFT<-lm(normd13C~char, data=taxon.2)
# lmFT<-lm(normd13C~time+temp, data=taxon.2)
# 
# 
# taxon.2_215<-taxon.2[taxon.2$temp!="215",]
# FT.testFC<-aov(normd13C~TT, data=taxon.2_215) ##### no diff between any samples p =0.76
# posthoc<-TukeyHSD(x=BAR.testFC, 'TT')
# 
# lmFT<-lm(normd13C~char, data=taxon.2_215)
# 
# 
# 
# taxon.3<-dataFC[dataFC$Species==unique(dataFC$Species)[3],]
# Oat.testFC<-aov(normd13C~TT, data=taxon.3) ##### sig diff between all samples
# posthoc<-TukeyHSD(x=Oat.testFC, 'TT')
# lmOat<-lm(normd13C~time+temp, data=taxon.3)
# 
# taxon.3_215<-taxon.3[taxon.3$temp!="215",]
# Oat.testFC<-aov(normd13C~TT, data=taxon.3_215) ##### no diff between any samples p =0.76
# posthoc<-TukeyHSD(x=Oat.testFC, 'TT')
# 
# lmOat<-lm(normd13C~char, data=taxon.3_215)
# 
# 
# taxon.4<-dataFC[dataFC$Species==unique(dataFC$Species)[4],]
# 
# Rye.testFC<-aov(normd13C~TT, data=taxon.4) ##### sig diff between all samples
# posthoc<-TukeyHSD(x=Rye.testFC, 'TT')
# lmRye<-lm(normd13C~char, data=taxon.4)
# 
# Rye.testFC<-aov(normd13C~TT, data=taxon.4) ##### sig diff between all samples
# posthoc<-TukeyHSD(x=Rye.testFC, 'TT')
# lmRye<-lm(normd13C~char, data=taxon.4)
# lmRye<-lm(normd13C~time+temp, data=taxon.4)
# 
# taxon.4_215<-taxon.4[taxon.4$temp!="215",]
# rye.testFC<-aov(normd13C~TT, data=taxon.4_215) ##### no diff between any samples p =0.76
# posthoc<-TukeyHSD(x=rye.testFC, 'TT')
# 
# lmrye<-lm(normd13C~char, data=taxon.4_215)


######Compare LM for the just charred material for each taxon:
#####

just.charred <- data3[data3$char=="charred",]
TT <- paste(just.charred$temp, just.charred$time, sep= "")
just.charred <- data.frame(just.charred, TT)

charred.lm2 <- lm(normd13C ~ Species + temp + time, data=just.charred) ######### used in paper - significant p value for temp
summary(charred.lm2)


no300<-just.charred[just.charred$temp!="300",]
no300.lm2 <- lm(normd13C ~ Species + temp + time, data=no300) 
summary(no300.lm2)


no215<-just.charred[just.charred$temp!="215",]
no215.lm2 <- lm(normd13C ~ Species + temp + time, data=no215) 
summary(no215.lm2)


no215300<-no215[no215$temp!="300",]
no215300.lm2 <- lm(normd13C ~ Species + temp + time, data=no215300)
summary(no215300.lm2)

# 
#####
# taxon.1 <- just.charred[just.charred$Species==unique(just.charred$Species)[1],]
# BH.test <- aov(normd13C ~ TT, data=taxon.1)  ##No SIG DIFFS AMONG OAT p = 0.0947 
# summary(BH.test)
# posthoc<-TukeyHSD(x=BH.test, 'TT')
# 
# taxon.1$TT<-factor(taxon.1$TT, levels=c("2154", "2158","21524", "2304", "2308", "23024", "2454", "2458", "24524", "2604", "2608", "26024", "3004", "3008", "30024"))
# plot(taxon.1$normd13C~taxon.1$TT, ylab="d13C, HB", xlab="temp/time", ylim=c(-27 ,-24.5))
# par(new=TRUE)
# beeswarm(taxon.1$normd13C~taxon.1$TT, ylab="d13C, HB", xlab="temp/time",ylim=c(-27, -24.5))
# 
# 
# taxon.2 <- just.charred[just.charred$Species==unique(just.charred$Species)[2],]
# BW.test <- aov(normd13C ~ TT, data=taxon.2)  ##SIG DIFFS AMONG RYE! p = 0.000166??? - issue with the sample from 245 24hrs - re run as they look like a misslabel - way to low
# summary(BW.test)
# taxon.2$TT<-factor(taxon.2$TT, levels=c("2154", "2158","21524", "2304", "2308", "23024", "2454", "2458", "24524", "2604", "2608", "26024", "3004", "3008", "30024"))
# plot(taxon.2$normd13C~taxon.2$TT,ylab="d13C, BW", xlab="temp/time", ylim=c(-27.5,-26.5))
# par(new=TRUE)
# beeswarm(taxon.2$normd13C~taxon.2$TT, ylab="d13C, BW", xlab="temp/time",ylim=c(-27.5,-26.5 ))
# 
# 
# posthoc<-TukeyHSD(x=BW.test, 'TT')
# 
# taxon.3 <- just.charred[just.charred$Species==unique(just.charred$Species)[3],]
# OAT.test <- aov(normd13C ~ TT, data=taxon.3)  ##NO SIG DIFFS AMONG SPT p = 0.129
# summary(OAT.test)
# taxon.3$TT<-factor(taxon.3$TT, levels=c("2154", "2158","21524", "2304", "2308", "23024", "2454", "2458", "24524", "2604", "2608", "26024", "3004", "3008", "30024"))
# 
# plot(taxon.3$normd13C~taxon.3$TT,ylab="d13C, oat", xlab="temp/time",ylim=c(-28.5, -26.5))
# par(new=TRUE)
# beeswarm(taxon.3$normd13C~taxon.3$TT, ylab="d13C, oat", xlab="temp/time",ylim=c(-28.5, -26.5))
# 
# taxon.4 <- just.charred[just.charred$Species==unique(just.charred$Species)[4],]
# RYE.test <- aov(normd13C ~ TT, data=taxon.4)  ##NO SIG DIFFS AMONG SPT p = 0.129
# summary(RYE.test)
# taxon.4$TT<-factor(taxon.4$TT, levels=c("2154", "2158","21524", "2304", "2308", "23024", "2454", "2458", "24524", "2604", "2608", "26024", "3004", "3008", "30024"))
# 
# plot(taxon.4$normd13C~taxon.4$TT,ylab="d13C, rye", xlab="temp/time", ylim=c( -26.5,-24.5))
# 
# par(new=TRUE)
# beeswarm(taxon.4$normd13C~taxon.4$TT, ylab="d13C, rye", xlab="temp/time",ylim=c( -26.5,-24.5))
# 

##############################
data2<-data3

summary <- data.frame(ddply(data2, c("Species", "temp", "time"), nrow))
summary2 <- ddply(data2, c("Species", "temp", "time"), function(x) c(d15N=mean(x$normd15N), sd=sd(x$normd15N), pcN=mean(x$pcN), d13C=mean(x$normd13C), sd=sd(x$normd13C), pcN=mean(x$pcC)))
#summary2 <- ddply(data, c("Species", "temp", "time"), function(x) c( d13C=mean(x$normd13C), sd=sd(x$normd13C), pcC=mean(x$pcC), CN=mean(x$CN_Crun)))

#####

# ##Plot how the model fits the data  #Cmodelfit.pdf
# taxon.names <- c("HB","BW","Oat", "Rye")
# 
# par(mfrow=c(1,1))
# par(mar=c(5,9,3,9))
# plot(as.numeric(data2$char)-1, data2$normd13C, pch=1, col=as.numeric(data2$Species), xlim=c(-0.2,1.2 ), ylim=c(-28.5, -23.8), cex=1.3, axes=F, xlab="", ylab=expression(paste(delta^{13}, "C (\u2030)")))
# axis(2)
# axis(1, labels=c("Charred", "Fresh"), at=c(0, 1))
# first <- coef(lm1)[[1]]
# avalues <- c(first, coef(lm1)[[2]]+first, coef(lm1)[[3]]+first, coef(lm1)[[4]]+first)
# for (i in 1:length(unique(data2$Species))){
# 	ilist <- c(1,2,3,4)
# abline(a=avalues[i], b=-0.106, col=ilist[i])
# }
# legend("topright", taxon.names, fill=as.numeric(unique(data2$Species)), cex=0.8, box.lty=0)

# 
# ##Colour - coloured by temperature = Cbytaxoncolour.pdf
# mycols1 <- c(rainbow(12))
# mycols2 <- c(mycols1[12], mycols1[1:11])
# mycols <- c("black", rev(mycols2))
# palette(mycols)
# par(oma=c(3,6,3,3))
# taxon.names <- c("Barley","Oat", "Rye", "Spelt")
# par(mfrow=c(1,5))
# par(mar=c(3,0,1,3))
# batch <- data2[data2$Species==unique(data2$Species)[1],]
# TT <- paste(batch$temp, batch$time, sep= "")
# plot(as.numeric(batch$char)-1, batch$normd13C, pch=1, col=as.numeric(as.factor(TT)), xlim=c(-0.2,1.2 ), axes=F, xlab="", ylim=c(min(data2$normd13C), max(data2$normd13C)), main=taxon.names[1], ylab="", cex=1.2)
# axis(1, labels=c("Charred", "Fresh"), at=c(0, 1))
# axis(2)
# #axis(4,labels=F, line=-11.5)
# mtext(expression(paste(delta^{13}, "C (\u2030)")), side=2, line=2, cex=0.8)
# abline(lm(batch$normd13C ~ batch$char), col="grey", lwd=0.8)
# for (i in 2){
# #par(mar=c(3,0,1,3))
# batch <- data2[data2$Species==unique(data2$Species)[i],]
# TT <- paste(batch$temp, batch$time, sep= "")
# plot(as.numeric(batch$char)-1, batch$normd13C, pch=1, col=as.numeric(as.factor(TT)), xlim=c(-0.2,1.2 ), axes=F, xlab="", ylim=c(min(data2$normd13C), max(data2$normd13C)), main=taxon.names[i], ylab="", cex=1.2)
# axis(1, labels=c("Charred", "Fresh"), at=c(0, 1))
# axis(2)
# #axis(4,labels=F, line=-11)
# abline(lm(batch$normd13C ~ batch$char), lwd=0.8, col="grey")
# }
# for (i in 3){
# par(mar=c(3,0,1,3))
# batch <- data2[data2$Species==unique(data2$Species)[i],]
# TT <- paste(batch$temp, batch$time, sep= "")
# plot(as.numeric(batch$char)-1, batch$normd13C, pch=1, col=as.numeric(as.factor(TT)), xlim=c(-0.2,1.2 ), axes=F, xlab="", ylim=c(min(data2$normd13C), max(data2$normd13C)), main=taxon.names[i], ylab="", cex=1.2)
# axis(1, labels=c("Charred", "Fresh"), at=c(0, 1))
# axis(2)
# mtext(expression(paste(delta^{13}, "C (\u2030)")), side=2, cex=0.8, line=2)
# #axis(4,labels=F, line=-11)
# abline(lm(batch$normd13C ~ batch$char), lwd=0.8, col="grey")
# }
# 
# par(mar=c(9,5,5,0))
# mycols <- c(rev(mycols2), "black")
# palette(mycols)
# plot(c(rep(1:3, 4), 4), c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), 1), pch=22, bg=mycols, cex=3, axes=F, xlim=c(0.5, 4.5), ylim=c(0.5, 4.5), xlab="", ylab="")
# mtext("Charring temperatures", side=3, line=2, cex=0.7)
# mtext("and times", side=3, line=1, cex=0.7)
# mtext(expression(paste("Temperature [",degree,"C]")), side=2, line=1.5, cex=0.7)
# mtext("Time (hours)       ", side=1, line=1, cex=0.7)
# axis(2, labels=c("215","230", "245", "260"), at=c(1:4), las=2, lwd=0, line=-1, tick=T)
# axis(1, labels=c("4", "8", "24"), at=c(1:3), lwd=0, line=-1)
# axis(1, labels=c("Fresh"), at=4, las=2, tick=F, line=-1)


#######
Ccharoff <- rep(0, nrow(data2))
for (i in 1:length(unique(data2$Species))){
offset.start <- summary2$d13C[summary2$Species==unique(data2$Species)[i] & summary2$time==0]
offsetted.values <- (data2$normd13C[data2$Species==unique(data2$Species)[i]])-offset.start
Ccharoff[data2$Species==unique(data2$Species)[i]] <- offsetted.values}
data3 <- data.frame(data2, Ccharoff)
#data3 <- data.frame(data, Ncharoff,Ccharoff)

summary3 <- ddply(data3, c("Species", "temp", "time"), function(x) c( d13C=mean(x$normd13C), sd=sd(x$normd13C), pcC=mean(x$pcC), CN=mean(x$CN_Crun), charoff=mean(x$Ccharoff)))

####Figure 1 in paper
png(filename="FigureCcharoff.png",
    units="cm",
    width=14, 
    height=13, 
    res=300,
    pointsize=8)

par(oma=c(4,6,4,2))
taxon.names <- c("Barley",  "BW", "Oat","Rye")
data.sort <- data3[order(data3$Species, data3$temp, data3$time),]
batch <- data.sort[data.sort$Species != "NB" & data.sort$Species != "HBA"& data.sort$Species != "PBA", ]
summary4 <- ddply(data3, c("Species"), function(x) c( d13C=mean(x$normd13C), sd=sd(x$normd13C), pcC=mean(x$pcC), charoff=mean(x$Ccharoff)))

pch.list <- rep(0, length(batch$time))
pch.list[batch$time==0] <- 8
pch.list[batch$time==4] <- 21
pch.list[batch$time==8] <- 22
pch.list[batch$time==24] <- 23
col.list <- rep(0, length(batch$temp))
col.list[batch$temp==0] <- "black"
col.list[batch$temp==215] <- "white"
col.list[batch$temp==230] <- "lightgray"
col.list[batch$temp==245] <- "darkgray"
col.list[batch$temp==260] <- "black"
col.list[batch$temp==300]<-"red"
par(mfrow=c(1,5))
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[1],]
plot(batch2$Ccharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-0.9, 1.25), xlim=c(-1, 50), main=taxon.names[1])
#legend("topleft", c( "Uncharred","215", "230", "245", "260",  "4h", "8h", "24h"), pch=c(8, 22,22,22,22,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray", "black", "white", "white", "white"), bg="white", pt.cex=c(1.2,2, 2, 2, 2, 1.2,1.2,1.2), cex=1.2, box.lty=0)
axis(2)
box()
abline(h=0)
#abline(h=0.63)
mtext(expression(paste(Delta^{13}, "C (\u2030) uncharred-charred")), side=2, line=2, cex=0.8)

for (i in 2){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$Ccharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-0.9, 1.25), xlim=c(-1, 50), main=taxon.names[i])
#axis(2)
box()
abline(h=0)
}
taxon.names <- c( "Oat","Rye")
batch <- data.sort[data.sort$Species != "Barley" & data.sort$Species != "BW", ]
pch.list <- rep(0, length(batch$time))
pch.list[batch$time==0] <- 8
pch.list[batch$time==4] <- 21
pch.list[batch$time==8] <- 22
pch.list[batch$time==24] <- 23
col.list <- rep(0, length(batch$temp))
col.list[batch$temp==0] <- "black"
col.list[batch$temp==215] <- "white"
col.list[batch$temp==230] <- "lightgray"
col.list[batch$temp==245] <- "darkgray"
col.list[batch$temp==260] <- "black"
col.list[batch$temp==300]<-"red"

for (i in 1){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$Ccharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-0.9, 1.25), xlim=c(-1, 50), main=taxon.names[i])
#axis(2)
box()
abline(h=0)
}

for (i in 2){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$Ccharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-0.9, 1.25), xlim=c(-1, 50), main=taxon.names[i])
#axis(2)
box()
abline(h=0)
}

plot(batch2$pcN, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 45), ylim=c(38, 85), type="n")
legend(1, 70, c("", "215 ˚C", "230 ˚C", "245 ˚C", "260 ˚C", "300 ˚C","", "Uncharred", "4h", "8h", "24h"), pch=c(NA, 22,22,22,22,22,NA, 8,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray", "black", "red",NA, NA, "white", "white", "white"), bg="white", xpd=T, pt.cex=c(2, 2,2,2,2,2,2,1,1,1,1), cex=1)
text(14, 69, "Temperature", xpd=T)
text(7.7, 62.2, "Time", xpd=T)
dev.off()


###### 
# 
# ##Plot offsets by the mean of each burn for each taxon  #Cwithinbatch.pdf ------<<<<<< Can't work out yet
# par(oma=c(10,5,10,3))
# taxon.names <- c("Oat",  "Rye", "Spelt")
# par(mfrow=c(3,1))
# par(mar=c(0,5,0,7))
# 
# for (i in 1:3){
# 	batch <- data2[data2$Species==unique(data2$S)[i], c( 7, 22)]  #selectnormd13C
# 	batch$Packno <- as.numeric(paste(batch$CharrNo))
# 	batch.sort <- batch[order(batch$Packno),]
# 	offsetted.value <- 0
# 	for (j in 1:16){
# 		small.batch <- batch.sort[batch.sort$Packno==j,]
# 		new.value <- small.batch$normd13C-mean(small.batch$normd13C)
# 		offsetted.value <- c(offsetted.value, new.value)
# 	}
# 	offsetted.value <- offsetted.value[2:length(offsetted.value)]
# 	plot(offsetted.value, batch.sort$Packno, xlim=c(-1, 1), axes=F, ylim=c(0, 16), ylab="", xlab="")
# 	axis(2, labels=c(expression(paste("4 h at 215 ",degree,"C ")),
# 	expression(paste("8 h at 215 ",degree,"C ")),
# 	expression(paste("24 h at 215 ",degree,"C ")),
# 	expression(paste("4 h at 230 ",degree,"C ")),
# 	expression(paste("8 h at 230 ",degree,"C ")),
# 	expression(paste("24 h at 230 ",degree,"C ")),
# 	expression(paste("4 h at 245 ",degree,"C ")),
# 	expression(paste("8 h at 245 ",degree,"C ")),
# 	expression(paste("24 h at 245 ",degree,"C ")),
# 	expression(paste("4 h at 260 ",degree,"C ")),
# 	expression(paste("8 h at 260 ",degree,"C ")),
# 	expression(paste("24 h at 260 ",degree,"C ")),
# 	expression(paste("4 h at 300 ",degree,"C ")),
# 	expression(paste("8 h at 300 ",degree,"C ")),
# 	expression(paste("24 h at 300 ",degree,"C ")),	
# 	"Fresh"), 
# 	at=1:16, las=2, cex.axis=0.9)
# 	abline(v=0, col="grey", lwd=1.5)
# 	mtext(taxon.names[i], side=4, las=2, line=0.5, cex=0.7)
# 	box()
# }
# axis(1)
# mtext(expression(paste(Delta^{13}, "C (\u2030)")), side=1, line=2.2, cex=0.7)
# mtext("A, B and C replicates from each burn", side=1, cex=0.6, line=2.8)
# mtext("normalized to the mean value of each burn", side=1, line=3.4, cex=0.6)
# mtext("Within-batch variability", side=3, line=22, cex=0.8)
# 
# for (i in 1:3){
# 	batch <- data2[data2$Species==unique(data2$S)[i], c( 12, 22)]  #selectnormd15N and CharrNo
# 	batch$Packno <- as.numeric(paste(batch$CharrNo))
# 	batch.sort <- batch[order(batch$Packno),]
# 	offsetted.value <- 0
# 	for (j in 1:16){
# 		small.batch <- batch.sort[batch.sort$Packno==j,]
# 		new.value <- small.batch$normd15N-mean(small.batch$normd15N)
# 		offsetted.value <- c(offsetted.value, new.value)
# 	}
# 	offsetted.value <- offsetted.value[2:length(offsetted.value)]
# 	plot(offsetted.value, batch.sort$Packno, xlim=c(-1.1, 1.1), axes=F, ylim=c(0, 16), ylab="", xlab="")
# 	axis(2, labels=c(expression(paste("4 h at 215 ",degree,"C ")),
# 	expression(paste("8 h at 215 ",degree,"C ")),
# 	expression(paste("24 h at 215 ",degree,"C ")),
# 	expression(paste("4 h at 230 ",degree,"C ")),
# 	expression(paste("8 h at 230 ",degree,"C ")),
# 	expression(paste("24 h at 230 ",degree,"C ")),
# 	expression(paste("4 h at 245 ",degree,"C ")),
# 	expression(paste("8 h at 245 ",degree,"C ")),
# 	expression(paste("24 h at 245 ",degree,"C ")),
# 	expression(paste("4 h at 260 ",degree,"C ")),
# 	expression(paste("8 h at 260 ",degree,"C ")),
# 	expression(paste("24 h at 260 ",degree,"C ")),
# 	expression(paste("4 h at 300 ",degree,"C ")),
# 	expression(paste("8 h at 300 ",degree,"C ")),
# 	expression(paste("24 h at 300 ",degree,"C ")),
# 		"Fresh") , at=1:16, las=2, cex.axis=0.7)
# 	abline(v=0, col="grey", lwd=1.5)
# 	mtext(taxon.names[i], side=4, las=2, line=0.5, cex=0.7)
# 	box()
# }# }
# # box()
# axis(1)
# mtext(expression(paste(Delta^{15}, "N (\u2030)")), side=1, line=2.2, cex=0.7)
# mtext("A, B and C replicates from each burn", side=1, cex=0.6, line=2.8)
# mtext("normalized to the mean value of each burn", side=1, line=3.4, cex=0.6)
# mtext("Within-batch variability", side=3, line=22, cex=0.8)
#####
##### charring %C

library(plyr)
library(nlme)
library(beeswarm)
library(lm.beta)

data2<-data

TT2 <- paste(data$Species,data2$temp, data2$time, sep= "")
data2<-data.frame(data2,TT2)


RYE<-data2[data2$ID=="RCSPARERYE24524A"|data2$ID=="RCSPARERYE24524B"|data2$ID=="RCSPARERYE24524C",]
data2<-data2[data2$TT!="rye24524",]
data2<-rbind(data2,RYE)

data2<-data2[data2$Species!="RyeX",]
data2<-data2[data2$Species!="Wheat",]
data2<-data2[data2$ID!="RYE2304C",]



data2$Species<-gsub('Free-T wheat','BW',data2$Species)

data2$Species<-gsub('HBH','Barley',data2$Species)

data2<-data2[data2$ID!="BAR0F",]
data2<-data2[data2$ID!="BAR0E",]
data2<-data2[data2$ID!="BAR0D",]

data2<- data2[data2$Species=="BW"|data2$Species=="Barley"|data2$Species=="rye"|data2$Species=="oat"|data2$Species=="Free-T wheat"|data2$Species=="spelt",]

data2$normd15N<-as.numeric(data2$normd15N)
data2$d15Nsd<-as.numeric(data2$d15Nsd)

par(oma=c(1,6,1,2))
taxon.names <- c("Barley", "BW", "Oat","Rye", "Spelt")
data.sort <- data2[order(data2$Species, data2$temp, data2$time),]
batch <- data.sort[data.sort$Species != "NB" & data.sort$Species != "HBA"& data.sort$Species != "PBA", ]
pch.list <- rep(0, length(batch$time))
pch.list[batch$time==0] <- 8
pch.list[batch$time==4] <- 21
pch.list[batch$time==8] <- 22
pch.list[batch$time==24] <- 23
col.list <- rep(0, length(batch$temp))
col.list[batch$temp==0] <- "black"
col.list[batch$temp==215] <- "white"
col.list[batch$temp==230] <- "lightgray"
col.list[batch$temp==245] <- "darkgray"
col.list[batch$temp==260] <- "black"
col.list[batch$temp==300] <- "red"

par(mfrow=c(1,6))
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[1],]
plot(batch2$pcC, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 52), main=taxon.names[1], ylim=c(30, 90))
axis(2)
box()

mtext("% C", side=2, line=2, cex=1)

batch2 <- batch[batch$Species==unique(batch$Species)[2],]
plot(batch2$pcC, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 52), main=taxon.names[2], ylim=c(30, 90))
axis(2, labels=F)
box()

batch3 <- data.sort[data.sort$Species != "Barley" & data.sort$Species != "BW", ]
pch.list <- rep(0, length(batch3$time))
pch.list[batch3$time==0] <- 8
pch.list[batch3$time==4] <- 21
pch.list[batch3$time==8] <- 22
pch.list[batch3$time==24] <- 23
col.list <- rep(0, length(batch3$temp))
col.list[batch3$temp==0] <- "black"
col.list[batch3$temp==215] <- "white"
col.list[batch3$temp==230] <- "lightgray"
col.list[batch3$temp==245] <- "darkgray"
col.list[batch3$temp==260] <- "black"
col.list[batch3$temp==300] <- "red"
for (i in 3:5){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$pcC, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 50), main=taxon.names[i], ylim=c(30, 90))
axis(2, labels=F)
box()
abline(h=0)
}
plot(batch2$pcC, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 41), ylim=c(30, 85), type="n")
legend(1, 87, c("", "215 ˚C", "230 ˚C", "245 ˚C", "260 ˚C", "300 ˚C" ,"", "Uncharred", "4h", "8h", "24h"), pch=c(NA, 22,22,22,22,22,NA, 8,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray", "black","red", NA, NA, "white", "white", "white"), bg="white", xpd=T, pt.cex=c(2, 2,2,2,2,2,1,1,1,1,1))
text(14, 86, "Temperature", xpd=T)
text(7, 78.5, "Time", xpd=T)
#####
#charring d13C
par(oma=c(1,6,1,2))
taxon.names <- c("Oat",  "Rye", "Spelt")
data.sort <- data2[order(data2$Species, data2$temp, data2$time),]
batch <- data.sort[data.sort$Species != "NB" & data.sort$Species != "HBA"& data.sort$Species != "PBA", ]
pch.list <- rep(0, length(batch$time))
pch.list[batch$time==0] <- 8
pch.list[batch$time==4] <- 21
pch.list[batch$time==8] <- 22
pch.list[batch$time==24] <- 23
col.list <- rep(0, length(batch$temp))
col.list[batch$temp==0] <- "black"
col.list[batch$temp==215] <- "white"
col.list[batch$temp==230] <- "lightgray"
col.list[batch$temp==245] <- "darkgray"
col.list[batch$temp==260] <- "gray45"
col.list[batch$temp==300] <- "black"

par(mfrow=c(1,4))
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[1],]
plot(batch2$normd13C, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 50), main=taxon.names[1], ylim=c(-29,-24))
axis(2)
box()

mtext(expression(paste(delta^{13}, "C (\u2030)")), side=2, line=2, cex=1)
for (i in 2:3){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$normd13C, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 50), main=taxon.names[i], ylim=c(-29,-24))
axis(2, labels=F)
box()
abline(h=0)
}
plot(batch2$pcC, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 41), ylim=c(38, 85), type="n")
legend(1, 70, c("", "215 ˚C", "230 ˚C", "245 ˚C", "260 ˚C", "300 ˚C" ,"", "Uncharred", "4h", "8h", "24h"), pch=c(NA, 22,22,22,22,22,NA, 8,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray","gray45", "black", NA, NA, "white", "white", "white"), bg="white", xpd=T, pt.cex=c(2, 2,2,2,2,2,1,1,1,1,1))
text(12.5, 68.5, "Temperature", xpd=T)
text(6, 62.5, "Time", xpd=T)

#charring CN_Crun
#library(plotrix)
par(oma=c(1,6,1,2))
taxon.names <- c("Oat",  "Rye", "Spelt")
data.sort <- data2[order(data2$Species, data2$temp, data2$time),]
batch <- data.sort[data.sort$Species != "NB" & data.sort$Species != "HBA"& data.sort$Species != "PBA", ]
#batch <- data.sort
pch.list <- rep(0, length(batch$time))
pch.list[batch$time==0] <- 8
pch.list[batch$time==4] <- 21
pch.list[batch$time==8] <- 22
pch.list[batch$time==24] <- 23
col.list <- rep(0, length(batch$temp))
col.list[batch$temp==0] <- "black"
col.list[batch$temp==215] <- "white"
col.list[batch$temp==230] <- "lightgray"
col.list[batch$temp==245] <- "darkgray"
col.list[batch$temp==260] <- "gray45"
col.list[batch$temp==300] <- "black"

par(mfrow=c(1,4))
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[1],]
plot(batch2$CN_Crun, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 50), main=taxon.names[1], ylim=c(15,60))
axis(2)
box()
#abline(h=0)
#abline(h=72.5)
mtext("C/N (molar)", side=2, line=2, cex=1)
for (i in 2:3){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$CN_Crun, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 50), main=taxon.names[i], ylim=c(15,60))
axis(2, labels=F)
box()
abline(h=0)
}
plot(batch2$CN_Crun, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 41), ylim=c(38, 85), type="n")
legend(1, 70, c("", "215 ˚C", "230 ˚C", "245 ˚C", "260 ˚C", "300 ˚C" ,"", "Uncharred", "4h", "8h", "24h"), pch=c(NA, 22,22,22,22,22,NA, 8,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray","gray45", "black", NA, NA, "white", "white", "white"), bg="white", xpd=T, pt.cex=c(2, 2,2,2,2,2,1,1,1,1,1))
text(11, 69, "Temperature", xpd=T)
text(6, 62, "Time", xpd=T)


#####
##################just charred
mycols1 <- c(rainbow(5))
mycols <- c(rev(mycols1))
palette(c("yellow","orange","red", "purple","blue"))
TT <- paste(batch$temp, batch$time, sep= "")
data2<-just.charred

taxon.names <- c("Barley","Wheat","Oat", "Rye")
par(mfrow=c(4,1))
par(oma=c(3,6,3,3))
par(mar=c(1,0,1,0))
batch <- data2[data2$Species==unique(data2$Species)[1],]
TT <- paste(batch$temp, batch$time, sep= "")

plot(batch$TT, batch$normd13C, axes=F, xlab="", ylim=c(min(batch$normd13C), max(batch$normd13C)), main=taxon.names[1], ylab="", cex=1.2,col=c("yellow","yellow","yellow","orange","orange","orange","red","red","red", "purple","purple","purple","blue","blue","blue"))
axis(2)
#axis(4,labels=F, line=-11.5)
mtext(expression(paste(delta^{13}, "C (\u2030)")), side=2, line=2, cex=0.8)

for (i in 2){
par(mar=c(1,0,1,0))
batch <- data2[data2$Species==unique(data2$Species)[i],]
TT <- paste(batch$temp, batch$time, sep= "")

plot(batch$TT, batch$normd13C, axes =F, xlab="", ylim=c(min(batch$normd13C), max(batch$normd13C)), main=taxon.names[2], ylab="", cex=1.2,col=c("yellow","yellow","yellow","orange","orange","orange","red","red","red", "purple","purple","purple","blue","blue","blue"))
mtext(expression(paste(delta^{13}, "C (\u2030)")), side=2, cex=0.8, line=2)
axis(2)
#axis(4,labels=F, line=-11)

}



for (i in 3){
par(mar=c(1,0,1,0))
batch <- data2[data2$Species==unique(data2$Species)[i],]
TT <- paste(batch$temp, batch$time, sep= "")
plot(batch$TT, batch$normd13C, axes =F, xlab="", ylim=c(min(batch$normd13C), max(batch$normd13C)), main=taxon.names[3], ylab="", cex=1.2, col=c("yellow","yellow","yellow","orange","orange","orange","red","red","red", "purple","purple","purple","blue","blue","blue"))
axis(2)
mtext(expression(paste(delta^{13}, "C (\u2030)")), side=2, cex=0.8, line=2)
#axis(4,labels=F, line=-11)

}
for (i in 4){
par(mar=c(1,0,1,0))
batch <- data2[data2$Species==unique(data2$Species)[i],]
TT <- paste(batch$temp, batch$time, sep= "")
batch$TT<-factor(batch$TT, levels=c("2154", "2158","21524", "2304", "2308", "23024", "2454", "2458", "24524", "2604", "2608", "26024", "3004", "3008", "30024"))
plot(batch$TT, batch$normd13C, axes =F, xlab="", ylim=c(min(batch$normd13C), max(batch$normd13C)), main=taxon.names[4], ylab="", cex=1.2,col=c("yellow","yellow","yellow","orange","orange","orange","red","red","red", "purple","purple","purple","blue","blue","blue"))

axis(1, at= c(1:15), label=c("4hrs", "8hrs","24hrs", "4hrs", "8hrs", "24hrs", "4hrs", "8hrs", "24hrs", "4hrs", "8hrs", "24hrs", "4hrs", "8hrs", "24hrs" ), cex=0.8)

axis(2)
mtext(expression(paste(delta^{13}, "C (\u2030)")), side=2, cex=0.8, line=2)
#axis(4,labels=F, line=-11)

}

palette(c("yellow","orange","red", "purple","blue"))

par(mar=c(10,2,10,0))

 plot(c(100, 200), c(0, 1000), axes=F, xlab="", ylab="", col="white")
 rect(100, 0, 120, 160, col = "blue") 
 rect(100, 200, 120, 360, col = "purple")
 rect(100, 400, 120, 560, col = "red")
 rect(100, 600, 120, 760, col = "orange")
  rect(100, 800, 120, 960, col = "yellow")
  text(140, 880, expression(paste("215",degree,"C")), cex=3)
   text(140, 680, expression(paste("230",degree,"C")), cex=3)
 text(140, 480, expression(paste("245",degree,"C")), cex=3)
   text(140, 280, expression(paste("260",degree,"C")), cex=3)
 text(140, 80, expression(paste("300",degree,"C")), cex=3)





plot(c(rep(1:1, 5), 4), c(rep(1, 3), rep(2, 3),rep(5,3), 1), pch=22, bg=c("yellow","orange","red", "purple","blue"), cex=3, axes=F, xlim=c(0.5, 5), ylim=c(0.5, 5), xlab="", ylab="")
mtext("Charring temperatures", side=3, line=2, cex=0.7)
mtext("and times", side=3, line=1, cex=0.7)
mtext(expression(paste("Temperature [",degree,"C]")), side=2, line=1.5, cex=0.7)
mtext("Time (hours)       ", side=1, line=1, cex=0.7)
axis(2, labels=c("215","230", "245", "260", "300"), at=c(1:5), las=2, lwd=0, line=-1, tick=T)
axis(1, labels=c("4", "8", "24"), at=c(1:3), lwd=0, line=-1)
axis(1, labels=c("Fresh"), at=4, las=2, tick=F, line=-1)


#############################################

########. N  versions of all the figures ###############################

library(plyr)
library(nlme)
library(beeswarm)
library(lm.beta)

data2<-data

data2$Species<-droplevels(data2$Species)

TT2 <- paste(data$Species,data2$temp, data2$time, sep= "")
data2<-data.frame(data2,TT2)

###Remove bad Rye samples 
RYE<-data2[data2$ID=="RCSPARERYE24524A"|data2$ID=="RCSPARERYE24524B"|data2$ID=="RCSPARERYE24524C",]
data2<-data2[data2$TT!="rye24524",]
data2<-rbind(data2,RYE)

data2<-data2[data2$Species!="RyeX",]
data2<-data2[data2$Species!="Wheat",]
data2<-data2[data2$ID!="RYE2304C",]

### Only using rye, oat, BW and HB
data2<- data2[data2$Species=="BW"|data2$Species=="HBH"|data2$Species=="rye"|data2$Species=="oat"|data2$Species=="Free-T wheat"|data2$Species=="Barley",]

data2$Species<-gsub('Free-T wheat','BW',data2$Species)

data2$Species<-gsub('HBH','Barley',data2$Species)

data2$normd15N<-as.numeric(data2$normd15N)
data2$d15Nsd<-as.numeric(data2$d15Nsd)

data2<-data2[data2$ID!="BAR0F",]
data2<-data2[data2$ID!="BAR0E",]
data2<-data2[data2$ID!="BAR0D",]
data3<- data2[data2$Species=="BW"|data2$Species=="HBH"|data2$Species=="rye"|data2$Species=="oat"|data2$Species=="Free-T wheat"|data2$Species=="Barley"|data2$Species=="spelt",]

####Table 6. lm models to calculate charring offset 
## 215-300 degress C
lm1 <- lm(normd15N28thoct ~ Species, data = data2)
summary(lm1)
lm2<-lm(normd15N28thoct~char+Species,data=data2)
summary(lm2)
confint(lm2)#### and this one for CI

### 215 to 260 degrees C
no300data2<-data2[data2$temp!="300",]
lm1 <- lm(normd15N28thoct ~ Species, data = no300data2)
summary(lm1)
lm2 <- lm(normd15N28thoct ~ char + Species, data=no300data2)### data in paper
summary(lm2)
confint(lm2)#### and this one for CI

### 230-300 degrees C
no215data2<-data2[data2$temp!="215",]
lm1 <- lm(normd15N28thoct ~ Species, data = no215data2)
summary(lm1)
lm2 <- lm(normd15N28thoct ~ char + Species, data=no215data2)### data in paper
summary(lm2)
confint(lm2)#### and this one for CI

### 230-260 degrees C

no300215data2<-no300data2[no300data2$temp!="215",]
lm1 <- lm(normd15N28thoct ~ Species, data = no300215data2)
summary(lm1)
lm2 <- lm(normd15N28thoct ~ char + Species, data=no300215data2)### data in paper
summary(lm2)
confint(lm2)#### and this one for CI

######
##############################

data2$normd15N28thoct<-as.numeric(data2$normd15N28thoct)
summary <- data.frame(ddply(data2, c("Species", "temp", "time"), nrow))
summary2 <- ddply(data2, c("Species", "temp", "time"), function(x) c(d15N=mean(x$normd15N28thoct), sd=sd(x$normd15N28thoct), pcN=mean(x$pcN), d13C=mean(x$normd13C), sd=sd(x$normd13C), pcN=mean(x$pcC)))

Ncharoff <- rep(0, nrow(data2))
for (i in 1:length(unique(data2$Species))){
offset.start <- summary2$d15N[summary2$Species==unique(data2$Species)[i] & summary2$time==0]
offsetted.values <- (data2$normd15N[data2$Species==unique(data2$Species)[i]])-offset.start
Ncharoff[data2$Species==unique(data2$Species)[i]] <- offsetted.values}

data3 <- data.frame(data2, Ncharoff)

summary3 <- ddply(data3, c("Species", "temp", "time"), function(x) c( d15N=mean(x$normd15N28thoct), sd=sd(x$normd15N28thoct), pcN=mean(x$pcN), CN=mean(x$CN_Crun), charoff=mean(x$Ncharoff)))
#Figure 2
png(filename="FigureNcharoff.png",
    units="cm",
    width=14, 
    height=13, 
    res=300,
    pointsize=8)

par(oma=c(5,6,5,2))
taxon.names <- c("Barley","BW")
data.sort <- data3[order(data3$Species, data3$temp, data3$time),]
batch <- data.sort[data.sort$Species != "NB" & data.sort$Species != "HBA"& data.sort$Species != "PBA", ]
summary4 <- ddply(data3, c("Species"), function(x) c( d15N=mean(x$normd15N), sd=sd(x$normd15N), pcN=mean(x$pcN), CN=mean(x$CN_Crun), charoff=mean(x$Ncharoff)))

pch.list <- rep(0, length(batch$time))
pch.list[batch$time==0] <- 8
pch.list[batch$time==4] <- 21
pch.list[batch$time==8] <- 22
pch.list[batch$time==24] <- 23
col.list <- rep(0, length(batch$temp))
col.list[batch$temp==0] <- "black"
col.list[batch$temp==215] <- "white"
col.list[batch$temp==230] <- "lightgray"
col.list[batch$temp==245] <- "darkgray"
col.list[batch$temp==260] <- "black"
col.list[batch$temp==300]<-"red"
par(mfrow=c(1,5))
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[1],]
plot(batch2$Ncharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-1.5,2.7), xlim=c(-1, 50), main=taxon.names[1])
#legend("topleft", c( "Uncharred","215", "230", "245", "260",  "4h", "8h", "24h"), pch=c(8, 22,22,22,22,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray", "black", "white", "white", "white"), bg="white", pt.cex=c(1.2,2, 2, 2, 2, 1.2,1.2,1.2), cex=1.2, box.lty=0)
axis(2)
box()
abline(h=0)
#abline(h=0.63)
mtext(expression(paste(Delta^{15}, "N (\u2030) uncharred-charred")), side=2, line=2, cex=0.8)

for (i in 2){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$Ncharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-1.5, 2.7), xlim=c(-1, 50), main=taxon.names[i])
#axis(2)
box()
abline(h=0)
}


taxon.names <- c( "Oat","Rye")
batch <- data.sort[data.sort$Species != "Barley" & data.sort$Species != "BW", ]
pch.list <- rep(0, length(batch$time))
pch.list[batch$time==0] <- 8
pch.list[batch$time==4] <- 21
pch.list[batch$time==8] <- 22
pch.list[batch$time==24] <- 23
col.list <- rep(0, length(batch$temp))
col.list[batch$temp==0] <- "black"
col.list[batch$temp==215] <- "white"
col.list[batch$temp==230] <- "lightgray"
col.list[batch$temp==245] <- "darkgray"
col.list[batch$temp==260] <- "black"
col.list[batch$temp==300]<-"red"


for (i in 1){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$Ncharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-1.5, 2.7), xlim=c(-1, 50), main=taxon.names[i])
#axis(2)
box()
abline(h=0)
}

for (i in 2){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$Ncharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-1.5, 2.7), xlim=c(-1, 50), main=taxon.names[i])
#axis(2)
box()
abline(h=0)
}


plot(batch2$pcN, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 45), ylim=c(38, 85), type="n")
legend(1, 70, c("", "215 ˚C", "230 ˚C", "245 ˚C", "260 ˚C", "300 ˚C","", "Uncharred", "4h", "8h", "24h"), pch=c(NA, 22,22,22,22,22,NA, 8,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray", "black", "red",NA, NA, "white", "white", "white"), bg="white", xpd=T, pt.cex=c(2, 2,2,2,2,2,2,1,1,1,1), cex=1)
text(14, 68.7, "Temperature", xpd=T)
text(8, 62, "Time", xpd=T)


dev.off()

#charring %N

library(plyr)
library(nlme)
library(beeswarm)
library(lm.beta)

data2<-data

data2$Species<-droplevels(data2$Species)

TT2 <- paste(data$Species,data2$temp, data2$time, sep= "")
data2<-data.frame(data2,TT2)

###Remove bad Rye samples 
RYE<-data2[data2$ID=="RCSPARERYE24524A"|data2$ID=="RCSPARERYE24524B"|data2$ID=="RCSPARERYE24524C",]
data2<-data2[data2$TT!="rye24524",]
data2<-rbind(data2,RYE)

data2<-data2[data2$Species!="RyeX",]
data2<-data2[data2$Species!="Wheat",]
data2<-data2[data2$ID!="RYE2304C",]

### Only using rye, oat, BW and HB

data2$Species<-gsub('Free-T wheat','BW',data2$Species)

data2$Species<-gsub('HBH','Barley',data2$Species)

data2$normd15N<-as.numeric(data2$normd15N)
data2$d15Nsd<-as.numeric(data2$d15Nsd)

data2<-data2[data2$ID!="BAR0F",]
data2<-data2[data2$ID!="BAR0E",]
data2<-data2[data2$ID!="BAR0D",]
data3<- data2[data2$Species=="BW"|data2$Species=="HBH"|data2$Species=="rye"|data2$Species=="oat"|data2$Species=="Free-T wheat"|data2$Species=="Barley"|data2$Species=="spelt",]

par(oma=c(1,6,1,2))
taxon.names <- c("Barley", "BW","Oat","Rye", "Spelt")
data.sort <- data3[order(data3$Species, data3$temp, data3$time),]
batch <- data.sort[data.sort$Species != "NB" & data.sort$Species != "HBA"& data.sort$Species != "PBA", ]
#batch <- data.sort
pch.list <- rep(0, length(batch$time))
pch.list[batch$time==0] <- 8
pch.list[batch$time==4] <- 21
pch.list[batch$time==8] <- 22
pch.list[batch$time==24] <- 23
col.list <- rep(0, length(batch$temp))
col.list[batch$temp==0] <- "black"
col.list[batch$temp==215] <- "white"
col.list[batch$temp==230] <- "lightgray"
col.list[batch$temp==245] <- "darkgray"
col.list[batch$temp==260] <- "black"
col.list[batch$temp==300] <- "red"

par(mfrow=c(1,6))
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[1],]
plot(batch2$pcN, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 55), main=taxon.names[1], ylim=c(0, 6))
axis(2)
box()
mtext("% N", side=2, line=2, cex=1)
batch2 <- batch[batch$Species==unique(batch$Species)[2],]
plot(batch2$pcN, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 55), main=taxon.names[2], ylim=c(0, 6))
box()

batch3 <- data.sort[data.sort$Species != "Barley" & data.sort$Species != "BW", ]
pch.list <- rep(0, length(batch3$time))
pch.list[batch3$time==0] <- 8
pch.list[batch3$time==4] <- 21
pch.list[batch3$time==8] <- 22
pch.list[batch3$time==24] <- 23
col.list <- rep(0, length(batch3$temp))
col.list[batch3$temp==0] <- "black"
col.list[batch3$temp==215] <- "white"
col.list[batch3$temp==230] <- "lightgray"
col.list[batch3$temp==245] <- "darkgray"
col.list[batch3$temp==260] <- "black"
col.list[batch3$temp==300] <- "red"

for (i in 3:5){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$pcN, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 50), main=taxon.names[i], ylim=c(0, 6))
axis(2, labels=F)
box()

}
plot(batch2$pcN, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 41), ylim=c(0, 6), type="n")
legend(1, 6, c("", "215 ˚C", "230 ˚C", "245 ˚C", "260 ˚C", "300 ˚C" ,"", "Uncharred", "4h", "8h", "24h"), pch=c(NA, 22,22,22,22,22,NA, 8,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray", "black","red", NA, NA, "white", "white", "white"), bg="white", xpd=T, pt.cex=c(2, 2,2,2,2,2,1,1,1,1,1))
text(14, 5.9, "Temperature", xpd=T)
text(6, 5.1, "Time", xpd=T)
# 
# 
# #charring new CN ration
# #library(plotrix)
# par(oma=c(1,6,1,2))
# taxon.names <- c("Barley", "BW","Oat", "Rye","Spelt")
# data.sort <- data2[order(data2$Species, data2$temp, data2$time),]
# batch <- data.sort[data.sort$Species != "NB" & data.sort$Species != "HBA"& data.sort$Species != "PBA", ]
# #batch <- data.sort
# pch.list <- rep(0, length(batch$time))
# pch.list[batch$time==0] <- 8
# pch.list[batch$time==4] <- 21
# pch.list[batch$time==8] <- 22
# pch.list[batch$time==24] <- 23
# col.list <- rep(0, length(batch$temp))
# col.list[batch$temp==0] <- "black"
# col.list[batch$temp==215] <- "white"
# col.list[batch$temp==230] <- "lightgray"
# col.list[batch$temp==245] <- "darkgray"
# col.list[batch$temp==260] <- "gray45"
# col.list[batch$temp==300] <- "black"
# 
# par(mfrow=c(1,6))
# par(mar=c(3,0,1,0))
# batch2 <- batch[batch$Species==unique(batch$Species)[1],]
# plot(batch2$newCN, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 50), main=taxon.names[1], ylim=c(5,55))
# axis(2)
# box()
# #abline(h=0)
# #abline(h=72.5)
# mtext("C/N (molar)", side=2, line=2, cex=1)
# for (i in 2:3){
# par(mar=c(3,0,1,0))
# batch2 <- batch[batch$Species==unique(batch$Species)[i],]
# plot(batch2$newCN, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 50), main=taxon.names[i], ylim=c(5,55))
# axis(2, labels=F)
# box()
# abline(h=0)
# }
# plot(batch2$CN_Crun, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, xlim=c(-1, 41), ylim=c(38, 85), type="n")
# legend(1, 70, c("", "215 ˚C", "230 ˚C", "245 ˚C", "260 ˚C", "300 ˚C" ,"", "Uncharred", "4h", "8h", "24h"), pch=c(NA, 22,22,22,22,22,NA, 8,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray","gray45", "black", NA, NA, "white", "white", "white"), bg="white", xpd=T, pt.cex=c(2, 2,2,2,2,2,1,1,1,1,1))
# text(11, 69, "Temperature", xpd=T)
# text(6, 62, "Time", xpd=T

#### Mass loss graphs

library(plyr)
library(nlme)
library(beeswarm)
library(lm.beta)

par(oma=c(0,0,0,0))
par(mar=c(5,5,5,10))
data.sort <- pcwloss[order(pcwloss$Species, pcwloss$temp, pcwloss$time),]
batch<-data.sort
colnames(batch)[1]<-"ID"


pch.list <- rep(0, length(batch$time))
pch.list[batch$time==0] <- 8
pch.list[batch$time==4] <- 21
pch.list[batch$time==8] <- 22
pch.list[batch$time==24] <- 23
col.list <- rep(0, length(batch$temp))
col.list[batch$temp==0] <- "black"
col.list[batch$temp==215] <- "white"
col.list[batch$temp==230] <- "lightgray"
col.list[batch$temp==245] <- "darkgray"
col.list[batch$temp==260] <- "black"
col.list[batch$temp==300] <- "red"

plot(batch$pcwtloss, bg=c(col.list), pch=c(pch.list), ylab="% mass loss", axes=F, xlab="", cex=1.2)
axis(1, at=c(1.5, 17, 16.5*2, 16.5*3, 16.5*4), labels=c("Barley", "BW", "Oat", "Rye", "Spelt"), cex.axis=0.85, crt=80)
box()
axis(2)

legend(85.5, 70.7, c("", "215 ˚C", "230 ˚C", "245 ˚C", "260 ˚C", "300 ˚C","", "Uncharred", "4h", "8h", "24h"), pch=c(NA, 22,22,22,22,22,NA, 8,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray", "black", "red", NA, NA, "white", "white", "white"), bg="white", xpd=T, pt.cex=c(2, 2,2,2,2,2,2,1,1,1,1,1))
text(91.5, 68, "Temperature", xpd=T)
text(88, 52, "Time", xpd=T)


