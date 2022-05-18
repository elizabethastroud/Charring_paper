library(plyr)
library(nlme)
library(beeswarm)
library(lm.beta)
data<-read.csv("Charringdataclean15.csv")
setwd("~/OneDrive - Nexus365/Feedsax/Charring results/Charring experiment/Experimental_charring/")
data2<-data

TT2 <- paste(data2$Species,data2$temp, data2$time, sep= "")
data2<-data.frame(data2,TT2)

###work out if barley and wheat can be compared between Nitsch et all and new data

bar_whe<-data2[data2$TT2=="Barley00"|data2$TT2=="BW00",]
bar_whe$normd15N<-as.numeric(bar_whe$normd15N)

bar<-bar_whe[bar_whe$TT2=="Barley00",]
whe<-bar_whe[bar_whe$TT2=="BW00",]

boxplot(whe$normd13C~whe$CharrNo)
t.test(whe$normd13C~whe$CharrNo) #p = 0.227
boxplot(whe$normd15N~whe$CharrNo)
t.test(whe$normd15N~whe$CharrNo) #p =0.9178

boxplot(bar$normd13C~bar$CharrNo)
t.test(bar$normd13C~bar$CharrNo)# p=0.2965
boxplot(bar$normd15N~bar$CharrNo)
t.test(bar$normd15N~bar$CharrNo) # p =0.2848

###### Kragten spreadsheet values - all data

mean(data2$d15Nsd, na.rm= TRUE)
mean(data2$d13Csd, na.rm= TRUE)

### Only using rye, oat, BW and HB
data3<- data2[data2$Species=="BW"|data2$Species=="HBH"|data2$Species=="rye"|data2$Species=="oat"|data2$Species=="Free-T wheat"|data2$Species=="Barley",]
data3$normd15N<-as.numeric(data3$normd15N)
data3$d15Nsd<-as.numeric(data3$d15Nsd)
###### Kragten spreadsheet values - just wheat, barely, rye and oat
mean(data3$d15Nsd, na.rm= TRUE)
mean(data3$d13Csd, na.rm= TRUE)
summary <- data.frame(ddply(data3, c("Species", "temp", "time"), nrow))

#Szpak values carbon
library(dplyr)
library(multiway)
RawStandards<-read.csv("RawStandards.csv")
RepCar<-read.csv("RepCar.csv")
all.standards<- RawStandards %>% 
  group_by(RunfileC,ID)%>%
  summarise(Number=n(), d13Cmean=mean(normd13C),  d13Csd=sd(normd13C)) %>%
  as.data.frame()

all.standards$srm<-(all.standards$Number-1)*(all.standards$d13Csd^2)
dfsrm<-sum(all.standards$Number)-nrow(all.standards)
Ssrm<-sqrt(sum(all.standards$srm)/dfsrm)

checkC<-subset(all.standards, all.standards$Standard=="P2"|all.standards$Standard=="SALANINE")

CheckS.1<--28.19#### P2
CheckS.2<--27.11##### Alanine
CheckS.1sd<-0.14#### P2
Check.2sd<-0.124

y="P2" #### change if using different check standards 
fun1<-function(x,y) if(x==y) {CheckS.1} else {CheckS.2}
checkC$known<-mapply(fun1, checkC$Standard, y)

y="P2"
fun1<-function(x,y) if(x==y) {CheckS.1sd} else {CheckS.2sd}
checkC$knownsd<-mapply(fun1, checkC$Standard, y)

checkC$Diff_measured_known<-checkC$d13Cmean-checkC$known

RMSbias<-sqrt(sumsq(checkC$Diff_measured_known)/nrow(checkC))
u_cref<-sqrt(sumsq(checkC$knownsd)/nrow(checkC))

x<-list(RMSbias, u_cref)
u_bias<-sqrt(sumsq(x))


RepCar$Sd<-apply(subset(RepCar,select = c("normd13C_DulpA","normd13C_DulpB")),1,sd)
RepCar$Mean<-apply(subset(RepCar,select = c("normd13C_DulpA","normd13C_DulpB")),1,mean)
RepCar$Number<-2
RepCar$RepSsrm<-1*(RepCar$Sd^2)
dfrep<-sum(RepCar$Number)-nrow(RepCar)
Srep<-sqrt((sum(RepCar$RepSsrm))/dfrep)

uRw<-sqrt((Ssrm^2)+(Srep^2)/2)
y<-list(u_bias, uRw)
Uc<-sqrt(sumsq(y))

#Szpak values Nitrogen
library(dplyr)
library(multiway)
RawStandardsN<-read.csv("RawStandardsN.csv")
RepNit<-read.csv("RepNit.csv")
all.standardsN<- RawStandardsN %>% 
  group_by(RunfileN,ID)%>%
  summarise(Number=n(), d15Nmean=mean(normd15N),  d15Nsd=sd(normd15N)) %>%
  as.data.frame()

all.standardsN$srm<-(all.standardsN$Number-1)*(all.standardsN$d15Nsd^2)
dfsrmN<-sum(all.standardsN$Number)-nrow(all.standardsN)
SsrmN<-sqrt(sum(all.standardsN$srm)/dfsrmN)


checkN<-subset(all.standardsN, all.standardsN$ID=="P2"|all.standardsN$ID=="SALANINE"|all.standardsN$ID=="LEU")
checkN<-checkN[-18,]# remove P2 from 200827 as used as calibration standard
CheckS.1<--1.57#### P2
CheckS.2<--1.57##### Alanine
CheckS.3<- 6.36###Leucine

y="LEU"
fun1<-function(x,y) if(x==y) {CheckS.3} else {CheckS.2}         
checkN$known<-mapply(fun1, checkN$ID,  y)

CheckS.1sd<-0.19#### P2
CheckS.2sd<-0.21
CheckS.3sd<-0.14

y="P2"

fun1<-function(x,y) if(x==y) {CheckS.1sd} else {CheckS.2sd}

checkN$knownsd<-mapply(fun1, checkN$ID, y)

checkN$knownsd[checkN$ID=="LEU"]<-CheckS.3sd

checkN$Diff_measured_known<-checkN$d15Nmean-checkN$known

RMSbias<-sqrt(sumsq(checkN$Diff_measured_known)/nrow(checkN))
u_cref<-sqrt(sumsq(checkN$knownsd)/nrow(checkN))

x<-list(RMSbias, u_cref)
u_bias<-sqrt(sumsq(x))

RepNit$Sd<-apply(subset(RepNit,select = c("normd15N_DulpA","normd15N_DulpB")),1,sd)

RepNit$Mean<-apply(subset(RepNit,select = c("normd15N_DulpA","normd15N_DulpB")),1,mean)

RepNit$Number<-2

RepNit$RepSsrm<-1*(RepNit$Sd^2)
dfrep<-sum(RepNit$Number)-nrow(RepNit)
Srep<-sqrt((sum(RepNit$RepSsrm))/dfrep)

uRw<-sqrt((SsrmN^2)+(Srep^2)/2)
y<-list(u_bias, uRw)
Uc<-sqrt(sumsq(y))


###########################################################
#      C versions of all the figures and statistics
###########################################################

 summarycharred<- ddply(data3, c("Species","char"), function(x) c(d15N=mean(x$normd15N), sd=sd(x$normd15N), max=max(x$normd15N), min=min(x$normd15N), d13C=mean(x$normd13C), sd=sd(x$normd13C), max=max(x$normd13C), min=min(x$normd13C)))
# range -27.9 to -24.9 uncharred and -28.3 to -24.6 charred

######Compare LM for the just charred material for each taxon: reported in table 3


just.charred <- data3[data3$char=="charred",]
TT <- paste(just.charred$temp, just.charred$time, sep= "")
just.charred <- data.frame(just.charred, TT)

# 215- 300 degrees
charred.lm2 <- lm(normd13C ~ Species + temp + time, data=just.charred) ######### used in paper - significant p value for temp
summary(charred.lm2)
#Temp est = 2.903e-03, p = 7.96e-05
#Time est = 2.533e-03, p = 0.298

# 215-260
no300<-just.charred[just.charred$temp!="300",]
no300.lm2 <- lm(normd13C ~ Species + temp + time, data=no300) 
summary(no300.lm2)

#Temp est = 0.003, p = 0.0363
#Time est = 0.003, p = 0.2835

#230-300
no215<-just.charred[just.charred$temp!="215",]
no215.lm2 <- lm(normd13C ~ Species + temp + time, data=no215) 
summary(no215.lm2)
#Temp est = 1.583e-03, p = 0.056
#Time est =2.268e-03, p = 0.362

# 230-260
no215300<-no215[no215$temp!="300",]
no215300.lm2 <- lm(normd13C ~ Species + temp + time, data=no215300)
summary(no215300.lm2)

#Temp est = -0.002242, p = 0.271
#time est = 0.002822, p =0.328

## Compare random slopes and intercepts models to just multiple regression These are the values reported in table 5
###All temperatures 215-300
lm1 <- lm(normd13C ~ Species, data = data3)
summary(lm1)
#adj R2=0.8667
#p value <2.2e-16

lm2 <- lm(normd13C ~ char + Species, data=data3)#### and this one for CI - this was are saying predict d13c from charring and species - so species has a direct effect 
summary(lm2)
#adj R2 0.8678
#p value <2.2e-16
#charred-fresh p = 0.106
# est(beta) = -0.11705
confint(lm2)
#charred fresh -0.259, 0.025059


### without 300 data
no300data3<-data3[data3$temp!="300",]
lm1 <- lm(normd13C ~ Species, data = no300data3)
summary(lm1)
#adj R2 0.8716
#p value <2.2e-16

lm2 <- lm(normd13C ~ char + Species, data=no300data3)### data in X paper
summary(lm2)
#adj R2 0.8719
#p value <2.2e-16
#charred-fresh p =0.26 
# est(beta) = -0.0825
confint(lm2)#### and this one for CI
#charred fresh -0.22656, 0.0616

### without 215 data
no215data3<-data3[data3$temp!="215",]
lm1 <- lm(normd13C ~ Species, data = no215data3)
summary(lm1)
#adj R2 0.8869
#p value <2.2e-16

lm2 <- lm(normd13C ~ char + Species, data=no215data3)### data in paper
summary(lm2)
#adj R2 0.8902
#p value <2.2e-16
#charred-fresh p =0.0171
# est(beta) = -0.157
confint(lm2)#### and this one for CI
#charred- fresh -0.2858, -0.02831935
### without 215 and 300 data
no215300data3<-no215data3[no215data3$temp!="300",]
lm1 <- lm(normd13C ~ Species, data = no215300data3)
summary(lm1)
#adj R2 0.8931
#p value <2.2e-16

lm2 <- lm(normd13C ~ char + Species, data=no215300data3)### data in paper
summary(lm2)
#adj R2 0.8953
#p value <2.2e-16
#charred-fresh p = 0.0621
# est(beta) = -0.12451
confint(lm2)
#charred- fresh -0.255398, 0.0063819
##############################

summary <- data.frame(ddply(data3, c("Species", "temp", "time"), nrow))
summary2 <- ddply(data3, c("Species", "temp", "time"), function(x) c(d15N=mean(x$normd15N), sd=sd(x$normd15N), pcN=mean(x$pcN), d13C=mean(x$normd13C), sd=sd(x$normd13C), pcN=mean(x$pcC)))

#######
Ccharoff <- rep(0, nrow(data3))
for (i in 1:length(unique(data3$Species))){
offset.start <- summary2$d13C[summary2$Species==unique(data3$Species)[i] & summary2$time==0]
offsetted.values <- (data3$normd13C[data3$Species==unique(data3$Species)[i]])-offset.start
Ccharoff[data3$Species==unique(data3$Species)[i]] <- offsetted.values}
datan <- data.frame(data3, Ccharoff)


summary3 <- ddply(datan, c("Species", "temp", "time"), function(x) c( d13C=mean(x$normd13C), sd=sd(x$normd13C), pcC=mean(x$pcC), CN=mean(x$CN_Crun), charoff=mean(x$Ccharoff)))

####Figure 1 in paper
png(filename="FigureC1.png",
    units="cm",
    width=14, 
    height=13, 
    res=300,
    pointsize=8)

par(oma=c(4,6,4,2))
taxon.names <- c("Barley",  "BW", "Oat","Rye")
data.sort <- datan[order(datan$Species, datan$temp, datan$time),]
batch <- data.sort[data.sort$Species != "NB" & data.sort$Species != "HBA"& data.sort$Species != "PBA", ]
summary4 <- ddply(datan, c("Species"), function(x) c( d13C=mean(x$normd13C), sd=sd(x$normd13C), pcC=mean(x$pcC), charoff=mean(x$Ccharoff)))

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
##### charring %C

library(plyr)
library(nlme)
library(beeswarm)
library(lm.beta)
data<-read.csv("~/OneDrive - Nexus365/Feedsax/Charring results/Charring experiment/Experimental_charring/Charringdataclean15.csv")

TT2 <- paste(data$Species,data$temp, data$time, sep= "")
data2<-data.frame(data,TT2)


data2<- data2[data2$Species=="BW"|data2$Species=="Barley"|data2$Species=="rye"|data2$Species=="oat"|data2$Species=="Free-T wheat"|data2$Species=="spelt",]

data2$normd15N<-as.numeric(data2$normd15N)
data2$d15Nsd<-as.numeric(data2$d15Nsd)
png(filename="FigureCpercent.png",
    units="cm",
    width=14, 
    height=13, 
    res=300,
    pointsize=8)

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
text(7, 79, "Time", xpd=T)
dev.off()

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


###########################################################
#      N  versions of all the figures and statistics
###########################################################
library(plyr)
library(nlme)
library(beeswarm)
library(lm.beta)

# data<-read_csv("Charringdataclean15.csv")

TT2 <- paste(data$Species,data$temp, data$time, sep= "")
data2<-data.frame(data,TT2)

### Only using rye, oat, BW and HB
data3<- data2[data2$Species=="BW"|data2$Species=="HBH"|data2$Species=="rye"|data2$Species=="oat"|data2$Species=="Free-T wheat"|data2$Species=="Barley",]
data3$normd15N<-as.numeric(data3$normd15N)
data3$d15Nsd<-as.numeric(data3$d15Nsd)


######Compare LM for the just charred material for each taxon: reported in table5
#####

just.charred <- data3[data3$char=="charred",]
TT <- paste(just.charred$temp, just.charred$time, sep= "")
just.charred <- data.frame(just.charred, TT)

# 215- 300 degrees
charred.lm2 <- lm(normd15N ~ Species + temp + time, data=just.charred) ######### used in paper - significant p value for temp
summary(charred.lm2)

#temp, est = -0.003747, p value = 0.019514*
#time, est = 0.013651, p value = 0.011778*

# 215-260
no300<-just.charred[just.charred$temp!="300",]
no300.lm2 <- lm(normd15N ~ Species + temp + time, data=no300) 
summary(no300.lm2)
#temp, est = -0.00997, p value = 0.00152*
#time, est = 0.01408, p value = 0.0200*


#230-300
no215<-just.charred[just.charred$temp!="215",]
no215.lm2 <- lm(normd15N ~ Species + temp + time, data=no215) 
summary(no215.lm2)

#temp, est =-0.003119 , p value = 0.10878
#time, est = 0.013017, p value = 0.02719*

# 230-260
no215300<-no215[no215$temp!="300",]
no215300.lm2 <- lm(normd15N ~ Species + temp + time, data=no215300)
summary(no215300.lm2)
#temp, est = -0.014838, p value = 0.002590*
#time, est = 0.013379, p value = 0.052131.


####Table 6. lm models to calculate charring offset 
## 215-300 degrees C
lm1 <- lm(normd15N ~ Species, data = data3)
summary(lm1)
#adj R2 0.8331
# p value <2.2e-16
lm2<-lm(normd15N~char+Species, data= data3)
summary(lm2)
#adj R2 0.8359
# p value <2.2e-16
#est (beta) -0.32549
# p value charred-fresh 0.0403*
confint(lm2)#### and this one for CI
# -0.6363860, -0.014589

### 215 to 260 degrees C
no300data3<-data3[data3$temp!="300",]
lm1 <- lm(normd15N ~ Species, data = no300data3)
summary(lm1)
#adj R2 0.8508
# p value <2.2e-16

lm2 <- lm(normd15N ~ char + Species, data=no300data3)### data in paper
summary(lm2)
#adj R2 0.8436
# p value <2.2e-16
#est (beta) -0.3198
# p value charred-fresh 0.051.
confint(lm2)#### and this one for CI
#-0.6410482, 0.00137829

### 230-300 degrees C
no215data3<-data3[data3$temp!="215",]
lm1 <- lm(normd15N ~ Species, data = no215data3)
summary(lm1)
#adj R2 0.8245
# p value <2.2e-16

lm2 <- lm(normd15N ~ char + Species, data=no215data3)### data in paper
summary(lm2)
#adj R2 0.828
# p value <2.2e-16
# est (beta) -0.3151
# p value = 0.041
confint(lm2)#### and this one for CI
# -0.6172327 to -0.01299759
### 230-260 degrees C

no300215data3<-no300data3[no300data3$temp!="215",]
lm1 <- lm(normd15N ~ Species, data = no300215data3)
summary(lm1)
#adj R2 0.8268
# p value <2.2e-16

lm2 <- lm(normd15N ~ char + Species, data=no300215data3)### data in paper
summary(lm2)
#adj R2 0.83
# p value <2.2e-16
# est (beta) -0.3038
# p value = 0.0629.

confint(lm2)#### and this one for CI

#-0.6241998 to 0.016618

######Residual standard error table 7

dataFC<- data3
TT <- paste(dataFC$temp, dataFC$time, sep= "")
dataFC <- data.frame(dataFC, TT)
taxon.1<-dataFC[dataFC$Species==unique(dataFC$Species)[1],]
lmBAR<-lm(normd13C~char, data=taxon.1)
summary(lmBAR)
#0.2468
taxon.2<-dataFC[dataFC$Species==unique(dataFC$Species)[2],]
lmFT<-lm(normd13C~char, data=taxon.2)
summary(lmFT)
#0.185

taxon.3<-dataFC[dataFC$Species==unique(dataFC$Species)[3],]
lmOat<-lm(normd13C~char, data=taxon.3)
summary(lmOat)
#0.378

taxon.4<-dataFC[dataFC$Species==unique(dataFC$Species)[4],]
lmRye<-lm(normd13C~char, data=taxon.4)
summary(lmRye)
#0.3293
lmBARN<-lm(normd15N~char, data=taxon.1)
summary(lmBARN)
#0.4791
lmFTN<-lm(normd15N~char, data=taxon.2)
summary(lmFTN)
#0.4281

lmOatN<-lm(normd15N~char, data=taxon.3)
summary(lmOatN)
#0.5293
lmRyeN<-lm(normd15N~char, data=taxon.4)
summary(lmRyeN)
#0.9468

##############################

summary <- data.frame(ddply(data3, c("Species", "temp", "time"), nrow))
summary2 <- ddply(data3, c("Species", "temp", "time"), function(x) c(d15N=mean(x$normd15N), sd=sd(x$normd15N), pcN=mean(x$pcN), d13C=mean(x$normd13C), sd=sd(x$normd13C), pcN=mean(x$pcC)))

Ncharoff <- rep(0, nrow(data3))
for (i in 1:length(unique(data3$Species))){
offset.start <- summary2$d15N[summary2$Species==unique(data3$Species)[i] & summary2$time==0]
offsetted.values <- (data3$normd15N[data3$Species==unique(data3$Species)[i]])-offset.start
Ncharoff[data3$Species==unique(data3$Species)[i]] <- offsetted.values}

data4 <- data.frame(data3, Ncharoff)

summary3 <- ddply(data4, c("Species", "temp", "time"), function(x) c( d15N=mean(x$normd15N), sd=sd(x$normd15N), pcN=mean(x$pcN), CN=mean(x$CN_Crun), charoff=mean(x$Ncharoff)))
#Figure 2
png(filename="FigureNcharoffcleannew.png",
    units="cm",
    width=14, 
    height=13, 
    res=300,
    pointsize=8)

par(oma=c(5,6,5,2))
taxon.names <- c("Barley","BW")
data.sort <- data4[order(data4$Species, data4$temp, data4$time),]
batch <- data.sort[data.sort$Species != "NB" & data.sort$Species != "HBA"& data.sort$Species != "PBA", ]
summary4 <- ddply(data4, c("Species"), function(x) c( d15N=mean(x$normd15N), sd=sd(x$normd15N), pcN=mean(x$pcN), CN=mean(x$CN_Crun), charoff=mean(x$Ncharoff)))

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
plot(batch2$Ncharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-2,2.7), xlim=c(-1, 50), main=taxon.names[1])
#legend("topleft", c( "Uncharred","215", "230", "245", "260",  "4h", "8h", "24h"), pch=c(8, 22,22,22,22,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray", "black", "white", "white", "white"), bg="white", pt.cex=c(1.2,2, 2, 2, 2, 1.2,1.2,1.2), cex=1.2, box.lty=0)
axis(2)
box()
abline(h=0)
#abline(h=0.63)
mtext(expression(paste(Delta^{15}, "N (\u2030) uncharred-charred")), side=2, line=2, cex=0.8)

for (i in 2){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$Ncharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-2, 2.7), xlim=c(-1, 50), main=taxon.names[i])
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
plot(batch2$Ncharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-2, 2.7), xlim=c(-1, 50), main=taxon.names[i])
#axis(2)
box()
abline(h=0)
}

for (i in 2){
par(mar=c(3,0,1,0))
batch2 <- batch[batch$Species==unique(batch$Species)[i],]
plot(batch2$Ncharoff, bg=c(col.list), pch=c(pch.list), ylab="", axes=F, xlab="", cex=1.5, ylim=c(-2, 2.7), xlim=c(-1, 50), main=taxon.names[i])
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

### Only using rye, oat, BW and HB

data2$Species<-gsub('Free-T wheat','BW',data2$Species)

data2$Species<-gsub('HBH','Barley',data2$Species)

data2$normd15N<-as.numeric(data2$normd15N)
data2$d15Nsd<-as.numeric(data2$d15Nsd)

data3<- data2[data2$Species=="BW"|data2$Species=="HBH"|data2$Species=="rye"|data2$Species=="oat"|data2$Species=="Free-T wheat"|data2$Species=="Barley"|data2$Species=="spelt",]
png(filename="FigureNpercent.png",
    units="cm",
    width=14, 
    height=13, 
    res=300,
    pointsize=8)

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
legend(1, 6.2, c("", "215 ˚C", "230 ˚C", "245 ˚C", "260 ˚C", "300 ˚C" ,"", "Uncharred", "4h", "8h", "24h"), pch=c(NA, 22,22,22,22,22,NA, 8,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray", "black","red", NA, NA, "white", "white", "white"), bg="white", xpd=T, pt.cex=c(2, 2,2,2,2,2,1,1,1,1,1))
text(14, 6.1, "Temperature", xpd=T)
text(6, 5.3, "Time", xpd=T)

dev.off()

#### Mass loss graphs

pcwloss <- read.csv("pcwloss.csv")
png(filename="percentloss.png",
    units="cm",
    width=14, 
    height=13, 
    res=300,
    pointsize=8)

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

legend(85, 70.7, c("", "215 ˚C", "230 ˚C", "245 ˚C", "260 ˚C", "300 ˚C","", "Uncharred", "4h", "8h", "24h"), pch=c(NA, 22,22,22,22,22,NA, 8,21,22,23), col="black", pt.bg=c(NA, "white", "lightgray", "darkgray", "black", "red", NA, NA, "white", "white", "white"), bg="white", xpd=T, pt.cex=c(2, 2,2,2,2,2,2,1,1,1,1,1))
text(93, 68, "Temperature", xpd=T)
text(88.7, 52, "Time", xpd=T)

dev.off()

