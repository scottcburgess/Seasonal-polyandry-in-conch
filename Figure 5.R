### This code produces Figure 5, and associated analyses, in: 
### Hooks AP, Burgess SC. Variation in polyandry, reproductive output, 
# and within-brood genetic diversity in a marine snail population 
# across seasons and years. Marine Ecology Progress Series
# Code written by Alex Hooks, hooksap@gmail.com
# with input from Scott Burgess, sburgess@bio.fsu.edu

# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"

# Load required libraries
library(glmmTMB)
library(DHARMa)

# Import data
fieldDat_mother <- read.csv("Data A for Figure 5.csv", header=T)
dat<- read.csv("Data B for Figure 5.csv", header=T)
eff.sire<- read.csv("Data for Figure 3 and 4.csv", header=T,stringsAsFactors=T)

# Change time.season and year to a factor
fieldDat_mother$time.season<-as.factor(fieldDat_mother$time.season)
fieldDat_mother$year<-as.factor(fieldDat_mother$year)


# Did the number of sires per brood explain the 
# estimated total number of hatchlings within a brood?
# not plotted, just reported under 3.3 Reproductive output, in the paper

range(eff.sire$capsule.num) # Mothers laid between 6 to 21 egg capsules on an egg string 
range(eff.sire$total.output) # and produced 355 to 2190 estimated total number of hatchlings within an egg string

m1 <- lm(total.output ~ poly(number.sire,2) * time.season * year, data=eff.sire)
m1a <- lm(total.output ~ number.sire * time.season * year, data=eff.sire)
anova(m1,m1a,test="F") # no curvi-linear

m2 <- lm(total.output ~ number.sire * time.season * year, data=eff.sire)
m3 <- lm(total.output ~ number.sire + time.season + year, data=eff.sire)
anova(m2,m3) # no interactions

m4 <- lm(total.output ~ number.sire + time.season, data=eff.sire)
m5 <- lm(total.output ~ number.sire + year, data=eff.sire)
m6 <- lm(total.output ~ time.season + year, data=eff.sire)
anova(m3,m4) # year did not affect total.output (F= 0.082, df=1,20, p=0.778)
anova(m3,m5) # time.season did not affect total.output (F= 0.0176, df=1,20, p=0.900)
anova(m3,m6) # number.sire (F=3.21, df=1,20, p=0.088)

# However, the number of sires per brood did not explain the estimated 
# total number of hatchlings within a brood (ANOVA: F1,20=3.21, p=0.088, Fig 3b). 

# Check model assumptions
simulationOutput <- simulateResiduals(fittedModel = m3, plot = F)
plot(simulationOutput) # looks fine

# Get the estimated parameters and 95% confidence intervals
cbind.data.frame(round(confint(m3),2),fit=round(coef(m3),2))
# The 95% confidence interval for the effect of the number of sires 
# on the total number of hatchlings was -93.78 – 7.11.




# Number of egg capsules per mother
m1<-glmmTMB(capsule.num~time.season*year, family="nbinom1",data= fieldDat_mother)
m2<-glmmTMB(capsule.num~time.season+year, family="nbinom1", data= fieldDat_mother)
m3<-glmmTMB(capsule.num~time.season, family="nbinom1", data= fieldDat_mother)
m4<-glmmTMB(capsule.num~year, family="nbinom1", data= fieldDat_mother)
m5<-glmmTMB(capsule.num~1, family="nbinom1", data= fieldDat_mother)
anova(m1,m2,test="Chisq") #NS (X2=0.1786, df=1, p= 0.6726)
anova(m2,m3,test="Chisq") #NS (X2=0.2197, df=1, p=0.6393), no year affect
anova(m2,m4,test="Chisq") #NS (X2=0.1469, df=1, p=0.7015), no time.season affect
anova(m3,m5,test="Chisq") #NS
anova(m4,m5,test="Chisq") #NS
a.year.p.value <- round(anova(m2,m3,test="Chisq")$'Pr(>Chisq)'[2],3) # populates p values onto figures
a.time.season.p.value <- round(anova(m2,m4,test="Chisq")$'Pr(>Chisq)'[2],3) # for plotting

# Check model assumptions
simulationOutput <- simulateResiduals(fittedModel = m2, plot = F)
plot(simulationOutput) # looks fine

# Get predicted values for plotting
m2<-glmmTMB(capsule.num~time.season+year, family="nbinom1", data= fieldDat_mother)
seasonvec<- unique(fieldDat_mother$time.season)
yearvec<- unique(fieldDat_mother$year)
newdata.ec <- expand.grid(time.season=seasonvec, year=yearvec)
pred <- predict(m2,newdata.ec,se.fit = T,type="response")
newdata.ec <- data.frame(
  newdata.ec, 
  fit = pred$fit, 
  plo = pred$fit - (2 * pred$se.fit), 
  phi = pred$fit + (2 * pred$se.fit)) 


# Number of hatchlings per egg capsule
# removed egg capsules that died in this analysis 

# Using raw averages (offspring.num not rounded to nearest whole number)
m1 <- lm(offspring.num ~ time.season * year, data=fieldDat_mother)
m2 <- lm(offspring.num ~ time.season + year, data=fieldDat_mother)
m3 <- lm(offspring.num ~ time.season, data=fieldDat_mother)
m4 <- lm(offspring.num ~ year, data=fieldDat_mother)
m5 <- lm(offspring.num ~ 1, data=fieldDat_mother)
anova(m1,m2) #NS no interaction
anova(m2,m3) #NS no affect of year (F=0.2026, df=1,52, p=0.655)
anova(m2,m4) #S time.season on offspring number (F=6.5383, df=1,52, p=0.014)
anova(m3,m5) #S
anova(m4,m5) #NS
b.year.p.value <- round(data.frame(anova(m2,m3))$Pr[2],3) #populates p values onto figures
b.time.season.p.value <- round(data.frame(anova(m2,m4))$Pr[2],3) # for plotting

# Check model assumptions
simulationOutput <- simulateResiduals(fittedModel = m2, plot = F)
plot(simulationOutput) # looks fine

# On average, there were 72.27 (61.09 – 83.46, 95% CI) 
# hatchlings per egg capsule early in the reproductive season 
# and 90.67 (81.53 – 99.80, 95% CI) hatchlings per capsule 
# late in the reproductive season
round(predict(m3,list(time.season=c("E","L")),interval="confidence"),2)


# Get predicted values for plotting
seasonvec<- unique(fieldDat_mother $time.season)
yearvec<- unique(fieldDat_mother $year)
newdata.hat <- expand.grid(time.season=seasonvec, year=yearvec)
pred <- predict(m2,newdata.hat,se.fit = T,type="response")
newdata.hat <- data.frame(
  newdata.hat, 
  fit = pred$fit, 
  plo = pred$fit - (2 * pred$se.fit), 
  phi = pred$fit + (2 * pred$se.fit)) 


# Egg capsule viability 2020 
dat2020 <- dat[dat$year=="2020",] # 2020 data only (egg capsules in maintained in field, field survival), multiple egg capsules per mother

m1 <- glmmTMB(survive ~ time.season * capsule.size, data = dat2020, family="binomial")
m2 <- glmmTMB(survive ~ time.season + capsule.size, data = dat2020, family="binomial")
m3 <- glmmTMB(survive ~ time.season, data = dat2020, family="binomial")
m4 <- glmmTMB(survive ~ capsule.size, data = dat2020, family="binomial")
m5 <- glmmTMB(survive ~ 1, data = dat2020, family="binomial")
anova(m1,m2,test="Chisq") # NS no interaction
anova(m2,m3,test="Chisq") # NS no affect of capsule size (X2=-1.2965, df=-1, p=0.255)
anova(m2,m4,test="Chisq") # S time.season affect on survival (X2=-25.691,df=-1, p<0.001)
anova(m3,m5,test="Chisq") # S
anova(m4,m5,test="Chisq") # NS

c.capsule.size.p.value <- round(anova(m2,m3,test="Chisq")$'Pr(>Chisq)'[2],3) 

# The odds of egg capsule survival late in the season was
# 18% (9% – 38%, 95% CI) of that early in the reproductive season
exp(confint(m2))[2,]

# Get predicted values for plotting
seasonvec <- sort(unique(dat2020$time.season))
capsule.sizevec <- seq(min(dat2020$capsule.size),max(dat2020$capsule.size),length=100)
newdata.s <- expand.grid(time.season=seasonvec, capsule.size= capsule.sizevec)
pred <- as.data.frame(predict(m2,newdata.s,se.fit=T))
newdata.s <- data.frame(
  newdata.s,
  fit = pred$fit, 
  plo = pred$fit - 2 * pred$se.fit, 
  phi = pred$fit + 2 * pred$se.fit) 




# Total number of hatchlings per mother
# (includes failed egg capsules as zeros)
dat0<-dat
# make egg capsules that died have zero hatchlings (instead of NA)
dat0[is.na(dat0)] <- 0 
fieldDat0<- with(dat0,  aggregate(cbind(offspring.num, capsule.size), by = list(mother.id, year, time.season, capsule.num), FUN = median))
names(fieldDat0) <- c("mother.id", "year", "time.season", "capsule.num","offspring.num","capsule.size")
fieldDat0$total.output <- fieldDat0$offspring.num * fieldDat0$capsule.num 
fieldDat0$total.output<-round(fieldDat0$total.output,0) #rounding total.output for hurdle model
head(fieldDat0)


# split by year due to distribution differences (i.e., zeros in 2020)
# 2018
m1<-glmmTMB(total.output ~ time.season, family="nbinom1",data=fieldDat0[fieldDat0$year=="2018",])
m2<-glmmTMB(total.output ~ 1, family="nbinom1",data=fieldDat0[fieldDat0$year=="2018",])
anova(m1,m2) # NS (X<0.001,df=1, p=0.995), no affect of time.season
d.total.output18.p.value <- data.frame(p=round(anova(m1,m2)$'Pr(>Chisq)'[2],3))

# Check model assumptions
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput) 

# Get predicted values for plotting
seasonvec<- unique(fieldDat0[fieldDat0$year=="2018",]$time.season)
newdata.18 <- expand.grid(time.season=seasonvec)
pred <- predict(m1,newdata.18,se.fit = T,type="response")
newdata.18 <- data.frame(
  newdata.18, 
  fit = pred$fit, 
  plo = pred$fit - (2 * pred$se.fit), 
  phi = pred$fit + (2 * pred$se.fit)) 


# 2020
m1 <- glmmTMB(total.output ~ time.season,family="nbinom1", data=fieldDat0[fieldDat0$year=="2020",])
# Check model assumptions
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput) # no good

m1 <- glmmTMB(total.output ~ time.season,ziformula = ~time.season,family="truncated_nbinom2", data=fieldDat0[fieldDat0$year=="2020",])
# Check model assumptions
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput) # good  

m2 <- glmmTMB(total.output ~ 1,ziformula = ~time.season,family="truncated_nbinom2", data=fieldDat0[fieldDat0$year=="2020",])
anova(m1,m2,test="Chisq") # NS (X2=0.004, df=1, p=0.9496), no affect time.season
d.total.output2020.p.value <- data.frame(p=round(anova(m1,m2)$'Pr(>Chisq)'[2],3))


# Get predicted values for plotting
seasonvec<- unique(fieldDat0[fieldDat0$year=="2020",]$time.season)
newdata.20 <- expand.grid(time.season=seasonvec)
pred <- predict(m1,newdata.20,se.fit = T,type="response")
newdata.20 <- data.frame(
  newdata.20, 
  fit = pred$fit, 
  plo = pred$fit - (2 * pred$se.fit), 
  phi = pred$fit + (2 * pred$se.fit)) 






# Make plot - Figure 5
offset <- 0.15

quartz(width=8.5,height=7)
par(mfrow=c(2,2), oma=c(5,3,1,1.5), mar=c(2,4,3,3))
cols <- c(adjustcolor("dodgerblue",alpha.f=0.5),adjustcolor("tomato",alpha.f=0.5))

# panel a
plot(c(0.5,2.5),c(0,22),type="n",las=1,cex=2, cex.axis=1.5,ylab="",xlab="",bty="l",xaxt="n")
axis(side=1,line=0,at=c(1,2),labels=c("Early\n(April)","Late\n(July)"),cex.axis=1.5,padj=0.8)
mtext("Number of egg capsules\n per mother",side=2,line=3,cex=1.2)
mtext("a)",side=3,adj=0,cex=1.5,line=1)
mtext(eval(paste("time of season p =",a.time.season.p.value)),side=1,adj=.9,cex=1,line=-2.25)
mtext(eval(paste("year p =",a.year.p.value)),side=1,adj=.91,cex=1,line=-1.25)

#Early
set.seed(1);with(fieldDat_mother[fieldDat_mother$year=="2018" & fieldDat_mother$time.season=="E",],points(jitter(as.numeric(as.factor(time.season))-offset,0.3), capsule.num,col=cols[1],pch=19,cex=.75))
with(newdata.ec[newdata.ec$time.season=="E" & newdata.ec$year =="2018",], segments(c(1-offset),plo,c(1-offset),phi,lwd=2.3,col="dodgerblue4"))
with(newdata.ec[newdata.ec$time.season=="E" & newdata.ec$year =="2018",], points(c(1-offset),fit,pch=21,cex=2,lwd=2,bg=cols[1],col="dodgerblue4"))

set.seed(1);with(fieldDat_mother[fieldDat_mother$year=="2020" & fieldDat_mother$time.season=="E",],points(jitter(as.numeric(as.factor(time.season))+offset,0.3), capsule.num,col=cols[1],pch=17,cex=.75))
with(newdata.ec[newdata.ec$time.season=="E" & newdata.ec$year =="2020",], segments(c(1+offset),plo,c(1+offset),phi,lwd=2.3,col="dodgerblue4"))
with(newdata.ec[newdata.ec$time.season=="E" & newdata.ec$year =="2020",], points(c(1+offset),fit,pch=24,cex=2,lwd=2,bg=cols[1],col="dodgerblue4"))

#Late
set.seed(1);with(fieldDat_mother[fieldDat_mother$year=="2018" & fieldDat_mother$time.season=="L",],points(jitter(as.numeric(as.factor(time.season))-offset,0.3), capsule.num,col=cols[2],pch=19,cex=.75))
with(newdata.ec[newdata.ec$time.season=="L" & newdata.ec$year =="2018",], segments(c(2-offset),plo,c(2-offset),phi,lwd=2.3,col="tomato4"))
with(newdata.ec[newdata.ec$time.season=="L" & newdata.ec$year =="2018",], points(c(2-offset),fit,pch=21,cex=2,lwd=2,bg=cols[2],col="tomato4"))

set.seed(1);with(fieldDat_mother[fieldDat_mother$year=="2020" & fieldDat_mother$time.season=="L",],points(jitter(as.numeric(as.factor(time.season))+offset,0.3), capsule.num,col=cols[2],pch=17,cex=.75))
with(newdata.ec[newdata.ec$time.season=="L" & newdata.ec$year =="2020",], segments(c(2+offset),plo,c(2+offset),phi,lwd=2.3,col="tomato4"))
with(newdata.ec[newdata.ec$time.season=="L" & newdata.ec$year =="2020",], points(c(2+offset),fit,pch=24,cex=2,lwd=2,bg=cols[2],col="tomato4"))


# panel b
plot(c(0.5,2.5),c(0,175),type="n",las=1,cex=2, cex.axis=1.5,ylab="",xlab="",bty="l",xaxt="n")
axis(side=1,line=0,at=c(1,2),labels=c("Early\n(April)","Late\n(July)"),cex.axis=1.5,padj=0.8)
mtext("Number of hatchlings\n per egg capsule",side=2,line=3.5,cex=1.2)
mtext("b)",side=3,adj=0,cex=1.5,line=1)
mtext(eval(paste("time of season p =",b.time.season.p.value)),side=1,adj=.9,cex=1,line=-2.25)
mtext(eval(paste("year p =",b.year.p.value)),side=1,adj=.91,cex=1,line=-1.25)


#Early
set.seed(1);with(fieldDat_mother[fieldDat_mother$year=="2018" & fieldDat_mother$time.season=="E",],points(jitter(as.numeric(as.factor(time.season))-offset,0.3), offspring.num,col=cols[1],pch=19,cex=.75))
with(newdata.hat[newdata.hat$time.season=="E" & newdata.hat$year =="2018",], segments(c(1-offset),plo,c(1-offset),phi,lwd=2.3,col="dodgerblue4"))
with(newdata.hat[newdata.hat$time.season=="E" & newdata.hat$year =="2018",], points(c(1-offset),fit,pch=21,cex=2,lwd=2,bg=cols[1],col="dodgerblue4"))

set.seed(1);with(fieldDat_mother[fieldDat_mother$year=="2020" & fieldDat_mother$time.season=="E",],points(jitter(as.numeric(as.factor(time.season))+offset,0.3), offspring.num,col=cols[1],pch=17,cex=.75))
with(newdata.hat[newdata.hat$time.season=="E" & newdata.hat$year =="2020",], segments(c(1+offset),plo,c(1+offset),phi,lwd=2.3,col="dodgerblue4"))
with(newdata.hat[newdata.hat$time.season=="E" & newdata.hat$year =="2020",], points(c(1+offset),fit,pch=24,cex=2,lwd=2,bg=cols[1],col="dodgerblue4"))

#Late
set.seed(1);with(fieldDat_mother[fieldDat_mother$year=="2018" & fieldDat_mother$time.season=="L",],points(jitter(as.numeric(as.factor(time.season))-offset,0.3), offspring.num,col=cols[2],pch=19,cex=.75))
with(newdata.hat[newdata.hat$time.season=="L" & newdata.hat$year =="2018",], segments(c(2-offset),plo,c(2-offset),phi,lwd=2.3,col="tomato4"))
with(newdata.hat[newdata.hat$time.season=="L" & newdata.hat$year =="2018",], points(c(2-offset),fit,pch=21,cex=2,lwd=2,bg=cols[2],col="tomato4"))

set.seed(1);with(fieldDat_mother[fieldDat_mother$year=="2020" & fieldDat_mother$time.season=="L",],points(jitter(as.numeric(as.factor(time.season))+offset,0.3), offspring.num,col=cols[2],pch=17,cex=.75))
with(newdata.hat[newdata.hat$time.season=="L" & newdata.hat$year =="2020",], segments(c(2+offset),plo,c(2+offset),phi,lwd=2.3,col="tomato4"))
with(newdata.hat[newdata.hat$time.season=="L" & newdata.hat$year =="2020",], points(c(2+offset),fit,pch=24,cex=2,lwd=2,bg=cols[2],col="tomato4"))

par(new=T,xpd=NA)
legend("topright",inset = c(.35, -.1),legend=("2018"),pch=21,cex=1.2,lwd=2,bty="n",pt.cex = 1.5,x.intersp=.5)
legend('topright',inset = c(.04, -.1),legend=("2020"),pch=24,cex=1.2,lwd=2,bty="n",pt.cex = 1.5,x.intersp=.3)



# panel c
plot(c(min(capsule.sizevec), max(capsule.sizevec)), c(0,1.1),type="n",las=1,cex=2, cex.axis=1.5,ylab="",xlab="",bty="l",xaxt="n")
axis(side=1,line=0,at= seq(30,55,5),cex.axis=1.5)
mtext("Egg capsule viability (2020)",side=2,line=4,cex=1.2)
mtext("c)",side=3,adj=0,cex=1.5,line=1)
mtext("time of season p<0.001",side=1,adj=.9,cex=1,line=-2.25)
mtext(eval(paste("capsule.size p =",c.capsule.size.p.value)),side=1,adj=.91,cex=1,line=-1.25)
mtext("Capsule diameter (mm)", side=1, line=3, cex=1.2)

with(newdata.s[newdata.s$time.season=="E",],
	polygon(c(capsule.size,rev(capsule.size)),c(plogis(plo),rev(plogis(phi))),col=cols[1],border=F))
with(newdata.s[newdata.s$time.season=="E",], lines(capsule.size,plogis(fit),lwd=2, col="dodgerblue4"))
text(40,1,"Early",cex=1.2)

with(newdata.s[newdata.s$time.season=="L",],
	polygon(c(capsule.size,rev(capsule.size)),c(plogis(plo),rev(plogis(phi))),col=cols[2],border=F))
with(newdata.s[newdata.s$time.season=="L",], lines(capsule.size,plogis(fit),lwd=2,col="tomato4"))
text(40,0.4,"Late",cex=1.2)


# panel d
plot(c(0.5,2.5),c(-800,3500),type="n",las=1,cex=2, cex.axis=1.5,ylab="",xlab="",bty="l",xaxt="n")
axis(side=1,line=0,at=c(1,2),labels=c("Early\n(April)","Late\n(July)"),cex.axis=1.5,padj=0.8) 
mtext("Estimated total number \n of hatchlings per mother",side=2,line=4,cex=1.2)
mtext("d)",side=3,adj=0,cex=1.5,line=1)
mtext(eval(paste("time of season (2018) p =",d.total.output18.p.value)),side=1,adj=.9,cex=1,line=-2.25)
mtext(eval(paste("time of season (2020) p =",d.total.output2020.p.value)),side=1,adj=.91,cex=1,line=-1.25)

#Early
set.seed(1);with(fieldDat0[fieldDat0$year=="2018" & fieldDat0$time.season=="E",],points(jitter(as.numeric(as.factor(time.season))-offset,0.3), total.output,col=cols[1],pch=19,cex=.75))
with(newdata.18[newdata.18$time.season=="E",], segments(c(1-offset),(plo),c(1-offset),(phi),lwd=2.3,col="dodgerblue4"))
with(newdata.18[newdata.18$time.season=="E" ,], points(c(1-offset),(fit),pch=21,cex=2,lwd=2,bg=cols[1],col="dodgerblue4"))

set.seed(1);with(fieldDat0[fieldDat0$year=="2020" & fieldDat0$time.season=="E",],points(jitter(as.numeric(as.factor(time.season))+offset,0.3), total.output,col=cols[1],pch=17,cex=.75))
with(newdata.20[newdata.20$time.season=="E",] , segments(c(1+offset),(plo),c(1+offset),(phi),lwd=2.3,col="dodgerblue4"))
with(newdata.20[newdata.20$time.season=="E",], points(c(1+offset),(fit),pch=24,cex=2,lwd=2,bg=cols[1],col="dodgerblue4"))

#Late
set.seed(1);with(fieldDat0[fieldDat0$year=="2018" & fieldDat0$time.season=="L",],points(jitter(rep(2,length(time.season))-offset,0.3),total.output,col=cols[2],pch=19,cex=.75))
with(newdata.18[newdata.18$time.season=="L",], segments(c(2-offset),(plo),c(2-offset),(phi),lwd=2.3,col="tomato4"))
with(newdata.18[newdata.18$time.season=="L",], points(c(2-offset),(fit),pch=21,cex=2,lwd=2,bg=cols[2],col="tomato4"))

set.seed(1);with(fieldDat0[fieldDat0$year=="2020" & fieldDat0$time.season=="L",],points(jitter(rep(2,length(time.season))+offset,0.3), total.output,col=cols[2],pch=19,cex=.75))
with(newdata.20[newdata.20$time.season=="L",], segments(c(2+offset),plo,c(2+offset),phi,lwd=2.3,col="tomato4"))
with(newdata.20[newdata.20$time.season=="L",], points(c(2+offset),fit,pch=24,cex=2,lwd=2,bg=cols[2],col="tomato4"))

