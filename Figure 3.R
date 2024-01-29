### This code produces Figure 3, and associated analyses, in: 
### Hooks AP, Burgess SC. Variation in polyandry, reproductive output, 
# and within-brood genetic diversity in a marine snail population 
# across seasons and years. Marine Ecology Progress Series
# Code written by Alex Hooks, hooksap@gmail.com
# with input from Scott Burgess, sburgess@bio.fsu.edu

# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"

# Load required libraries
library(glmmTMB) 

# Import data
eff.sire<- read.csv("Data for Figure 3 and 4.csv", header=T,stringsAsFactors=T)

# Change year to a factor
eff.sire$year <- factor(eff.sire$year)
summary(eff.sire)
head(eff.sire)


# Did the number of sires vary between years or seasons?
m1 <- glmmTMB(number.sire ~ year * time.season, family='nbinom1',data=eff.sire)
summary(m1)
m2 <- glmmTMB(number.sire ~ year + time.season, family='nbinom1', data=eff.sire)
anova(m1,m2,test="Chisq") # no interaction
m3 <- glmmTMB(number.sire ~ year, family='nbinom1', data=eff.sire)
anova(m2,m3,test="Chisq") # no effect of time.season (X2=2.2695, df=1, p=0.1319)
m4 <- glmmTMB(number.sire ~ time.season, family='nbinom1', data=eff.sire)
anova(m2,m4,test="Chisq") # no effect of year (X2=2.2761, df=1, p=0.1314)
m5 <- glmmTMB(number.sire ~ 1, family='nbinom1', data=eff.sire)
anova(m3,m5,test="Chisq") # no effect of year
anova(m4,m5,test="Chisq") # no effect of time.season
a.year.p.value <- round(anova(m2,m4,test="Chisq")$'Pr(>Chisq)'[2],3) # populates p values onto figures
a.time.season.p.value <- round(anova(m2,m3,test="Chisq")$'Pr(>Chisq)'[2],3) # for plotting

# On average, there were 9.13 (95% CI 7.64 - 10.90) sires per egg string per mother
exp(confint(m5)) # 7.64 - 10.89
# One mother had only 2 sires per egg string, while one mother had 19 sires per egg string 
range(eff.sire$number.sire) #range 2-19
mean(eff.sire$number.sire) #mean 9.125

m2 <- glmmTMB(number.sire ~ year + time.season, family='nbinom1', data=eff.sire)
seasonvec<- unique(eff.sire$time.season)
yearvec<- unique(eff.sire$year)
newdata.s <- expand.grid(time.season=seasonvec, year=yearvec)
pred <- predict(m2,newdata.s,se.fit = T,type="response")
newdata.s <- data.frame(
  newdata.s, 
  fit = pred$fit, 
  plo = pred$fit - (2 * pred$se.fit), 
  phi = pred$fit + (2 * pred$se.fit)) 




# Make plot - Figure 3
quartz(width=9,height=5)
par(mfrow=c(1,2),mar=c(3,4,1,2),oma=c(2,2,2,2))

#pane a histogram
with(eff.sire,hist(number.sire,breaks=c(0,2,4,6,8,10,12,14,16,18,20),main="",las=1,ylab="",xlab="",xaxt="n",cex.axis=1.5))
axis(side=1,at=seq(0,20,by=4),cex.axis=1.5)
mtext(side=2,"Number of mothers",line=2.5,cex=1.5)
mtext(side=1,"Number of sires at hatching\nwithin a brood",line=4,cex=1.5)
mtext("a)",side=3,adj=0,cex=1.5,line=1)


#pane b
offset <- 0.15

cols <- c(adjustcolor("dodgerblue",alpha.f=0.5),adjustcolor("tomato",alpha.f=0.5))

plot(c(0.5,2.5),c(0,20),type="n",las=1,cex=2, cex.axis=1.5,ylab="",xlab="",bty="l",xaxt="n")
axis(side=1,line=0,at=c(1,2),labels=c("Early\n(April)","Late\n(July)"),cex.axis=1.5,padj=0.8)
mtext(side=2,"Number of sires at hatching\nwithin a brood",line=2.5,cex=1.5)
mtext("b)",side=3,adj=0,cex=1.5,line=1)
mtext(eval(paste("time of season p =",a.time.season.p.value)),side=1,adj=.9,cex=1,line=-2.25)
mtext(eval(paste("year p =",a.year.p.value)),side=1,adj=.91,cex=1,line=-1.25)

#Early
set.seed(1);with(eff.sire[eff.sire$year=="2018" & eff.sire$time.season=="E",],points(jitter(as.numeric(as.factor(time.season))-offset,0.3), number.sire,col=cols[1],pch=19,cex=.75))
with(newdata.s[newdata.s$time.season=="E" & newdata.s$year =="2018",], segments(c(1-offset),plo,c(1-offset),phi,lwd=2.3,col="dodgerblue4"))
with(newdata.s[newdata.s$time.season=="E" & newdata.s$year =="2018",], points(c(1-offset),fit,pch=21,cex=2,lwd=2,bg=cols[1],col="dodgerblue4"))

set.seed(1);with(eff.sire[eff.sire$year=="2020" & eff.sire$time.season=="E",],points(jitter(as.numeric(as.factor(time.season))+offset,0.3), number.sire,col=cols[1],pch=17,cex=.75))
with(newdata.s[newdata.s$time.season=="E" & newdata.s$year =="2020",], segments(c(1+offset),plo,c(1+offset),phi,lwd=2.3,col="dodgerblue4"))
with(newdata.s[newdata.s$time.season=="E" & newdata.s$year =="2020",], points(c(1+offset),fit,pch=24,cex=2,lwd=2,bg=cols[1],col="dodgerblue4"))

#Late
set.seed(1);with(eff.sire[eff.sire$year=="2018" & eff.sire$time.season=="L",],points(jitter(as.numeric(as.factor(time.season))-offset,0.3), number.sire,col=cols[2],pch=19,cex=.75))
with(newdata.s[newdata.s$time.season=="L" & newdata.s$year =="2018",], segments(c(2-offset),plo,c(2-offset),phi,lwd=2.3,col="tomato4"))
with(newdata.s[newdata.s$time.season=="L" & newdata.s$year =="2018",], points(c(2-offset),fit,pch=21,cex=2,lwd=2,bg=cols[2],col="tomato4"))

set.seed(1);with(eff.sire[eff.sire$year=="2020" & eff.sire$time.season=="L",],points(jitter(as.numeric(as.factor(time.season))+offset,0.3), number.sire,col=cols[2],pch=17,cex=.75))
with(newdata.s[newdata.s$time.season=="L" & newdata.s$year =="2020",], segments(c(2+offset),plo,c(2+offset),phi,lwd=2.3,col="tomato4"))
with(newdata.s[newdata.s$time.season=="L" & newdata.s$year =="2020",], points(c(2+offset),fit,pch=24,cex=2,lwd=2,bg=cols[2],col="tomato4"))


par(new=T,xpd=NA)
legend("topright",inset = c(.4, -.1),legend=("2018"),pch=21,cex=1.5,lwd=2,bty="n",pt.cex = 1.5,x.intersp=.4)
legend('topright',inset = c(.04, -.1),legend=("2020"),pch=24,cex=1.5,lwd=2,bty="n",pt.cex = 1.5,x.intersp=.3)

