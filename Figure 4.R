### This code produces Figure 4, and associated analyses, in: 
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
eff.sire<- read.csv("Data for Figure 3 and 4.csv", header=T,stringsAsFactors=T)

# Change year to a factor
eff.sire$year <- factor(eff.sire$year)
summary(eff.sire)
head(eff.sire)


# Set up vectors for calculating predicted values for plotting
timevec <- c(rep("E",100),rep("L",100))
eff.sire.vec <- c(
  seq(with(eff.sire[eff.sire$time.season=="E",], min(number.sire)),
      with(eff.sire[eff.sire$time.season=="E",], max(number.sire)),
      length=100),
  seq(with(eff.sire[eff.sire$time.season=="L",], min(number.sire)),
      with(eff.sire[eff.sire$time.season=="L",], max(number.sire)),
      length=100)
)
year.vec <- c(rep("2018",length(timevec)),rep("2020",length(timevec)))



# Allelic richness
# Is there evidence for curvilinearity? 
m1 <- lm(Ar ~ poly(number.sire,2) * time.season * year, data=eff.sire)
m1a <- lm(Ar ~ number.sire * time.season * year, data=eff.sire)
anova(m1,m1a) # no curvilinearity

m2 <- update(m1a,.~.-number.sire:time.season:year)
anova(m1a,m2) # There is NOT a three-way interaction between number, time.season, and year
# F(1,16)= 4.24, p=0.056) 

m3 <- lm(Ar ~ number.sire + time.season + year, data=eff.sire)
anova(m1a,m3) # No evidence for any 2 or 3 way interactions 

m4 <- lm(Ar ~ time.season + year, data=eff.sire)
anova(m3,m4) # The number of sires increased Ar (F(1,20)=6.93, p = 0.016)


# Check assumptions
simulationOutput <- simulateResiduals(fittedModel = m3, plot = F)
plot(simulationOutput) # looks fine

# Remove the two mothers with 15+ sires to check their influence
m1 <- lm(Ar ~ poly(number.sire,2) * time.season * year, data=eff.sire[eff.sire$number.sire<15,])
m1a <- lm(Ar ~ number.sire * time.season * year, data=eff.sire[eff.sire$number.sire<15,])
anova(m1,m1a) # Still no evidence for curvi-linear relationship, so its not caused by >15 mothers.
m3 <- lm(Ar ~ number.sire + time.season + year, data=eff.sire[eff.sire$number.sire<15,])
m3a <- lm(Ar ~ poly(number.sire,2) + time.season + year, data=eff.sire[eff.sire$number.sire<15,])
anova(m3,m3a) # Still no evidence for curvi-linear relationship, so its not caused by >15 mothers.

# Move forward with this model:
m3 <- lm(Ar ~ number.sire + time.season + year, data=eff.sire)
# Extract predicted values for plotting
newdata.ar <- cbind.data.frame(time.season = timevec, number.sire = eff.sire.vec, year=year.vec)
p <- predict(m3,newdata.ar,se.fit=T)
newdata.ar$fit <- p$fit
newdata.ar$se <- p$se.fit
newdata.ar$plo <- newdata.ar$fit - 2*newdata.ar$se
newdata.ar$phi <- newdata.ar$fit + 2*newdata.ar $se



# Expected heterozygosity - Hexp
m1 <- lm(Hexp ~ poly(number.sire,2) * time.season * year, data=eff.sire)
m1a <- lm(Hexp ~ number.sire * time.season * year, data=eff.sire)
anova(m1,m1a) # no evidence for curvi-linear relationship

m2 <- update(m1a,.~.-number.sire:time.season:year)
anova(m1a,m2) # There is NOT a three-way interaction between number, time.season, and year

m3 <- lm(Hexp ~ number.sire + time.season + year, data=eff.sire)
anova(m1a,m3) # None of the two way interactions were significant

m4 <- lm(Hexp ~ time.season + year, data=eff.sire)
anova(m3,m4) # The number of sires increased Hexp (F(1,20)=7.90, p = 0.012)


# library('AICcmodavg')
# aictab(list(m1=m1,m1a=m1a,m2=m2,m3=m3))

# Check assumptions
simulationOutput <- simulateResiduals(fittedModel = m3, plot = F)
plot(simulationOutput) # looks fine


# Move forward with this model:
m3 <- glmmTMB(Hexp ~ number.sire + time.season + year, data=eff.sire)
# Extract predicted values for plotting
newdata.he <- cbind.data.frame(time.season = timevec, number.sire = eff.sire.vec, year=year.vec)
p <- predict(m3,newdata.ar,se.fit=T)
newdata.he$fit <- p$fit
newdata.he$se <- p$se.fit
newdata.he$plo <- newdata.he$fit - 2*newdata.he$se
newdata.he$phi <- newdata.he$fit + 2*newdata.he $se




# Make plot - Figure 4
quartz(width=9,height=7)
par(mfrow=c(2,2),oma=c(4.5,4.5,1.5,1.5))
cols <- c(adjustcolor("dodgerblue",alpha.f=0.5),adjustcolor("tomato",alpha.f=0.5))

# Panel a, 2018
par(mar=c(3,4,4,1.5))
plot(c(0,20),c(0,10),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l",axes=F)
axis(side=1,at=seq(0,20,by=4),cex.axis=1.5)
axis(side=2,at=seq(0,10, by=2),las=1,cex.axis=1.5)
mtext("Allelic richness",side=2,line=4,cex=1.5,padj=0.4)
mtext("a)",side=3,adj=0,cex=1.5, line=0)
mtext("2018",side=3,adj=.5,cex=1.5,line=2)

y1 <- eff.sire[eff.sire$time.season=="E" & eff.sire$year=="2018",]
with(y1, points(number.sire, Ar,pch=19,col=cols[1]))
y2 <- eff.sire[eff.sire$time.season=="L" & eff.sire$year=="2018",]
with(y2, points(number.sire, Ar,pch=17,col=cols[2]))

xmin <- range(y1$number.sire)
y <- newdata.ar[newdata.ar$year=="2018" &
                  newdata.ar$time.season=="E" &
                  newdata.ar$number.sire > xmin[1] &
                  newdata.ar$number.sire < xmin[2], ]
with(y, polygon(c(number.sire, rev(number.sire)), c(plo,rev(phi)),lwd=2,border=F, col=cols[1]))
with(y, lines(number.sire,fit,lwd=2,col="dodgerblue"))

xmin <- range(y2$number.sire)
y <- newdata.ar[newdata.ar$year=="2018" &
                  newdata.ar$time.season=="L" &
                  newdata.ar$number.sire > xmin[1] &
                  newdata.ar$number.sire < xmin[2], ]
with(y,polygon(c(number.sire, rev(number.sire)), c(plo,rev(phi)),lwd=2,border=F, col=cols[2]))
with(y, lines(number.sire,fit,lwd=2,col="tomato"))


# Panel a, 2020
par(mar=c(3,2,4,1.5))
plot(c(0,20),c(0,10),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l",axes=F)
axis(side=1,at=seq(0,20,by=4),cex.axis=1.5)
axis(side=2,at=seq(0,10, by=2),las=1,cex.axis=1.5)
#mtext("Number of sires",side=1,line=3,cex=1.2)
mtext("2020",side=3,adj=.5,cex=1.5,line=2)

y1 <- eff.sire[eff.sire$time.season=="E" & eff.sire$year=="2020",]
with(y1, points(number.sire, Ar,pch=19,cex=1.5,col=cols[1]))
y2 <- eff.sire[eff.sire$time.season=="L" & eff.sire$year=="2020",]
with(y2, points(number.sire, Ar,pch=17,cex=1.5,col=cols[2]))

xmin <- range(y1$number.sire)
y <- newdata.ar[newdata.ar$year=="2020" &
                  newdata.ar$time.season=="E" &
                  newdata.ar$number.sire > xmin[1] &
                  newdata.ar$number.sire < xmin[2], ]
with(y, polygon(c(number.sire, rev(number.sire)), c(plo,rev(phi)),lwd=2,border=F, col=cols[1]))
with(y, lines(number.sire,fit,lwd=2,col="dodgerblue"))

xmin <- range(y2$number.sire)
y <- newdata.ar[newdata.ar$year=="2020" &
                  newdata.ar$time.season=="L" &
                  newdata.ar$number.sire > xmin[1] &
                  newdata.ar$number.sire < xmin[2], ]
with(y, polygon(c(number.sire, rev(number.sire)), c(plo,rev(phi)),lwd=2,border=F, col=cols[2]))
with(y, lines(number.sire,fit,lwd=2,col="tomato"))

# Add legend
par(new=T)
legend("topleft", legend=c("Early (April)","Late (July)"),pch=c(19,17),col=c("dodgerblue","tomato"),cex=1.4,bty="n")

# Panel b, 2018
par(mar=c(3,4,3,1.5))
plot(c(0,20),c(0.6,0.8),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l",axes=F)
axis(side=1,at=seq(0,20,by=4),cex.axis=1.5)
axis(side=2,at=seq(0,1, by=0.05),las=1,cex.axis=1.5)
mtext("Expected Heterozygosity",side=2,line=4.5,cex=1.5,padj=0.4)
mtext("Number of sires",side=1,line=3,cex=1.5)
mtext("b)",side=3,adj=0,cex=1.5, line=0)

y1 <- eff.sire[eff.sire$time.season=="E" & eff.sire$year=="2018",]
with(y1, points(number.sire, Hexp,pch=19,col=cols[1]))
y2 <- eff.sire[eff.sire$time.season=="L" & eff.sire$year=="2018",]
with(y2, points(number.sire, Hexp,pch=17,col=cols[2]))

xmin <- range(y1$number.sire)
y <- newdata.he[newdata.he$year=="2018" &
                  newdata.he$time.season=="E" &
                  newdata.he$number.sire > xmin[1] &
                  newdata.he$number.sire < xmin[2], ]
with(y,polygon(c(number.sire, rev(number.sire)), c(plo,rev(phi)),lwd=2,border=F, col=cols[1]))
with(y, lines(number.sire,fit,lwd=2,col="dodgerblue"))

xmin <- range(y2$number.sire)
y <- newdata.he[newdata.he$year=="2018" &
                  newdata.he$time.season=="L" &
                  newdata.he$number.sire > xmin[1] &
                  newdata.he$number.sire < xmin[2], ]
with(y, polygon(c(number.sire, rev(number.sire)), c(plo,rev(phi)),lwd=2,border=F, col=cols[2]))
with(y, lines(number.sire,fit,lwd=2,col="tomato"))


# Pabel b, 2020
par(mar=c(3,2,3,1.5))
plot(c(0,20),c(0.6,0.8),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l",axes=F)
axis(side=1,at=seq(0,20,by=4),cex.axis=1.5)
axis(side=2,at=seq(0,1, by=0.05),las=1,cex.axis=1.5)
mtext("Number of sires",side=1,line=3,cex=1.5)
#mtext("2020",side=3,adj=.5,cex=1.5,line=-.5)

y1 <- eff.sire[eff.sire$time.season=="E" & eff.sire$year=="2020",]
with(y1, points(number.sire, Hexp,pch=19,col=cols[1]))
y2 <- eff.sire[eff.sire$time.season=="L" & eff.sire$year=="2020",]
with(y2, points(number.sire, Hexp,pch=17,col=cols[2]))

xmin <- range(y1$number.sire)
y <- newdata.he[newdata.he$year=="2020" &
                  newdata.he$time.season=="E" &
                  newdata.he$number.sire > xmin[1] &
                  newdata.he$number.sire < xmin[2], ]
with(y,polygon(c(number.sire, rev(number.sire)), c(plo,rev(phi)),lwd=2,border=F, col=cols[1]))
with(y, lines(number.sire,fit,lwd=2,col="dodgerblue"))

xmin <- range(y2$number.sire)
y <- newdata.he[newdata.he$year=="2020" &
                  newdata.he$time.season=="L" &
                  newdata.he$number.sire > xmin[1] &
                  newdata.he$number.sire < xmin[2], ]
with(y, polygon(c(number.sire, rev(number.sire)), c(plo,rev(phi)),lwd=2,border=F, col=cols[2]))
with(y, lines(number.sire,fit,lwd=2,col="tomato"))

