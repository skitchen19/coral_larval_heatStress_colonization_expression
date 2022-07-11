library(lme4)
library(car)
library(MuMIn)
library(rms)
library(nlme)
library(multcomp)

col<-read.csv("C:/Users/Sheila's Comp/Documents/GitHub/coral_larval_heatStress_colonization_expression/Colonization/col_success_algal_density_combined.csv", head=T)
attach(col)
names(col)

col$temp<-as.factor(col$temperature)
#col$day<-as.numeric(col$Day)
col$day<-as.factor(col$Day)
col$tube<-as.factor(col$replicate_tube)
col$int <- interaction(col$temperature, col$Day)

mod1=glmer(infection~(1|tube)+day*temp, family=binomial, data=col)

summary(mod1)
r.squaredGLMM(mod1)

exp(1.34)

mod2=glmer(infection~(1|tube)+day+temp, family=binomial, data=col)
summary(mod2)
resid2=residuals(mod2)

anova(mod2, mod1, test="Chisq")


#generate frequency/proportion table
with(col, table(infection.1., Day, temperature))->tabCount
new<-ftable(tabCount, col.vars="infection.1.")
prop.table(new,1)*100 #compares across the row

#plot frequencies
per<-read.csv("C:/Users/Sheila's Comp/Documents/GitHub/coral_larval_heatStress_colonization_expression/Colonization/colonization_avg.csv", head=T)
attach(per)
names(per)

par(mar=c(5,8,5,2),bty="l",lwd=2)
plot(day,average,type="n", cex.lab=1.7,font.lab=2,ylim=c(0,100), xlim=c(0,14.25), yaxs="i", ylab=" ", xlab="Days", las=1, cex.axis=1.2,cex=1.2, font.axis=2)
title(ylab = "Average % colonized", cex.lab = 1.7, font.lab=2,
      line = 4.5)


segments(day[1:10], average[1:10]-stder[1:10],day, average[1:10]+stder[1:10], lwd=2, col="black")
epsilon = 0.1
segments(day[1:10]-epsilon,average[1:10]-stder[1:10],day[1:10]+epsilon,average[1:10]-stder[1:10],lwd=2, col="black")
segments(day[1:10]-epsilon,average[1:10]+stder[1:10],day[1:10]+epsilon,average[1:10]+stder[1:10],lwd=2, col="black")

segments(day[11:20], average[11:20]-stder[11:20],day[11:20], average[11:20]+stder[11:20], lwd=2, col="sienna1")
epsilon = 0.1
segments(day[11:20]-epsilon,average[11:20]-stder[11:20],day[11:20]+epsilon,average[11:20]-stder[11:20],lwd=2, col="sienna1")
segments(day[11:20]-epsilon,average[11:20]+stder[11:20],day[11:20]+epsilon,average[11:20]+stder[11:20],lwd=2, col="sienna1")

lines(day[1:10], jitter(average[1:10]), col="black", bg="white", pch=19, type="o", lty=2,lwd=2, cex=1.5)
lines(day[11:20],jitter(average[11:20]),col="sienna1", bg="white", pch=19, type="o", lwd=2, cex=1.5, lty=2)

legend("bottomright", inset=.05, pch=c(19), pt.bg=c("white","white"), col=c("black", "sienna1"),
  	c("27 °C","32 °C"), horiz=TRUE, lty=c(2,2),lwd=2, cex=1.3,bty="n",text.font=2)

ggplot(col3, aes(x=day, y=average, color=temp)) + 
    geom_errorbar(aes(ymin=average-stder, ymax=average+stder), width=.1) +
    geom_line()
    geom_point()

###############
#Algal density#
library(lme4)
library(lattice)
library(ggplot2)
library(MASS)
library(dplyr)

col<-col %>% filter(!day == 5)    

col$temp<-as.factor(col$temperature)
col$day<-as.numeric(col$Day)
#col$day<-as.factor(col$Day)
col$count<-as.numeric(col$Algal_count)

ggplot(col,aes(x=col$day,y=col$count,colour=col$temp)) +
  geom_point() + 
  theme_bw() + theme(panel.margin=unit(0,"lines")) +
  scale_color_manual(values=c("#3B9AB2","#F21A00"))

mod3<- glm(count~temp*day, family="quasipoisson", data=col)
mod4<- glm(count~temp+day, family="quasipoisson", data=col)

summary(mod3)
summary(mod4)

anova(mod4,mod3, test="Chisq")

r.squaredGLMM(mod3)

sum(resid(mod3, type = "pearson")^2)/mod3$df.res
sum(resid(mod4, type = "pearson")^2)/mod4$df.res

aov1<- glmmPQL(count~temp+day,random=~1|tube, family="quasipoisson", data=col)
summary(aov1)

aov2<- glmmPQL(count~temp*day,random=~1|tube, family="quasipoisson", data=col)
summary(aov2)

##To break axis##
library(plotrix)

ag<-read.csv("C:/Users/Sheila's Comp/Documents/GitHub/coral_larval_heatStress_colonization_expression/Colonization/algal_avg.csv", head=T)
attach(ag)

par(mar=c(5,8,5,2),bty="l",lwd=2)
plot(dy,avg, type="n", cex.lab=1.7,font.axis=2, font.lab=2, xlim=c(0,14), xlab="Days", ylab=" ",cex.axis=1.2, yaxt="n",ylim=c(0,700))
title(ylab = "Average # of symbionts per larva", cex.lab = 1.7, font.lab=2,
      line = 4.5)


ygap <- ifelse(avg > 500, avg-500, avg)
yat <- pretty(ygap)
yat <- yat[yat!=500]
ylab <- ifelse(yat > 500, yat+500, yat)
axis(2,at=yat, labels=ylab,las=1, font.axis=2, cex.axis=1.2)
abline(h=seq(500,510,1), col="white")
axis.break(2,500,style="slash") 

segments(dy[1:9], ygap[1:9]-stder[1:9],dy[1:9], ygap[1:9]+stder[1:9], lwd=2, col="black")
epsilon = 0.1
segments(dy[1:9]-epsilon,ygap[1:9]-stder[1:9],dy[1:9]+epsilon,ygap[1:9]-stder[1:9],lwd=2, col="black")
segments(dy[1:9]-epsilon,ygap[1:9]+stder[1:9],dy[1:9]+epsilon,ygap[1:9]+stder[1:9],lwd=2, col="black")
lines(dy[1:9], ygap[1:9], col="black",bg="white", pch=19, type="o", lwd=2, cex=1.2, lty=2)

segments(dy[10:18], ygap[10:18]-stder[10:18],dy[10:18], ygap[10:18]+stder[10:18], lwd=2, col="sienna1")
segments(dy[10:18]-epsilon,avg[10:18]-stder[10:18],dy[10:18]+epsilon,avg[10:18]-stder[10:18],lwd=2, col="sienna1")
segments(dy[10:18]-epsilon,avg[10:18]+stder[10:18],dy[10:18]+epsilon,avg[10:18]+stder[10:18],lwd=2, col="sienna1")
lines(dy[10:18],ygap[10:18],col="sienna1",bg="white", pch=19, type="o", lwd=2, cex=1.2, lty=2)


legend("topleft", inset=.05,pch=19, pt.bg=c("white", "white"), col=c("black", "sienna1"),
  	c("27 °C","32 °C"), horiz=TRUE, cex=1.3, lty=c(2,2),lwd=2, bty="n", text.font=2)