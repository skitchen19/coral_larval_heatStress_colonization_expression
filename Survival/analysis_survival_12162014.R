####################################################
## Analysis of larva survival data  ################
## Duo Jiang #######################################
## Dec. 16, 2014 ###################################
####################################################

library(survival)
library(interval)
library(lme4)

#variables: larva.num, time1, time2, event, temp, sym, left, right
data=read.table(file.choose(), header=T)
data=read.table("survival.txt", header=T)
attach(data)

##################################################
## Interval-censored proportional hazard model ###
## tube effect incorporated as stratification ####
##################################################

#### 1.Weibull model #######
sreg.weibull.both<-survreg(Surv(time1,time2,type="interval2")~temp+sym+strata(replicate))
summary(sreg.weibull.both)
sreg.weibull.inter<-survreg(Surv(time1,time2,type="interval2")~temp*sym+strata(replicate))
summary(sreg.weibull.inter)
anova(sreg.weibull.both, sreg.weibull.inter)
## fold change and 95% confidence intervals ##
exp(cbind(-sreg.weibull.both$coef,-confint(sreg.weibull.both))["temp",])
exp(cbind(sreg.weibull.both$coef,confint(sreg.weibull.both))["sym",])
ConvertWeibull(sreg.weibull.both, conf.level = 0.95)

## plot fitted survival probabilities for different treat groups ##
plot(0:14, (0:14)/14*0.3+0.7,type='n', font=14, font.main=15, xaxs="i", yaxs="i",xlab="Days",cex.lab=1.2, ylab="Survival probability", 
	main="Temperature and symbiosis affect\n how survival probability decreases with time")
pct<-1:999/1000
ptime=numeric(999)
for (i in 1:24){
ptime<- ptime + predict(sreg.weibull.both,newdata=data.frame(temp=0, sym=0, replicate=i),type='quantile',p=pct)
}
ptime = ptime/24
lines(ptime,1-pct,col="black",lty=1, lwd=2)
ptime=numeric(999)
for (i in 1:24){
ptime<- ptime + predict(sreg.weibull.both,newdata=data.frame(temp=0, sym=1, replicate=i),type='quantile',p=pct)
}
ptime = ptime/24
lines(ptime,1-pct,col="darkorange",lty=1, lwd=2)
ptime=numeric(999)
for (i in 1:24){
ptime<- ptime + predict(sreg.weibull.both,newdata=data.frame(temp=1, sym=0, replicate=i),type='quantile',p=pct)
}
ptime = ptime/24
lines(ptime,1-pct,col="black",lty=2, lwd=2)
ptime=numeric(999)
for (i in 1:24){
ptime<- ptime + predict(sreg.weibull.both,newdata=data.frame(temp=1, sym=1, replicate=i),type='quantile',p=pct)
}
ptime = ptime/24
lines(ptime,1-pct,col="darkorange",lty=2, lwd=2)
legend(0, 0.77, "27, - sym", col="black", lty=1, cex=0.9, bty='n')
legend(0, 0.755, "32, - sym", col="black", lty=2, cex=0.9, bty='n')
legend(0, 0.74, "27, + sym", col="darkorange", lty=1, cex=0.9, bty='n')
legend(0, 0.725, "32, + sym",col="darkorange", lty=2, cex=0.9, bty='n')

## plot survival probabilities for individual tubes in the four treatment groups ##
pct<-1:999/1000
ptime=numeric(999)
par(mfrow=c(2,2))

plot(0:14, (0:14)/14*0.4+0.6,type='n',xlab="",font=14, font.main=15, xaxs="i", yaxs="i",cex.lab=1.2, ylab="Survival probability", 
	main="27, - sym")
for(i in 1:6){
ptime<-predict(sreg.weibull.both,newdata=data.frame(temp=0, sym=0, replicate=i),type='quantile',p=pct)
lines(ptime,1-pct,col="black",lty=1, lwd=2)
}

plot(0:14, (0:14)/14*0.4+0.6,type='n',xlab=" ", font=14, font.main=15, xaxs="i", yaxs="i",cex.lab=1.2, ylab="Survival probability", 
	main="27, + sym")
for(i in 7:12){
ptime<-predict(sreg.weibull.both,newdata=data.frame(temp=0, sym=1, replicate=i),type='quantile',p=pct)
lines(ptime,1-pct,col="darkorange",lty=1, lwd=2)
}

plot(0:14, (0:14)/14*0.4+0.6,type='n',xlab="Days",font=14, font.main=15, xaxs="i", yaxs="i",cex.lab=1.2, ylab="Survival probability", 
	main="32, - sym")
for(i in 13:18){
ptime<-predict(sreg.weibull.both,newdata=data.frame(temp=1, sym=0, replicate=i),type='quantile',p=pct)
lines(ptime,1-pct,col="black",lty=2, lwd=2)
}

plot(0:14, (0:14)/14*0.4+0.6,type='n',xlab="Days", font=14, font.main=15, xaxs="i", yaxs="i",cex.lab=1.2,ylab="Survival probability", 
	main="32, + sym")
for(i in 19:24){
ptime<-predict(sreg.weibull.both,newdata=data.frame(temp=1, sym=1, replicate=i),type='quantile',p=pct)
lines(ptime,1-pct,col="darkorange",lty=2, lwd=2)
}


#### As an alternative to a Weibull model, we can use the log-logistic model, #######
#### which gives very similar results ##############
sreg.loglogistic.both<-survreg(Surv(time1,time2,type="interval2")~temp+sym+strata(replicate), dist="loglogistic")
summary(sreg.loglogistic.both)
sreg.loglogistic.inter<-survreg(Surv(time1,time2,type="interval2")~temp*sym+strata(replicate), dist="loglogistic")
summary(sreg.loglogistic.inter)
anova(sreg.loglogistic.both, sreg.loglogistic.inter)


#####################################
## Survive beyond Day 14 or not? ####
#####################################
## creating binary response ##
missing = numeric(length(time1))
surv14 = numeric(length(time1))+1
for (i in 1:length(time1)){
	if(!is.na(time2[i])){
		if(time2[i]<=14)surv14[i]=0
	}
	if(time1[i]<14 && is.na(time2[i])) missing[i]=1
}

### glmm: tube effects incorporated as random effects ###
data.2=data.frame(surv14.2=surv14[missing==0], temp.2=temp[missing==0], sym.2=sym[missing==0], replicate.2=replicate[missing==0])
fit.glmer.inter=glmer(surv14.2~temp.2*sym.2+(1|(replicate.2)), family=binomial(logit),data=data.2)
summary(fit.glmer.inter)# model including interaction
fit.glmer.both=glmer(surv14.2~temp.2+sym.2+(1|(replicate.2)), family=binomial(logit),data=data.2)
summary(fit.glmer.both) # model excluding interaction (the better fitting model)
anova(fit.glmer.both, fit.glmer.inter) # interaction is not significant

## effects of temperature and symbiosis on odds of survival
# temperature effect
(temp.effect=exp(coef(summary(fit.glmer.both))[2,1]))
# 95% confidence interval
c(temp.effect*exp(-1.96*coef(summary(fit.glmer.both))[2,2]),temp.effect*exp(1.96*coef(summary(fit.glmer.both))[2,2]))
# symbiosis effect
(sym.effect=exp(coef(summary(fit.glmer.both))[3,1]))
# 95% confidence interval
c(sym.effect*exp(-1.96*coef(summary(fit.glmer.both))[3,2]),sym.effect*exp(1.96*coef(summary(fit.glmer.both))[3,2]))


############################### ##################################
## Can a larva survive beyond a given day (0.5, 1, 2, 3, 5, 7)? ##
##################################################################
c = 3 # specify the day used as the cutoff 
missing = numeric(length(time1))
survc = numeric(length(time1))+1
for (i in 1:length(time1)){
	if(!is.na(time2[i])){
		if(time2[i]<=c)survc[i]=0
	}
	if(!is.na(time1[i])){
		if(time1[i]<c && is.na(time2[i])) missing[i]=1
	}
}
data.3=data.frame(survc.3=survc[missing==0], temp.3=temp[missing==0], sym.3=sym[missing==0], replicate.3=replicate[missing==0])
fit.glmer.inter.3=glmer(survc.3~temp.3*sym.3+(1|(replicate.3)), family=binomial(logit),data=data.3)
summary(fit.glmer.inter.3)
fit.glmer.both.3=glmer(survc.3~temp.3+sym.3+(1|(replicate.3)), family=binomial(logit),data=data.3)
summary(fit.glmer.both.3)
anova(fit.glmer.both.3, fit.glmer.inter.3) # interaction is not significant for any given cutoff


