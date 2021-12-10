# metamorphosis analysis
met <- read.csv("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Supplmental Figures/metamorph_input.csv",
                       check.names=FALSE, header=T, 
                       na.strings=c("", "NA"), stringsAsFactors=FALSE, sep=",",quote="")
head(met)

library(lme4)
library(car)
library(MuMIn)
library(rms)
library(nlme)
library(multcomp)

met$temp<-as.factor(met$temp)
met$tube<-as.factor(met$replicate)
met$symbiont<-as.factor(met$symbiont)
met$int <- interaction(met$temp, met$symbiont)
met$polyp<-as.factor(met$polyp)

mod3<-glm(polyp ~temp+symbiont,data=met, family="binomial")
anova(mod3)
summary(mod3)
confint(mod3) # 95% CI for the coefficients
exp(coef(mod3)) # exponentiated coefficients
exp(confint(mod3)) # 95% CI for exponentiated coefficients
predict(mod3, type="response") # predicted values
cdplot(polyp ~temp, data=met) 

mod4<-glm(polyp ~temp*symbiont,data=met, family="binomial")
anova(mod4)
summary(mod4)

anova(mod3, mod4, test="Chisq")
drop1(mod4, test = "Chisq")

