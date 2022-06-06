setwd("~/Activité Professionnelle/CNRS 2018-2020/Publications/2 - Mismatch Energétique")

rm(list=ls())

library(AICcmodavg)
library(bbmle)
library(cowplot)
library(dataframes2xls)
library(data.table)
library(deSolve)  
library(dplyr)
library(emdbook)
library(FME)
library(forcats)
library(ggplot2)
library(grid)
library(gridExtra)
library(lme4)
library(lmtest)
library(nlme)
library(plyr)
library(quantmod)
library(reshape2)
library(zoo)

###############################################################
###############################################################
##### EXTRACTION OF THE PARAMETERS FROM EXPERIMENTAL DATA #####
###############################################################
###############################################################

# Dataframe with respiration rate and ingestion rate (µgC/jour)
Tab=read.table("Data_Mismatch.txt", h=T, dec=".") 
names(Tab)
summary(Tab)
Tab=na.omit(Tab)

# Considerate individuals as a qualitative variable
Tab$Indiv=as.factor(Tab$Indiv)

# Suppress rows with negative values
Tab=Tab[Tab$Respi>0,]
min(Tab$Respi)
Tab=Tab[Tab$Nutri>0,]
min(Tab$Nutri)

# Suppress outlier individuals
Tab=Tab[!Tab$Indiv=="12",]
Tab=Tab[!Tab$Indiv=="30",]
Tab=Tab[!Tab$Indiv=="76",]
Tab=Tab[!Tab$Indiv=="78",]

# Check for the normality of data
shapiro.test(log(Tab$Respi))
qqnorm(log(Tab$Respi))
qqline(log(Tab$Respi))
hist(log(Tab$Respi))

shapiro.test(log(Tab$Nutri))
qqnorm(log(Tab$Nutri))
qqline(log(Tab$Nutri))
hist(log(Tab$Nutri))

# Check for homocesdasticity of data
bartlett.test(Tab$Respi~Tab$Temp)
bartlett.test(Tab$Nutri~Tab$Temp)

# Convert temperature into the Boltzmann term (°K)
InvTem=1/((Tab$Temp+273.15)*8.62*10^(-5))
InvTemC=(InvTem-mean(InvTem))


#################################
### Metabolic rate parameters ###
#################################

# Comparison of linear and quadratic models
ModelLM1=lm(log(Respi)~log(Masse)+InvTemC, data=Tab)
ModelPOLY1=glm(log(Respi)~log(Masse)+ poly(InvTemC, degree=2, raw=F), data=Tab)
anova(ModelLM1,ModelPOLY1)
lrtest(ModelLM1,ModelPOLY1)
AIC(ModelLM1,ModelPOLY1)

# Linear model metabolic rate parameters
ModelMLinear=lm(log(Respi)~log(Masse)+InvTemC, data=Tab)
summary(ModelMLinear)
confint(ModelMLinear)
anova(ModelMLinear)
coef(ModelMLinear)

# Quadratic model metabolic rate parameters
ModelMQuadra=glm(log(Respi)~log(Masse)+poly(InvTemC, degree=2, raw=T), data=Tab)
summary(ModelMQuadra)
confint(ModelMQuadra)
anova(ModelMQuadra)
coef(ModelMQuadra)

# Check for the normality of residuals
shapiro.test(ModelMQuadra$residuals)
qqnorm(ModelMQuadra$residuals)
qqline(ModelMQuadra$residuals)
hist(ModelMQuadra$residuals)

# Exponential plot
PlotM=ggplot(data=Tab, aes(x=Temp, y=log(Respi))) +
  geom_smooth(method="lm", formula=y~poly(x,2), color="black", linetype="solid", size=1, se=F) +     
  geom_point(stat="identity", color="turquoise3", size=3) + ylim(1.0,5.0) + xlim(5,21) +
  ylab(expression('Ln Routine metabolic rate'~'('*µg~C~d^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  annotate("text", x=8.50, y=5.00, label=expression('y = - 0.21'*x^2*~'- 0.67x + 2.42'), size=6) +
  annotate("text", x=8.50, y=4.70, label=expression(italic(P) *~'< 0.001'), size=6) +
  annotate("text", x=8.50, y=4.40, label=expression(italic(R^2) *~'= 0.83'), size=6)


# Functions for metabolic rate (µgC/day)
Boltz=8.62*10^(-5)
Temp=seq(0,30,0.1)
MeanInvTem=mean(InvTem)        # Mean inverse temperature
MeanMass=mean(Tab$Masse)       # Mean Gammarus body mass

# With linear function
Meta=function(Temp,Mass){exp(2.34492)*Mass^0.58414*exp(-0.70152*(1/((Temp+273.15)*Boltz)-MeanInvTem))}
Meta(20,MeanMass)

# With quadratic function
MetaQuadra=function(Temp,Mass){exp(2.41599)*Mass^0.62308*exp(-0.66731*(1/((Temp+273.15)*Boltz)-MeanInvTem))*exp(-0.21153*(1/((Temp+273.15)*Boltz)-MeanInvTem)^2)}
MetaQuadra(20,MeanMass)




#################################
### Ingestion rate parameters ###
#################################

# Comparison of linear and quadratic models
ModelLM2=lm(log(Nutri)~log(Masse)+InvTemC, data=Tab)
ModelPOLY2=glm(log(Nutri)~log(Masse)+poly(InvTemC, degree=2, raw=F), data=Tab)
anova(ModelLM2,ModelPOLY2)
lrtest(ModelLM2,ModelPOLY2)
AIC(ModelLM2,ModelPOLY2)

# Linear model ingestion rate parameters
ModelILinear=lm(log(Nutri)~log(Masse)+InvTemC, data=Tab)
summary(ModelILinear)
confint(ModelILinear)
anova(ModelILinear)
coef(ModelILinear)

# Quadratic model ingestion rate parameters
ModelIQuadra=glm(log(Nutri)~log(Masse)+poly(InvTemC, degree=2, raw=T), data=Tab)
summary(ModelIQuadra)
confint(ModelIQuadra)
anova(ModelIQuadra)
coef(ModelIQuadra)

# Check for the normality of residuals
shapiro.test(ModelIQuadra$residuals)
qqnorm(ModelIQuadra$residuals)
qqline(ModelIQuadra$residuals)
hist(ModelIQuadra$residuals)

# Exponential plot
PlotI=ggplot(data=Tab, aes(x=Temp, y=log(Nutri))) +
  geom_smooth(method="lm", formula=y~poly(x,2), color="black", linetype="solid", size=1, se=F) + 
  geom_point(stat="identity", color="tomato1", size=3) + ylim(4.0,8.0) + xlim(5,21) +
  ylab(expression('Ln Ingestion rate'~'('*µg~C~d^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  annotate("text", x=8.50, y=8.00, label=expression('y = - 0.19'*x^2*~'- 0.32x + 5.27'), size=6) +
  annotate("text", x=8.50, y=7.70, label=expression(italic(P) *~'< 0.001'), size=6) +
  annotate("text", x=8.50, y=7.40, label=expression(italic(R^2) *~'= 0.42'), size=6)


# Functions for ingestion rate (µgC/day)
Boltz=8.62*10^(-5)
Temp=seq(0,30,0.1)
MeanInvTem=mean(InvTem)        # Mean inverse temperature
MeanMass=mean(Tab$Masse)       # Mean Gammarus body mass

# With linear function
Ing=function(Temp,Mass){exp(5.12999)*Mass^0.83128*exp(-0.28227*(1/((Temp+273.15)*Boltz)-MeanInvTem))}
Ing(20,MeanMass)

# With quadratic function
IngQuadra=function(Temp,Mass){exp(5.26814)*Mass^0.81654*exp(-0.31876*(1/((Temp+273.15)*Boltz)-MeanInvTem))*exp(-0.18909*(1/((Temp+273.15)*Boltz)-MeanInvTem)^2)}
IngQuadra(20,MeanMass)




###################################################
### Panel for metabolic rate and ingestion rate ###
###################################################

tiff('Metabolic and Ingestion Rates.tiff', units="in", width=15, height=8, res=300)
Panel=plot_grid(PlotM, PlotI, align="h", vjust=1, nrow=1, ncol=2)
Xaxis=textGrob(expression('Temperature (°C)'), gp=gpar(fontface="bold", fontsize=20))
grid.arrange(arrangeGrob(Panel, bottom=Xaxis))
dev.off()



###############################
### Assimilation efficiency ###
###############################

# Logistic model assimilation efficiency
AssimQuadra=function(Temp){(exp(-0.84730)*exp(0.16400*((Temp+273.15)-285.65)/(Boltz*285.65*(Temp+273.15))))/(1+(exp(-0.84730)*exp(0.16400*((Temp+273.15)-285.65)/(Boltz*285.65*(Temp+273.15)))))}
AssimQuadra(20)

# Calculate assimilation efficiency 
Tab$Assim=AssimQuadra(Tab$Temp)

# Check for the normality of data
shapiro.test(Tab$Assim)
qqnorm(Tab$Assim)
qqline(Tab$Assim)
hist(Tab$Assim)

# Linear model energetic efficiency parameters
ModelA=lm(Assim~Temp, data=Tab)
summary(ModelA)
confint(ModelA)
anova(ModelA)
coef(ModelA)




############################
### Energetic efficiency ###
############################

# Calculate energetic efficiency
Tab$Energ=(Tab$Nutri/Tab$Respi)*Tab$Assim
Tab=Tab[Tab$Energ<30,]

# Check for the normality of data
shapiro.test(Tab$Energ)
qqnorm(Tab$Energ)
qqline(Tab$Energ)
hist(Tab$Energ)

# Linear model energetic efficiency parameters
ModelE=lm(Energ~Temp, data=Tab)
summary(ModelE)
confint(ModelE)
anova(ModelE)
coef(ModelE)

# Energetic efficiency plot
tiff('Energetic Efficiencies.tiff', units="in", width=8, height=8, res=300)
ggplot(data=Tab, aes(x=Temp, y=Energ)) +
  geom_smooth(method="lm", formula=y~x, color="black", linetype="solid", size=1, se=F) +
  geom_point(stat="identity", color="chartreuse3", size=3) + ylim(0,20) + xlim(5,22) +
  ylab(expression('Energetic efficiency')) + xlab(expression('Temperature (°C)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=20)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  annotate("text", x=8.50, y=20.0, label=expression('y = -0.26x + 11.21'), size=6) +
  annotate("text", x=8.50, y=18.5, label=expression(italic(P) *~'< 0.001'), size=6) +
  annotate("text", x=8.50, y=17.1, label=expression(italic(R^2) *~'= 0.19'), size=6)
dev.off()




##########################################################
##########################################################
##### PARAMETRIZATION OF THE CONSUMER-RESOURCE MODEL #####
##########################################################
##########################################################

#############################
### Attack rate parameter ###
#############################

# Duration of the feeding experiment (days)
Duration=2

# Mean Gammarus body mass (mg)
MeanMass=mean(Tab$Masse)

# Mean initial leaf discs mass: 10.25 ± 0.68 mg
LeafMass=10.25*1000

# Converted from mg to in µgC with a factor: 0.45
LeafMass=LeafMass*0.45

# Calculation of K-value and attack rate (quadratic model)
K=-log(1-IngQuadra(Temp,MeanMass)*2/LeafMass)
Attack=K/2

# Attack rate function based on exponential decay (quadratic model)
AttackQuadra=function(Temp,Mass){-log(1-IngQuadra(Temp,Mass)*2/LeafMass)/2}
AttackQuadra(20,MeanMass,2)




###############################
### Handling time parameter ###
###############################

# Handling time function based on exponential decay (quadratic model)
HandlingQuadra=function(Temp,Mass){1/(IngQuadra(Temp,Mass)/1000)}
HandlingQuadra(20,MeanMass)




#########################################
### Assimilation efficiency parameter ###
#########################################

# Assimilation efficiency function based on exponential decay (quadratic model)
AssimQuadra=function(Temp){(exp(-0.84730)*exp(0.16400*((Temp+273.15)-285.65)/(Boltz*285.65*(Temp+273.15))))/(1+(exp(-0.84730)*exp(0.16400*((Temp+273.15)-285.65)/(Boltz*285.65*(Temp+273.15)))))}
AssimQuadra(20)




######################################
### Energetic efficiency parameter ###
######################################

EnergQuadra=function(Temp,MeanMass){
  LeafMass=10.25*0.45*1000                      # Leaf litter mass of the feeding (in µgC)
  Metabo=MetaQuadra(Temp,MeanMass)              # Gammarus metabolic rate (in µgC/day)
  a=AttackQuadra(Temp,MeanMass)                 # Gammarus attack rate (in µgC/day)
  h=HandlingQuadra(Temp,MeanMass)               # Gammarus handling time (in 1/day)
  A=AssimQuadra(Temp)                           # Gammarus assimilation efficiency
  Inges=(a*LeafMass)/(1+a*h*LeafMass)           # Gammarus consumption rate (in µgC/day)
  Energ=A*Inges/Metabo                          # Gammarus energetic efficiency
  return(Energ)
}
EnergQuadra(20,MeanMass)




###############################
### Leaf decomposition rate ###
###############################

# Reference temperature
TRef=10

# Function for the leaf litter microbial decomposition rate
DecompLeaf=function(Temp){0.00956*exp(-0.37000*(1/((Temp+273.15)*Boltz)-TRef))} 

# Plot of leaf litter decomposition rate
plot(DecompLeaf(Temp)*10^7~Temp, ylab=expression('Leaf decomposition rate '~'('*10^7*mg~C~day^-1*')'), xlab="Temperature (°C)",
     las=1, cex.lab=1, cex.axis=1, pch=16, cex=0.1, lwd=1, col="black")




#############################
### Leaf respiration rate ###
#############################

# Reference temperature
TRef=10

# Function for the leaf litter respiration rate
RespLeaf=function(Temp){2.5/0.45*10^-3*exp(-0.65000*(1/((Temp+273.15)*Boltz)-TRef))} 
RespLeaf(20)

# Plot of leaf litter respiration rate
plot(RespLeaf(Temp)*10^10~Temp, ylab=expression('Leaf respiration rate '~'('*10^10*mg~C~day^-1*')'), xlab="Temperature (°C)",
     las=1, cex.lab=1, cex.axis=1, pch=16, cex=0.1, lwd=1, col="black")




###############################
### Consumer-resource model ###
###############################

GammLeafModel=function(Temp,GammMass,Leaf,Gamm){
  Nutri=function(Times,state,parameters){
    with(as.list(c(state,parameters)), {
      Fr=a*L/(1+a*h*L)                    # Holling type II functional response 
      dL=-G*Fr-Dl*L                       # Biomass changes of leaf litter stock
      dG=G*(A*Fr-M)                       # Biomass changes of Gammarus population 
      list(c(dL,dG))
    })
  }
  
# Define the years   
y1=365; y2=2*y1; y3=3*y1; y4=4*y1; y5=5*y1; y6=6*y1

# Time points to trigger litter fall
FallTime=c(seq(1,15), seq(y1+1,y1+15),seq(y2+1,y2+15),seq(y3+1,y3+15),seq(y4+1,y4+15),seq(y5+1,y5+15),seq(y6+1,y6+15))

# Leaf litter fall function
FallFunc=function(time,y,parms) {
  with(as.list(c(y,parms)), {
    return(c(L+Leaf,G))
  })
}

# Model parameters
Parameters=c(
  M=MetaQuadra(Temp,GammMass)/1000,              # Gammarus metabolic rate (in mgC/day)
  a=AttackQuadra(Temp,GammMass),                 # Gammarus attack rate (in mgC/day)
  h=HandlingQuadra(Temp,GammMass),               # Gammarus handling time (in 1/day)
  A=AssimQuadra(Temp),                           # Gammarus assimilation efficiency
  Dl=DecompLeaf(Temp),                           # Leaf microbial decomposition (in 1/day)
  Rl=RespLeaf(Temp))                             # Leaf respiration (in mgC/mgleaf/day)


# Time and starting conditions
Times=seq(0,365*7,by=1)        # Time in days for 7 years
State=c(L=Leaf,G=Gamm)         # Starting biomasses (in g/m2)
  

# Model output
Out=as.data.frame(ode(time=Times, func=Nutri, y=State, events=list(func=FallFunc, time=FallTime), parms=Parameters))
return(Out)
} 




####################################
### TEMPERATURE-SIZE RULE MODELS ###
####################################

Cf=6.5                                 # Average conversion factor from dry to fresh mass

# Average TSR response
TSRA=function(Temp,Mass){
  Slope=-3.90-0.53*log10(Mass)         # Slope of change in mass per carbon          
  Dx=log(1+Slope/100)                  # Proportion of change in mass per C
  Const=exp(log(Mass)-12.5*Dx)         # Constant of change in mass at 12.5°C
  DryMass=Const*exp(Dx*(Temp))         # Dry body mass (mg)
  FreshMass=DryMass/Cf                 # Fresh body mass (mg)
  return(DryMass)
}

# Maximum TSR response
TSRM=function(Temp,Mass){
  Slope=-8.0                           # Slope of change in mass per carbon          
  Dx=log(1+Slope/100)                  # Proportion of change in mass per C
  Const=exp(log(Mass)-12.5*Dx)         # Constant of change in mass at 12.5°C
  DryMass=Const*exp(Dx*(Temp))         # Dry body mass (mg)
  FreshMass=DryMass/Cf                 # Fresh body mass (mg)
  return(DryMass)
}




############################################################
############################################################
##### TESTING DIFFERENT LITTER DECOMPOSITION SCENARIOS #####
############################################################
############################################################

#####################################
### SCENARIO 0: STANDARD SCENARIO ###
#####################################

# Leaf fall = 300 gC/m2/an = 300 000 mgC/m2/an
# Gammarus density = 30 mgDM/m2 = 15 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
tiff('Population Dynamics SD.tiff', units="in", width=15, height=8, res=1000)
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5SD=GammLeafModel(Temp=5,GammMass=4.26,Leaf=300000/Days,Gamm=15)
Test5SD=Test5SD %>% mutate(L=if_else(L<10^-3, 0, L))
Test5SD=Test5SD %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5SD[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1), lwd=c(1.5,2))
mtext("5°C", cex=1.5, side=3, line=-3, at=300)

Test10SD=GammLeafModel(Temp=10,GammMass=4.26,Leaf=300000/Days,Gamm=15)
Test10SD=Test10SD %>% mutate(L=if_else(L<10^-3, 0, L))
Test10SD=Test10SD %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10SD[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10°C", cex=1.5, side=3, line=-3, at=300)

Test15SD=GammLeafModel(Temp=15,GammMass=4.26,Leaf=300000/Days,Gamm=15)
Test15SD=Test15SD %>% mutate(L=if_else(L<10^-3, 0, L))
Test15SD=Test15SD %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15SD[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15°C", cex=1.5, side=3, line=-3, at=300)

Test20SD=GammLeafModel(Temp=20,GammMass=4.26,Leaf=300000/Days,Gamm=15)
Test20SD=Test20SD %>% mutate(L=if_else(L<10^-3, 0, L))
Test20SD=Test20SD %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20SD[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20°C", cex=1.5, side=3, line=-3, at=300)

Test25SD=GammLeafModel(Temp=25,GammMass=4.26,Leaf=300000/Days,Gamm=15)
Test25SD=Test25SD %>% mutate(L=if_else(L<10^-3, 0, L))
Test25SD=Test25SD %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25SD[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25°C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))
dev.off()


# Calculate population metabolism and leaf ingestion
Test5SD$M=(MetaQuadra(5,4.26)/1000)*Test5SD[,3]
Test5SD$I=(AttackQuadra(5,4.26)*Test5SD[,2]/(1+AttackQuadra(5,4.26)*1/(IngQuadra(5,4.26)/1000)*Test5SD[,2]))*0.30*Test5SD[,3]
Test5SD=Test5SD %>% mutate(M=if_else(M<0, 0, M))
Test5SD=Test5SD %>% mutate(I=if_else(I<0, 0, I))
Test10SD$M=(MetaQuadra(10,4.26)/1000)*Test10SD[,3]
Test10SD$I=(AttackQuadra(10,4.26)*Test10SD[,2]/(1+AttackQuadra(10,4.26)*1/(IngQuadra(10,4.26)/1000)*Test10SD[,2]))*0.30*Test10SD[,3]
Test10SD=Test10SD %>% mutate(M=if_else(M<0, 0, M))
Test10SD=Test10SD %>% mutate(I=if_else(I<0, 0, I))
Test15SD$M=(MetaQuadra(15,4.26)/1000)*Test15SD[,3]
Test15SD$I=(AttackQuadra(15,4.26)*Test15SD[,2]/(1+AttackQuadra(15,4.26)*1/(IngQuadra(15,4.26)/1000)*Test15SD[,2]))*0.30*Test15SD[,3]
Test15SD=Test15SD %>% mutate(M=if_else(M<0, 0, M))
Test15SD=Test15SD %>% mutate(I=if_else(I<0, 0, I))
Test20SD$M=(MetaQuadra(20,4.26)/1000)*Test20SD[,3]
Test20SD$I=(AttackQuadra(20,4.26)*Test20SD[,2]/(1+AttackQuadra(20,4.26)*1/(IngQuadra(20,4.26)/1000)*Test20SD[,2]))*0.30*Test20SD[,3]
Test20SD=Test20SD %>% mutate(M=if_else(M<0, 0, M))
Test20SD=Test20SD %>% mutate(I=if_else(I<0, 0, I))
Test25SD$M=(MetaQuadra(25,4.26)/1000)*Test25SD[,3]
Test25SD$I=(AttackQuadra(25,4.26)*Test25SD[,2]/(1+AttackQuadra(25,4.26)*1/(IngQuadra(25,4.26)/1000)*Test25SD[,2]))*0.30*Test25SD[,3]
Test25SD=Test25SD %>% mutate(M=if_else(M<0, 0, M))
Test25SD=Test25SD %>% mutate(I=if_else(I<0, 0, I))

# Combine the dataframes
TestSD=rbind(Test5SD,Test10SD,Test15SD,Test20SD,Test25SD)
TestSD$Temp=rep(c("5","10","15","20","25"), each=2556)
TestSD$Temp=factor(TestSD$Temp, levels=unique(TestSD$Temp))
colnames(TestSD)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestSD=subset(TestSD, !Time==2555); TestSD$Year=rep(rep(1:(nrow(TestSD)/(365*5)), each=365),5)
CutSD=as.data.frame(setDT(TestSD)[TestSD[, tail(.I, -16), by=Year]$V1]); CutSD=CutSD[!CutSD$Year=="1",]
MeanL=as.data.frame(setDT(CutSD)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutSD)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLSD=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGSD=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestSD=subset(TestSD, !Time==2555); TestSD$Year=rep(rep(1:(nrow(TestSD)/(365*5)), each=365),5)
CutSD=as.data.frame(setDT(TestSD)[TestSD[, tail(.I, -16), by=Year]$V1]); CutSD=CutSD[!CutSD$Year=="1",]
ThLSD=CutSD[CutSD$L>60000,]; ThGSD=CutSD[CutSD$G>5000,]
TimeLSD=as.data.frame(setDT(ThLSD)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGSD=as.data.frame(setDT(ThGSD)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLSD=as.data.frame(setDT(TimeLSD)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGSD=as.data.frame(setDT(TimeGSD)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLSD=as.data.frame(setDT(TestSD)[, .(MaxLSD=findPeaks(L), MinLSD=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGSD=as.data.frame(setDT(TestSD)[, .(MaxGSD=findPeaks(G), MinGSD=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5SD=c(CycleLSD[1,2]:CycleLSD[1,3],CycleLSD[2,2]:CycleLSD[2,3],CycleLSD[3,2]:CycleLSD[3,3],CycleLSD[4,2]:CycleLSD[4,3],CycleLSD[5,2]:CycleLSD[5,3],CycleLSD[6,2]:CycleLSD[6,3],CycleLSD[7,2]:CycleLSD[7,3])+t0
CycleL5SD=c(rep("A",length(CycleLSD[1,2]:CycleLSD[1,3])),rep("B",length(CycleLSD[2,2]:CycleLSD[2,3])),rep("C",length(CycleLSD[3,2]:CycleLSD[3,3])),rep("D",length(CycleLSD[4,2]:CycleLSD[4,3])),rep("E",length(CycleLSD[5,2]:CycleLSD[5,3])),rep("F",length(CycleLSD[6,2]:CycleLSD[6,3])),rep("G",length(CycleLSD[7,2]:CycleLSD[7,3])))
CutL10SD=c(CycleLSD[8,2]:CycleLSD[8,3],CycleLSD[9,2]:CycleLSD[9,3],CycleLSD[10,2]:CycleLSD[10,3],CycleLSD[11,2]:CycleLSD[11,3],CycleLSD[12,2]:CycleLSD[12,3],CycleLSD[13,2]:CycleLSD[13,3],CycleLSD[14,2]:CycleLSD[14,3])+t1
CycleL10SD=c(rep("A",length(CycleLSD[8,2]:CycleLSD[8,3])),rep("B",length(CycleLSD[9,2]:CycleLSD[9,3])),rep("C",length(CycleLSD[10,2]:CycleLSD[10,3])),rep("D",length(CycleLSD[11,2]:CycleLSD[11,3])),rep("E",length(CycleLSD[12,2]:CycleLSD[12,3])),rep("F",length(CycleLSD[13,2]:CycleLSD[13,3])),rep("G",length(CycleLSD[14,2]:CycleLSD[14,3])))
CutL15SD=c(CycleLSD[15,2]:CycleLSD[15,3],CycleLSD[16,2]:CycleLSD[16,3],CycleLSD[17,2]:CycleLSD[17,3],CycleLSD[18,2]:CycleLSD[18,3],CycleLSD[19,2]:CycleLSD[19,3],CycleLSD[20,2]:CycleLSD[20,3],CycleLSD[21,2]:CycleLSD[21,3])+t2
CycleL15SD=c(rep("A",length(CycleLSD[15,2]:CycleLSD[15,3])),rep("B",length(CycleLSD[16,2]:CycleLSD[16,3])),rep("C",length(CycleLSD[17,2]:CycleLSD[17,3])),rep("D",length(CycleLSD[18,2]:CycleLSD[18,3])),rep("E",length(CycleLSD[19,2]:CycleLSD[19,3])),rep("F",length(CycleLSD[20,2]:CycleLSD[20,3])),rep("G",length(CycleLSD[21,2]:CycleLSD[21,3])))
CutL20SD=c(CycleLSD[22,2]:CycleLSD[22,3],CycleLSD[23,2]:CycleLSD[23,3],CycleLSD[24,2]:CycleLSD[24,3],CycleLSD[25,2]:CycleLSD[25,3],CycleLSD[26,2]:CycleLSD[26,3],CycleLSD[27,2]:CycleLSD[27,3],CycleLSD[28,2]:CycleLSD[28,3])+t3
CycleL20SD=c(rep("A",length(CycleLSD[22,2]:CycleLSD[22,3])),rep("B",length(CycleLSD[23,2]:CycleLSD[23,3])),rep("C",length(CycleLSD[24,2]:CycleLSD[24,3])),rep("D",length(CycleLSD[25,2]:CycleLSD[25,3])),rep("E",length(CycleLSD[26,2]:CycleLSD[26,3])),rep("F",length(CycleLSD[27,2]:CycleLSD[27,3])),rep("G",length(CycleLSD[28,2]:CycleLSD[28,3])))
CutL25SD=c(CycleLSD[29,2]:CycleLSD[29,3],CycleLSD[30,2]:CycleLSD[30,3],CycleLSD[31,2]:CycleLSD[31,3],CycleLSD[32,2]:CycleLSD[32,3],CycleLSD[33,2]:CycleLSD[33,3],CycleLSD[34,2]:CycleLSD[34,3],CycleLSD[35,2]:CycleLSD[35,3])+t4
CycleL25SD=c(rep("A",length(CycleLSD[29,2]:CycleLSD[29,3])),rep("B",length(CycleLSD[30,2]:CycleLSD[30,3])),rep("C",length(CycleLSD[31,2]:CycleLSD[31,3])),rep("D",length(CycleLSD[32,2]:CycleLSD[32,3])),rep("E",length(CycleLSD[33,2]:CycleLSD[33,3])),rep("F",length(CycleLSD[34,2]:CycleLSD[34,3])),rep("G",length(CycleLSD[35,2]:CycleLSD[35,3])))

CutLSD=TestSD[c(CutL5SD,CutL10SD,CutL15SD,CutL20SD,CutL25SD),]
CutLSD$Cycle=c(CycleL5SD,CycleL10SD,CycleL15SD,CycleL20SD,CycleL25SD)

CutG5SD=c(CycleGSD[1,2]:CycleGSD[1,3],CycleGSD[2,2]:CycleGSD[2,3],CycleGSD[3,2]:CycleGSD[3,3],CycleGSD[4,2]:CycleGSD[4,3],CycleGSD[5,2]:CycleGSD[5,3],CycleGSD[6,2]:CycleGSD[6,3],CycleGSD[7,2]:CycleGSD[7,3])+t0
CycleG5SD=c(rep("A",length(CycleGSD[1,2]:CycleGSD[1,3])),rep("B",length(CycleGSD[2,2]:CycleGSD[2,3])),rep("C",length(CycleGSD[3,2]:CycleGSD[3,3])),rep("D",length(CycleGSD[4,2]:CycleGSD[4,3])),rep("E",length(CycleGSD[5,2]:CycleGSD[5,3])),rep("F",length(CycleGSD[6,2]:CycleGSD[6,3])),rep("G",length(CycleGSD[7,2]:CycleGSD[7,3])))
CutG10SD=c(CycleGSD[8,2]:CycleGSD[8,3],CycleGSD[9,2]:CycleGSD[9,3],CycleGSD[10,2]:CycleGSD[10,3],CycleGSD[11,2]:CycleGSD[11,3],CycleGSD[12,2]:CycleGSD[12,3],CycleGSD[13,2]:CycleGSD[13,3],CycleGSD[14,2]:CycleGSD[14,3])+t1
CycleG10SD=c(rep("A",length(CycleGSD[8,2]:CycleGSD[8,3])),rep("B",length(CycleGSD[9,2]:CycleGSD[9,3])),rep("C",length(CycleGSD[10,2]:CycleGSD[10,3])),rep("D",length(CycleGSD[11,2]:CycleGSD[11,3])),rep("E",length(CycleGSD[12,2]:CycleGSD[12,3])),rep("F",length(CycleGSD[13,2]:CycleGSD[13,3])),rep("G",length(CycleGSD[14,2]:CycleGSD[14,3])))
CutG15SD=c(CycleGSD[15,2]:CycleGSD[15,3],CycleGSD[16,2]:CycleGSD[16,3],CycleGSD[17,2]:CycleGSD[17,3],CycleGSD[18,2]:CycleGSD[18,3],CycleGSD[19,2]:CycleGSD[19,3],CycleGSD[20,2]:CycleGSD[20,3],CycleGSD[21,2]:CycleGSD[21,3])+t2
CycleG15SD=c(rep("A",length(CycleGSD[15,2]:CycleGSD[15,3])),rep("B",length(CycleGSD[16,2]:CycleGSD[16,3])),rep("C",length(CycleGSD[17,2]:CycleGSD[17,3])),rep("D",length(CycleGSD[18,2]:CycleGSD[18,3])),rep("E",length(CycleGSD[19,2]:CycleGSD[19,3])),rep("F",length(CycleGSD[20,2]:CycleGSD[20,3])),rep("G",length(CycleGSD[21,2]:CycleGSD[21,3])))
CutG20SD=c(CycleGSD[22,2]:CycleGSD[22,3],CycleGSD[23,2]:CycleGSD[23,3],CycleGSD[24,2]:CycleGSD[24,3],CycleGSD[25,2]:CycleGSD[25,3],CycleGSD[26,2]:CycleGSD[26,3],CycleGSD[27,2]:CycleGSD[27,3],CycleGSD[28,2]:CycleGSD[28,3])+t3
CycleG20SD=c(rep("A",length(CycleGSD[22,2]:CycleGSD[22,3])),rep("B",length(CycleGSD[23,2]:CycleGSD[23,3])),rep("C",length(CycleGSD[24,2]:CycleGSD[24,3])),rep("D",length(CycleGSD[25,2]:CycleGSD[25,3])),rep("E",length(CycleGSD[26,2]:CycleGSD[26,3])),rep("F",length(CycleGSD[27,2]:CycleGSD[27,3])),rep("G",length(CycleGSD[28,2]:CycleGSD[28,3])))
CutG25SD=c(CycleGSD[29,2]:CycleGSD[29,3],CycleGSD[30,2]:CycleGSD[30,3],CycleGSD[31,2]:CycleGSD[31,3],CycleGSD[32,2]:CycleGSD[32,3],CycleGSD[33,2]:CycleGSD[33,3],CycleGSD[34,2]:CycleGSD[34,3],CycleGSD[35,2]:CycleGSD[35,3])+t4
CycleG25SD=c(rep("A",length(CycleGSD[29,2]:CycleGSD[29,3])),rep("B",length(CycleGSD[30,2]:CycleGSD[30,3])),rep("C",length(CycleGSD[31,2]:CycleGSD[31,3])),rep("D",length(CycleGSD[32,2]:CycleGSD[32,3])),rep("E",length(CycleGSD[33,2]:CycleGSD[33,3])),rep("F",length(CycleGSD[34,2]:CycleGSD[34,3])),rep("G",length(CycleGSD[35,2]:CycleGSD[35,3])))

CutGSD=TestSD[c(CutG5SD,CutG10SD,CutG15SD,CutG20SD,CutG25SD),]
CutGSD$Cycle=c(CycleG5SD,CycleG10SD,CycleG15SD,CycleG20SD,CycleG25SD)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLSD=CutLSD[!CutLSD$Cycle=="A",]; CutGSD=CutGSD[!CutGSD$Cycle=="A",]
MinLSD=setDT(CutLSD)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLSD=setDT(CutLSD)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGSD=setDT(CutGSD)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGSD=setDT(CutGSD)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLSD=rbind(MinLSD[,c(3,1,2,4)],MaxLSD[,c(3,1,2,4)])
BiomGSD=rbind(MinGSD[,c(3,1,2,5)],MaxGSD[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLSD=as.data.frame(BiomLSD %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLSD=SlopLSD[order(as.numeric(SlopLSD$Temp)),]; SlopLSD$Slope=unlist(SlopLSD$Slope)
BiomLSD=cbind(SlopLSD,Max=BiomLSD[BiomLSD$L>100000]$L)
MaxiLSD=setDT(BiomLSD)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLSD=setDT(BiomLSD)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGSD=as.data.frame(BiomGSD %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGSD=SlopGSD[order(as.numeric(SlopGSD$Temp)),]; SlopGSD$Slope=unlist(SlopGSD$Slope)
BiomGSD=cbind(SlopGSD,Max=BiomGSD[BiomGSD$G>10000]$G)
MaxiGSD=setDT(BiomGSD)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGSD=setDT(BiomGSD)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL0=data.frame(MeanLSD[,1],MeanLSD[,2],MeanLSD[,3],TimeLSD[,2],TimeLSD[,3],MaxiLSD[,2],MaxiLSD[,3],MiniLSD[,2],MiniLSD[,3])
colnames(DataL0)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL","MeanMaxiL","SdMaxiL","MeanSlopL","SdSlopL")
DataL0$Scenario=rep("SD",5)

DataG0=data.frame(MeanGSD[,1],MeanGSD[,2],MeanGSD[,3],TimeGSD[,2],TimeGSD[,3],MaxiGSD[,2],MaxiGSD[,3],MiniGSD[,2],MiniGSD[,3])
colnames(DataG0)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG","MeanMaxiG","SdMaxiG","MeanSlopG","SdSlopG")
DataG0$Scenario=rep("SD",5)


########################################
### SCENARIO 1: AVERAGE TSR RESPONSE ###
########################################

# Leaf fall = 300 gC/m2/an = 300 000 mgC/m2/an
# Gammarus density = 30 mgDM/m2 = 15 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
tiff('Population Dynamics TSRA.tiff', units="in", width=15, height=8, res=1000)
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5TSRA=GammLeafModel(Temp=5,GammMass=TSRA(5,4.26),Leaf=300000/Days,Gamm=15)
Test5TSRA=Test5TSRA %>% mutate(L=if_else(L<10^-3, 0, L))
Test5TSRA=Test5TSRA %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5TSRA[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5 °C", cex=1.5, side=3, line=-3, at=300)

Test10TSRA=GammLeafModel(Temp=10,GammMass=TSRA(10,4.26),Leaf=300000/Days,Gamm=15)
Test10TSRA=Test10TSRA %>% mutate(L=if_else(L<10^-3, 0, L))
Test10TSRA=Test10TSRA %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10TSRA[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10 °C", cex=1.5, side=3, line=-3, at=300)

Test15TSRA=GammLeafModel(Temp=15,GammMass=TSRA(15,4.26),Leaf=300000/Days,Gamm=15)
Test15TSRA=Test15TSRA %>% mutate(L=if_else(L<10^-3, 0, L))
Test15TSRA=Test15TSRA %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15TSRA[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15 °C", cex=1.5, side=3, line=-3, at=300)

Test20TSRA=GammLeafModel(Temp=20,GammMass=TSRA(20,4.26),Leaf=300000/Days,Gamm=15)
Test20TSRA=Test20TSRA %>% mutate(L=if_else(L<10^-3, 0, L))
Test20TSRA=Test20TSRA %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20TSRA[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20 °C", cex=1.5, side=3, line=-3, at=300)

Test25TSRA=GammLeafModel(Temp=25,GammMass=TSRA(25,4.26),Leaf=300000/Days,Gamm=15)
Test25TSRA=Test25TSRA %>% mutate(L=if_else(L<10^-3, 0, L))
Test25TSRA=Test25TSRA %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25TSRA[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25 °C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))
dev.off()


# Calculate population metabolism and leaf ingestion
Test5TSRA$M=(MetaQuadra(5,TSRA(5,4.26))/1000)*Test5TSRA[,3]
Test5TSRA$I=(AttackQuadra(5,TSRA(5,4.26))*Test5TSRA[,2]/(1+AttackQuadra(5,TSRA(5,4.26))*1/(IngQuadra(5,TSRA(5,4.26))/1000)*Test5TSRA[,2]))*0.30*Test5TSRA[,3]
Test10TSRA$M=(MetaQuadra(10,TSRA(5,4.26))/1000)*Test10TSRA[,3]
Test10TSRA$I=(AttackQuadra(10,TSRA(5,4.26))*Test10TSRA[,2]/(1+AttackQuadra(10,TSRA(5,4.26))*1/(IngQuadra(10,TSRA(5,4.26))/1000)*Test10TSRA[,2]))*0.30*Test10TSRA[,3]
Test15TSRA$M=(MetaQuadra(15,TSRA(5,4.26))/1000)*Test15TSRA[,3]
Test15TSRA$I=(AttackQuadra(15,TSRA(5,4.26))*Test15TSRA[,2]/(1+AttackQuadra(15,TSRA(5,4.26))*1/(IngQuadra(15,TSRA(5,4.26))/1000)*Test15TSRA[,2]))*0.30*Test15TSRA[,3]
Test20TSRA$M=(MetaQuadra(20,TSRA(5,4.26))/1000)*Test20TSRA[,3]
Test20TSRA$I=(AttackQuadra(20,TSRA(5,4.26))*Test20TSRA[,2]/(1+AttackQuadra(20,TSRA(5,4.26))*1/(IngQuadra(20,TSRA(5,4.26))/1000)*Test20TSRA[,2]))*0.30*Test20TSRA[,3]
Test25TSRA$M=(MetaQuadra(25,TSRA(5,4.26))/1000)*Test25TSRA[,3]
Test25TSRA$I=(AttackQuadra(25,TSRA(5,4.26))*Test25TSRA[,2]/(1+AttackQuadra(25,TSRA(5,4.26))*1/(IngQuadra(25,TSRA(5,4.26))/1000)*Test25TSRA[,2]))*0.30*Test25TSRA[,3]

# Combine the dataframes
TestTSRA=rbind(Test5TSRA,Test10TSRA,Test15TSRA,Test20TSRA,Test25TSRA)
TestTSRA$Temp=rep(c("5","10","15","20","25"), each=2556)
TestTSRA$Temp=factor(TestTSRA$Temp, levels=unique(TestTSRA$Temp))
colnames(TestTSRA)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestTSRA=subset(TestTSRA, !Time==2555); TestTSRA$Year=rep(rep(1:(nrow(TestTSRA)/(365*5)), each=365),5)
CutTSRA=as.data.frame(setDT(TestTSRA)[TestTSRA[, tail(.I, -16), by=Year]$V1]); CutTSRA=CutTSRA[!CutTSRA$Year=="1",]
MeanL=as.data.frame(setDT(CutTSRA)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutTSRA)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLTSRA=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGTSRA=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestTSRA=subset(TestTSRA, !Time==2555); TestTSRA$Year=rep(rep(1:(nrow(TestTSRA)/(365*5)), each=365),5)
CutTSRA=as.data.frame(setDT(TestTSRA)[TestTSRA[, tail(.I, -16), by=Year]$V1]); CutTSRA=CutTSRA[!CutTSRA$Year=="1",]
ThLTSRA=CutTSRA[CutTSRA$L>60000,]; ThGTSRA=CutTSRA[CutTSRA$G>5000,]
TimeLTSRA=as.data.frame(setDT(ThLTSRA)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGTSRA=as.data.frame(setDT(ThGTSRA)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLTSRA=as.data.frame(setDT(TimeLTSRA)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGTSRA=as.data.frame(setDT(TimeGTSRA)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLTSRA=as.data.frame(setDT(TestTSRA)[, .(MaxLTSRA=findPeaks(L), MinLTSRA=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGTSRA=as.data.frame(setDT(TestTSRA)[, .(MaxGTSRA=findPeaks(G), MinGTSRA=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5TSRA=c(CycleLTSRA[1,2]:CycleLTSRA[1,3],CycleLTSRA[2,2]:CycleLTSRA[2,3],CycleLTSRA[3,2]:CycleLTSRA[3,3],CycleLTSRA[4,2]:CycleLTSRA[4,3],CycleLTSRA[5,2]:CycleLTSRA[5,3],CycleLTSRA[6,2]:CycleLTSRA[6,3],CycleLTSRA[7,2]:CycleLTSRA[7,3])+t0
CycleL5TSRA=c(rep("A",length(CycleLTSRA[1,2]:CycleLTSRA[1,3])),rep("B",length(CycleLTSRA[2,2]:CycleLTSRA[2,3])),rep("C",length(CycleLTSRA[3,2]:CycleLTSRA[3,3])),rep("D",length(CycleLTSRA[4,2]:CycleLTSRA[4,3])),rep("E",length(CycleLTSRA[5,2]:CycleLTSRA[5,3])),rep("F",length(CycleLTSRA[6,2]:CycleLTSRA[6,3])),rep("G",length(CycleLTSRA[7,2]:CycleLTSRA[7,3])))
CutL10TSRA=c(CycleLTSRA[8,2]:CycleLTSRA[8,3],CycleLTSRA[9,2]:CycleLTSRA[9,3],CycleLTSRA[10,2]:CycleLTSRA[10,3],CycleLTSRA[11,2]:CycleLTSRA[11,3],CycleLTSRA[12,2]:CycleLTSRA[12,3],CycleLTSRA[13,2]:CycleLTSRA[13,3],CycleLTSRA[14,2]:CycleLTSRA[14,3])+t1
CycleL10TSRA=c(rep("A",length(CycleLTSRA[8,2]:CycleLTSRA[8,3])),rep("B",length(CycleLTSRA[9,2]:CycleLTSRA[9,3])),rep("C",length(CycleLTSRA[10,2]:CycleLTSRA[10,3])),rep("D",length(CycleLTSRA[11,2]:CycleLTSRA[11,3])),rep("E",length(CycleLTSRA[12,2]:CycleLTSRA[12,3])),rep("F",length(CycleLTSRA[13,2]:CycleLTSRA[13,3])),rep("G",length(CycleLTSRA[14,2]:CycleLTSRA[14,3])))
CutL15TSRA=c(CycleLTSRA[15,2]:CycleLTSRA[15,3],CycleLTSRA[16,2]:CycleLTSRA[16,3],CycleLTSRA[17,2]:CycleLTSRA[17,3],CycleLTSRA[18,2]:CycleLTSRA[18,3],CycleLTSRA[19,2]:CycleLTSRA[19,3],CycleLTSRA[20,2]:CycleLTSRA[20,3],CycleLTSRA[21,2]:CycleLTSRA[21,3])+t2
CycleL15TSRA=c(rep("A",length(CycleLTSRA[15,2]:CycleLTSRA[15,3])),rep("B",length(CycleLTSRA[16,2]:CycleLTSRA[16,3])),rep("C",length(CycleLTSRA[17,2]:CycleLTSRA[17,3])),rep("D",length(CycleLTSRA[18,2]:CycleLTSRA[18,3])),rep("E",length(CycleLTSRA[19,2]:CycleLTSRA[19,3])),rep("F",length(CycleLTSRA[20,2]:CycleLTSRA[20,3])),rep("G",length(CycleLTSRA[21,2]:CycleLTSRA[21,3])))
CutL20TSRA=c(CycleLTSRA[22,2]:CycleLTSRA[22,3],CycleLTSRA[23,2]:CycleLTSRA[23,3],CycleLTSRA[24,2]:CycleLTSRA[24,3],CycleLTSRA[25,2]:CycleLTSRA[25,3],CycleLTSRA[26,2]:CycleLTSRA[26,3],CycleLTSRA[27,2]:CycleLTSRA[27,3],CycleLTSRA[28,2]:CycleLTSRA[28,3])+t3
CycleL20TSRA=c(rep("A",length(CycleLTSRA[22,2]:CycleLTSRA[22,3])),rep("B",length(CycleLTSRA[23,2]:CycleLTSRA[23,3])),rep("C",length(CycleLTSRA[24,2]:CycleLTSRA[24,3])),rep("D",length(CycleLTSRA[25,2]:CycleLTSRA[25,3])),rep("E",length(CycleLTSRA[26,2]:CycleLTSRA[26,3])),rep("F",length(CycleLTSRA[27,2]:CycleLTSRA[27,3])),rep("G",length(CycleLTSRA[28,2]:CycleLTSRA[28,3])))
CutL25TSRA=c(CycleLTSRA[29,2]:CycleLTSRA[29,3],CycleLTSRA[30,2]:CycleLTSRA[30,3],CycleLTSRA[31,2]:CycleLTSRA[31,3],CycleLTSRA[32,2]:CycleLTSRA[32,3],CycleLTSRA[33,2]:CycleLTSRA[33,3],CycleLTSRA[34,2]:CycleLTSRA[34,3],CycleLTSRA[35,2]:CycleLTSRA[35,3])+t4
CycleL25TSRA=c(rep("A",length(CycleLTSRA[29,2]:CycleLTSRA[29,3])),rep("B",length(CycleLTSRA[30,2]:CycleLTSRA[30,3])),rep("C",length(CycleLTSRA[31,2]:CycleLTSRA[31,3])),rep("D",length(CycleLTSRA[32,2]:CycleLTSRA[32,3])),rep("E",length(CycleLTSRA[33,2]:CycleLTSRA[33,3])),rep("F",length(CycleLTSRA[34,2]:CycleLTSRA[34,3])),rep("G",length(CycleLTSRA[35,2]:CycleLTSRA[35,3])))

CutLTSRA=TestTSRA[c(CutL5TSRA,CutL10TSRA,CutL15TSRA,CutL20TSRA,CutL25TSRA),]
CutLTSRA$Cycle=c(CycleL5TSRA,CycleL10TSRA,CycleL15TSRA,CycleL20TSRA,CycleL25TSRA)

CutG5TSRA=c(CycleGTSRA[1,2]:CycleGTSRA[1,3],CycleGTSRA[2,2]:CycleGTSRA[2,3],CycleGTSRA[3,2]:CycleGTSRA[3,3],CycleGTSRA[4,2]:CycleGTSRA[4,3],CycleGTSRA[5,2]:CycleGTSRA[5,3],CycleGTSRA[6,2]:CycleGTSRA[6,3],CycleGTSRA[7,2]:CycleGTSRA[7,3])+t0
CycleG5TSRA=c(rep("A",length(CycleGTSRA[1,2]:CycleGTSRA[1,3])),rep("B",length(CycleGTSRA[2,2]:CycleGTSRA[2,3])),rep("C",length(CycleGTSRA[3,2]:CycleGTSRA[3,3])),rep("D",length(CycleGTSRA[4,2]:CycleGTSRA[4,3])),rep("E",length(CycleGTSRA[5,2]:CycleGTSRA[5,3])),rep("F",length(CycleGTSRA[6,2]:CycleGTSRA[6,3])),rep("G",length(CycleGTSRA[7,2]:CycleGTSRA[7,3])))
CutG10TSRA=c(CycleGTSRA[8,2]:CycleGTSRA[8,3],CycleGTSRA[9,2]:CycleGTSRA[9,3],CycleGTSRA[10,2]:CycleGTSRA[10,3],CycleGTSRA[11,2]:CycleGTSRA[11,3],CycleGTSRA[12,2]:CycleGTSRA[12,3],CycleGTSRA[13,2]:CycleGTSRA[13,3],CycleGTSRA[14,2]:CycleGTSRA[14,3])+t1
CycleG10TSRA=c(rep("A",length(CycleGTSRA[8,2]:CycleGTSRA[8,3])),rep("B",length(CycleGTSRA[9,2]:CycleGTSRA[9,3])),rep("C",length(CycleGTSRA[10,2]:CycleGTSRA[10,3])),rep("D",length(CycleGTSRA[11,2]:CycleGTSRA[11,3])),rep("E",length(CycleGTSRA[12,2]:CycleGTSRA[12,3])),rep("F",length(CycleGTSRA[13,2]:CycleGTSRA[13,3])),rep("G",length(CycleGTSRA[14,2]:CycleGTSRA[14,3])))
CutG15TSRA=c(CycleGTSRA[15,2]:CycleGTSRA[15,3],CycleGTSRA[16,2]:CycleGTSRA[16,3],CycleGTSRA[17,2]:CycleGTSRA[17,3],CycleGTSRA[18,2]:CycleGTSRA[18,3],CycleGTSRA[19,2]:CycleGTSRA[19,3],CycleGTSRA[20,2]:CycleGTSRA[20,3],CycleGTSRA[21,2]:CycleGTSRA[21,3])+t2
CycleG15TSRA=c(rep("A",length(CycleGTSRA[15,2]:CycleGTSRA[15,3])),rep("B",length(CycleGTSRA[16,2]:CycleGTSRA[16,3])),rep("C",length(CycleGTSRA[17,2]:CycleGTSRA[17,3])),rep("D",length(CycleGTSRA[18,2]:CycleGTSRA[18,3])),rep("E",length(CycleGTSRA[19,2]:CycleGTSRA[19,3])),rep("F",length(CycleGTSRA[20,2]:CycleGTSRA[20,3])),rep("G",length(CycleGTSRA[21,2]:CycleGTSRA[21,3])))
CutG20TSRA=c(CycleGTSRA[22,2]:CycleGTSRA[22,3],CycleGTSRA[23,2]:CycleGTSRA[23,3],CycleGTSRA[24,2]:CycleGTSRA[24,3],CycleGTSRA[25,2]:CycleGTSRA[25,3],CycleGTSRA[26,2]:CycleGTSRA[26,3],CycleGTSRA[27,2]:CycleGTSRA[27,3],CycleGTSRA[28,2]:CycleGTSRA[28,3])+t3
CycleG20TSRA=c(rep("A",length(CycleGTSRA[22,2]:CycleGTSRA[22,3])),rep("B",length(CycleGTSRA[23,2]:CycleGTSRA[23,3])),rep("C",length(CycleGTSRA[24,2]:CycleGTSRA[24,3])),rep("D",length(CycleGTSRA[25,2]:CycleGTSRA[25,3])),rep("E",length(CycleGTSRA[26,2]:CycleGTSRA[26,3])),rep("F",length(CycleGTSRA[27,2]:CycleGTSRA[27,3])),rep("G",length(CycleGTSRA[28,2]:CycleGTSRA[28,3])))
CutG25TSRA=c(CycleGTSRA[29,2]:CycleGTSRA[29,3],CycleGTSRA[30,2]:CycleGTSRA[30,3],CycleGTSRA[31,2]:CycleGTSRA[31,3],CycleGTSRA[32,2]:CycleGTSRA[32,3],CycleGTSRA[33,2]:CycleGTSRA[33,3],CycleGTSRA[34,2]:CycleGTSRA[34,3],CycleGTSRA[35,2]:CycleGTSRA[35,3])+t4
CycleG25TSRA=c(rep("A",length(CycleGTSRA[29,2]:CycleGTSRA[29,3])),rep("B",length(CycleGTSRA[30,2]:CycleGTSRA[30,3])),rep("C",length(CycleGTSRA[31,2]:CycleGTSRA[31,3])),rep("D",length(CycleGTSRA[32,2]:CycleGTSRA[32,3])),rep("E",length(CycleGTSRA[33,2]:CycleGTSRA[33,3])),rep("F",length(CycleGTSRA[34,2]:CycleGTSRA[34,3])),rep("G",length(CycleGTSRA[35,2]:CycleGTSRA[35,3])))

CutGTSRA=TestTSRA[c(CutG5TSRA,CutG10TSRA,CutG15TSRA,CutG20TSRA,CutG25TSRA),]
CutGTSRA$Cycle=c(CycleG5TSRA,CycleG10TSRA,CycleG15TSRA,CycleG20TSRA,CycleG25TSRA)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLTSRA=CutLTSRA[!CutLTSRA$Cycle=="A",]; CutGTSRA=CutGTSRA[!CutGTSRA$Cycle=="A",]
MinLTSRA=setDT(CutLTSRA)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLTSRA=setDT(CutLTSRA)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGTSRA=setDT(CutGTSRA)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGTSRA=setDT(CutGTSRA)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLTSRA=rbind(MinLTSRA[,c(3,1,2,4)],MaxLTSRA[,c(3,1,2,4)])
BiomGTSRA=rbind(MinGTSRA[,c(3,1,2,5)],MaxGTSRA[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLTSRA=as.data.frame(BiomLTSRA %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLTSRA=SlopLTSRA[order(as.numeric(SlopLTSRA$Temp)),]; SlopLTSRA$Slope=unlist(SlopLTSRA$Slope)
BiomLTSRA=cbind(SlopLTSRA,Max=BiomLTSRA[BiomLTSRA$L>100000]$L)
MaxiLTSRA=setDT(BiomLTSRA)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLTSRA=setDT(BiomLTSRA)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGTSRA=as.data.frame(BiomGTSRA %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGTSRA=SlopGTSRA[order(as.numeric(SlopGTSRA$Temp)),]; SlopGTSRA$Slope=unlist(SlopGTSRA$Slope)
BiomGTSRA=cbind(SlopGTSRA,Max=BiomGTSRA[BiomGTSRA$G>10000]$G)
MaxiGTSRA=setDT(BiomGTSRA)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGTSRA=setDT(BiomGTSRA)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL1=data.frame(MeanLTSRA[,1],MeanLTSRA[,2],MeanLTSRA[,3],TimeLTSRA[,2],TimeLTSRA[,3],MaxiLTSRA[,2],MaxiLTSRA[,3],MiniLTSRA[,2],MiniLTSRA[,3])
colnames(DataL1)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL","MeanMaxiL","SdMaxiL","MeanSlopL","SdSlopL")
DataL1$Scenario=rep("TSRA",5)

DataG1=data.frame(MeanGTSRA[,1],MeanGTSRA[,2],MeanGTSRA[,3],TimeGTSRA[,2],TimeGTSRA[,3],MaxiGTSRA[,2],MaxiGTSRA[,3],MiniGTSRA[,2],MiniGTSRA[,3])
colnames(DataG1)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG","MeanMaxiG","SdMaxiG","MeanSlopG","SdSlopG")
DataG1$Scenario=rep("TSRA",5)


########################################
### SCENARIO 2: MAXIMUM TSR RESPONSE ###
########################################

# Leaf fall = 300 gC/m2/an = 300 000 mgC/m2/an
# Gammarus density = 30 mgDM/m2 = 15 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
tiff('Population Dynamics TSRM.tiff', units="in", width=15, height=8, res=1000)
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5TSRM=GammLeafModel(Temp=5,GammMass=TSRM(5,4.26),Leaf=300000/Days,Gamm=15)
Test5TSRM=Test5TSRM %>% mutate(L=if_else(L<10^-3, 0, L))
Test5TSRM=Test5TSRM %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5TSRM[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5 °C", cex=1.5, side=3, line=-3, at=300)

Test10TSRM=GammLeafModel(Temp=10,GammMass=TSRM(10,4.26),Leaf=300000/Days,Gamm=15)
Test10TSRM=Test10TSRM %>% mutate(L=if_else(L<10^-3, 0, L))
Test10TSRM=Test10TSRM %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10TSRM[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10 °C", cex=1.5, side=3, line=-3, at=300)

Test15TSRM=GammLeafModel(Temp=15,GammMass=TSRM(15,4.26),Leaf=300000/Days,Gamm=15)
Test15TSRM=Test15TSRM %>% mutate(L=if_else(L<10^-3, 0, L))
Test15TSRM=Test15TSRM %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15TSRM[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15 °C", cex=1.5, side=3, line=-3, at=300)

Test20TSRM=GammLeafModel(Temp=20,GammMass=TSRM(20,4.26),Leaf=300000/Days,Gamm=15)
Test20TSRM=Test20TSRM %>% mutate(L=if_else(L<10^-3, 0, L))
Test20TSRM=Test20TSRM %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20TSRM[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20 °C", cex=1.5, side=3, line=-3, at=300)

Test25TSRM=GammLeafModel(Temp=25,GammMass=TSRM(25,4.26),Leaf=300000/Days,Gamm=15)
Test25TSRM=Test25TSRM %>% mutate(L=if_else(L<10^-3, 0, L))
Test25TSRM=Test25TSRM %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25TSRM[,-1]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25 °C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))
dev.off()


# Calculate population metabolism and leaf ingestion
Test5TSRM$M=(MetaQuadra(5,TSRM(5,4.26))/1000)*Test5TSRM[,3]
Test5TSRM$I=(AttackQuadra(5,TSRM(5,4.26))*Test5TSRM[,2]/(1+AttackQuadra(5,TSRM(5,4.26))*1/(IngQuadra(5,TSRM(5,4.26))/1000)*Test5TSRM[,2]))*0.30*Test5TSRM[,3]
Test10TSRM$M=(MetaQuadra(10,TSRM(5,4.26))/1000)*Test10TSRM[,3]
Test10TSRM$I=(AttackQuadra(10,TSRM(5,4.26))*Test10TSRM[,2]/(1+AttackQuadra(10,TSRM(5,4.26))*1/(IngQuadra(10,TSRM(5,4.26))/1000)*Test10TSRM[,2]))*0.30*Test10TSRM[,3]
Test15TSRM$M=(MetaQuadra(15,TSRM(5,4.26))/1000)*Test15TSRM[,3]
Test15TSRM$I=(AttackQuadra(15,TSRM(5,4.26))*Test15TSRM[,2]/(1+AttackQuadra(15,TSRM(5,4.26))*1/(IngQuadra(15,TSRM(5,4.26))/1000)*Test15TSRM[,2]))*0.30*Test15TSRM[,3]
Test20TSRM$M=(MetaQuadra(20,TSRM(5,4.26))/1000)*Test20TSRM[,3]
Test20TSRM$I=(AttackQuadra(20,TSRM(5,4.26))*Test20TSRM[,2]/(1+AttackQuadra(20,TSRM(5,4.26))*1/(IngQuadra(20,TSRM(5,4.26))/1000)*Test20TSRM[,2]))*0.30*Test20TSRM[,3]
Test25TSRM$M=(MetaQuadra(25,TSRM(5,4.26))/1000)*Test25TSRM[,3]
Test25TSRM$I=(AttackQuadra(25,TSRM(5,4.26))*Test25TSRM[,2]/(1+AttackQuadra(25,TSRM(5,4.26))*1/(IngQuadra(25,TSRM(5,4.26))/1000)*Test25TSRM[,2]))*0.30*Test25TSRM[,3]

# Combine the dataframes
TestTSRM=rbind(Test5TSRM,Test10TSRM,Test15TSRM,Test20TSRM,Test25TSRM)
TestTSRM$Temp=rep(c("5","10","15","20","25"), each=2556)
TestTSRM$Temp=factor(TestTSRM$Temp, levels=unique(TestTSRM$Temp))
colnames(TestTSRM)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestTSRM=subset(TestTSRM, !Time==2555); TestTSRM$Year=rep(rep(1:(nrow(TestTSRM)/(365*5)), each=365),5)
CutTSRM=as.data.frame(setDT(TestTSRM)[TestTSRM[, tail(.I, -16), by=Year]$V1]); CutTSRM=CutTSRM[!CutTSRM$Year=="1",]
MeanL=as.data.frame(setDT(CutTSRM)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutTSRM)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLTSRM=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGTSRM=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestTSRM=subset(TestTSRM, !Time==2555); TestTSRM$Year=rep(rep(1:(nrow(TestTSRM)/(365*5)), each=365),5)
CutTSRM=as.data.frame(setDT(TestTSRM)[TestTSRM[, tail(.I, -16), by=Year]$V1]); CutTSRM=CutTSRM[!CutTSRM$Year=="1",]
ThLTSRM=CutTSRM[CutTSRM$L>60000,]; ThGTSRM=CutTSRM[CutTSRM$G>5000,]
TimeLTSRM=as.data.frame(setDT(ThLTSRM)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGTSRM=as.data.frame(setDT(ThGTSRM)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLTSRM=as.data.frame(setDT(TimeLTSRM)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGTSRM=as.data.frame(setDT(TimeGTSRM)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLTSRM=as.data.frame(setDT(TestTSRM)[, .(MaxLTSRM=findPeaks(L), MinLTSRM=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGTSRM=as.data.frame(setDT(TestTSRM)[, .(MaxGTSRM=findPeaks(G), MinGTSRM=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5TSRM=c(CycleLTSRM[1,2]:CycleLTSRM[1,3],CycleLTSRM[2,2]:CycleLTSRM[2,3],CycleLTSRM[3,2]:CycleLTSRM[3,3],CycleLTSRM[4,2]:CycleLTSRM[4,3],CycleLTSRM[5,2]:CycleLTSRM[5,3],CycleLTSRM[6,2]:CycleLTSRM[6,3],CycleLTSRM[7,2]:CycleLTSRM[7,3])+t0
CycleL5TSRM=c(rep("A",length(CycleLTSRM[1,2]:CycleLTSRM[1,3])),rep("B",length(CycleLTSRM[2,2]:CycleLTSRM[2,3])),rep("C",length(CycleLTSRM[3,2]:CycleLTSRM[3,3])),rep("D",length(CycleLTSRM[4,2]:CycleLTSRM[4,3])),rep("E",length(CycleLTSRM[5,2]:CycleLTSRM[5,3])),rep("F",length(CycleLTSRM[6,2]:CycleLTSRM[6,3])),rep("G",length(CycleLTSRM[7,2]:CycleLTSRM[7,3])))
CutL10TSRM=c(CycleLTSRM[8,2]:CycleLTSRM[8,3],CycleLTSRM[9,2]:CycleLTSRM[9,3],CycleLTSRM[10,2]:CycleLTSRM[10,3],CycleLTSRM[11,2]:CycleLTSRM[11,3],CycleLTSRM[12,2]:CycleLTSRM[12,3],CycleLTSRM[13,2]:CycleLTSRM[13,3],CycleLTSRM[14,2]:CycleLTSRM[14,3])+t1
CycleL10TSRM=c(rep("A",length(CycleLTSRM[8,2]:CycleLTSRM[8,3])),rep("B",length(CycleLTSRM[9,2]:CycleLTSRM[9,3])),rep("C",length(CycleLTSRM[10,2]:CycleLTSRM[10,3])),rep("D",length(CycleLTSRM[11,2]:CycleLTSRM[11,3])),rep("E",length(CycleLTSRM[12,2]:CycleLTSRM[12,3])),rep("F",length(CycleLTSRM[13,2]:CycleLTSRM[13,3])),rep("G",length(CycleLTSRM[14,2]:CycleLTSRM[14,3])))
CutL15TSRM=c(CycleLTSRM[15,2]:CycleLTSRM[15,3],CycleLTSRM[16,2]:CycleLTSRM[16,3],CycleLTSRM[17,2]:CycleLTSRM[17,3],CycleLTSRM[18,2]:CycleLTSRM[18,3],CycleLTSRM[19,2]:CycleLTSRM[19,3],CycleLTSRM[20,2]:CycleLTSRM[20,3],CycleLTSRM[21,2]:CycleLTSRM[21,3])+t2
CycleL15TSRM=c(rep("A",length(CycleLTSRM[15,2]:CycleLTSRM[15,3])),rep("B",length(CycleLTSRM[16,2]:CycleLTSRM[16,3])),rep("C",length(CycleLTSRM[17,2]:CycleLTSRM[17,3])),rep("D",length(CycleLTSRM[18,2]:CycleLTSRM[18,3])),rep("E",length(CycleLTSRM[19,2]:CycleLTSRM[19,3])),rep("F",length(CycleLTSRM[20,2]:CycleLTSRM[20,3])),rep("G",length(CycleLTSRM[21,2]:CycleLTSRM[21,3])))
CutL20TSRM=c(CycleLTSRM[22,2]:CycleLTSRM[22,3],CycleLTSRM[23,2]:CycleLTSRM[23,3],CycleLTSRM[24,2]:CycleLTSRM[24,3],CycleLTSRM[25,2]:CycleLTSRM[25,3],CycleLTSRM[26,2]:CycleLTSRM[26,3],CycleLTSRM[27,2]:CycleLTSRM[27,3],CycleLTSRM[28,2]:CycleLTSRM[28,3])+t3
CycleL20TSRM=c(rep("A",length(CycleLTSRM[22,2]:CycleLTSRM[22,3])),rep("B",length(CycleLTSRM[23,2]:CycleLTSRM[23,3])),rep("C",length(CycleLTSRM[24,2]:CycleLTSRM[24,3])),rep("D",length(CycleLTSRM[25,2]:CycleLTSRM[25,3])),rep("E",length(CycleLTSRM[26,2]:CycleLTSRM[26,3])),rep("F",length(CycleLTSRM[27,2]:CycleLTSRM[27,3])),rep("G",length(CycleLTSRM[28,2]:CycleLTSRM[28,3])))
CutL25TSRM=c(CycleLTSRM[29,2]:CycleLTSRM[29,3],CycleLTSRM[30,2]:CycleLTSRM[30,3],CycleLTSRM[31,2]:CycleLTSRM[31,3],CycleLTSRM[32,2]:CycleLTSRM[32,3],CycleLTSRM[33,2]:CycleLTSRM[33,3],CycleLTSRM[34,2]:CycleLTSRM[34,3],CycleLTSRM[35,2]:CycleLTSRM[35,3])+t4
CycleL25TSRM=c(rep("A",length(CycleLTSRM[29,2]:CycleLTSRM[29,3])),rep("B",length(CycleLTSRM[30,2]:CycleLTSRM[30,3])),rep("C",length(CycleLTSRM[31,2]:CycleLTSRM[31,3])),rep("D",length(CycleLTSRM[32,2]:CycleLTSRM[32,3])),rep("E",length(CycleLTSRM[33,2]:CycleLTSRM[33,3])),rep("F",length(CycleLTSRM[34,2]:CycleLTSRM[34,3])),rep("G",length(CycleLTSRM[35,2]:CycleLTSRM[35,3])))

CutLTSRM=TestTSRM[c(CutL5TSRM,CutL10TSRM,CutL15TSRM,CutL20TSRM,CutL25TSRM),]
CutLTSRM$Cycle=c(CycleL5TSRM,CycleL10TSRM,CycleL15TSRM,CycleL20TSRM,CycleL25TSRM)

CutG5TSRM=c(CycleGTSRM[1,2]:CycleGTSRM[1,3],CycleGTSRM[2,2]:CycleGTSRM[2,3],CycleGTSRM[3,2]:CycleGTSRM[3,3],CycleGTSRM[4,2]:CycleGTSRM[4,3],CycleGTSRM[5,2]:CycleGTSRM[5,3],CycleGTSRM[6,2]:CycleGTSRM[6,3],CycleGTSRM[7,2]:CycleGTSRM[7,3])+t0
CycleG5TSRM=c(rep("A",length(CycleGTSRM[1,2]:CycleGTSRM[1,3])),rep("B",length(CycleGTSRM[2,2]:CycleGTSRM[2,3])),rep("C",length(CycleGTSRM[3,2]:CycleGTSRM[3,3])),rep("D",length(CycleGTSRM[4,2]:CycleGTSRM[4,3])),rep("E",length(CycleGTSRM[5,2]:CycleGTSRM[5,3])),rep("F",length(CycleGTSRM[6,2]:CycleGTSRM[6,3])),rep("G",length(CycleGTSRM[7,2]:CycleGTSRM[7,3])))
CutG10TSRM=c(CycleGTSRM[8,2]:CycleGTSRM[8,3],CycleGTSRM[9,2]:CycleGTSRM[9,3],CycleGTSRM[10,2]:CycleGTSRM[10,3],CycleGTSRM[11,2]:CycleGTSRM[11,3],CycleGTSRM[12,2]:CycleGTSRM[12,3],CycleGTSRM[13,2]:CycleGTSRM[13,3],CycleGTSRM[14,2]:CycleGTSRM[14,3])+t1
CycleG10TSRM=c(rep("A",length(CycleGTSRM[8,2]:CycleGTSRM[8,3])),rep("B",length(CycleGTSRM[9,2]:CycleGTSRM[9,3])),rep("C",length(CycleGTSRM[10,2]:CycleGTSRM[10,3])),rep("D",length(CycleGTSRM[11,2]:CycleGTSRM[11,3])),rep("E",length(CycleGTSRM[12,2]:CycleGTSRM[12,3])),rep("F",length(CycleGTSRM[13,2]:CycleGTSRM[13,3])),rep("G",length(CycleGTSRM[14,2]:CycleGTSRM[14,3])))
CutG15TSRM=c(CycleGTSRM[15,2]:CycleGTSRM[15,3],CycleGTSRM[16,2]:CycleGTSRM[16,3],CycleGTSRM[17,2]:CycleGTSRM[17,3],CycleGTSRM[18,2]:CycleGTSRM[18,3],CycleGTSRM[19,2]:CycleGTSRM[19,3],CycleGTSRM[20,2]:CycleGTSRM[20,3],CycleGTSRM[21,2]:CycleGTSRM[21,3])+t2
CycleG15TSRM=c(rep("A",length(CycleGTSRM[15,2]:CycleGTSRM[15,3])),rep("B",length(CycleGTSRM[16,2]:CycleGTSRM[16,3])),rep("C",length(CycleGTSRM[17,2]:CycleGTSRM[17,3])),rep("D",length(CycleGTSRM[18,2]:CycleGTSRM[18,3])),rep("E",length(CycleGTSRM[19,2]:CycleGTSRM[19,3])),rep("F",length(CycleGTSRM[20,2]:CycleGTSRM[20,3])),rep("G",length(CycleGTSRM[21,2]:CycleGTSRM[21,3])))
CutG20TSRM=c(CycleGTSRM[22,2]:CycleGTSRM[22,3],CycleGTSRM[23,2]:CycleGTSRM[23,3],CycleGTSRM[24,2]:CycleGTSRM[24,3],CycleGTSRM[25,2]:CycleGTSRM[25,3],CycleGTSRM[26,2]:CycleGTSRM[26,3],CycleGTSRM[27,2]:CycleGTSRM[27,3],CycleGTSRM[28,2]:CycleGTSRM[28,3])+t3
CycleG20TSRM=c(rep("A",length(CycleGTSRM[22,2]:CycleGTSRM[22,3])),rep("B",length(CycleGTSRM[23,2]:CycleGTSRM[23,3])),rep("C",length(CycleGTSRM[24,2]:CycleGTSRM[24,3])),rep("D",length(CycleGTSRM[25,2]:CycleGTSRM[25,3])),rep("E",length(CycleGTSRM[26,2]:CycleGTSRM[26,3])),rep("F",length(CycleGTSRM[27,2]:CycleGTSRM[27,3])),rep("G",length(CycleGTSRM[28,2]:CycleGTSRM[28,3])))
CutG25TSRM=c(CycleGTSRM[29,2]:CycleGTSRM[29,3],CycleGTSRM[30,2]:CycleGTSRM[30,3],CycleGTSRM[31,2]:CycleGTSRM[31,3],CycleGTSRM[32,2]:CycleGTSRM[32,3],CycleGTSRM[33,2]:CycleGTSRM[33,3],CycleGTSRM[34,2]:CycleGTSRM[34,3],CycleGTSRM[35,2]:CycleGTSRM[35,3])+t4
CycleG25TSRM=c(rep("A",length(CycleGTSRM[29,2]:CycleGTSRM[29,3])),rep("B",length(CycleGTSRM[30,2]:CycleGTSRM[30,3])),rep("C",length(CycleGTSRM[31,2]:CycleGTSRM[31,3])),rep("D",length(CycleGTSRM[32,2]:CycleGTSRM[32,3])),rep("E",length(CycleGTSRM[33,2]:CycleGTSRM[33,3])),rep("F",length(CycleGTSRM[34,2]:CycleGTSRM[34,3])),rep("G",length(CycleGTSRM[35,2]:CycleGTSRM[35,3])))

CutGTSRM=TestTSRM[c(CutG5TSRM,CutG10TSRM,CutG15TSRM,CutG20TSRM,CutG25TSRM),]
CutGTSRM$Cycle=c(CycleG5TSRM,CycleG10TSRM,CycleG15TSRM,CycleG20TSRM,CycleG25TSRM)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLTSRM=CutLTSRM[!CutLTSRM$Cycle=="A",]; CutGTSRM=CutGTSRM[!CutGTSRM$Cycle=="A",]
MinLTSRM=setDT(CutLTSRM)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLTSRM=setDT(CutLTSRM)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGTSRM=setDT(CutGTSRM)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGTSRM=setDT(CutGTSRM)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLTSRM=rbind(MinLTSRM[,c(3,1,2,4)],MaxLTSRM[,c(3,1,2,4)])
BiomGTSRM=rbind(MinGTSRM[,c(3,1,2,5)],MaxGTSRM[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLTSRM=as.data.frame(BiomLTSRM %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLTSRM=SlopLTSRM[order(as.numeric(SlopLTSRM$Temp)),]; SlopLTSRM$Slope=unlist(SlopLTSRM$Slope)
BiomLTSRM=cbind(SlopLTSRM,Max=BiomLTSRM[BiomLTSRM$L>100000]$L)
MaxiLTSRM=setDT(BiomLTSRM)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLTSRM=setDT(BiomLTSRM)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGTSRM=as.data.frame(BiomGTSRM %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGTSRM=SlopGTSRM[order(as.numeric(SlopGTSRM$Temp)),]; SlopGTSRM$Slope=unlist(SlopGTSRM$Slope)
BiomGTSRM=cbind(SlopGTSRM,Max=BiomGTSRM[BiomGTSRM$G>10000]$G)
MaxiGTSRM=setDT(BiomGTSRM)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGTSRM=setDT(BiomGTSRM)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL2=data.frame(MeanLTSRM[,1],MeanLTSRM[,2],MeanLTSRM[,3],TimeLTSRM[,2],TimeLTSRM[,3],MaxiLTSRM[,2],MaxiLTSRM[,3],MiniLTSRM[,2],MiniLTSRM[,3])
colnames(DataL2)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL","MeanMaxiL","SdMaxiL","MeanSlopL","SdSlopL")
DataL2$Scenario=rep("TSRM",5)

DataG2=data.frame(MeanGTSRM[,1],MeanGTSRM[,2],MeanGTSRM[,3],TimeGTSRM[,2],TimeGTSRM[,3],MaxiGTSRM[,2],MaxiGTSRM[,3],MiniGTSRM[,2],MiniGTSRM[,3])
colnames(DataG2)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG","MeanMaxiG","SdMaxiG","MeanSlopG","SdSlopG")
DataG2$Scenario=rep("TSRM",5)




############################################################
### LEAF LITTER AND GAMMARUS BIOMASSES SCENARIOS OUTPUTS ###
############################################################

# Bind dataframes of all scenarios
DataLP=rbind(DataL0,DataL1,DataL2)
DataGP=rbind(DataG0,DataG1,DataG2)
DataLP$Temperature=as.numeric(as.character(DataLP$Temperature))
DataGP$Temperature=as.numeric(as.character(DataGP$Temperature))

# Create a predicted dataframe
Temperature=expand.grid(Temperature=seq(5,25,by=0.5))
Scenario=c(rep("SD",41),rep("TSRA",41),rep("TSRM",41))
PredL=data.frame(Temperature,Scenario)
PredG=data.frame(Temperature,Scenario)

# Predict on short temperature intervals
PredL$PredBiomL=data.frame(predict(lmList(MeanBiomL ~ poly(Temperature,2)|Scenario, data=DataLP),PredL))[[1]]
PredL$PredTimeL=data.frame(predict(lmList(MeanTimeL ~ poly(Temperature,2)|Scenario, data=DataLP),PredL))[[1]]
PredL$PredMaxiL=data.frame(predict(lmList(MeanMaxiL ~ poly(Temperature,2)|Scenario, data=DataLP),PredL))[[1]]
PredL$PredSlopL=data.frame(predict(lmList(MeanSlopL ~ poly(Temperature,2)|Scenario, data=DataLP),PredL))[[1]]

PredG$PredBiomG=data.frame(predict(lmList(MeanBiomG ~ poly(Temperature,2)|Scenario, data=DataGP),PredG))[[1]]
PredG$PredTimeG=data.frame(predict(lmList(MeanTimeG ~ poly(Temperature,2)|Scenario, data=DataGP),PredG))[[1]]
PredG$PredMaxiG=data.frame(predict(lmList(MeanMaxiG ~ poly(Temperature,2)|Scenario, data=DataGP),PredG))[[1]]
PredG$PredSlopG=data.frame(predict(lmList(MeanSlopG ~ poly(Temperature,2)|Scenario, data=DataGP),PredG))[[1]]


####################################
### Panel of population dynamics ###
####################################

tiff('Population Dynamics.tiff', units="in", width=15, height=15, res=1000)
par(mfrow=c(3,3))
par(oma=c(4,4,4,4))
par(mar=c(3,3,0,0))

matplot(Test5SD[,-c(1,4,5)]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5 °C", cex=1.5, side=3, line=1, at=1277.5)
matplot(Test15SD[,-c(1,4,5)]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15 °C", cex=1.5, side=3, line=1, at=1277.5)
matplot(Test25SD[,-c(1,4,5)]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25 °C", cex=1.5, side=3, line=1, at=1277.5)
mtext(bquote(TSR['R']), cex=1.5, side=4, line=2, at=2)
matplot(Test5TSRA[,-c(1,4,5)]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
matplot(Test15TSRA[,-c(1,4,5)]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
matplot(Test25TSRA[,-c(1,4,5)]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext(bquote(TSR['A']), cex=1.5, side=4, line=2, at=2)
matplot(Test5TSRM[,-c(1,4,5)]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
matplot(Test15TSRM[,-c(1,4,5)]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
matplot(Test25TSRM[,-c(1,4,5)]/10^5, type="l", ylim=c(0,6), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext(bquote(TSR['M']), cex=1.5, side=4, line=2, at=2)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))
dev.off()


##########################
### Plotting biomasses ###
##########################

# Plot leaf litter mean annual biomasses
PlotBiomL=ggplot(DataLP, aes(x=Temperature, y=MeanBiomL/10^4, group=Scenario)) +
  geom_point(aes(color=Scenario), size=4, pch=16, position=position_dodge(0.5)) +
  geom_line(data=PredL, aes(x=Temperature, y=PredBiomL/10^4, color=Scenario, linetype=Scenario), size=1.0) +
  geom_errorbar(aes(ymin=(MeanBiomL-SdBiomL)/10^4, ymax=(MeanBiomL+SdBiomL)/10^4), color="grey50", size=1.0, width=1.0, position=position_dodge(0.5)) +
  ylab(expression('Mean leaf litter biomass'~'('*10^4~mg~C~m^-2*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,25,by=5), limits=c(4,26)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  scale_color_manual(values=c(TSRA="steelblue2", TSRM="tomato2", SD="black")) +
  scale_linetype_manual(values=c(TSRA="solid", TSRM="solid", SD="dotted")) +
  theme(legend.position="none")

# Plot Gammarus mean annual biomasses
PlotBiomG=ggplot(DataGP, aes(x=Temperature, y=MeanBiomG/10^4, group=Scenario)) +
  geom_point(aes(color=Scenario), size=4, pch=16, position=position_dodge(0.5)) +
  geom_line(data=PredG, aes(x=Temperature, y=PredBiomG/10^4, color=Scenario, linetype=Scenario), size=1.0, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=(MeanBiomG-SdBiomG)/10^4, ymax=(MeanBiomG+SdBiomG)/10^4), color="grey50", size=1.0, width=1.0, position=position_dodge(0.5)) +
  ylab(expression('Mean'*~italic(Gammarus)*~'biomass'~'('*10^4~mg~C~m^-2*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,3,by=1), limits=c(0,3)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,25,by=5), limits=c(4,26)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  scale_color_manual(values=c(TSRA="steelblue2", TSRM="tomato2", SD="black")) +
  scale_linetype_manual(values=c(TSRA="solid", TSRM="solid", SD="dotted")) +
  theme(legend.position="none")

tiff('Mean Biomass Scenarios.tiff', units="in", width=15, height=8, res=1000)
Panel=plot_grid(PlotBiomL, PlotBiomG, align="h", vjust=1, nrow=1, ncol=2)
Xaxis=textGrob(expression('Temperature (°C)'), gp=gpar(fontface="bold", fontsize=20))
grid.arrange(arrangeGrob(Panel, bottom=Xaxis))
dev.off()


##################################
### Plotting persistence times ###
##################################

# Plot leaf litter mean annual persistence times
PlotTimeL=ggplot(DataLP, aes(x=Temperature, y=MeanTimeL, group=Scenario)) +
  geom_point(aes(color=Scenario), size=4, pch=16, position=position_dodge(0.5)) +
  geom_line(data=PredL, aes(x=Temperature, y=PredTimeL, color=Scenario, linetype=Scenario), size=1.0, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=MeanTimeL-SdTimeL, ymax=MeanTimeL+SdTimeL), color="grey50", size=1.0, width=1.0, position=position_dodge(0.5)) +
  ylab(expression('Leaf litter persistence time (d)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,100,by=25), limits=c(0,100)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,25,by=5), limits=c(4,26)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  scale_color_manual(values=c(TSRA="steelblue2", TSRM="tomato2", SD="black")) +
  scale_linetype_manual(values=c(TSRA="solid", TSRM="solid", SD="dotted")) +
  theme(legend.position="none")

# Plot Gammarus mean annual persistence times
PlotTimeG=ggplot(DataGP, aes(x=Temperature, y=MeanTimeG, group=Scenario)) +
  geom_point(aes(color=Scenario), size=4, pch=16, position=position_dodge(0.5)) +
  geom_line(data=PredG, aes(x=Temperature, y=PredTimeG, color=Scenario, linetype=Scenario), size=1.0, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=MeanTimeG-SdTimeG, ymax=MeanTimeG+SdTimeG), color="grey50", size=1.0, width=1.0, position=position_dodge(0.5)) +
  ylab(expression(italic(Gammarus)*~'persistence time (d)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,300,by=100), limits=c(0,300)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,25,by=5), limits=c(4,26)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  scale_color_manual(values=c(TSRA="steelblue2", TSRM="tomato2", SD="black")) +
  scale_linetype_manual(values=c(TSRA="solid", TSRM="solid", SD="dotted")) +
  theme(legend.position="none")

tiff('Persistence Time Scenarios.tiff', units="in", width=15, height=8, res=1000)
Panel=plot_grid(PlotTimeL, PlotTimeG, align="h", vjust=1, nrow=1, ncol=2)
Xaxis=textGrob(expression('Temperature (°C)'), gp=gpar(fontface="bold", fontsize=20))
grid.arrange(arrangeGrob(Panel, bottom=Xaxis))
dev.off()


##########################################
### Plotting biomass decreasing slopes ###
##########################################

# Plot leaf litter mean annual biomass slopes
PlotSlopL=ggplot(DataLP, aes(x=Temperature, y=abs(MeanSlopL)/10^4, group=Scenario)) +
  geom_point(aes(color=Scenario), position=position_dodge(0.5), size=4, pch=16) +
  geom_line(data=PredL, aes(x=Temperature, y=abs(PredSlopL)/10^4, color=Scenario, linetype=Scenario), size=1.0, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=(abs(MeanSlopL)-abs(SdSlopL))/10^4, ymax=(abs(MeanSlopL)+abs(SdSlopL))/10^4), color="grey50", size=1.0, width=1.0, position=position_dodge(0.5)) +
  ylab(expression('Mean leaf litter biomass slope'~'('*10^4~mg~C~m^-2~day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.5,by=0.5), limits=c(0,1.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,25,by=5), limits=c(4,26)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  scale_color_manual(values=c(TSRA="steelblue2", TSRM="tomato2", SD="black")) +
  scale_linetype_manual(values=c(TSRA="solid", TSRM="solid", SD="dotted")) +
  theme(legend.position="none")

# Plot Gammarus mean annual biomass slopes
PlotSlopG=ggplot(DataGP, aes(x=Temperature, y=abs(MeanSlopG)/10^3, group=Scenario)) +
  geom_point(aes(color=Scenario), position=position_dodge(0.5), size=4, pch=16) +
  geom_line(data=PredG, aes(x=Temperature, y=abs(PredSlopG)/10^3, color=Scenario, linetype=Scenario), size=1.0, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=(abs(MeanSlopG)-abs(SdSlopG))/10^3, ymax=(abs(MeanSlopG)+abs(SdSlopG))/10^3), color="grey50", size=1.0, width=1.0, position=position_dodge(0.5)) +
  ylab(expression('Mean'*~italic(Gammarus)*~'biomass slope'~'('*10^3~mg~C~m^-2~day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.5,by=0.5), limits=c(0,1.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,25,by=5), limits=c(4,26)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  scale_color_manual(values=c(TSRA="steelblue2", TSRM="tomato2", SD="black")) +
  scale_linetype_manual(values=c(TSRA="solid", TSRM="solid", SD="dotted")) +
  theme(legend.position="none")

tiff('Mean Biomass Slope Scenarios.tiff', units="in", width=15, height=8, res=1000)
Panel=plot_grid(PlotSlopL, PlotSlopG, align="h", vjust=1, nrow=1, ncol=2)
Xaxis=textGrob(expression('Temperature (°C)'), gp=gpar(fontface="bold", fontsize=20))
grid.arrange(arrangeGrob(Panel, bottom=Xaxis))
dev.off()




####################################################################################################
####################################################################################################




#####################################################################################
### SUPPLEMENTARY ANALYSIS 1: TESTING LEAF LITTER AND GAMMARUS BIOMASS VARIATIONS ###
#####################################################################################

#################################################
### SCENARIO 3: LEAF LITTER BIOMASS REDUCTION ###
#################################################

#################################################
### SUBSCENARIO 1: GAMMARUS BIOMASS REDUCTION ###
#################################################

# Leaf fall = 100 gC/m2/an = 100 000 mgC/m2/an
# Gammarus density = 10 mgDM/m2 = 5 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5LRGR=GammLeafModel(Temp=5,GammMass=TSRA(5,4.26),Leaf=100000/Days,Gamm=5)
Test5LRGR=Test5LRGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test5LRGR=Test5LRGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5LRGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5°C", cex=1.5, side=3, line=-3, at=300)

Test10LRGR=GammLeafModel(Temp=10,GammMass=TSRA(10,4.26),Leaf=100000/Days,Gamm=5)
Test10LRGR=Test10LRGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test10LRGR=Test10LRGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10LRGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10°C", cex=1.5, side=3, line=-3, at=300)

Test15LRGR=GammLeafModel(Temp=15,GammMass=TSRA(15,4.26),Leaf=100000/Days,Gamm=5)
Test15LRGR=Test15LRGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test15LRGR=Test15LRGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15LRGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15°C", cex=1.5, side=3, line=-3, at=300)

Test20LRGR=GammLeafModel(Temp=20,GammMass=TSRA(20,4.26),Leaf=100000/Days,Gamm=5)
Test20LRGR=Test20LRGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test20LRGR=Test20LRGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20LRGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20°C", cex=1.5, side=3, line=-3, at=300)

Test25LRGR=GammLeafModel(Temp=25,GammMass=TSRA(25,4.26),Leaf=100000/Days,Gamm=5)
Test25LRGR=Test25LRGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test25LRGR=Test25LRGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25LRGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25°C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))


# Calculate population metabolism and leaf ingestion
Test5LRGR$M=(MetaQuadra(5,TSRM(5,4.26))/1000)*Test5LRGR[,3]
Test5LRGR$I=(AttackQuadra(5,TSRM(5,4.26))*Test5LRGR[,2]/(1+AttackQuadra(5,TSRM(5,4.26))*1/(IngQuadra(5,TSRM(5,4.26))/1000)*Test5LRGR[,2]))*0.30*Test5LRGR[,3]
Test10LRGR$M=(MetaQuadra(10,TSRM(5,4.26))/1000)*Test10LRGR[,3]
Test10LRGR$I=(AttackQuadra(10,TSRM(5,4.26))*Test10LRGR[,2]/(1+AttackQuadra(10,TSRM(5,4.26))*1/(IngQuadra(10,TSRM(5,4.26))/1000)*Test10LRGR[,2]))*0.30*Test10LRGR[,3]
Test15LRGR$M=(MetaQuadra(15,TSRM(5,4.26))/1000)*Test15LRGR[,3]
Test15LRGR$I=(AttackQuadra(15,TSRM(5,4.26))*Test15LRGR[,2]/(1+AttackQuadra(15,TSRM(5,4.26))*1/(IngQuadra(15,TSRM(5,4.26))/1000)*Test15LRGR[,2]))*0.30*Test15LRGR[,3]
Test20LRGR$M=(MetaQuadra(20,TSRM(5,4.26))/1000)*Test20LRGR[,3]
Test20LRGR$I=(AttackQuadra(20,TSRM(5,4.26))*Test20LRGR[,2]/(1+AttackQuadra(20,TSRM(5,4.26))*1/(IngQuadra(20,TSRM(5,4.26))/1000)*Test20LRGR[,2]))*0.30*Test20LRGR[,3]
Test25LRGR$M=(MetaQuadra(25,TSRM(5,4.26))/1000)*Test25LRGR[,3]
Test25LRGR$I=(AttackQuadra(25,TSRM(5,4.26))*Test25LRGR[,2]/(1+AttackQuadra(25,TSRM(5,4.26))*1/(IngQuadra(25,TSRM(5,4.26))/1000)*Test25LRGR[,2]))*0.30*Test25LRGR[,3]

# Combine the dataframes
TestLRGR=rbind(Test5LRGR,Test10LRGR,Test15LRGR,Test20LRGR,Test25LRGR)
TestLRGR$Temp=rep(c("5","10","15","20","25"), each=2556)
TestLRGR$Temp=factor(TestLRGR$Temp, levels=unique(TestLRGR$Temp))
colnames(TestLRGR)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestLRGR=subset(TestLRGR, !Time==2555); TestLRGR$Year=rep(rep(1:(nrow(TestLRGR)/(365*5)), each=365),5)
CutLRGR=as.data.frame(setDT(TestLRGR)[TestLRGR[, tail(.I, -16), by=Year]$V1]); CutLRGR=CutLRGR[!CutLRGR$Year=="1",]
MeanL=as.data.frame(setDT(CutLRGR)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutLRGR)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLLRGR=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGLRGR=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestLRGR=subset(TestLRGR, !Time==2555); TestLRGR$Year=rep(rep(1:(nrow(TestLRGR)/(365*5)), each=365),5)
CutLRGR=as.data.frame(setDT(TestLRGR)[TestLRGR[, tail(.I, -16), by=Year]$V1]); CutLRGR=CutLRGR[!CutLRGR$Year=="1",]
ThLLRGR=CutLRGR[CutLRGR$L>60000,]; ThGLRGR=CutLRGR[CutLRGR$G>5000,]
TimeLLRGR=as.data.frame(setDT(ThLLRGR)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGLRGR=as.data.frame(setDT(ThGLRGR)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLLRGR=as.data.frame(setDT(TimeLLRGR)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGLRGR=as.data.frame(setDT(TimeGLRGR)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLLRGR=as.data.frame(setDT(TestLRGR)[, .(MaxLLRGR=findPeaks(L), MinLLRGR=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGLRGR=as.data.frame(setDT(TestLRGR)[, .(MaxGLRGR=findPeaks(G), MinGLRGR=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5LRGR=c(CycleLLRGR[1,2]:CycleLLRGR[1,3],CycleLLRGR[2,2]:CycleLLRGR[2,3],CycleLLRGR[3,2]:CycleLLRGR[3,3],CycleLLRGR[4,2]:CycleLLRGR[4,3],CycleLLRGR[5,2]:CycleLLRGR[5,3],CycleLLRGR[6,2]:CycleLLRGR[6,3],CycleLLRGR[7,2]:CycleLLRGR[7,3])+t0
CycleL5LRGR=c(rep("A",length(CycleLLRGR[1,2]:CycleLLRGR[1,3])),rep("B",length(CycleLLRGR[2,2]:CycleLLRGR[2,3])),rep("C",length(CycleLLRGR[3,2]:CycleLLRGR[3,3])),rep("D",length(CycleLLRGR[4,2]:CycleLLRGR[4,3])),rep("E",length(CycleLLRGR[5,2]:CycleLLRGR[5,3])),rep("F",length(CycleLLRGR[6,2]:CycleLLRGR[6,3])),rep("G",length(CycleLLRGR[7,2]:CycleLLRGR[7,3])))
CutL10LRGR=c(CycleLLRGR[8,2]:CycleLLRGR[8,3],CycleLLRGR[9,2]:CycleLLRGR[9,3],CycleLLRGR[10,2]:CycleLLRGR[10,3],CycleLLRGR[11,2]:CycleLLRGR[11,3],CycleLLRGR[12,2]:CycleLLRGR[12,3],CycleLLRGR[13,2]:CycleLLRGR[13,3],CycleLLRGR[14,2]:CycleLLRGR[14,3])+t1
CycleL10LRGR=c(rep("A",length(CycleLLRGR[8,2]:CycleLLRGR[8,3])),rep("B",length(CycleLLRGR[9,2]:CycleLLRGR[9,3])),rep("C",length(CycleLLRGR[10,2]:CycleLLRGR[10,3])),rep("D",length(CycleLLRGR[11,2]:CycleLLRGR[11,3])),rep("E",length(CycleLLRGR[12,2]:CycleLLRGR[12,3])),rep("F",length(CycleLLRGR[13,2]:CycleLLRGR[13,3])),rep("G",length(CycleLLRGR[14,2]:CycleLLRGR[14,3])))
CutL15LRGR=c(CycleLLRGR[15,2]:CycleLLRGR[15,3],CycleLLRGR[16,2]:CycleLLRGR[16,3],CycleLLRGR[17,2]:CycleLLRGR[17,3],CycleLLRGR[18,2]:CycleLLRGR[18,3],CycleLLRGR[19,2]:CycleLLRGR[19,3],CycleLLRGR[20,2]:CycleLLRGR[20,3],CycleLLRGR[21,2]:CycleLLRGR[21,3])+t2
CycleL15LRGR=c(rep("A",length(CycleLLRGR[15,2]:CycleLLRGR[15,3])),rep("B",length(CycleLLRGR[16,2]:CycleLLRGR[16,3])),rep("C",length(CycleLLRGR[17,2]:CycleLLRGR[17,3])),rep("D",length(CycleLLRGR[18,2]:CycleLLRGR[18,3])),rep("E",length(CycleLLRGR[19,2]:CycleLLRGR[19,3])),rep("F",length(CycleLLRGR[20,2]:CycleLLRGR[20,3])),rep("G",length(CycleLLRGR[21,2]:CycleLLRGR[21,3])))
CutL20LRGR=c(CycleLLRGR[22,2]:CycleLLRGR[22,3],CycleLLRGR[23,2]:CycleLLRGR[23,3],CycleLLRGR[24,2]:CycleLLRGR[24,3],CycleLLRGR[25,2]:CycleLLRGR[25,3],CycleLLRGR[26,2]:CycleLLRGR[26,3],CycleLLRGR[27,2]:CycleLLRGR[27,3],CycleLLRGR[28,2]:CycleLLRGR[28,3])+t3
CycleL20LRGR=c(rep("A",length(CycleLLRGR[22,2]:CycleLLRGR[22,3])),rep("B",length(CycleLLRGR[23,2]:CycleLLRGR[23,3])),rep("C",length(CycleLLRGR[24,2]:CycleLLRGR[24,3])),rep("D",length(CycleLLRGR[25,2]:CycleLLRGR[25,3])),rep("E",length(CycleLLRGR[26,2]:CycleLLRGR[26,3])),rep("F",length(CycleLLRGR[27,2]:CycleLLRGR[27,3])),rep("G",length(CycleLLRGR[28,2]:CycleLLRGR[28,3])))
CutL25LRGR=c(CycleLLRGR[29,2]:CycleLLRGR[29,3],CycleLLRGR[30,2]:CycleLLRGR[30,3],CycleLLRGR[31,2]:CycleLLRGR[31,3],CycleLLRGR[32,2]:CycleLLRGR[32,3],CycleLLRGR[33,2]:CycleLLRGR[33,3],CycleLLRGR[34,2]:CycleLLRGR[34,3],CycleLLRGR[35,2]:CycleLLRGR[35,3])+t4
CycleL25LRGR=c(rep("A",length(CycleLLRGR[29,2]:CycleLLRGR[29,3])),rep("B",length(CycleLLRGR[30,2]:CycleLLRGR[30,3])),rep("C",length(CycleLLRGR[31,2]:CycleLLRGR[31,3])),rep("D",length(CycleLLRGR[32,2]:CycleLLRGR[32,3])),rep("E",length(CycleLLRGR[33,2]:CycleLLRGR[33,3])),rep("F",length(CycleLLRGR[34,2]:CycleLLRGR[34,3])),rep("G",length(CycleLLRGR[35,2]:CycleLLRGR[35,3])))

CutLLRGR=TestLRGR[c(CutL5LRGR,CutL10LRGR,CutL15LRGR,CutL20LRGR,CutL25LRGR),]
CutLLRGR$Cycle=c(CycleL5LRGR,CycleL10LRGR,CycleL15LRGR,CycleL20LRGR,CycleL25LRGR)

CutG5LRGR=c(CycleGLRGR[1,2]:CycleGLRGR[1,3],CycleGLRGR[2,2]:CycleGLRGR[2,3],CycleGLRGR[3,2]:CycleGLRGR[3,3],CycleGLRGR[4,2]:CycleGLRGR[4,3],CycleGLRGR[5,2]:CycleGLRGR[5,3],CycleGLRGR[6,2]:CycleGLRGR[6,3],CycleGLRGR[7,2]:CycleGLRGR[7,3])+t0
CycleG5LRGR=c(rep("A",length(CycleGLRGR[1,2]:CycleGLRGR[1,3])),rep("B",length(CycleGLRGR[2,2]:CycleGLRGR[2,3])),rep("C",length(CycleGLRGR[3,2]:CycleGLRGR[3,3])),rep("D",length(CycleGLRGR[4,2]:CycleGLRGR[4,3])),rep("E",length(CycleGLRGR[5,2]:CycleGLRGR[5,3])),rep("F",length(CycleGLRGR[6,2]:CycleGLRGR[6,3])),rep("G",length(CycleGLRGR[7,2]:CycleGLRGR[7,3])))
CutG10LRGR=c(CycleGLRGR[8,2]:CycleGLRGR[8,3],CycleGLRGR[9,2]:CycleGLRGR[9,3],CycleGLRGR[10,2]:CycleGLRGR[10,3],CycleGLRGR[11,2]:CycleGLRGR[11,3],CycleGLRGR[12,2]:CycleGLRGR[12,3],CycleGLRGR[13,2]:CycleGLRGR[13,3],CycleGLRGR[14,2]:CycleGLRGR[14,3])+t1
CycleG10LRGR=c(rep("A",length(CycleGLRGR[8,2]:CycleGLRGR[8,3])),rep("B",length(CycleGLRGR[9,2]:CycleGLRGR[9,3])),rep("C",length(CycleGLRGR[10,2]:CycleGLRGR[10,3])),rep("D",length(CycleGLRGR[11,2]:CycleGLRGR[11,3])),rep("E",length(CycleGLRGR[12,2]:CycleGLRGR[12,3])),rep("F",length(CycleGLRGR[13,2]:CycleGLRGR[13,3])),rep("G",length(CycleGLRGR[14,2]:CycleGLRGR[14,3])))
CutG15LRGR=c(CycleGLRGR[15,2]:CycleGLRGR[15,3],CycleGLRGR[16,2]:CycleGLRGR[16,3],CycleGLRGR[17,2]:CycleGLRGR[17,3],CycleGLRGR[18,2]:CycleGLRGR[18,3],CycleGLRGR[19,2]:CycleGLRGR[19,3],CycleGLRGR[20,2]:CycleGLRGR[20,3],CycleGLRGR[21,2]:CycleGLRGR[21,3])+t2
CycleG15LRGR=c(rep("A",length(CycleGLRGR[15,2]:CycleGLRGR[15,3])),rep("B",length(CycleGLRGR[16,2]:CycleGLRGR[16,3])),rep("C",length(CycleGLRGR[17,2]:CycleGLRGR[17,3])),rep("D",length(CycleGLRGR[18,2]:CycleGLRGR[18,3])),rep("E",length(CycleGLRGR[19,2]:CycleGLRGR[19,3])),rep("F",length(CycleGLRGR[20,2]:CycleGLRGR[20,3])),rep("G",length(CycleGLRGR[21,2]:CycleGLRGR[21,3])))
CutG20LRGR=c(CycleGLRGR[22,2]:CycleGLRGR[22,3],CycleGLRGR[23,2]:CycleGLRGR[23,3],CycleGLRGR[24,2]:CycleGLRGR[24,3],CycleGLRGR[25,2]:CycleGLRGR[25,3],CycleGLRGR[26,2]:CycleGLRGR[26,3],CycleGLRGR[27,2]:CycleGLRGR[27,3],CycleGLRGR[28,2]:CycleGLRGR[28,3])+t3
CycleG20LRGR=c(rep("A",length(CycleGLRGR[22,2]:CycleGLRGR[22,3])),rep("B",length(CycleGLRGR[23,2]:CycleGLRGR[23,3])),rep("C",length(CycleGLRGR[24,2]:CycleGLRGR[24,3])),rep("D",length(CycleGLRGR[25,2]:CycleGLRGR[25,3])),rep("E",length(CycleGLRGR[26,2]:CycleGLRGR[26,3])),rep("F",length(CycleGLRGR[27,2]:CycleGLRGR[27,3])),rep("G",length(CycleGLRGR[28,2]:CycleGLRGR[28,3])))
CutG25LRGR=c(CycleGLRGR[29,2]:CycleGLRGR[29,3],CycleGLRGR[30,2]:CycleGLRGR[30,3],CycleGLRGR[31,2]:CycleGLRGR[31,3],CycleGLRGR[32,2]:CycleGLRGR[32,3],CycleGLRGR[33,2]:CycleGLRGR[33,3],CycleGLRGR[34,2]:CycleGLRGR[34,3],CycleGLRGR[35,2]:CycleGLRGR[35,3])+t4
CycleG25LRGR=c(rep("A",length(CycleGLRGR[29,2]:CycleGLRGR[29,3])),rep("B",length(CycleGLRGR[30,2]:CycleGLRGR[30,3])),rep("C",length(CycleGLRGR[31,2]:CycleGLRGR[31,3])),rep("D",length(CycleGLRGR[32,2]:CycleGLRGR[32,3])),rep("E",length(CycleGLRGR[33,2]:CycleGLRGR[33,3])),rep("F",length(CycleGLRGR[34,2]:CycleGLRGR[34,3])),rep("G",length(CycleGLRGR[35,2]:CycleGLRGR[35,3])))

CutGLRGR=TestLRGR[c(CutG5LRGR,CutG10LRGR,CutG15LRGR,CutG20LRGR,CutG25LRGR),]
CutGLRGR$Cycle=c(CycleG5LRGR,CycleG10LRGR,CycleG15LRGR,CycleG20LRGR,CycleG25LRGR)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLLRGR=CutLLRGR[!CutLLRGR$Cycle=="A",]; CutGLRGR=CutGLRGR[!CutGLRGR$Cycle=="A",]
MinLLRGR=setDT(CutLLRGR)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLLRGR=setDT(CutLLRGR)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGLRGR=setDT(CutGLRGR)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGLRGR=setDT(CutGLRGR)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLLRGR=rbind(MinLLRGR[,c(3,1,2,4)],MaxLLRGR[,c(3,1,2,4)])
BiomGLRGR=rbind(MinGLRGR[,c(3,1,2,5)],MaxGLRGR[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLLRGR=as.data.frame(BiomLLRGR %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLLRGR=SlopLLRGR[order(as.numeric(SlopLLRGR$Temp)),]; SlopLLRGR$Slope=unlist(SlopLLRGR$Slope)
BiomLLRGR=cbind(SlopLLRGR,Max=BiomLLRGR[BiomLLRGR$L>100000]$L)
MaxiLLRGR=setDT(BiomLLRGR)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLLRGR=setDT(BiomLLRGR)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGLRGR=as.data.frame(BiomGLRGR %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGLRGR=SlopGLRGR[order(as.numeric(SlopGLRGR$Temp)),]; SlopGLRGR$Slope=unlist(SlopGLRGR$Slope)
BiomGLRGR=cbind(SlopGLRGR,Max=BiomGLRGR[BiomGLRGR$G>10000]$G)
MaxiGLRGR=setDT(BiomGLRGR)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGLRGR=setDT(BiomGLRGR)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL31=data.frame(MeanLLRGR[,1],MeanLLRGR[,2],MeanLLRGR[,3],TimeLLRGR[,2],TimeLLRGR[,3])
colnames(DataL31)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL")
DataL31$Scenario=rep("LRGR",5)
DataL31$Subscenario=rep("GR",5)

DataG31=data.frame(MeanGLRGR[,1],MeanGLRGR[,2],MeanGLRGR[,3],TimeGLRGR[,2],TimeGLRGR[,3])
colnames(DataG31)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG")
DataG31$Scenario=rep("LRGR",5)
DataG31$Subscenario=rep("GR",5)


################################################
### SUBSCENARIO 2: GAMMARUS BIOMASS CONSTANT ###
################################################

# Leaf fall = 100 gC/m2/an = 100 000 mgC/m2/an
# Gammarus density = 30 mgDM/m2 = 15 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5LRGC=GammLeafModel(Temp=5,GammMass=TSRA(5,4.26),Leaf=100000/Days,Gamm=15)
Test5LRGC=Test5LRGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test5LRGC=Test5LRGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5LRGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5°C", cex=1.5, side=3, line=-3, at=300)

Test10LRGC=GammLeafModel(Temp=10,GammMass=TSRA(10,4.26),Leaf=100000/Days,Gamm=15)
Test10LRGC=Test10LRGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test10LRGC=Test10LRGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10LRGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10°C", cex=1.5, side=3, line=-3, at=300)

Test15LRGC=GammLeafModel(Temp=15,GammMass=TSRA(15,4.26),Leaf=100000/Days,Gamm=15)
Test15LRGC=Test15LRGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test15LRGC=Test15LRGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15LRGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15°C", cex=1.5, side=3, line=-3, at=300)

Test20LRGC=GammLeafModel(Temp=20,GammMass=TSRA(20,4.26),Leaf=100000/Days,Gamm=15)
Test20LRGC=Test20LRGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test20LRGC=Test20LRGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20LRGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20°C", cex=1.5, side=3, line=-3, at=300)

Test25LRGC=GammLeafModel(Temp=25,GammMass=TSRA(25,4.26),Leaf=100000/Days,Gamm=15)
Test25LRGC=Test25LRGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test25LRGC=Test25LRGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25LRGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25°C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))


# Calculate population metabolism and leaf ingestion
Test5LRGC$M=(MetaQuadra(5,TSRM(5,4.26))/1000)*Test5LRGC[,3]
Test5LRGC$I=(AttackQuadra(5,TSRM(5,4.26))*Test5LRGC[,2]/(1+AttackQuadra(5,TSRM(5,4.26))*1/(IngQuadra(5,TSRM(5,4.26))/1000)*Test5LRGC[,2]))*0.30*Test5LRGC[,3]
Test10LRGC$M=(MetaQuadra(10,TSRM(5,4.26))/1000)*Test10LRGC[,3]
Test10LRGC$I=(AttackQuadra(10,TSRM(5,4.26))*Test10LRGC[,2]/(1+AttackQuadra(10,TSRM(5,4.26))*1/(IngQuadra(10,TSRM(5,4.26))/1000)*Test10LRGC[,2]))*0.30*Test10LRGC[,3]
Test15LRGC$M=(MetaQuadra(15,TSRM(5,4.26))/1000)*Test15LRGC[,3]
Test15LRGC$I=(AttackQuadra(15,TSRM(5,4.26))*Test15LRGC[,2]/(1+AttackQuadra(15,TSRM(5,4.26))*1/(IngQuadra(15,TSRM(5,4.26))/1000)*Test15LRGC[,2]))*0.30*Test15LRGC[,3]
Test20LRGC$M=(MetaQuadra(20,TSRM(5,4.26))/1000)*Test20LRGC[,3]
Test20LRGC$I=(AttackQuadra(20,TSRM(5,4.26))*Test20LRGC[,2]/(1+AttackQuadra(20,TSRM(5,4.26))*1/(IngQuadra(20,TSRM(5,4.26))/1000)*Test20LRGC[,2]))*0.30*Test20LRGC[,3]
Test25LRGC$M=(MetaQuadra(25,TSRM(5,4.26))/1000)*Test25LRGC[,3]
Test25LRGC$I=(AttackQuadra(25,TSRM(5,4.26))*Test25LRGC[,2]/(1+AttackQuadra(25,TSRM(5,4.26))*1/(IngQuadra(25,TSRM(5,4.26))/1000)*Test25LRGC[,2]))*0.30*Test25LRGC[,3]

# Combine the dataframes
TestLRGC=rbind(Test5LRGC,Test10LRGC,Test15LRGC,Test20LRGC,Test25LRGC)
TestLRGC$Temp=rep(c("5","10","15","20","25"), each=2556)
TestLRGC$Temp=factor(TestLRGC$Temp, levels=unique(TestLRGC$Temp))
colnames(TestLRGC)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestLRGC=subset(TestLRGC, !Time==2555); TestLRGC$Year=rep(rep(1:(nrow(TestLRGC)/(365*5)), each=365),5)
CutLRGC=as.data.frame(setDT(TestLRGC)[TestLRGC[, tail(.I, -16), by=Year]$V1]); CutLRGC=CutLRGC[!CutLRGC$Year=="1",]
MeanL=as.data.frame(setDT(CutLRGC)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutLRGC)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLLRGC=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGLRGC=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestLRGC=subset(TestLRGC, !Time==2555); TestLRGC$Year=rep(rep(1:(nrow(TestLRGC)/(365*5)), each=365),5)
CutLRGC=as.data.frame(setDT(TestLRGC)[TestLRGC[, tail(.I, -16), by=Year]$V1]); CutLRGC=CutLRGC[!CutLRGC$Year=="1",]
ThLLRGC=CutLRGC[CutLRGC$L>60000,]; ThGLRGC=CutLRGC[CutLRGC$G>5000,]
TimeLLRGC=as.data.frame(setDT(ThLLRGC)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGLRGC=as.data.frame(setDT(ThGLRGC)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLLRGC=as.data.frame(setDT(TimeLLRGC)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGLRGC=as.data.frame(setDT(TimeGLRGC)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLLRGC=as.data.frame(setDT(TestLRGC)[, .(MaxLLRGC=findPeaks(L), MinLLRGC=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGLRGC=as.data.frame(setDT(TestLRGC)[, .(MaxGLRGC=findPeaks(G), MinGLRGC=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5LRGC=c(CycleLLRGC[1,2]:CycleLLRGC[1,3],CycleLLRGC[2,2]:CycleLLRGC[2,3],CycleLLRGC[3,2]:CycleLLRGC[3,3],CycleLLRGC[4,2]:CycleLLRGC[4,3],CycleLLRGC[5,2]:CycleLLRGC[5,3],CycleLLRGC[6,2]:CycleLLRGC[6,3],CycleLLRGC[7,2]:CycleLLRGC[7,3])+t0
CycleL5LRGC=c(rep("A",length(CycleLLRGC[1,2]:CycleLLRGC[1,3])),rep("B",length(CycleLLRGC[2,2]:CycleLLRGC[2,3])),rep("C",length(CycleLLRGC[3,2]:CycleLLRGC[3,3])),rep("D",length(CycleLLRGC[4,2]:CycleLLRGC[4,3])),rep("E",length(CycleLLRGC[5,2]:CycleLLRGC[5,3])),rep("F",length(CycleLLRGC[6,2]:CycleLLRGC[6,3])),rep("G",length(CycleLLRGC[7,2]:CycleLLRGC[7,3])))
CutL10LRGC=c(CycleLLRGC[8,2]:CycleLLRGC[8,3],CycleLLRGC[9,2]:CycleLLRGC[9,3],CycleLLRGC[10,2]:CycleLLRGC[10,3],CycleLLRGC[11,2]:CycleLLRGC[11,3],CycleLLRGC[12,2]:CycleLLRGC[12,3],CycleLLRGC[13,2]:CycleLLRGC[13,3],CycleLLRGC[14,2]:CycleLLRGC[14,3])+t1
CycleL10LRGC=c(rep("A",length(CycleLLRGC[8,2]:CycleLLRGC[8,3])),rep("B",length(CycleLLRGC[9,2]:CycleLLRGC[9,3])),rep("C",length(CycleLLRGC[10,2]:CycleLLRGC[10,3])),rep("D",length(CycleLLRGC[11,2]:CycleLLRGC[11,3])),rep("E",length(CycleLLRGC[12,2]:CycleLLRGC[12,3])),rep("F",length(CycleLLRGC[13,2]:CycleLLRGC[13,3])),rep("G",length(CycleLLRGC[14,2]:CycleLLRGC[14,3])))
CutL15LRGC=c(CycleLLRGC[15,2]:CycleLLRGC[15,3],CycleLLRGC[16,2]:CycleLLRGC[16,3],CycleLLRGC[17,2]:CycleLLRGC[17,3],CycleLLRGC[18,2]:CycleLLRGC[18,3],CycleLLRGC[19,2]:CycleLLRGC[19,3],CycleLLRGC[20,2]:CycleLLRGC[20,3],CycleLLRGC[21,2]:CycleLLRGC[21,3])+t2
CycleL15LRGC=c(rep("A",length(CycleLLRGC[15,2]:CycleLLRGC[15,3])),rep("B",length(CycleLLRGC[16,2]:CycleLLRGC[16,3])),rep("C",length(CycleLLRGC[17,2]:CycleLLRGC[17,3])),rep("D",length(CycleLLRGC[18,2]:CycleLLRGC[18,3])),rep("E",length(CycleLLRGC[19,2]:CycleLLRGC[19,3])),rep("F",length(CycleLLRGC[20,2]:CycleLLRGC[20,3])),rep("G",length(CycleLLRGC[21,2]:CycleLLRGC[21,3])))
CutL20LRGC=c(CycleLLRGC[22,2]:CycleLLRGC[22,3],CycleLLRGC[23,2]:CycleLLRGC[23,3],CycleLLRGC[24,2]:CycleLLRGC[24,3],CycleLLRGC[25,2]:CycleLLRGC[25,3],CycleLLRGC[26,2]:CycleLLRGC[26,3],CycleLLRGC[27,2]:CycleLLRGC[27,3],CycleLLRGC[28,2]:CycleLLRGC[28,3])+t3
CycleL20LRGC=c(rep("A",length(CycleLLRGC[22,2]:CycleLLRGC[22,3])),rep("B",length(CycleLLRGC[23,2]:CycleLLRGC[23,3])),rep("C",length(CycleLLRGC[24,2]:CycleLLRGC[24,3])),rep("D",length(CycleLLRGC[25,2]:CycleLLRGC[25,3])),rep("E",length(CycleLLRGC[26,2]:CycleLLRGC[26,3])),rep("F",length(CycleLLRGC[27,2]:CycleLLRGC[27,3])),rep("G",length(CycleLLRGC[28,2]:CycleLLRGC[28,3])))
CutL25LRGC=c(CycleLLRGC[29,2]:CycleLLRGC[29,3],CycleLLRGC[30,2]:CycleLLRGC[30,3],CycleLLRGC[31,2]:CycleLLRGC[31,3],CycleLLRGC[32,2]:CycleLLRGC[32,3],CycleLLRGC[33,2]:CycleLLRGC[33,3],CycleLLRGC[34,2]:CycleLLRGC[34,3],CycleLLRGC[35,2]:CycleLLRGC[35,3])+t4
CycleL25LRGC=c(rep("A",length(CycleLLRGC[29,2]:CycleLLRGC[29,3])),rep("B",length(CycleLLRGC[30,2]:CycleLLRGC[30,3])),rep("C",length(CycleLLRGC[31,2]:CycleLLRGC[31,3])),rep("D",length(CycleLLRGC[32,2]:CycleLLRGC[32,3])),rep("E",length(CycleLLRGC[33,2]:CycleLLRGC[33,3])),rep("F",length(CycleLLRGC[34,2]:CycleLLRGC[34,3])),rep("G",length(CycleLLRGC[35,2]:CycleLLRGC[35,3])))

CutLLRGC=TestLRGC[c(CutL5LRGC,CutL10LRGC,CutL15LRGC,CutL20LRGC,CutL25LRGC),]
CutLLRGC$Cycle=c(CycleL5LRGC,CycleL10LRGC,CycleL15LRGC,CycleL20LRGC,CycleL25LRGC)

CutG5LRGC=c(CycleGLRGC[1,2]:CycleGLRGC[1,3],CycleGLRGC[2,2]:CycleGLRGC[2,3],CycleGLRGC[3,2]:CycleGLRGC[3,3],CycleGLRGC[4,2]:CycleGLRGC[4,3],CycleGLRGC[5,2]:CycleGLRGC[5,3],CycleGLRGC[6,2]:CycleGLRGC[6,3],CycleGLRGC[7,2]:CycleGLRGC[7,3])+t0
CycleG5LRGC=c(rep("A",length(CycleGLRGC[1,2]:CycleGLRGC[1,3])),rep("B",length(CycleGLRGC[2,2]:CycleGLRGC[2,3])),rep("C",length(CycleGLRGC[3,2]:CycleGLRGC[3,3])),rep("D",length(CycleGLRGC[4,2]:CycleGLRGC[4,3])),rep("E",length(CycleGLRGC[5,2]:CycleGLRGC[5,3])),rep("F",length(CycleGLRGC[6,2]:CycleGLRGC[6,3])),rep("G",length(CycleGLRGC[7,2]:CycleGLRGC[7,3])))
CutG10LRGC=c(CycleGLRGC[8,2]:CycleGLRGC[8,3],CycleGLRGC[9,2]:CycleGLRGC[9,3],CycleGLRGC[10,2]:CycleGLRGC[10,3],CycleGLRGC[11,2]:CycleGLRGC[11,3],CycleGLRGC[12,2]:CycleGLRGC[12,3],CycleGLRGC[13,2]:CycleGLRGC[13,3],CycleGLRGC[14,2]:CycleGLRGC[14,3])+t1
CycleG10LRGC=c(rep("A",length(CycleGLRGC[8,2]:CycleGLRGC[8,3])),rep("B",length(CycleGLRGC[9,2]:CycleGLRGC[9,3])),rep("C",length(CycleGLRGC[10,2]:CycleGLRGC[10,3])),rep("D",length(CycleGLRGC[11,2]:CycleGLRGC[11,3])),rep("E",length(CycleGLRGC[12,2]:CycleGLRGC[12,3])),rep("F",length(CycleGLRGC[13,2]:CycleGLRGC[13,3])),rep("G",length(CycleGLRGC[14,2]:CycleGLRGC[14,3])))
CutG15LRGC=c(CycleGLRGC[15,2]:CycleGLRGC[15,3],CycleGLRGC[16,2]:CycleGLRGC[16,3],CycleGLRGC[17,2]:CycleGLRGC[17,3],CycleGLRGC[18,2]:CycleGLRGC[18,3],CycleGLRGC[19,2]:CycleGLRGC[19,3],CycleGLRGC[20,2]:CycleGLRGC[20,3],CycleGLRGC[21,2]:CycleGLRGC[21,3])+t2
CycleG15LRGC=c(rep("A",length(CycleGLRGC[15,2]:CycleGLRGC[15,3])),rep("B",length(CycleGLRGC[16,2]:CycleGLRGC[16,3])),rep("C",length(CycleGLRGC[17,2]:CycleGLRGC[17,3])),rep("D",length(CycleGLRGC[18,2]:CycleGLRGC[18,3])),rep("E",length(CycleGLRGC[19,2]:CycleGLRGC[19,3])),rep("F",length(CycleGLRGC[20,2]:CycleGLRGC[20,3])),rep("G",length(CycleGLRGC[21,2]:CycleGLRGC[21,3])))
CutG20LRGC=c(CycleGLRGC[22,2]:CycleGLRGC[22,3],CycleGLRGC[23,2]:CycleGLRGC[23,3],CycleGLRGC[24,2]:CycleGLRGC[24,3],CycleGLRGC[25,2]:CycleGLRGC[25,3],CycleGLRGC[26,2]:CycleGLRGC[26,3],CycleGLRGC[27,2]:CycleGLRGC[27,3],CycleGLRGC[28,2]:CycleGLRGC[28,3])+t3
CycleG20LRGC=c(rep("A",length(CycleGLRGC[22,2]:CycleGLRGC[22,3])),rep("B",length(CycleGLRGC[23,2]:CycleGLRGC[23,3])),rep("C",length(CycleGLRGC[24,2]:CycleGLRGC[24,3])),rep("D",length(CycleGLRGC[25,2]:CycleGLRGC[25,3])),rep("E",length(CycleGLRGC[26,2]:CycleGLRGC[26,3])),rep("F",length(CycleGLRGC[27,2]:CycleGLRGC[27,3])),rep("G",length(CycleGLRGC[28,2]:CycleGLRGC[28,3])))
CutG25LRGC=c(CycleGLRGC[29,2]:CycleGLRGC[29,3],CycleGLRGC[30,2]:CycleGLRGC[30,3],CycleGLRGC[31,2]:CycleGLRGC[31,3],CycleGLRGC[32,2]:CycleGLRGC[32,3],CycleGLRGC[33,2]:CycleGLRGC[33,3],CycleGLRGC[34,2]:CycleGLRGC[34,3],CycleGLRGC[35,2]:CycleGLRGC[35,3])+t4
CycleG25LRGC=c(rep("A",length(CycleGLRGC[29,2]:CycleGLRGC[29,3])),rep("B",length(CycleGLRGC[30,2]:CycleGLRGC[30,3])),rep("C",length(CycleGLRGC[31,2]:CycleGLRGC[31,3])),rep("D",length(CycleGLRGC[32,2]:CycleGLRGC[32,3])),rep("E",length(CycleGLRGC[33,2]:CycleGLRGC[33,3])),rep("F",length(CycleGLRGC[34,2]:CycleGLRGC[34,3])),rep("G",length(CycleGLRGC[35,2]:CycleGLRGC[35,3])))

CutGLRGC=TestLRGC[c(CutG5LRGC,CutG10LRGC,CutG15LRGC,CutG20LRGC,CutG25LRGC),]
CutGLRGC$Cycle=c(CycleG5LRGC,CycleG10LRGC,CycleG15LRGC,CycleG20LRGC,CycleG25LRGC)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLLRGC=CutLLRGC[!CutLLRGC$Cycle=="A",]; CutGLRGC=CutGLRGC[!CutGLRGC$Cycle=="A",]
MinLLRGC=setDT(CutLLRGC)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLLRGC=setDT(CutLLRGC)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGLRGC=setDT(CutGLRGC)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGLRGC=setDT(CutGLRGC)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLLRGC=rbind(MinLLRGC[,c(3,1,2,4)],MaxLLRGC[,c(3,1,2,4)])
BiomGLRGC=rbind(MinGLRGC[,c(3,1,2,5)],MaxGLRGC[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLLRGC=as.data.frame(BiomLLRGC %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLLRGC=SlopLLRGC[order(as.numeric(SlopLLRGC$Temp)),]; SlopLLRGC$Slope=unlist(SlopLLRGC$Slope)
BiomLLRGC=cbind(SlopLLRGC,Max=BiomLLRGC[BiomLLRGC$L>100000]$L)
MaxiLLRGC=setDT(BiomLLRGC)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLLRGC=setDT(BiomLLRGC)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGLRGC=as.data.frame(BiomGLRGC %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGLRGC=SlopGLRGC[order(as.numeric(SlopGLRGC$Temp)),]; SlopGLRGC$Slope=unlist(SlopGLRGC$Slope)
BiomGLRGC=cbind(SlopGLRGC,Max=BiomGLRGC[BiomGLRGC$G>10000]$G)
MaxiGLRGC=setDT(BiomGLRGC)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGLRGC=setDT(BiomGLRGC)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL32=data.frame(MeanLLRGC[,1],MeanLLRGC[,2],MeanLLRGC[,3],TimeLLRGC[,2],TimeLLRGC[,3])
colnames(DataL32)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL")
DataL32$Scenario=rep("LRGC",5)
DataL32$Subscenario=rep("GC",5)

DataG32=data.frame(MeanGLRGC[,1],MeanGLRGC[,2],MeanGLRGC[,3],TimeGLRGC[,2],TimeGLRGC[,3])
colnames(DataG32)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG")
DataG32$Scenario=rep("LRGC",5)
DataG32$Subscenario=rep("GC",5)


################################################
### SUBSCENARIO 3: GAMMARUS BIOMASS INCREASE ###
################################################

# Leaf fall = 100 gC/m2/an = 100 000 mgC/m2/an
# Gammarus density = 90 mgDM/m2 = 45 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5LRGI=GammLeafModel(Temp=5,GammMass=TSRA(5,4.26),Leaf=100000/Days,Gamm=45)
Test5LRGI=Test5LRGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test5LRGI=Test5LRGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5LRGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5°C", cex=1.5, side=3, line=-3, at=300)

Test10LRGI=GammLeafModel(Temp=10,GammMass=TSRA(10,4.26),Leaf=100000/Days,Gamm=45)
Test10LRGI=Test10LRGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test10LRGI=Test10LRGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10LRGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10°C", cex=1.5, side=3, line=-3, at=300)

Test15LRGI=GammLeafModel(Temp=15,GammMass=TSRA(15,4.26),Leaf=100000/Days,Gamm=45)
Test15LRGI=Test15LRGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test15LRGI=Test15LRGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15LRGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15°C", cex=1.5, side=3, line=-3, at=300)

Test20LRGI=GammLeafModel(Temp=20,GammMass=TSRA(20,4.26),Leaf=100000/Days,Gamm=45)
Test20LRGI=Test20LRGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test20LRGI=Test20LRGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20LRGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20°C", cex=1.5, side=3, line=-3, at=300)

Test25LRGI=GammLeafModel(Temp=25,GammMass=TSRA(25,4.26),Leaf=100000/Days,Gamm=45)
Test25LRGI=Test25LRGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test25LRGI=Test25LRGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25LRGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25°C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))


# Calculate population metabolism and leaf ingestion
Test5LRGI$M=(MetaQuadra(5,TSRM(5,4.26))/1000)*Test5LRGI[,3]
Test5LRGI$I=(AttackQuadra(5,TSRM(5,4.26))*Test5LRGI[,2]/(1+AttackQuadra(5,TSRM(5,4.26))*1/(IngQuadra(5,TSRM(5,4.26))/1000)*Test5LRGI[,2]))*0.30*Test5LRGI[,3]
Test10LRGI$M=(MetaQuadra(10,TSRM(5,4.26))/1000)*Test10LRGI[,3]
Test10LRGI$I=(AttackQuadra(10,TSRM(5,4.26))*Test10LRGI[,2]/(1+AttackQuadra(10,TSRM(5,4.26))*1/(IngQuadra(10,TSRM(5,4.26))/1000)*Test10LRGI[,2]))*0.30*Test10LRGI[,3]
Test15LRGI$M=(MetaQuadra(15,TSRM(5,4.26))/1000)*Test15LRGI[,3]
Test15LRGI$I=(AttackQuadra(15,TSRM(5,4.26))*Test15LRGI[,2]/(1+AttackQuadra(15,TSRM(5,4.26))*1/(IngQuadra(15,TSRM(5,4.26))/1000)*Test15LRGI[,2]))*0.30*Test15LRGI[,3]
Test20LRGI$M=(MetaQuadra(20,TSRM(5,4.26))/1000)*Test20LRGI[,3]
Test20LRGI$I=(AttackQuadra(20,TSRM(5,4.26))*Test20LRGI[,2]/(1+AttackQuadra(20,TSRM(5,4.26))*1/(IngQuadra(20,TSRM(5,4.26))/1000)*Test20LRGI[,2]))*0.30*Test20LRGI[,3]
Test25LRGI$M=(MetaQuadra(25,TSRM(5,4.26))/1000)*Test25LRGI[,3]
Test25LRGI$I=(AttackQuadra(25,TSRM(5,4.26))*Test25LRGI[,2]/(1+AttackQuadra(25,TSRM(5,4.26))*1/(IngQuadra(25,TSRM(5,4.26))/1000)*Test25LRGI[,2]))*0.30*Test25LRGI[,3]

# Combine the dataframes
TestLRGI=rbind(Test5LRGI,Test10LRGI,Test15LRGI,Test20LRGI,Test25LRGI)
TestLRGI$Temp=rep(c("5","10","15","20","25"), each=2556)
TestLRGI$Temp=factor(TestLRGI$Temp, levels=unique(TestLRGI$Temp))
colnames(TestLRGI)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestLRGI=subset(TestLRGI, !Time==2555); TestLRGI$Year=rep(rep(1:(nrow(TestLRGI)/(365*5)), each=365),5)
CutLRGI=as.data.frame(setDT(TestLRGI)[TestLRGI[, tail(.I, -16), by=Year]$V1]); CutLRGI=CutLRGI[!CutLRGI$Year=="1",]
MeanL=as.data.frame(setDT(CutLRGI)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutLRGI)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLLRGI=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGLRGI=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestLRGI=subset(TestLRGI, !Time==2555); TestLRGI$Year=rep(rep(1:(nrow(TestLRGI)/(365*5)), each=365),5)
CutLRGI=as.data.frame(setDT(TestLRGI)[TestLRGI[, tail(.I, -16), by=Year]$V1]); CutLRGI=CutLRGI[!CutLRGI$Year=="1",]
ThLLRGI=CutLRGI[CutLRGI$L>60000,]; ThGLRGI=CutLRGI[CutLRGI$G>5000,]
TimeLLRGI=as.data.frame(setDT(ThLLRGI)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGLRGI=as.data.frame(setDT(ThGLRGI)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLLRGI=as.data.frame(setDT(TimeLLRGI)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGLRGI=as.data.frame(setDT(TimeGLRGI)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLLRGI=as.data.frame(setDT(TestLRGI)[, .(MaxLLRGI=findPeaks(L), MinLLRGI=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGLRGI=as.data.frame(setDT(TestLRGI)[, .(MaxGLRGI=findPeaks(G), MinGLRGI=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5LRGI=c(CycleLLRGI[1,2]:CycleLLRGI[1,3],CycleLLRGI[2,2]:CycleLLRGI[2,3],CycleLLRGI[3,2]:CycleLLRGI[3,3],CycleLLRGI[4,2]:CycleLLRGI[4,3],CycleLLRGI[5,2]:CycleLLRGI[5,3],CycleLLRGI[6,2]:CycleLLRGI[6,3],CycleLLRGI[7,2]:CycleLLRGI[7,3])+t0
CycleL5LRGI=c(rep("A",length(CycleLLRGI[1,2]:CycleLLRGI[1,3])),rep("B",length(CycleLLRGI[2,2]:CycleLLRGI[2,3])),rep("C",length(CycleLLRGI[3,2]:CycleLLRGI[3,3])),rep("D",length(CycleLLRGI[4,2]:CycleLLRGI[4,3])),rep("E",length(CycleLLRGI[5,2]:CycleLLRGI[5,3])),rep("F",length(CycleLLRGI[6,2]:CycleLLRGI[6,3])),rep("G",length(CycleLLRGI[7,2]:CycleLLRGI[7,3])))
CutL10LRGI=c(CycleLLRGI[8,2]:CycleLLRGI[8,3],CycleLLRGI[9,2]:CycleLLRGI[9,3],CycleLLRGI[10,2]:CycleLLRGI[10,3],CycleLLRGI[11,2]:CycleLLRGI[11,3],CycleLLRGI[12,2]:CycleLLRGI[12,3],CycleLLRGI[13,2]:CycleLLRGI[13,3],CycleLLRGI[14,2]:CycleLLRGI[14,3])+t1
CycleL10LRGI=c(rep("A",length(CycleLLRGI[8,2]:CycleLLRGI[8,3])),rep("B",length(CycleLLRGI[9,2]:CycleLLRGI[9,3])),rep("C",length(CycleLLRGI[10,2]:CycleLLRGI[10,3])),rep("D",length(CycleLLRGI[11,2]:CycleLLRGI[11,3])),rep("E",length(CycleLLRGI[12,2]:CycleLLRGI[12,3])),rep("F",length(CycleLLRGI[13,2]:CycleLLRGI[13,3])),rep("G",length(CycleLLRGI[14,2]:CycleLLRGI[14,3])))
CutL15LRGI=c(CycleLLRGI[15,2]:CycleLLRGI[15,3],CycleLLRGI[16,2]:CycleLLRGI[16,3],CycleLLRGI[17,2]:CycleLLRGI[17,3],CycleLLRGI[18,2]:CycleLLRGI[18,3],CycleLLRGI[19,2]:CycleLLRGI[19,3],CycleLLRGI[20,2]:CycleLLRGI[20,3],CycleLLRGI[21,2]:CycleLLRGI[21,3])+t2
CycleL15LRGI=c(rep("A",length(CycleLLRGI[15,2]:CycleLLRGI[15,3])),rep("B",length(CycleLLRGI[16,2]:CycleLLRGI[16,3])),rep("C",length(CycleLLRGI[17,2]:CycleLLRGI[17,3])),rep("D",length(CycleLLRGI[18,2]:CycleLLRGI[18,3])),rep("E",length(CycleLLRGI[19,2]:CycleLLRGI[19,3])),rep("F",length(CycleLLRGI[20,2]:CycleLLRGI[20,3])),rep("G",length(CycleLLRGI[21,2]:CycleLLRGI[21,3])))
CutL20LRGI=c(CycleLLRGI[22,2]:CycleLLRGI[22,3],CycleLLRGI[23,2]:CycleLLRGI[23,3],CycleLLRGI[24,2]:CycleLLRGI[24,3],CycleLLRGI[25,2]:CycleLLRGI[25,3],CycleLLRGI[26,2]:CycleLLRGI[26,3],CycleLLRGI[27,2]:CycleLLRGI[27,3],CycleLLRGI[28,2]:CycleLLRGI[28,3])+t3
CycleL20LRGI=c(rep("A",length(CycleLLRGI[22,2]:CycleLLRGI[22,3])),rep("B",length(CycleLLRGI[23,2]:CycleLLRGI[23,3])),rep("C",length(CycleLLRGI[24,2]:CycleLLRGI[24,3])),rep("D",length(CycleLLRGI[25,2]:CycleLLRGI[25,3])),rep("E",length(CycleLLRGI[26,2]:CycleLLRGI[26,3])),rep("F",length(CycleLLRGI[27,2]:CycleLLRGI[27,3])),rep("G",length(CycleLLRGI[28,2]:CycleLLRGI[28,3])))
CutL25LRGI=c(CycleLLRGI[29,2]:CycleLLRGI[29,3],CycleLLRGI[30,2]:CycleLLRGI[30,3],CycleLLRGI[31,2]:CycleLLRGI[31,3],CycleLLRGI[32,2]:CycleLLRGI[32,3],CycleLLRGI[33,2]:CycleLLRGI[33,3],CycleLLRGI[34,2]:CycleLLRGI[34,3],CycleLLRGI[35,2]:CycleLLRGI[35,3])+t4
CycleL25LRGI=c(rep("A",length(CycleLLRGI[29,2]:CycleLLRGI[29,3])),rep("B",length(CycleLLRGI[30,2]:CycleLLRGI[30,3])),rep("C",length(CycleLLRGI[31,2]:CycleLLRGI[31,3])),rep("D",length(CycleLLRGI[32,2]:CycleLLRGI[32,3])),rep("E",length(CycleLLRGI[33,2]:CycleLLRGI[33,3])),rep("F",length(CycleLLRGI[34,2]:CycleLLRGI[34,3])),rep("G",length(CycleLLRGI[35,2]:CycleLLRGI[35,3])))

CutLLRGI=TestLRGI[c(CutL5LRGI,CutL10LRGI,CutL15LRGI,CutL20LRGI,CutL25LRGI),]
CutLLRGI$Cycle=c(CycleL5LRGI,CycleL10LRGI,CycleL15LRGI,CycleL20LRGI,CycleL25LRGI)

CutG5LRGI=c(CycleGLRGI[1,2]:CycleGLRGI[1,3],CycleGLRGI[2,2]:CycleGLRGI[2,3],CycleGLRGI[3,2]:CycleGLRGI[3,3],CycleGLRGI[4,2]:CycleGLRGI[4,3],CycleGLRGI[5,2]:CycleGLRGI[5,3],CycleGLRGI[6,2]:CycleGLRGI[6,3],CycleGLRGI[7,2]:CycleGLRGI[7,3])+t0
CycleG5LRGI=c(rep("A",length(CycleGLRGI[1,2]:CycleGLRGI[1,3])),rep("B",length(CycleGLRGI[2,2]:CycleGLRGI[2,3])),rep("C",length(CycleGLRGI[3,2]:CycleGLRGI[3,3])),rep("D",length(CycleGLRGI[4,2]:CycleGLRGI[4,3])),rep("E",length(CycleGLRGI[5,2]:CycleGLRGI[5,3])),rep("F",length(CycleGLRGI[6,2]:CycleGLRGI[6,3])),rep("G",length(CycleGLRGI[7,2]:CycleGLRGI[7,3])))
CutG10LRGI=c(CycleGLRGI[8,2]:CycleGLRGI[8,3],CycleGLRGI[9,2]:CycleGLRGI[9,3],CycleGLRGI[10,2]:CycleGLRGI[10,3],CycleGLRGI[11,2]:CycleGLRGI[11,3],CycleGLRGI[12,2]:CycleGLRGI[12,3],CycleGLRGI[13,2]:CycleGLRGI[13,3],CycleGLRGI[14,2]:CycleGLRGI[14,3])+t1
CycleG10LRGI=c(rep("A",length(CycleGLRGI[8,2]:CycleGLRGI[8,3])),rep("B",length(CycleGLRGI[9,2]:CycleGLRGI[9,3])),rep("C",length(CycleGLRGI[10,2]:CycleGLRGI[10,3])),rep("D",length(CycleGLRGI[11,2]:CycleGLRGI[11,3])),rep("E",length(CycleGLRGI[12,2]:CycleGLRGI[12,3])),rep("F",length(CycleGLRGI[13,2]:CycleGLRGI[13,3])),rep("G",length(CycleGLRGI[14,2]:CycleGLRGI[14,3])))
CutG15LRGI=c(CycleGLRGI[15,2]:CycleGLRGI[15,3],CycleGLRGI[16,2]:CycleGLRGI[16,3],CycleGLRGI[17,2]:CycleGLRGI[17,3],CycleGLRGI[18,2]:CycleGLRGI[18,3],CycleGLRGI[19,2]:CycleGLRGI[19,3],CycleGLRGI[20,2]:CycleGLRGI[20,3],CycleGLRGI[21,2]:CycleGLRGI[21,3])+t2
CycleG15LRGI=c(rep("A",length(CycleGLRGI[15,2]:CycleGLRGI[15,3])),rep("B",length(CycleGLRGI[16,2]:CycleGLRGI[16,3])),rep("C",length(CycleGLRGI[17,2]:CycleGLRGI[17,3])),rep("D",length(CycleGLRGI[18,2]:CycleGLRGI[18,3])),rep("E",length(CycleGLRGI[19,2]:CycleGLRGI[19,3])),rep("F",length(CycleGLRGI[20,2]:CycleGLRGI[20,3])),rep("G",length(CycleGLRGI[21,2]:CycleGLRGI[21,3])))
CutG20LRGI=c(CycleGLRGI[22,2]:CycleGLRGI[22,3],CycleGLRGI[23,2]:CycleGLRGI[23,3],CycleGLRGI[24,2]:CycleGLRGI[24,3],CycleGLRGI[25,2]:CycleGLRGI[25,3],CycleGLRGI[26,2]:CycleGLRGI[26,3],CycleGLRGI[27,2]:CycleGLRGI[27,3],CycleGLRGI[28,2]:CycleGLRGI[28,3])+t3
CycleG20LRGI=c(rep("A",length(CycleGLRGI[22,2]:CycleGLRGI[22,3])),rep("B",length(CycleGLRGI[23,2]:CycleGLRGI[23,3])),rep("C",length(CycleGLRGI[24,2]:CycleGLRGI[24,3])),rep("D",length(CycleGLRGI[25,2]:CycleGLRGI[25,3])),rep("E",length(CycleGLRGI[26,2]:CycleGLRGI[26,3])),rep("F",length(CycleGLRGI[27,2]:CycleGLRGI[27,3])),rep("G",length(CycleGLRGI[28,2]:CycleGLRGI[28,3])))
CutG25LRGI=c(CycleGLRGI[29,2]:CycleGLRGI[29,3],CycleGLRGI[30,2]:CycleGLRGI[30,3],CycleGLRGI[31,2]:CycleGLRGI[31,3],CycleGLRGI[32,2]:CycleGLRGI[32,3],CycleGLRGI[33,2]:CycleGLRGI[33,3],CycleGLRGI[34,2]:CycleGLRGI[34,3],CycleGLRGI[35,2]:CycleGLRGI[35,3])+t4
CycleG25LRGI=c(rep("A",length(CycleGLRGI[29,2]:CycleGLRGI[29,3])),rep("B",length(CycleGLRGI[30,2]:CycleGLRGI[30,3])),rep("C",length(CycleGLRGI[31,2]:CycleGLRGI[31,3])),rep("D",length(CycleGLRGI[32,2]:CycleGLRGI[32,3])),rep("E",length(CycleGLRGI[33,2]:CycleGLRGI[33,3])),rep("F",length(CycleGLRGI[34,2]:CycleGLRGI[34,3])),rep("G",length(CycleGLRGI[35,2]:CycleGLRGI[35,3])))

CutGLRGI=TestLRGI[c(CutG5LRGI,CutG10LRGI,CutG15LRGI,CutG20LRGI,CutG25LRGI),]
CutGLRGI$Cycle=c(CycleG5LRGI,CycleG10LRGI,CycleG15LRGI,CycleG20LRGI,CycleG25LRGI)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLLRGI=CutLLRGI[!CutLLRGI$Cycle=="A",]; CutGLRGI=CutGLRGI[!CutGLRGI$Cycle=="A",]
MinLLRGI=setDT(CutLLRGI)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLLRGI=setDT(CutLLRGI)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGLRGI=setDT(CutGLRGI)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGLRGI=setDT(CutGLRGI)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLLRGI=rbind(MinLLRGI[,c(3,1,2,4)],MaxLLRGI[,c(3,1,2,4)])
BiomGLRGI=rbind(MinGLRGI[,c(3,1,2,5)],MaxGLRGI[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLLRGI=as.data.frame(BiomLLRGI %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLLRGI=SlopLLRGI[order(as.numeric(SlopLLRGI$Temp)),]; SlopLLRGI$Slope=unlist(SlopLLRGI$Slope)
BiomLLRGI=cbind(SlopLLRGI,Max=BiomLLRGI[BiomLLRGI$L>100000]$L)
MaxiLLRGI=setDT(BiomLLRGI)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLLRGI=setDT(BiomLLRGI)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGLRGI=as.data.frame(BiomGLRGI %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGLRGI=SlopGLRGI[order(as.numeric(SlopGLRGI$Temp)),]; SlopGLRGI$Slope=unlist(SlopGLRGI$Slope)
BiomGLRGI=cbind(SlopGLRGI,Max=BiomGLRGI[BiomGLRGI$G>10000]$G)
MaxiGLRGI=setDT(BiomGLRGI)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGLRGI=setDT(BiomGLRGI)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL33=data.frame(MeanLLRGI[,1],MeanLLRGI[,2],MeanLLRGI[,3],TimeLLRGI[,2],TimeLLRGI[,3])
colnames(DataL33)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL")
DataL33$Scenario=rep("LRGI",5)
DataL33$Subscenario=rep("GI",5)

DataG33=data.frame(MeanGLRGI[,1],MeanGLRGI[,2],MeanGLRGI[,3],TimeGLRGI[,2],TimeGLRGI[,3])
colnames(DataG33)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG")
DataG33$Scenario=rep("LRGI",5)
DataG33$Subscenario=rep("GI",5)




####################################################################################################




################################################
### SCENARIO 4: LEAF LITTER BIOMASS CONSTANT ###
################################################

#################################################
### SUBSCENARIO 1: GAMMARUS BIOMASS REDUCTION ###
#################################################

# Leaf fall = 300 gC/m2/an = 300 000 mgC/m2/an
# Gammarus density = 10 mgDM/m2 = 5 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5LCGR=GammLeafModel(Temp=5,GammMass=TSRA(5,4.26),Leaf=300000/Days,Gamm=5)
Test5LCGR=Test5LCGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test5LCGR=Test5LCGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5LCGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5°C", cex=1.5, side=3, line=-3, at=300)

Test10LCGR=GammLeafModel(Temp=10,GammMass=TSRA(10,4.26),Leaf=300000/Days,Gamm=5)
Test10LCGR=Test10LCGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test10LCGR=Test10LCGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10LCGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10°C", cex=1.5, side=3, line=-3, at=300)

Test15LCGR=GammLeafModel(Temp=15,GammMass=TSRA(15,4.26),Leaf=300000/Days,Gamm=5)
Test15LCGR=Test15LCGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test15LCGR=Test15LCGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15LCGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15°C", cex=1.5, side=3, line=-3, at=300)

Test20LCGR=GammLeafModel(Temp=20,GammMass=TSRA(20,4.26),Leaf=300000/Days,Gamm=5)
Test20LCGR=Test20LCGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test20LCGR=Test20LCGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20LCGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20°C", cex=1.5, side=3, line=-3, at=300)

Test25LCGR=GammLeafModel(Temp=25,GammMass=TSRA(25,4.26),Leaf=300000/Days,Gamm=5)
Test25LCGR=Test25LCGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test25LCGR=Test25LCGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25LCGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25°C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))


# Calculate population metabolism and leaf ingestion
Test5LCGR$M=(MetaQuadra(5,TSRM(5,4.26))/1000)*Test5LCGR[,3]
Test5LCGR$I=(AttackQuadra(5,TSRM(5,4.26))*Test5LCGR[,2]/(1+AttackQuadra(5,TSRM(5,4.26))*1/(IngQuadra(5,TSRM(5,4.26))/1000)*Test5LCGR[,2]))*0.30*Test5LCGR[,3]
Test10LCGR$M=(MetaQuadra(10,TSRM(5,4.26))/1000)*Test10LCGR[,3]
Test10LCGR$I=(AttackQuadra(10,TSRM(5,4.26))*Test10LCGR[,2]/(1+AttackQuadra(10,TSRM(5,4.26))*1/(IngQuadra(10,TSRM(5,4.26))/1000)*Test10LCGR[,2]))*0.30*Test10LCGR[,3]
Test15LCGR$M=(MetaQuadra(15,TSRM(5,4.26))/1000)*Test15LCGR[,3]
Test15LCGR$I=(AttackQuadra(15,TSRM(5,4.26))*Test15LCGR[,2]/(1+AttackQuadra(15,TSRM(5,4.26))*1/(IngQuadra(15,TSRM(5,4.26))/1000)*Test15LCGR[,2]))*0.30*Test15LCGR[,3]
Test20LCGR$M=(MetaQuadra(20,TSRM(5,4.26))/1000)*Test20LCGR[,3]
Test20LCGR$I=(AttackQuadra(20,TSRM(5,4.26))*Test20LCGR[,2]/(1+AttackQuadra(20,TSRM(5,4.26))*1/(IngQuadra(20,TSRM(5,4.26))/1000)*Test20LCGR[,2]))*0.30*Test20LCGR[,3]
Test25LCGR$M=(MetaQuadra(25,TSRM(5,4.26))/1000)*Test25LCGR[,3]
Test25LCGR$I=(AttackQuadra(25,TSRM(5,4.26))*Test25LCGR[,2]/(1+AttackQuadra(25,TSRM(5,4.26))*1/(IngQuadra(25,TSRM(5,4.26))/1000)*Test25LCGR[,2]))*0.30*Test25LCGR[,3]

# Combine the dataframes
TestLCGR=rbind(Test5LCGR,Test10LCGR,Test15LCGR,Test20LCGR,Test25LCGR)
TestLCGR$Temp=rep(c("5","10","15","20","25"), each=2556)
TestLCGR$Temp=factor(TestLCGR$Temp, levels=unique(TestLCGR$Temp))
colnames(TestLCGR)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestLCGR=subset(TestLCGR, !Time==2555); TestLCGR$Year=rep(rep(1:(nrow(TestLCGR)/(365*5)), each=365),5)
CutLCGR=as.data.frame(setDT(TestLCGR)[TestLCGR[, tail(.I, -16), by=Year]$V1]); CutLCGR=CutLCGR[!CutLCGR$Year=="1",]
MeanL=as.data.frame(setDT(CutLCGR)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutLCGR)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLLCGR=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGLCGR=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestLCGR=subset(TestLCGR, !Time==2555); TestLCGR$Year=rep(rep(1:(nrow(TestLCGR)/(365*5)), each=365),5)
CutLCGR=as.data.frame(setDT(TestLCGR)[TestLCGR[, tail(.I, -16), by=Year]$V1]); CutLCGR=CutLCGR[!CutLCGR$Year=="1",]
ThLLCGR=CutLCGR[CutLCGR$L>60000,]; ThGLCGR=CutLCGR[CutLCGR$G>5000,]
TimeLLCGR=as.data.frame(setDT(ThLLCGR)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGLCGR=as.data.frame(setDT(ThGLCGR)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLLCGR=as.data.frame(setDT(TimeLLCGR)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGLCGR=as.data.frame(setDT(TimeGLCGR)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLLCGR=as.data.frame(setDT(TestLCGR)[, .(MaxLLCGR=findPeaks(L), MinLLCGR=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGLCGR=as.data.frame(setDT(TestLCGR)[, .(MaxGLCGR=findPeaks(G), MinGLCGR=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5LCGR=c(CycleLLCGR[1,2]:CycleLLCGR[1,3],CycleLLCGR[2,2]:CycleLLCGR[2,3],CycleLLCGR[3,2]:CycleLLCGR[3,3],CycleLLCGR[4,2]:CycleLLCGR[4,3],CycleLLCGR[5,2]:CycleLLCGR[5,3],CycleLLCGR[6,2]:CycleLLCGR[6,3],CycleLLCGR[7,2]:CycleLLCGR[7,3])+t0
CycleL5LCGR=c(rep("A",length(CycleLLCGR[1,2]:CycleLLCGR[1,3])),rep("B",length(CycleLLCGR[2,2]:CycleLLCGR[2,3])),rep("C",length(CycleLLCGR[3,2]:CycleLLCGR[3,3])),rep("D",length(CycleLLCGR[4,2]:CycleLLCGR[4,3])),rep("E",length(CycleLLCGR[5,2]:CycleLLCGR[5,3])),rep("F",length(CycleLLCGR[6,2]:CycleLLCGR[6,3])),rep("G",length(CycleLLCGR[7,2]:CycleLLCGR[7,3])))
CutL10LCGR=c(CycleLLCGR[8,2]:CycleLLCGR[8,3],CycleLLCGR[9,2]:CycleLLCGR[9,3],CycleLLCGR[10,2]:CycleLLCGR[10,3],CycleLLCGR[11,2]:CycleLLCGR[11,3],CycleLLCGR[12,2]:CycleLLCGR[12,3],CycleLLCGR[13,2]:CycleLLCGR[13,3],CycleLLCGR[14,2]:CycleLLCGR[14,3])+t1
CycleL10LCGR=c(rep("A",length(CycleLLCGR[8,2]:CycleLLCGR[8,3])),rep("B",length(CycleLLCGR[9,2]:CycleLLCGR[9,3])),rep("C",length(CycleLLCGR[10,2]:CycleLLCGR[10,3])),rep("D",length(CycleLLCGR[11,2]:CycleLLCGR[11,3])),rep("E",length(CycleLLCGR[12,2]:CycleLLCGR[12,3])),rep("F",length(CycleLLCGR[13,2]:CycleLLCGR[13,3])),rep("G",length(CycleLLCGR[14,2]:CycleLLCGR[14,3])))
CutL15LCGR=c(CycleLLCGR[15,2]:CycleLLCGR[15,3],CycleLLCGR[16,2]:CycleLLCGR[16,3],CycleLLCGR[17,2]:CycleLLCGR[17,3],CycleLLCGR[18,2]:CycleLLCGR[18,3],CycleLLCGR[19,2]:CycleLLCGR[19,3],CycleLLCGR[20,2]:CycleLLCGR[20,3],CycleLLCGR[21,2]:CycleLLCGR[21,3])+t2
CycleL15LCGR=c(rep("A",length(CycleLLCGR[15,2]:CycleLLCGR[15,3])),rep("B",length(CycleLLCGR[16,2]:CycleLLCGR[16,3])),rep("C",length(CycleLLCGR[17,2]:CycleLLCGR[17,3])),rep("D",length(CycleLLCGR[18,2]:CycleLLCGR[18,3])),rep("E",length(CycleLLCGR[19,2]:CycleLLCGR[19,3])),rep("F",length(CycleLLCGR[20,2]:CycleLLCGR[20,3])),rep("G",length(CycleLLCGR[21,2]:CycleLLCGR[21,3])))
CutL20LCGR=c(CycleLLCGR[22,2]:CycleLLCGR[22,3],CycleLLCGR[23,2]:CycleLLCGR[23,3],CycleLLCGR[24,2]:CycleLLCGR[24,3],CycleLLCGR[25,2]:CycleLLCGR[25,3],CycleLLCGR[26,2]:CycleLLCGR[26,3],CycleLLCGR[27,2]:CycleLLCGR[27,3],CycleLLCGR[28,2]:CycleLLCGR[28,3])+t3
CycleL20LCGR=c(rep("A",length(CycleLLCGR[22,2]:CycleLLCGR[22,3])),rep("B",length(CycleLLCGR[23,2]:CycleLLCGR[23,3])),rep("C",length(CycleLLCGR[24,2]:CycleLLCGR[24,3])),rep("D",length(CycleLLCGR[25,2]:CycleLLCGR[25,3])),rep("E",length(CycleLLCGR[26,2]:CycleLLCGR[26,3])),rep("F",length(CycleLLCGR[27,2]:CycleLLCGR[27,3])),rep("G",length(CycleLLCGR[28,2]:CycleLLCGR[28,3])))
CutL25LCGR=c(CycleLLCGR[29,2]:CycleLLCGR[29,3],CycleLLCGR[30,2]:CycleLLCGR[30,3],CycleLLCGR[31,2]:CycleLLCGR[31,3],CycleLLCGR[32,2]:CycleLLCGR[32,3],CycleLLCGR[33,2]:CycleLLCGR[33,3],CycleLLCGR[34,2]:CycleLLCGR[34,3],CycleLLCGR[35,2]:CycleLLCGR[35,3])+t4
CycleL25LCGR=c(rep("A",length(CycleLLCGR[29,2]:CycleLLCGR[29,3])),rep("B",length(CycleLLCGR[30,2]:CycleLLCGR[30,3])),rep("C",length(CycleLLCGR[31,2]:CycleLLCGR[31,3])),rep("D",length(CycleLLCGR[32,2]:CycleLLCGR[32,3])),rep("E",length(CycleLLCGR[33,2]:CycleLLCGR[33,3])),rep("F",length(CycleLLCGR[34,2]:CycleLLCGR[34,3])),rep("G",length(CycleLLCGR[35,2]:CycleLLCGR[35,3])))

CutLLCGR=TestLCGR[c(CutL5LCGR,CutL10LCGR,CutL15LCGR,CutL20LCGR,CutL25LCGR),]
CutLLCGR$Cycle=c(CycleL5LCGR,CycleL10LCGR,CycleL15LCGR,CycleL20LCGR,CycleL25LCGR)

CutG5LCGR=c(CycleGLCGR[1,2]:CycleGLCGR[1,3],CycleGLCGR[2,2]:CycleGLCGR[2,3],CycleGLCGR[3,2]:CycleGLCGR[3,3],CycleGLCGR[4,2]:CycleGLCGR[4,3],CycleGLCGR[5,2]:CycleGLCGR[5,3],CycleGLCGR[6,2]:CycleGLCGR[6,3],CycleGLCGR[7,2]:CycleGLCGR[7,3])+t0
CycleG5LCGR=c(rep("A",length(CycleGLCGR[1,2]:CycleGLCGR[1,3])),rep("B",length(CycleGLCGR[2,2]:CycleGLCGR[2,3])),rep("C",length(CycleGLCGR[3,2]:CycleGLCGR[3,3])),rep("D",length(CycleGLCGR[4,2]:CycleGLCGR[4,3])),rep("E",length(CycleGLCGR[5,2]:CycleGLCGR[5,3])),rep("F",length(CycleGLCGR[6,2]:CycleGLCGR[6,3])),rep("G",length(CycleGLCGR[7,2]:CycleGLCGR[7,3])))
CutG10LCGR=c(CycleGLCGR[8,2]:CycleGLCGR[8,3],CycleGLCGR[9,2]:CycleGLCGR[9,3],CycleGLCGR[10,2]:CycleGLCGR[10,3],CycleGLCGR[11,2]:CycleGLCGR[11,3],CycleGLCGR[12,2]:CycleGLCGR[12,3],CycleGLCGR[13,2]:CycleGLCGR[13,3],CycleGLCGR[14,2]:CycleGLCGR[14,3])+t1
CycleG10LCGR=c(rep("A",length(CycleGLCGR[8,2]:CycleGLCGR[8,3])),rep("B",length(CycleGLCGR[9,2]:CycleGLCGR[9,3])),rep("C",length(CycleGLCGR[10,2]:CycleGLCGR[10,3])),rep("D",length(CycleGLCGR[11,2]:CycleGLCGR[11,3])),rep("E",length(CycleGLCGR[12,2]:CycleGLCGR[12,3])),rep("F",length(CycleGLCGR[13,2]:CycleGLCGR[13,3])),rep("G",length(CycleGLCGR[14,2]:CycleGLCGR[14,3])))
CutG15LCGR=c(CycleGLCGR[15,2]:CycleGLCGR[15,3],CycleGLCGR[16,2]:CycleGLCGR[16,3],CycleGLCGR[17,2]:CycleGLCGR[17,3],CycleGLCGR[18,2]:CycleGLCGR[18,3],CycleGLCGR[19,2]:CycleGLCGR[19,3],CycleGLCGR[20,2]:CycleGLCGR[20,3],CycleGLCGR[21,2]:CycleGLCGR[21,3])+t2
CycleG15LCGR=c(rep("A",length(CycleGLCGR[15,2]:CycleGLCGR[15,3])),rep("B",length(CycleGLCGR[16,2]:CycleGLCGR[16,3])),rep("C",length(CycleGLCGR[17,2]:CycleGLCGR[17,3])),rep("D",length(CycleGLCGR[18,2]:CycleGLCGR[18,3])),rep("E",length(CycleGLCGR[19,2]:CycleGLCGR[19,3])),rep("F",length(CycleGLCGR[20,2]:CycleGLCGR[20,3])),rep("G",length(CycleGLCGR[21,2]:CycleGLCGR[21,3])))
CutG20LCGR=c(CycleGLCGR[22,2]:CycleGLCGR[22,3],CycleGLCGR[23,2]:CycleGLCGR[23,3],CycleGLCGR[24,2]:CycleGLCGR[24,3],CycleGLCGR[25,2]:CycleGLCGR[25,3],CycleGLCGR[26,2]:CycleGLCGR[26,3],CycleGLCGR[27,2]:CycleGLCGR[27,3],CycleGLCGR[28,2]:CycleGLCGR[28,3])+t3
CycleG20LCGR=c(rep("A",length(CycleGLCGR[22,2]:CycleGLCGR[22,3])),rep("B",length(CycleGLCGR[23,2]:CycleGLCGR[23,3])),rep("C",length(CycleGLCGR[24,2]:CycleGLCGR[24,3])),rep("D",length(CycleGLCGR[25,2]:CycleGLCGR[25,3])),rep("E",length(CycleGLCGR[26,2]:CycleGLCGR[26,3])),rep("F",length(CycleGLCGR[27,2]:CycleGLCGR[27,3])),rep("G",length(CycleGLCGR[28,2]:CycleGLCGR[28,3])))
CutG25LCGR=c(CycleGLCGR[29,2]:CycleGLCGR[29,3],CycleGLCGR[30,2]:CycleGLCGR[30,3],CycleGLCGR[31,2]:CycleGLCGR[31,3],CycleGLCGR[32,2]:CycleGLCGR[32,3],CycleGLCGR[33,2]:CycleGLCGR[33,3],CycleGLCGR[34,2]:CycleGLCGR[34,3],CycleGLCGR[35,2]:CycleGLCGR[35,3])+t4
CycleG25LCGR=c(rep("A",length(CycleGLCGR[29,2]:CycleGLCGR[29,3])),rep("B",length(CycleGLCGR[30,2]:CycleGLCGR[30,3])),rep("C",length(CycleGLCGR[31,2]:CycleGLCGR[31,3])),rep("D",length(CycleGLCGR[32,2]:CycleGLCGR[32,3])),rep("E",length(CycleGLCGR[33,2]:CycleGLCGR[33,3])),rep("F",length(CycleGLCGR[34,2]:CycleGLCGR[34,3])),rep("G",length(CycleGLCGR[35,2]:CycleGLCGR[35,3])))

CutGLCGR=TestLCGR[c(CutG5LCGR,CutG10LCGR,CutG15LCGR,CutG20LCGR,CutG25LCGR),]
CutGLCGR$Cycle=c(CycleG5LCGR,CycleG10LCGR,CycleG15LCGR,CycleG20LCGR,CycleG25LCGR)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLLCGR=CutLLCGR[!CutLLCGR$Cycle=="A",]; CutGLCGR=CutGLCGR[!CutGLCGR$Cycle=="A",]
MinLLCGR=setDT(CutLLCGR)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLLCGR=setDT(CutLLCGR)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGLCGR=setDT(CutGLCGR)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGLCGR=setDT(CutGLCGR)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLLCGR=rbind(MinLLCGR[,c(3,1,2,4)],MaxLLCGR[,c(3,1,2,4)])
BiomGLCGR=rbind(MinGLCGR[,c(3,1,2,5)],MaxGLCGR[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLLCGR=as.data.frame(BiomLLCGR %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLLCGR=SlopLLCGR[order(as.numeric(SlopLLCGR$Temp)),]; SlopLLCGR$Slope=unlist(SlopLLCGR$Slope)
BiomLLCGR=cbind(SlopLLCGR,Max=BiomLLCGR[BiomLLCGR$L>100000]$L)
MaxiLLCGR=setDT(BiomLLCGR)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLLCGR=setDT(BiomLLCGR)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGLCGR=as.data.frame(BiomGLCGR %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGLCGR=SlopGLCGR[order(as.numeric(SlopGLCGR$Temp)),]; SlopGLCGR$Slope=unlist(SlopGLCGR$Slope)
BiomGLCGR=cbind(SlopGLCGR,Max=BiomGLCGR[BiomGLCGR$G>10000]$G)
MaxiGLCGR=setDT(BiomGLCGR)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGLCGR=setDT(BiomGLCGR)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL41=data.frame(MeanLLCGR[,1],MeanLLCGR[,2],MeanLLCGR[,3],TimeLLCGR[,2],TimeLLCGR[,3])
colnames(DataL41)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL")
DataL41$Scenario=rep("LCGR",5)
DataL41$Subscenario=rep("GR",5)

DataG41=data.frame(MeanGLCGR[,1],MeanGLCGR[,2],MeanGLCGR[,3],TimeGLCGR[,2],TimeGLCGR[,3])
colnames(DataG41)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG")
DataG41$Scenario=rep("LCGR",5)
DataG41$Subscenario=rep("GR",5)


################################################
### SUBSCENARIO 2: GAMMARUS BIOMASS CONSTANT ###
################################################

# Leaf fall = 300 gC/m2/an = 300 000 mgC/m2/an
# Gammarus density = 30 mgDM/m2 = 15 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5LCGC=GammLeafModel(Temp=5,GammMass=TSRA(5,4.26),Leaf=300000/Days,Gamm=15)
Test5LCGC=Test5LCGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test5LCGC=Test5LCGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5LCGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5°C", cex=1.5, side=3, line=-3, at=300)

Test10LCGC=GammLeafModel(Temp=10,GammMass=TSRA(10,4.26),Leaf=300000/Days,Gamm=15)
Test10LCGC=Test10LCGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test10LCGC=Test10LCGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10LCGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10°C", cex=1.5, side=3, line=-3, at=300)

Test15LCGC=GammLeafModel(Temp=15,GammMass=TSRA(15,4.26),Leaf=300000/Days,Gamm=15)
Test15LCGC=Test15LCGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test15LCGC=Test15LCGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15LCGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15°C", cex=1.5, side=3, line=-3, at=300)

Test20LCGC=GammLeafModel(Temp=20,GammMass=TSRA(20,4.26),Leaf=300000/Days,Gamm=15)
Test20LCGC=Test20LCGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test20LCGC=Test20LCGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20LCGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20°C", cex=1.5, side=3, line=-3, at=300)

Test25LCGC=GammLeafModel(Temp=25,GammMass=TSRA(25,4.26),Leaf=300000/Days,Gamm=15)
Test25LCGC=Test25LCGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test25LCGC=Test25LCGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25LCGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25°C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))


# Calculate population metabolism and leaf ingestion
Test5LCGC$M=(MetaQuadra(5,TSRM(5,4.26))/1000)*Test5LCGC[,3]
Test5LCGC$I=(AttackQuadra(5,TSRM(5,4.26))*Test5LCGC[,2]/(1+AttackQuadra(5,TSRM(5,4.26))*1/(IngQuadra(5,TSRM(5,4.26))/1000)*Test5LCGC[,2]))*0.30*Test5LCGC[,3]
Test10LCGC$M=(MetaQuadra(10,TSRM(5,4.26))/1000)*Test10LCGC[,3]
Test10LCGC$I=(AttackQuadra(10,TSRM(5,4.26))*Test10LCGC[,2]/(1+AttackQuadra(10,TSRM(5,4.26))*1/(IngQuadra(10,TSRM(5,4.26))/1000)*Test10LCGC[,2]))*0.30*Test10LCGC[,3]
Test15LCGC$M=(MetaQuadra(15,TSRM(5,4.26))/1000)*Test15LCGC[,3]
Test15LCGC$I=(AttackQuadra(15,TSRM(5,4.26))*Test15LCGC[,2]/(1+AttackQuadra(15,TSRM(5,4.26))*1/(IngQuadra(15,TSRM(5,4.26))/1000)*Test15LCGC[,2]))*0.30*Test15LCGC[,3]
Test20LCGC$M=(MetaQuadra(20,TSRM(5,4.26))/1000)*Test20LCGC[,3]
Test20LCGC$I=(AttackQuadra(20,TSRM(5,4.26))*Test20LCGC[,2]/(1+AttackQuadra(20,TSRM(5,4.26))*1/(IngQuadra(20,TSRM(5,4.26))/1000)*Test20LCGC[,2]))*0.30*Test20LCGC[,3]
Test25LCGC$M=(MetaQuadra(25,TSRM(5,4.26))/1000)*Test25LCGC[,3]
Test25LCGC$I=(AttackQuadra(25,TSRM(5,4.26))*Test25LCGC[,2]/(1+AttackQuadra(25,TSRM(5,4.26))*1/(IngQuadra(25,TSRM(5,4.26))/1000)*Test25LCGC[,2]))*0.30*Test25LCGC[,3]

# Combine the dataframes
TestLCGC=rbind(Test5LCGC,Test10LCGC,Test15LCGC,Test20LCGC,Test25LCGC)
TestLCGC$Temp=rep(c("5","10","15","20","25"), each=2556)
TestLCGC$Temp=factor(TestLCGC$Temp, levels=unique(TestLCGC$Temp))
colnames(TestLCGC)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestLCGC=subset(TestLCGC, !Time==2555); TestLCGC$Year=rep(rep(1:(nrow(TestLCGC)/(365*5)), each=365),5)
CutLCGC=as.data.frame(setDT(TestLCGC)[TestLCGC[, tail(.I, -16), by=Year]$V1]); CutLCGC=CutLCGC[!CutLCGC$Year=="1",]
MeanL=as.data.frame(setDT(CutLCGC)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutLCGC)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLLCGC=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGLCGC=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestLCGC=subset(TestLCGC, !Time==2555); TestLCGC$Year=rep(rep(1:(nrow(TestLCGC)/(365*5)), each=365),5)
CutLCGC=as.data.frame(setDT(TestLCGC)[TestLCGC[, tail(.I, -16), by=Year]$V1]); CutLCGC=CutLCGC[!CutLCGC$Year=="1",]
ThLLCGC=CutLCGC[CutLCGC$L>60000,]; ThGLCGC=CutLCGC[CutLCGC$G>5000,]
TimeLLCGC=as.data.frame(setDT(ThLLCGC)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGLCGC=as.data.frame(setDT(ThGLCGC)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLLCGC=as.data.frame(setDT(TimeLLCGC)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGLCGC=as.data.frame(setDT(TimeGLCGC)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLLCGC=as.data.frame(setDT(TestLCGC)[, .(MaxLLCGC=findPeaks(L), MinLLCGC=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGLCGC=as.data.frame(setDT(TestLCGC)[, .(MaxGLCGC=findPeaks(G), MinGLCGC=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5LCGC=c(CycleLLCGC[1,2]:CycleLLCGC[1,3],CycleLLCGC[2,2]:CycleLLCGC[2,3],CycleLLCGC[3,2]:CycleLLCGC[3,3],CycleLLCGC[4,2]:CycleLLCGC[4,3],CycleLLCGC[5,2]:CycleLLCGC[5,3],CycleLLCGC[6,2]:CycleLLCGC[6,3],CycleLLCGC[7,2]:CycleLLCGC[7,3])+t0
CycleL5LCGC=c(rep("A",length(CycleLLCGC[1,2]:CycleLLCGC[1,3])),rep("B",length(CycleLLCGC[2,2]:CycleLLCGC[2,3])),rep("C",length(CycleLLCGC[3,2]:CycleLLCGC[3,3])),rep("D",length(CycleLLCGC[4,2]:CycleLLCGC[4,3])),rep("E",length(CycleLLCGC[5,2]:CycleLLCGC[5,3])),rep("F",length(CycleLLCGC[6,2]:CycleLLCGC[6,3])),rep("G",length(CycleLLCGC[7,2]:CycleLLCGC[7,3])))
CutL10LCGC=c(CycleLLCGC[8,2]:CycleLLCGC[8,3],CycleLLCGC[9,2]:CycleLLCGC[9,3],CycleLLCGC[10,2]:CycleLLCGC[10,3],CycleLLCGC[11,2]:CycleLLCGC[11,3],CycleLLCGC[12,2]:CycleLLCGC[12,3],CycleLLCGC[13,2]:CycleLLCGC[13,3],CycleLLCGC[14,2]:CycleLLCGC[14,3])+t1
CycleL10LCGC=c(rep("A",length(CycleLLCGC[8,2]:CycleLLCGC[8,3])),rep("B",length(CycleLLCGC[9,2]:CycleLLCGC[9,3])),rep("C",length(CycleLLCGC[10,2]:CycleLLCGC[10,3])),rep("D",length(CycleLLCGC[11,2]:CycleLLCGC[11,3])),rep("E",length(CycleLLCGC[12,2]:CycleLLCGC[12,3])),rep("F",length(CycleLLCGC[13,2]:CycleLLCGC[13,3])),rep("G",length(CycleLLCGC[14,2]:CycleLLCGC[14,3])))
CutL15LCGC=c(CycleLLCGC[15,2]:CycleLLCGC[15,3],CycleLLCGC[16,2]:CycleLLCGC[16,3],CycleLLCGC[17,2]:CycleLLCGC[17,3],CycleLLCGC[18,2]:CycleLLCGC[18,3],CycleLLCGC[19,2]:CycleLLCGC[19,3],CycleLLCGC[20,2]:CycleLLCGC[20,3],CycleLLCGC[21,2]:CycleLLCGC[21,3])+t2
CycleL15LCGC=c(rep("A",length(CycleLLCGC[15,2]:CycleLLCGC[15,3])),rep("B",length(CycleLLCGC[16,2]:CycleLLCGC[16,3])),rep("C",length(CycleLLCGC[17,2]:CycleLLCGC[17,3])),rep("D",length(CycleLLCGC[18,2]:CycleLLCGC[18,3])),rep("E",length(CycleLLCGC[19,2]:CycleLLCGC[19,3])),rep("F",length(CycleLLCGC[20,2]:CycleLLCGC[20,3])),rep("G",length(CycleLLCGC[21,2]:CycleLLCGC[21,3])))
CutL20LCGC=c(CycleLLCGC[22,2]:CycleLLCGC[22,3],CycleLLCGC[23,2]:CycleLLCGC[23,3],CycleLLCGC[24,2]:CycleLLCGC[24,3],CycleLLCGC[25,2]:CycleLLCGC[25,3],CycleLLCGC[26,2]:CycleLLCGC[26,3],CycleLLCGC[27,2]:CycleLLCGC[27,3],CycleLLCGC[28,2]:CycleLLCGC[28,3])+t3
CycleL20LCGC=c(rep("A",length(CycleLLCGC[22,2]:CycleLLCGC[22,3])),rep("B",length(CycleLLCGC[23,2]:CycleLLCGC[23,3])),rep("C",length(CycleLLCGC[24,2]:CycleLLCGC[24,3])),rep("D",length(CycleLLCGC[25,2]:CycleLLCGC[25,3])),rep("E",length(CycleLLCGC[26,2]:CycleLLCGC[26,3])),rep("F",length(CycleLLCGC[27,2]:CycleLLCGC[27,3])),rep("G",length(CycleLLCGC[28,2]:CycleLLCGC[28,3])))
CutL25LCGC=c(CycleLLCGC[29,2]:CycleLLCGC[29,3],CycleLLCGC[30,2]:CycleLLCGC[30,3],CycleLLCGC[31,2]:CycleLLCGC[31,3],CycleLLCGC[32,2]:CycleLLCGC[32,3],CycleLLCGC[33,2]:CycleLLCGC[33,3],CycleLLCGC[34,2]:CycleLLCGC[34,3],CycleLLCGC[35,2]:CycleLLCGC[35,3])+t4
CycleL25LCGC=c(rep("A",length(CycleLLCGC[29,2]:CycleLLCGC[29,3])),rep("B",length(CycleLLCGC[30,2]:CycleLLCGC[30,3])),rep("C",length(CycleLLCGC[31,2]:CycleLLCGC[31,3])),rep("D",length(CycleLLCGC[32,2]:CycleLLCGC[32,3])),rep("E",length(CycleLLCGC[33,2]:CycleLLCGC[33,3])),rep("F",length(CycleLLCGC[34,2]:CycleLLCGC[34,3])),rep("G",length(CycleLLCGC[35,2]:CycleLLCGC[35,3])))

CutLLCGC=TestLCGC[c(CutL5LCGC,CutL10LCGC,CutL15LCGC,CutL20LCGC,CutL25LCGC),]
CutLLCGC$Cycle=c(CycleL5LCGC,CycleL10LCGC,CycleL15LCGC,CycleL20LCGC,CycleL25LCGC)

CutG5LCGC=c(CycleGLCGC[1,2]:CycleGLCGC[1,3],CycleGLCGC[2,2]:CycleGLCGC[2,3],CycleGLCGC[3,2]:CycleGLCGC[3,3],CycleGLCGC[4,2]:CycleGLCGC[4,3],CycleGLCGC[5,2]:CycleGLCGC[5,3],CycleGLCGC[6,2]:CycleGLCGC[6,3],CycleGLCGC[7,2]:CycleGLCGC[7,3])+t0
CycleG5LCGC=c(rep("A",length(CycleGLCGC[1,2]:CycleGLCGC[1,3])),rep("B",length(CycleGLCGC[2,2]:CycleGLCGC[2,3])),rep("C",length(CycleGLCGC[3,2]:CycleGLCGC[3,3])),rep("D",length(CycleGLCGC[4,2]:CycleGLCGC[4,3])),rep("E",length(CycleGLCGC[5,2]:CycleGLCGC[5,3])),rep("F",length(CycleGLCGC[6,2]:CycleGLCGC[6,3])),rep("G",length(CycleGLCGC[7,2]:CycleGLCGC[7,3])))
CutG10LCGC=c(CycleGLCGC[8,2]:CycleGLCGC[8,3],CycleGLCGC[9,2]:CycleGLCGC[9,3],CycleGLCGC[10,2]:CycleGLCGC[10,3],CycleGLCGC[11,2]:CycleGLCGC[11,3],CycleGLCGC[12,2]:CycleGLCGC[12,3],CycleGLCGC[13,2]:CycleGLCGC[13,3],CycleGLCGC[14,2]:CycleGLCGC[14,3])+t1
CycleG10LCGC=c(rep("A",length(CycleGLCGC[8,2]:CycleGLCGC[8,3])),rep("B",length(CycleGLCGC[9,2]:CycleGLCGC[9,3])),rep("C",length(CycleGLCGC[10,2]:CycleGLCGC[10,3])),rep("D",length(CycleGLCGC[11,2]:CycleGLCGC[11,3])),rep("E",length(CycleGLCGC[12,2]:CycleGLCGC[12,3])),rep("F",length(CycleGLCGC[13,2]:CycleGLCGC[13,3])),rep("G",length(CycleGLCGC[14,2]:CycleGLCGC[14,3])))
CutG15LCGC=c(CycleGLCGC[15,2]:CycleGLCGC[15,3],CycleGLCGC[16,2]:CycleGLCGC[16,3],CycleGLCGC[17,2]:CycleGLCGC[17,3],CycleGLCGC[18,2]:CycleGLCGC[18,3],CycleGLCGC[19,2]:CycleGLCGC[19,3],CycleGLCGC[20,2]:CycleGLCGC[20,3],CycleGLCGC[21,2]:CycleGLCGC[21,3])+t2
CycleG15LCGC=c(rep("A",length(CycleGLCGC[15,2]:CycleGLCGC[15,3])),rep("B",length(CycleGLCGC[16,2]:CycleGLCGC[16,3])),rep("C",length(CycleGLCGC[17,2]:CycleGLCGC[17,3])),rep("D",length(CycleGLCGC[18,2]:CycleGLCGC[18,3])),rep("E",length(CycleGLCGC[19,2]:CycleGLCGC[19,3])),rep("F",length(CycleGLCGC[20,2]:CycleGLCGC[20,3])),rep("G",length(CycleGLCGC[21,2]:CycleGLCGC[21,3])))
CutG20LCGC=c(CycleGLCGC[22,2]:CycleGLCGC[22,3],CycleGLCGC[23,2]:CycleGLCGC[23,3],CycleGLCGC[24,2]:CycleGLCGC[24,3],CycleGLCGC[25,2]:CycleGLCGC[25,3],CycleGLCGC[26,2]:CycleGLCGC[26,3],CycleGLCGC[27,2]:CycleGLCGC[27,3],CycleGLCGC[28,2]:CycleGLCGC[28,3])+t3
CycleG20LCGC=c(rep("A",length(CycleGLCGC[22,2]:CycleGLCGC[22,3])),rep("B",length(CycleGLCGC[23,2]:CycleGLCGC[23,3])),rep("C",length(CycleGLCGC[24,2]:CycleGLCGC[24,3])),rep("D",length(CycleGLCGC[25,2]:CycleGLCGC[25,3])),rep("E",length(CycleGLCGC[26,2]:CycleGLCGC[26,3])),rep("F",length(CycleGLCGC[27,2]:CycleGLCGC[27,3])),rep("G",length(CycleGLCGC[28,2]:CycleGLCGC[28,3])))
CutG25LCGC=c(CycleGLCGC[29,2]:CycleGLCGC[29,3],CycleGLCGC[30,2]:CycleGLCGC[30,3],CycleGLCGC[31,2]:CycleGLCGC[31,3],CycleGLCGC[32,2]:CycleGLCGC[32,3],CycleGLCGC[33,2]:CycleGLCGC[33,3],CycleGLCGC[34,2]:CycleGLCGC[34,3],CycleGLCGC[35,2]:CycleGLCGC[35,3])+t4
CycleG25LCGC=c(rep("A",length(CycleGLCGC[29,2]:CycleGLCGC[29,3])),rep("B",length(CycleGLCGC[30,2]:CycleGLCGC[30,3])),rep("C",length(CycleGLCGC[31,2]:CycleGLCGC[31,3])),rep("D",length(CycleGLCGC[32,2]:CycleGLCGC[32,3])),rep("E",length(CycleGLCGC[33,2]:CycleGLCGC[33,3])),rep("F",length(CycleGLCGC[34,2]:CycleGLCGC[34,3])),rep("G",length(CycleGLCGC[35,2]:CycleGLCGC[35,3])))

CutGLCGC=TestLCGC[c(CutG5LCGC,CutG10LCGC,CutG15LCGC,CutG20LCGC,CutG25LCGC),]
CutGLCGC$Cycle=c(CycleG5LCGC,CycleG10LCGC,CycleG15LCGC,CycleG20LCGC,CycleG25LCGC)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLLCGC=CutLLCGC[!CutLLCGC$Cycle=="A",]; CutGLCGC=CutGLCGC[!CutGLCGC$Cycle=="A",]
MinLLCGC=setDT(CutLLCGC)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLLCGC=setDT(CutLLCGC)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGLCGC=setDT(CutGLCGC)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGLCGC=setDT(CutGLCGC)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLLCGC=rbind(MinLLCGC[,c(3,1,2,4)],MaxLLCGC[,c(3,1,2,4)])
BiomGLCGC=rbind(MinGLCGC[,c(3,1,2,5)],MaxGLCGC[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLLCGC=as.data.frame(BiomLLCGC %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLLCGC=SlopLLCGC[order(as.numeric(SlopLLCGC$Temp)),]; SlopLLCGC$Slope=unlist(SlopLLCGC$Slope)
BiomLLCGC=cbind(SlopLLCGC,Max=BiomLLCGC[BiomLLCGC$L>100000]$L)
MaxiLLCGC=setDT(BiomLLCGC)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLLCGC=setDT(BiomLLCGC)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGLCGC=as.data.frame(BiomGLCGC %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGLCGC=SlopGLCGC[order(as.numeric(SlopGLCGC$Temp)),]; SlopGLCGC$Slope=unlist(SlopGLCGC$Slope)
BiomGLCGC=cbind(SlopGLCGC,Max=BiomGLCGC[BiomGLCGC$G>10000]$G)
MaxiGLCGC=setDT(BiomGLCGC)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGLCGC=setDT(BiomGLCGC)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL42=data.frame(MeanLLCGC[,1],MeanLLCGC[,2],MeanLLCGC[,3],TimeLLCGC[,2],TimeLLCGC[,3])
colnames(DataL42)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL")
DataL42$Scenario=rep("LCGC",5)
DataL42$Subscenario=rep("GC",5)

DataG42=data.frame(MeanGLCGC[,1],MeanGLCGC[,2],MeanGLCGC[,3],TimeGLCGC[,2],TimeGLCGC[,3])
colnames(DataG42)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG")
DataG42$Scenario=rep("LCGC",5)
DataG42$Subscenario=rep("GC",5)


################################################
### SUBSCENARIO 3: GAMMARUS BIOMASS INCREASE ###
################################################

# Leaf fall = 300 gC/m2/an = 300 000 mgC/m2/an
# Gammarus density = 90 mgDM/m2 = 45 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5LCGI=GammLeafModel(Temp=5,GammMass=TSRA(5,4.26),Leaf=300000/Days,Gamm=45)
Test5LCGI=Test5LCGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test5LCGI=Test5LCGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5LCGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5°C", cex=1.5, side=3, line=-3, at=300)

Test10LCGI=GammLeafModel(Temp=10,GammMass=TSRA(10,4.26),Leaf=300000/Days,Gamm=45)
Test10LCGI=Test10LCGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test10LCGI=Test10LCGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10LCGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10°C", cex=1.5, side=3, line=-3, at=300)

Test15LCGI=GammLeafModel(Temp=15,GammMass=TSRA(15,4.26),Leaf=300000/Days,Gamm=45)
Test15LCGI=Test15LCGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test15LCGI=Test15LCGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15LCGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15°C", cex=1.5, side=3, line=-3, at=300)

Test20LCGI=GammLeafModel(Temp=20,GammMass=TSRA(20,4.26),Leaf=300000/Days,Gamm=45)
Test20LCGI=Test20LCGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test20LCGI=Test20LCGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20LCGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20°C", cex=1.5, side=3, line=-3, at=300)

Test25LCGI=GammLeafModel(Temp=25,GammMass=TSRA(25,4.26),Leaf=300000/Days,Gamm=45)
Test25LCGI=Test25LCGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test25LCGI=Test25LCGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25LCGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25°C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))


# Calculate population metabolism and leaf ingestion
Test5LCGI$M=(MetaQuadra(5,TSRM(5,4.26))/1000)*Test5LCGI[,3]
Test5LCGI$I=(AttackQuadra(5,TSRM(5,4.26))*Test5LCGI[,2]/(1+AttackQuadra(5,TSRM(5,4.26))*1/(IngQuadra(5,TSRM(5,4.26))/1000)*Test5LCGI[,2]))*0.30*Test5LCGI[,3]
Test10LCGI$M=(MetaQuadra(10,TSRM(5,4.26))/1000)*Test10LCGI[,3]
Test10LCGI$I=(AttackQuadra(10,TSRM(5,4.26))*Test10LCGI[,2]/(1+AttackQuadra(10,TSRM(5,4.26))*1/(IngQuadra(10,TSRM(5,4.26))/1000)*Test10LCGI[,2]))*0.30*Test10LCGI[,3]
Test15LCGI$M=(MetaQuadra(15,TSRM(5,4.26))/1000)*Test15LCGI[,3]
Test15LCGI$I=(AttackQuadra(15,TSRM(5,4.26))*Test15LCGI[,2]/(1+AttackQuadra(15,TSRM(5,4.26))*1/(IngQuadra(15,TSRM(5,4.26))/1000)*Test15LCGI[,2]))*0.30*Test15LCGI[,3]
Test20LCGI$M=(MetaQuadra(20,TSRM(5,4.26))/1000)*Test20LCGI[,3]
Test20LCGI$I=(AttackQuadra(20,TSRM(5,4.26))*Test20LCGI[,2]/(1+AttackQuadra(20,TSRM(5,4.26))*1/(IngQuadra(20,TSRM(5,4.26))/1000)*Test20LCGI[,2]))*0.30*Test20LCGI[,3]
Test25LCGI$M=(MetaQuadra(25,TSRM(5,4.26))/1000)*Test25LCGI[,3]
Test25LCGI$I=(AttackQuadra(25,TSRM(5,4.26))*Test25LCGI[,2]/(1+AttackQuadra(25,TSRM(5,4.26))*1/(IngQuadra(25,TSRM(5,4.26))/1000)*Test25LCGI[,2]))*0.30*Test25LCGI[,3]

# Combine the dataframes
TestLCGI=rbind(Test5LCGI,Test10LCGI,Test15LCGI,Test20LCGI,Test25LCGI)
TestLCGI$Temp=rep(c("5","10","15","20","25"), each=2556)
TestLCGI$Temp=factor(TestLCGI$Temp, levels=unique(TestLCGI$Temp))
colnames(TestLCGI)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestLCGI=subset(TestLCGI, !Time==2555); TestLCGI$Year=rep(rep(1:(nrow(TestLCGI)/(365*5)), each=365),5)
CutLCGI=as.data.frame(setDT(TestLCGI)[TestLCGI[, tail(.I, -16), by=Year]$V1]); CutLCGI=CutLCGI[!CutLCGI$Year=="1",]
MeanL=as.data.frame(setDT(CutLCGI)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutLCGI)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLLCGI=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGLCGI=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestLCGI=subset(TestLCGI, !Time==2555); TestLCGI$Year=rep(rep(1:(nrow(TestLCGI)/(365*5)), each=365),5)
CutLCGI=as.data.frame(setDT(TestLCGI)[TestLCGI[, tail(.I, -16), by=Year]$V1]); CutLCGI=CutLCGI[!CutLCGI$Year=="1",]
ThLLCGI=CutLCGI[CutLCGI$L>60000,]; ThGLCGI=CutLCGI[CutLCGI$G>5000,]
TimeLLCGI=as.data.frame(setDT(ThLLCGI)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGLCGI=as.data.frame(setDT(ThGLCGI)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLLCGI=as.data.frame(setDT(TimeLLCGI)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGLCGI=as.data.frame(setDT(TimeGLCGI)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLLCGI=as.data.frame(setDT(TestLCGI)[, .(MaxLLCGI=findPeaks(L), MinLLCGI=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGLCGI=as.data.frame(setDT(TestLCGI)[, .(MaxGLCGI=findPeaks(G), MinGLCGI=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5LCGI=c(CycleLLCGI[1,2]:CycleLLCGI[1,3],CycleLLCGI[2,2]:CycleLLCGI[2,3],CycleLLCGI[3,2]:CycleLLCGI[3,3],CycleLLCGI[4,2]:CycleLLCGI[4,3],CycleLLCGI[5,2]:CycleLLCGI[5,3],CycleLLCGI[6,2]:CycleLLCGI[6,3],CycleLLCGI[7,2]:CycleLLCGI[7,3])+t0
CycleL5LCGI=c(rep("A",length(CycleLLCGI[1,2]:CycleLLCGI[1,3])),rep("B",length(CycleLLCGI[2,2]:CycleLLCGI[2,3])),rep("C",length(CycleLLCGI[3,2]:CycleLLCGI[3,3])),rep("D",length(CycleLLCGI[4,2]:CycleLLCGI[4,3])),rep("E",length(CycleLLCGI[5,2]:CycleLLCGI[5,3])),rep("F",length(CycleLLCGI[6,2]:CycleLLCGI[6,3])),rep("G",length(CycleLLCGI[7,2]:CycleLLCGI[7,3])))
CutL10LCGI=c(CycleLLCGI[8,2]:CycleLLCGI[8,3],CycleLLCGI[9,2]:CycleLLCGI[9,3],CycleLLCGI[10,2]:CycleLLCGI[10,3],CycleLLCGI[11,2]:CycleLLCGI[11,3],CycleLLCGI[12,2]:CycleLLCGI[12,3],CycleLLCGI[13,2]:CycleLLCGI[13,3],CycleLLCGI[14,2]:CycleLLCGI[14,3])+t1
CycleL10LCGI=c(rep("A",length(CycleLLCGI[8,2]:CycleLLCGI[8,3])),rep("B",length(CycleLLCGI[9,2]:CycleLLCGI[9,3])),rep("C",length(CycleLLCGI[10,2]:CycleLLCGI[10,3])),rep("D",length(CycleLLCGI[11,2]:CycleLLCGI[11,3])),rep("E",length(CycleLLCGI[12,2]:CycleLLCGI[12,3])),rep("F",length(CycleLLCGI[13,2]:CycleLLCGI[13,3])),rep("G",length(CycleLLCGI[14,2]:CycleLLCGI[14,3])))
CutL15LCGI=c(CycleLLCGI[15,2]:CycleLLCGI[15,3],CycleLLCGI[16,2]:CycleLLCGI[16,3],CycleLLCGI[17,2]:CycleLLCGI[17,3],CycleLLCGI[18,2]:CycleLLCGI[18,3],CycleLLCGI[19,2]:CycleLLCGI[19,3],CycleLLCGI[20,2]:CycleLLCGI[20,3],CycleLLCGI[21,2]:CycleLLCGI[21,3])+t2
CycleL15LCGI=c(rep("A",length(CycleLLCGI[15,2]:CycleLLCGI[15,3])),rep("B",length(CycleLLCGI[16,2]:CycleLLCGI[16,3])),rep("C",length(CycleLLCGI[17,2]:CycleLLCGI[17,3])),rep("D",length(CycleLLCGI[18,2]:CycleLLCGI[18,3])),rep("E",length(CycleLLCGI[19,2]:CycleLLCGI[19,3])),rep("F",length(CycleLLCGI[20,2]:CycleLLCGI[20,3])),rep("G",length(CycleLLCGI[21,2]:CycleLLCGI[21,3])))
CutL20LCGI=c(CycleLLCGI[22,2]:CycleLLCGI[22,3],CycleLLCGI[23,2]:CycleLLCGI[23,3],CycleLLCGI[24,2]:CycleLLCGI[24,3],CycleLLCGI[25,2]:CycleLLCGI[25,3],CycleLLCGI[26,2]:CycleLLCGI[26,3],CycleLLCGI[27,2]:CycleLLCGI[27,3],CycleLLCGI[28,2]:CycleLLCGI[28,3])+t3
CycleL20LCGI=c(rep("A",length(CycleLLCGI[22,2]:CycleLLCGI[22,3])),rep("B",length(CycleLLCGI[23,2]:CycleLLCGI[23,3])),rep("C",length(CycleLLCGI[24,2]:CycleLLCGI[24,3])),rep("D",length(CycleLLCGI[25,2]:CycleLLCGI[25,3])),rep("E",length(CycleLLCGI[26,2]:CycleLLCGI[26,3])),rep("F",length(CycleLLCGI[27,2]:CycleLLCGI[27,3])),rep("G",length(CycleLLCGI[28,2]:CycleLLCGI[28,3])))
CutL25LCGI=c(CycleLLCGI[29,2]:CycleLLCGI[29,3],CycleLLCGI[30,2]:CycleLLCGI[30,3],CycleLLCGI[31,2]:CycleLLCGI[31,3],CycleLLCGI[32,2]:CycleLLCGI[32,3],CycleLLCGI[33,2]:CycleLLCGI[33,3],CycleLLCGI[34,2]:CycleLLCGI[34,3],CycleLLCGI[35,2]:CycleLLCGI[35,3])+t4
CycleL25LCGI=c(rep("A",length(CycleLLCGI[29,2]:CycleLLCGI[29,3])),rep("B",length(CycleLLCGI[30,2]:CycleLLCGI[30,3])),rep("C",length(CycleLLCGI[31,2]:CycleLLCGI[31,3])),rep("D",length(CycleLLCGI[32,2]:CycleLLCGI[32,3])),rep("E",length(CycleLLCGI[33,2]:CycleLLCGI[33,3])),rep("F",length(CycleLLCGI[34,2]:CycleLLCGI[34,3])),rep("G",length(CycleLLCGI[35,2]:CycleLLCGI[35,3])))

CutLLCGI=TestLCGI[c(CutL5LCGI,CutL10LCGI,CutL15LCGI,CutL20LCGI,CutL25LCGI),]
CutLLCGI$Cycle=c(CycleL5LCGI,CycleL10LCGI,CycleL15LCGI,CycleL20LCGI,CycleL25LCGI)

CutG5LCGI=c(CycleGLCGI[1,2]:CycleGLCGI[1,3],CycleGLCGI[2,2]:CycleGLCGI[2,3],CycleGLCGI[3,2]:CycleGLCGI[3,3],CycleGLCGI[4,2]:CycleGLCGI[4,3],CycleGLCGI[5,2]:CycleGLCGI[5,3],CycleGLCGI[6,2]:CycleGLCGI[6,3],CycleGLCGI[7,2]:CycleGLCGI[7,3])+t0
CycleG5LCGI=c(rep("A",length(CycleGLCGI[1,2]:CycleGLCGI[1,3])),rep("B",length(CycleGLCGI[2,2]:CycleGLCGI[2,3])),rep("C",length(CycleGLCGI[3,2]:CycleGLCGI[3,3])),rep("D",length(CycleGLCGI[4,2]:CycleGLCGI[4,3])),rep("E",length(CycleGLCGI[5,2]:CycleGLCGI[5,3])),rep("F",length(CycleGLCGI[6,2]:CycleGLCGI[6,3])),rep("G",length(CycleGLCGI[7,2]:CycleGLCGI[7,3])))
CutG10LCGI=c(CycleGLCGI[8,2]:CycleGLCGI[8,3],CycleGLCGI[9,2]:CycleGLCGI[9,3],CycleGLCGI[10,2]:CycleGLCGI[10,3],CycleGLCGI[11,2]:CycleGLCGI[11,3],CycleGLCGI[12,2]:CycleGLCGI[12,3],CycleGLCGI[13,2]:CycleGLCGI[13,3],CycleGLCGI[14,2]:CycleGLCGI[14,3])+t1
CycleG10LCGI=c(rep("A",length(CycleGLCGI[8,2]:CycleGLCGI[8,3])),rep("B",length(CycleGLCGI[9,2]:CycleGLCGI[9,3])),rep("C",length(CycleGLCGI[10,2]:CycleGLCGI[10,3])),rep("D",length(CycleGLCGI[11,2]:CycleGLCGI[11,3])),rep("E",length(CycleGLCGI[12,2]:CycleGLCGI[12,3])),rep("F",length(CycleGLCGI[13,2]:CycleGLCGI[13,3])),rep("G",length(CycleGLCGI[14,2]:CycleGLCGI[14,3])))
CutG15LCGI=c(CycleGLCGI[15,2]:CycleGLCGI[15,3],CycleGLCGI[16,2]:CycleGLCGI[16,3],CycleGLCGI[17,2]:CycleGLCGI[17,3],CycleGLCGI[18,2]:CycleGLCGI[18,3],CycleGLCGI[19,2]:CycleGLCGI[19,3],CycleGLCGI[20,2]:CycleGLCGI[20,3],CycleGLCGI[21,2]:CycleGLCGI[21,3])+t2
CycleG15LCGI=c(rep("A",length(CycleGLCGI[15,2]:CycleGLCGI[15,3])),rep("B",length(CycleGLCGI[16,2]:CycleGLCGI[16,3])),rep("C",length(CycleGLCGI[17,2]:CycleGLCGI[17,3])),rep("D",length(CycleGLCGI[18,2]:CycleGLCGI[18,3])),rep("E",length(CycleGLCGI[19,2]:CycleGLCGI[19,3])),rep("F",length(CycleGLCGI[20,2]:CycleGLCGI[20,3])),rep("G",length(CycleGLCGI[21,2]:CycleGLCGI[21,3])))
CutG20LCGI=c(CycleGLCGI[22,2]:CycleGLCGI[22,3],CycleGLCGI[23,2]:CycleGLCGI[23,3],CycleGLCGI[24,2]:CycleGLCGI[24,3],CycleGLCGI[25,2]:CycleGLCGI[25,3],CycleGLCGI[26,2]:CycleGLCGI[26,3],CycleGLCGI[27,2]:CycleGLCGI[27,3],CycleGLCGI[28,2]:CycleGLCGI[28,3])+t3
CycleG20LCGI=c(rep("A",length(CycleGLCGI[22,2]:CycleGLCGI[22,3])),rep("B",length(CycleGLCGI[23,2]:CycleGLCGI[23,3])),rep("C",length(CycleGLCGI[24,2]:CycleGLCGI[24,3])),rep("D",length(CycleGLCGI[25,2]:CycleGLCGI[25,3])),rep("E",length(CycleGLCGI[26,2]:CycleGLCGI[26,3])),rep("F",length(CycleGLCGI[27,2]:CycleGLCGI[27,3])),rep("G",length(CycleGLCGI[28,2]:CycleGLCGI[28,3])))
CutG25LCGI=c(CycleGLCGI[29,2]:CycleGLCGI[29,3],CycleGLCGI[30,2]:CycleGLCGI[30,3],CycleGLCGI[31,2]:CycleGLCGI[31,3],CycleGLCGI[32,2]:CycleGLCGI[32,3],CycleGLCGI[33,2]:CycleGLCGI[33,3],CycleGLCGI[34,2]:CycleGLCGI[34,3],CycleGLCGI[35,2]:CycleGLCGI[35,3])+t4
CycleG25LCGI=c(rep("A",length(CycleGLCGI[29,2]:CycleGLCGI[29,3])),rep("B",length(CycleGLCGI[30,2]:CycleGLCGI[30,3])),rep("C",length(CycleGLCGI[31,2]:CycleGLCGI[31,3])),rep("D",length(CycleGLCGI[32,2]:CycleGLCGI[32,3])),rep("E",length(CycleGLCGI[33,2]:CycleGLCGI[33,3])),rep("F",length(CycleGLCGI[34,2]:CycleGLCGI[34,3])),rep("G",length(CycleGLCGI[35,2]:CycleGLCGI[35,3])))

CutGLCGI=TestLCGI[c(CutG5LCGI,CutG10LCGI,CutG15LCGI,CutG20LCGI,CutG25LCGI),]
CutGLCGI$Cycle=c(CycleG5LCGI,CycleG10LCGI,CycleG15LCGI,CycleG20LCGI,CycleG25LCGI)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLLCGI=CutLLCGI[!CutLLCGI$Cycle=="A",]; CutGLCGI=CutGLCGI[!CutGLCGI$Cycle=="A",]
MinLLCGI=setDT(CutLLCGI)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLLCGI=setDT(CutLLCGI)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGLCGI=setDT(CutGLCGI)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGLCGI=setDT(CutGLCGI)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLLCGI=rbind(MinLLCGI[,c(3,1,2,4)],MaxLLCGI[,c(3,1,2,4)])
BiomGLCGI=rbind(MinGLCGI[,c(3,1,2,5)],MaxGLCGI[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLLCGI=as.data.frame(BiomLLCGI %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLLCGI=SlopLLCGI[order(as.numeric(SlopLLCGI$Temp)),]; SlopLLCGI$Slope=unlist(SlopLLCGI$Slope)
BiomLLCGI=cbind(SlopLLCGI,Max=BiomLLCGI[BiomLLCGI$L>100000]$L)
MaxiLLCGI=setDT(BiomLLCGI)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLLCGI=setDT(BiomLLCGI)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGLCGI=as.data.frame(BiomGLCGI %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGLCGI=SlopGLCGI[order(as.numeric(SlopGLCGI$Temp)),]; SlopGLCGI$Slope=unlist(SlopGLCGI$Slope)
BiomGLCGI=cbind(SlopGLCGI,Max=BiomGLCGI[BiomGLCGI$G>10000]$G)
MaxiGLCGI=setDT(BiomGLCGI)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGLCGI=setDT(BiomGLCGI)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL43=data.frame(MeanLLCGI[,1],MeanLLCGI[,2],MeanLLCGI[,3],TimeLLCGI[,2],TimeLLCGI[,3])
colnames(DataL43)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL")
DataL43$Scenario=rep("LCGI",5)
DataL43$Subscenario=rep("GI",5)

DataG43=data.frame(MeanGLCGI[,1],MeanGLCGI[,2],MeanGLCGI[,3],TimeGLCGI[,2],TimeGLCGI[,3])
colnames(DataG43)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG")
DataG43$Scenario=rep("LCGI",5)
DataG43$Subscenario=rep("GI",5)




####################################################################################################




################################################
### SCENARIO 5: LEAF LITTER BIOMASS INCREASE ###
################################################

#################################################
### SUBSCENARIO 1: GAMMARUS BIOMASS REDUCTION ###
#################################################

# Leaf fall = 900 gC/m2/an = 900 000 mgC/m2/an
# Gammarus density = 10 mgDM/m2 = 5 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5LIGR=GammLeafModel(Temp=5,GammMass=TSRA(5,4.26),Leaf=900000/Days,Gamm=5)
Test5LIGR=Test5LIGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test5LIGR=Test5LIGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5LIGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5°C", cex=1.5, side=3, line=-3, at=300)

Test10LIGR=GammLeafModel(Temp=10,GammMass=TSRA(10,4.26),Leaf=900000/Days,Gamm=5)
Test10LIGR=Test10LIGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test10LIGR=Test10LIGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10LIGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10°C", cex=1.5, side=3, line=-3, at=300)

Test15LIGR=GammLeafModel(Temp=15,GammMass=TSRA(15,4.26),Leaf=900000/Days,Gamm=5)
Test15LIGR=Test15LIGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test15LIGR=Test15LIGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15LIGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15°C", cex=1.5, side=3, line=-3, at=300)

Test20LIGR=GammLeafModel(Temp=20,GammMass=TSRA(20,4.26),Leaf=900000/Days,Gamm=5)
Test20LIGR=Test20LIGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test20LIGR=Test20LIGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20LIGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20°C", cex=1.5, side=3, line=-3, at=300)

Test25LIGR=GammLeafModel(Temp=25,GammMass=TSRA(25,4.26),Leaf=900000/Days,Gamm=5)
Test25LIGR=Test25LIGR %>% mutate(L=if_else(L<10^-3, 0, L))
Test25LIGR=Test25LIGR %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25LIGR[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25°C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))


# Calculate population metabolism and leaf ingestion
Test5LIGR$M=(MetaQuadra(5,TSRM(5,4.26))/1000)*Test5LIGR[,3]
Test5LIGR$I=(AttackQuadra(5,TSRM(5,4.26))*Test5LIGR[,2]/(1+AttackQuadra(5,TSRM(5,4.26))*1/(IngQuadra(5,TSRM(5,4.26))/1000)*Test5LIGR[,2]))*0.30*Test5LIGR[,3]
Test10LIGR$M=(MetaQuadra(10,TSRM(5,4.26))/1000)*Test10LIGR[,3]
Test10LIGR$I=(AttackQuadra(10,TSRM(5,4.26))*Test10LIGR[,2]/(1+AttackQuadra(10,TSRM(5,4.26))*1/(IngQuadra(10,TSRM(5,4.26))/1000)*Test10LIGR[,2]))*0.30*Test10LIGR[,3]
Test15LIGR$M=(MetaQuadra(15,TSRM(5,4.26))/1000)*Test15LIGR[,3]
Test15LIGR$I=(AttackQuadra(15,TSRM(5,4.26))*Test15LIGR[,2]/(1+AttackQuadra(15,TSRM(5,4.26))*1/(IngQuadra(15,TSRM(5,4.26))/1000)*Test15LIGR[,2]))*0.30*Test15LIGR[,3]
Test20LIGR$M=(MetaQuadra(20,TSRM(5,4.26))/1000)*Test20LIGR[,3]
Test20LIGR$I=(AttackQuadra(20,TSRM(5,4.26))*Test20LIGR[,2]/(1+AttackQuadra(20,TSRM(5,4.26))*1/(IngQuadra(20,TSRM(5,4.26))/1000)*Test20LIGR[,2]))*0.30*Test20LIGR[,3]
Test25LIGR$M=(MetaQuadra(25,TSRM(5,4.26))/1000)*Test25LIGR[,3]
Test25LIGR$I=(AttackQuadra(25,TSRM(5,4.26))*Test25LIGR[,2]/(1+AttackQuadra(25,TSRM(5,4.26))*1/(IngQuadra(25,TSRM(5,4.26))/1000)*Test25LIGR[,2]))*0.30*Test25LIGR[,3]

# Combine the dataframes
TestLIGR=rbind(Test5LIGR,Test10LIGR,Test15LIGR,Test20LIGR,Test25LIGR)
TestLIGR$Temp=rep(c("5","10","15","20","25"), each=2556)
TestLIGR$Temp=factor(TestLIGR$Temp, levels=unique(TestLIGR$Temp))
colnames(TestLIGR)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestLIGR=subset(TestLIGR, !Time==2555); TestLIGR$Year=rep(rep(1:(nrow(TestLIGR)/(365*5)), each=365),5)
CutLIGR=as.data.frame(setDT(TestLIGR)[TestLIGR[, tail(.I, -16), by=Year]$V1]); CutLIGR=CutLIGR[!CutLIGR$Year=="1",]
MeanL=as.data.frame(setDT(CutLIGR)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutLIGR)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLLIGR=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGLIGR=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestLIGR=subset(TestLIGR, !Time==2555); TestLIGR$Year=rep(rep(1:(nrow(TestLIGR)/(365*5)), each=365),5)
CutLIGR=as.data.frame(setDT(TestLIGR)[TestLIGR[, tail(.I, -16), by=Year]$V1]); CutLIGR=CutLIGR[!CutLIGR$Year=="1",]
ThLLIGR=CutLIGR[CutLIGR$L>60000,]; ThGLIGR=CutLIGR[CutLIGR$G>5000,]
TimeLLIGR=as.data.frame(setDT(ThLLIGR)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGLIGR=as.data.frame(setDT(ThGLIGR)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLLIGR=as.data.frame(setDT(TimeLLIGR)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGLIGR=as.data.frame(setDT(TimeGLIGR)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLLIGR=as.data.frame(setDT(TestLIGR)[, .(MaxLLIGR=findPeaks(L), MinLLIGR=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGLIGR=as.data.frame(setDT(TestLIGR)[, .(MaxGLIGR=findPeaks(G), MinGLIGR=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5LIGR=c(CycleLLIGR[1,2]:CycleLLIGR[1,3],CycleLLIGR[2,2]:CycleLLIGR[2,3],CycleLLIGR[3,2]:CycleLLIGR[3,3],CycleLLIGR[4,2]:CycleLLIGR[4,3],CycleLLIGR[5,2]:CycleLLIGR[5,3],CycleLLIGR[6,2]:CycleLLIGR[6,3],CycleLLIGR[7,2]:CycleLLIGR[7,3])+t0
CycleL5LIGR=c(rep("A",length(CycleLLIGR[1,2]:CycleLLIGR[1,3])),rep("B",length(CycleLLIGR[2,2]:CycleLLIGR[2,3])),rep("C",length(CycleLLIGR[3,2]:CycleLLIGR[3,3])),rep("D",length(CycleLLIGR[4,2]:CycleLLIGR[4,3])),rep("E",length(CycleLLIGR[5,2]:CycleLLIGR[5,3])),rep("F",length(CycleLLIGR[6,2]:CycleLLIGR[6,3])),rep("G",length(CycleLLIGR[7,2]:CycleLLIGR[7,3])))
CutL10LIGR=c(CycleLLIGR[8,2]:CycleLLIGR[8,3],CycleLLIGR[9,2]:CycleLLIGR[9,3],CycleLLIGR[10,2]:CycleLLIGR[10,3],CycleLLIGR[11,2]:CycleLLIGR[11,3],CycleLLIGR[12,2]:CycleLLIGR[12,3],CycleLLIGR[13,2]:CycleLLIGR[13,3],CycleLLIGR[14,2]:CycleLLIGR[14,3])+t1
CycleL10LIGR=c(rep("A",length(CycleLLIGR[8,2]:CycleLLIGR[8,3])),rep("B",length(CycleLLIGR[9,2]:CycleLLIGR[9,3])),rep("C",length(CycleLLIGR[10,2]:CycleLLIGR[10,3])),rep("D",length(CycleLLIGR[11,2]:CycleLLIGR[11,3])),rep("E",length(CycleLLIGR[12,2]:CycleLLIGR[12,3])),rep("F",length(CycleLLIGR[13,2]:CycleLLIGR[13,3])),rep("G",length(CycleLLIGR[14,2]:CycleLLIGR[14,3])))
CutL15LIGR=c(CycleLLIGR[15,2]:CycleLLIGR[15,3],CycleLLIGR[16,2]:CycleLLIGR[16,3],CycleLLIGR[17,2]:CycleLLIGR[17,3],CycleLLIGR[18,2]:CycleLLIGR[18,3],CycleLLIGR[19,2]:CycleLLIGR[19,3],CycleLLIGR[20,2]:CycleLLIGR[20,3],CycleLLIGR[21,2]:CycleLLIGR[21,3])+t2
CycleL15LIGR=c(rep("A",length(CycleLLIGR[15,2]:CycleLLIGR[15,3])),rep("B",length(CycleLLIGR[16,2]:CycleLLIGR[16,3])),rep("C",length(CycleLLIGR[17,2]:CycleLLIGR[17,3])),rep("D",length(CycleLLIGR[18,2]:CycleLLIGR[18,3])),rep("E",length(CycleLLIGR[19,2]:CycleLLIGR[19,3])),rep("F",length(CycleLLIGR[20,2]:CycleLLIGR[20,3])),rep("G",length(CycleLLIGR[21,2]:CycleLLIGR[21,3])))
CutL20LIGR=c(CycleLLIGR[22,2]:CycleLLIGR[22,3],CycleLLIGR[23,2]:CycleLLIGR[23,3],CycleLLIGR[24,2]:CycleLLIGR[24,3],CycleLLIGR[25,2]:CycleLLIGR[25,3],CycleLLIGR[26,2]:CycleLLIGR[26,3],CycleLLIGR[27,2]:CycleLLIGR[27,3],CycleLLIGR[28,2]:CycleLLIGR[28,3])+t3
CycleL20LIGR=c(rep("A",length(CycleLLIGR[22,2]:CycleLLIGR[22,3])),rep("B",length(CycleLLIGR[23,2]:CycleLLIGR[23,3])),rep("C",length(CycleLLIGR[24,2]:CycleLLIGR[24,3])),rep("D",length(CycleLLIGR[25,2]:CycleLLIGR[25,3])),rep("E",length(CycleLLIGR[26,2]:CycleLLIGR[26,3])),rep("F",length(CycleLLIGR[27,2]:CycleLLIGR[27,3])),rep("G",length(CycleLLIGR[28,2]:CycleLLIGR[28,3])))
CutL25LIGR=c(CycleLLIGR[29,2]:CycleLLIGR[29,3],CycleLLIGR[30,2]:CycleLLIGR[30,3],CycleLLIGR[31,2]:CycleLLIGR[31,3],CycleLLIGR[32,2]:CycleLLIGR[32,3],CycleLLIGR[33,2]:CycleLLIGR[33,3],CycleLLIGR[34,2]:CycleLLIGR[34,3],CycleLLIGR[35,2]:CycleLLIGR[35,3])+t4
CycleL25LIGR=c(rep("A",length(CycleLLIGR[29,2]:CycleLLIGR[29,3])),rep("B",length(CycleLLIGR[30,2]:CycleLLIGR[30,3])),rep("C",length(CycleLLIGR[31,2]:CycleLLIGR[31,3])),rep("D",length(CycleLLIGR[32,2]:CycleLLIGR[32,3])),rep("E",length(CycleLLIGR[33,2]:CycleLLIGR[33,3])),rep("F",length(CycleLLIGR[34,2]:CycleLLIGR[34,3])),rep("G",length(CycleLLIGR[35,2]:CycleLLIGR[35,3])))

CutLLIGR=TestLIGR[c(CutL5LIGR,CutL10LIGR,CutL15LIGR,CutL20LIGR,CutL25LIGR),]
CutLLIGR$Cycle=c(CycleL5LIGR,CycleL10LIGR,CycleL15LIGR,CycleL20LIGR,CycleL25LIGR)

CutG5LIGR=c(CycleGLIGR[1,2]:CycleGLIGR[1,3],CycleGLIGR[2,2]:CycleGLIGR[2,3],CycleGLIGR[3,2]:CycleGLIGR[3,3],CycleGLIGR[4,2]:CycleGLIGR[4,3],CycleGLIGR[5,2]:CycleGLIGR[5,3],CycleGLIGR[6,2]:CycleGLIGR[6,3],CycleGLIGR[7,2]:CycleGLIGR[7,3])+t0
CycleG5LIGR=c(rep("A",length(CycleGLIGR[1,2]:CycleGLIGR[1,3])),rep("B",length(CycleGLIGR[2,2]:CycleGLIGR[2,3])),rep("C",length(CycleGLIGR[3,2]:CycleGLIGR[3,3])),rep("D",length(CycleGLIGR[4,2]:CycleGLIGR[4,3])),rep("E",length(CycleGLIGR[5,2]:CycleGLIGR[5,3])),rep("F",length(CycleGLIGR[6,2]:CycleGLIGR[6,3])),rep("G",length(CycleGLIGR[7,2]:CycleGLIGR[7,3])))
CutG10LIGR=c(CycleGLIGR[8,2]:CycleGLIGR[8,3],CycleGLIGR[9,2]:CycleGLIGR[9,3],CycleGLIGR[10,2]:CycleGLIGR[10,3],CycleGLIGR[11,2]:CycleGLIGR[11,3],CycleGLIGR[12,2]:CycleGLIGR[12,3],CycleGLIGR[13,2]:CycleGLIGR[13,3],CycleGLIGR[14,2]:CycleGLIGR[14,3])+t1
CycleG10LIGR=c(rep("A",length(CycleGLIGR[8,2]:CycleGLIGR[8,3])),rep("B",length(CycleGLIGR[9,2]:CycleGLIGR[9,3])),rep("C",length(CycleGLIGR[10,2]:CycleGLIGR[10,3])),rep("D",length(CycleGLIGR[11,2]:CycleGLIGR[11,3])),rep("E",length(CycleGLIGR[12,2]:CycleGLIGR[12,3])),rep("F",length(CycleGLIGR[13,2]:CycleGLIGR[13,3])),rep("G",length(CycleGLIGR[14,2]:CycleGLIGR[14,3])))
CutG15LIGR=c(CycleGLIGR[15,2]:CycleGLIGR[15,3],CycleGLIGR[16,2]:CycleGLIGR[16,3],CycleGLIGR[17,2]:CycleGLIGR[17,3],CycleGLIGR[18,2]:CycleGLIGR[18,3],CycleGLIGR[19,2]:CycleGLIGR[19,3],CycleGLIGR[20,2]:CycleGLIGR[20,3],CycleGLIGR[21,2]:CycleGLIGR[21,3])+t2
CycleG15LIGR=c(rep("A",length(CycleGLIGR[15,2]:CycleGLIGR[15,3])),rep("B",length(CycleGLIGR[16,2]:CycleGLIGR[16,3])),rep("C",length(CycleGLIGR[17,2]:CycleGLIGR[17,3])),rep("D",length(CycleGLIGR[18,2]:CycleGLIGR[18,3])),rep("E",length(CycleGLIGR[19,2]:CycleGLIGR[19,3])),rep("F",length(CycleGLIGR[20,2]:CycleGLIGR[20,3])),rep("G",length(CycleGLIGR[21,2]:CycleGLIGR[21,3])))
CutG20LIGR=c(CycleGLIGR[22,2]:CycleGLIGR[22,3],CycleGLIGR[23,2]:CycleGLIGR[23,3],CycleGLIGR[24,2]:CycleGLIGR[24,3],CycleGLIGR[25,2]:CycleGLIGR[25,3],CycleGLIGR[26,2]:CycleGLIGR[26,3],CycleGLIGR[27,2]:CycleGLIGR[27,3],CycleGLIGR[28,2]:CycleGLIGR[28,3])+t3
CycleG20LIGR=c(rep("A",length(CycleGLIGR[22,2]:CycleGLIGR[22,3])),rep("B",length(CycleGLIGR[23,2]:CycleGLIGR[23,3])),rep("C",length(CycleGLIGR[24,2]:CycleGLIGR[24,3])),rep("D",length(CycleGLIGR[25,2]:CycleGLIGR[25,3])),rep("E",length(CycleGLIGR[26,2]:CycleGLIGR[26,3])),rep("F",length(CycleGLIGR[27,2]:CycleGLIGR[27,3])),rep("G",length(CycleGLIGR[28,2]:CycleGLIGR[28,3])))
CutG25LIGR=c(CycleGLIGR[29,2]:CycleGLIGR[29,3],CycleGLIGR[30,2]:CycleGLIGR[30,3],CycleGLIGR[31,2]:CycleGLIGR[31,3],CycleGLIGR[32,2]:CycleGLIGR[32,3],CycleGLIGR[33,2]:CycleGLIGR[33,3],CycleGLIGR[34,2]:CycleGLIGR[34,3],CycleGLIGR[35,2]:CycleGLIGR[35,3])+t4
CycleG25LIGR=c(rep("A",length(CycleGLIGR[29,2]:CycleGLIGR[29,3])),rep("B",length(CycleGLIGR[30,2]:CycleGLIGR[30,3])),rep("C",length(CycleGLIGR[31,2]:CycleGLIGR[31,3])),rep("D",length(CycleGLIGR[32,2]:CycleGLIGR[32,3])),rep("E",length(CycleGLIGR[33,2]:CycleGLIGR[33,3])),rep("F",length(CycleGLIGR[34,2]:CycleGLIGR[34,3])),rep("G",length(CycleGLIGR[35,2]:CycleGLIGR[35,3])))

CutGLIGR=TestLIGR[c(CutG5LIGR,CutG10LIGR,CutG15LIGR,CutG20LIGR,CutG25LIGR),]
CutGLIGR$Cycle=c(CycleG5LIGR,CycleG10LIGR,CycleG15LIGR,CycleG20LIGR,CycleG25LIGR)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLLIGR=CutLLIGR[!CutLLIGR$Cycle=="A",]; CutGLIGR=CutGLIGR[!CutGLIGR$Cycle=="A",]
MinLLIGR=setDT(CutLLIGR)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLLIGR=setDT(CutLLIGR)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGLIGR=setDT(CutGLIGR)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGLIGR=setDT(CutGLIGR)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLLIGR=rbind(MinLLIGR[,c(3,1,2,4)],MaxLLIGR[,c(3,1,2,4)])
BiomGLIGR=rbind(MinGLIGR[,c(3,1,2,5)],MaxGLIGR[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLLIGR=as.data.frame(BiomLLIGR %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLLIGR=SlopLLIGR[order(as.numeric(SlopLLIGR$Temp)),]; SlopLLIGR$Slope=unlist(SlopLLIGR$Slope)
BiomLLIGR=cbind(SlopLLIGR,Max=BiomLLIGR[BiomLLIGR$L>100000]$L)
MaxiLLIGR=setDT(BiomLLIGR)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLLIGR=setDT(BiomLLIGR)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGLIGR=as.data.frame(BiomGLIGR %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGLIGR=SlopGLIGR[order(as.numeric(SlopGLIGR$Temp)),]; SlopGLIGR$Slope=unlist(SlopGLIGR$Slope)
BiomGLIGR=cbind(SlopGLIGR,Max=BiomGLIGR[BiomGLIGR$G>10000]$G)
MaxiGLIGR=setDT(BiomGLIGR)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGLIGR=setDT(BiomGLIGR)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL51=data.frame(MeanLLIGR[,1],MeanLLIGR[,2],MeanLLIGR[,3],TimeLLIGR[,2],TimeLLIGR[,3])
colnames(DataL51)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL")
DataL51$Scenario=rep("LIGR",5)
DataL51$Subscenario=rep("GR",5)

DataG51=data.frame(MeanGLIGR[,1],MeanGLIGR[,2],MeanGLIGR[,3],TimeGLIGR[,2],TimeGLIGR[,3])
colnames(DataG51)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG")
DataG51$Scenario=rep("LIGR",5)
DataG51$Subscenario=rep("GR",5)


################################################
### SUBSCENARIO 2: GAMMARUS BIOMASS CONSTANT ###
################################################

# Leaf fall = 900 gC/m2/an = 900 000 mgC/m2/an
# Gammarus density = 30 mgDM/m2 = 15 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5LIGC=GammLeafModel(Temp=5,GammMass=TSRA(5,4.26),Leaf=900000/Days,Gamm=15)
Test5LIGC=Test5LIGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test5LIGC=Test5LIGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5LIGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5°C", cex=1.5, side=3, line=-3, at=300)

Test10LIGC=GammLeafModel(Temp=10,GammMass=TSRA(10,4.26),Leaf=900000/Days,Gamm=15)
Test10LIGC=Test10LIGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test10LIGC=Test10LIGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10LIGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10°C", cex=1.5, side=3, line=-3, at=300)

Test15LIGC=GammLeafModel(Temp=15,GammMass=TSRA(15,4.26),Leaf=900000/Days,Gamm=15)
Test15LIGC=Test15LIGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test15LIGC=Test15LIGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15LIGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15°C", cex=1.5, side=3, line=-3, at=300)

Test20LIGC=GammLeafModel(Temp=20,GammMass=TSRA(20,4.26),Leaf=900000/Days,Gamm=15)
Test20LIGC=Test20LIGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test20LIGC=Test20LIGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20LIGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20°C", cex=1.5, side=3, line=-3, at=300)

Test25LIGC=GammLeafModel(Temp=25,GammMass=TSRA(25,4.26),Leaf=900000/Days,Gamm=15)
Test25LIGC=Test25LIGC %>% mutate(L=if_else(L<10^-3, 0, L))
Test25LIGC=Test25LIGC %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25LIGC[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25°C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))


# Calculate population metabolism and leaf ingestion
Test5LIGC$M=(MetaQuadra(5,TSRM(5,4.26))/1000)*Test5LIGC[,3]
Test5LIGC$I=(AttackQuadra(5,TSRM(5,4.26))*Test5LIGC[,2]/(1+AttackQuadra(5,TSRM(5,4.26))*1/(IngQuadra(5,TSRM(5,4.26))/1000)*Test5LIGC[,2]))*0.30*Test5LIGC[,3]
Test10LIGC$M=(MetaQuadra(10,TSRM(5,4.26))/1000)*Test10LIGC[,3]
Test10LIGC$I=(AttackQuadra(10,TSRM(5,4.26))*Test10LIGC[,2]/(1+AttackQuadra(10,TSRM(5,4.26))*1/(IngQuadra(10,TSRM(5,4.26))/1000)*Test10LIGC[,2]))*0.30*Test10LIGC[,3]
Test15LIGC$M=(MetaQuadra(15,TSRM(5,4.26))/1000)*Test15LIGC[,3]
Test15LIGC$I=(AttackQuadra(15,TSRM(5,4.26))*Test15LIGC[,2]/(1+AttackQuadra(15,TSRM(5,4.26))*1/(IngQuadra(15,TSRM(5,4.26))/1000)*Test15LIGC[,2]))*0.30*Test15LIGC[,3]
Test20LIGC$M=(MetaQuadra(20,TSRM(5,4.26))/1000)*Test20LIGC[,3]
Test20LIGC$I=(AttackQuadra(20,TSRM(5,4.26))*Test20LIGC[,2]/(1+AttackQuadra(20,TSRM(5,4.26))*1/(IngQuadra(20,TSRM(5,4.26))/1000)*Test20LIGC[,2]))*0.30*Test20LIGC[,3]
Test25LIGC$M=(MetaQuadra(25,TSRM(5,4.26))/1000)*Test25LIGC[,3]
Test25LIGC$I=(AttackQuadra(25,TSRM(5,4.26))*Test25LIGC[,2]/(1+AttackQuadra(25,TSRM(5,4.26))*1/(IngQuadra(25,TSRM(5,4.26))/1000)*Test25LIGC[,2]))*0.30*Test25LIGC[,3]

# Combine the dataframes
TestLIGC=rbind(Test5LIGC,Test10LIGC,Test15LIGC,Test20LIGC,Test25LIGC)
TestLIGC$Temp=rep(c("5","10","15","20","25"), each=2556)
TestLIGC$Temp=factor(TestLIGC$Temp, levels=unique(TestLIGC$Temp))
colnames(TestLIGC)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestLIGC=subset(TestLIGC, !Time==2555); TestLIGC$Year=rep(rep(1:(nrow(TestLIGC)/(365*5)), each=365),5)
CutLIGC=as.data.frame(setDT(TestLIGC)[TestLIGC[, tail(.I, -16), by=Year]$V1]); CutLIGC=CutLIGC[!CutLIGC$Year=="1",]
MeanL=as.data.frame(setDT(CutLIGC)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutLIGC)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLLIGC=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGLIGC=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestLIGC=subset(TestLIGC, !Time==2555); TestLIGC$Year=rep(rep(1:(nrow(TestLIGC)/(365*5)), each=365),5)
CutLIGC=as.data.frame(setDT(TestLIGC)[TestLIGC[, tail(.I, -16), by=Year]$V1]); CutLIGC=CutLIGC[!CutLIGC$Year=="1",]
ThLLIGC=CutLIGC[CutLIGC$L>60000,]; ThGLIGC=CutLIGC[CutLIGC$G>5000,]
TimeLLIGC=as.data.frame(setDT(ThLLIGC)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGLIGC=as.data.frame(setDT(ThGLIGC)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLLIGC=as.data.frame(setDT(TimeLLIGC)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGLIGC=as.data.frame(setDT(TimeGLIGC)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLLIGC=as.data.frame(setDT(TestLIGC)[, .(MaxLLIGC=findPeaks(L), MinLLIGC=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGLIGC=as.data.frame(setDT(TestLIGC)[, .(MaxGLIGC=findPeaks(G), MinGLIGC=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5LIGC=c(CycleLLIGC[1,2]:CycleLLIGC[1,3],CycleLLIGC[2,2]:CycleLLIGC[2,3],CycleLLIGC[3,2]:CycleLLIGC[3,3],CycleLLIGC[4,2]:CycleLLIGC[4,3],CycleLLIGC[5,2]:CycleLLIGC[5,3],CycleLLIGC[6,2]:CycleLLIGC[6,3],CycleLLIGC[7,2]:CycleLLIGC[7,3])+t0
CycleL5LIGC=c(rep("A",length(CycleLLIGC[1,2]:CycleLLIGC[1,3])),rep("B",length(CycleLLIGC[2,2]:CycleLLIGC[2,3])),rep("C",length(CycleLLIGC[3,2]:CycleLLIGC[3,3])),rep("D",length(CycleLLIGC[4,2]:CycleLLIGC[4,3])),rep("E",length(CycleLLIGC[5,2]:CycleLLIGC[5,3])),rep("F",length(CycleLLIGC[6,2]:CycleLLIGC[6,3])),rep("G",length(CycleLLIGC[7,2]:CycleLLIGC[7,3])))
CutL10LIGC=c(CycleLLIGC[8,2]:CycleLLIGC[8,3],CycleLLIGC[9,2]:CycleLLIGC[9,3],CycleLLIGC[10,2]:CycleLLIGC[10,3],CycleLLIGC[11,2]:CycleLLIGC[11,3],CycleLLIGC[12,2]:CycleLLIGC[12,3],CycleLLIGC[13,2]:CycleLLIGC[13,3],CycleLLIGC[14,2]:CycleLLIGC[14,3])+t1
CycleL10LIGC=c(rep("A",length(CycleLLIGC[8,2]:CycleLLIGC[8,3])),rep("B",length(CycleLLIGC[9,2]:CycleLLIGC[9,3])),rep("C",length(CycleLLIGC[10,2]:CycleLLIGC[10,3])),rep("D",length(CycleLLIGC[11,2]:CycleLLIGC[11,3])),rep("E",length(CycleLLIGC[12,2]:CycleLLIGC[12,3])),rep("F",length(CycleLLIGC[13,2]:CycleLLIGC[13,3])),rep("G",length(CycleLLIGC[14,2]:CycleLLIGC[14,3])))
CutL15LIGC=c(CycleLLIGC[15,2]:CycleLLIGC[15,3],CycleLLIGC[16,2]:CycleLLIGC[16,3],CycleLLIGC[17,2]:CycleLLIGC[17,3],CycleLLIGC[18,2]:CycleLLIGC[18,3],CycleLLIGC[19,2]:CycleLLIGC[19,3],CycleLLIGC[20,2]:CycleLLIGC[20,3],CycleLLIGC[21,2]:CycleLLIGC[21,3])+t2
CycleL15LIGC=c(rep("A",length(CycleLLIGC[15,2]:CycleLLIGC[15,3])),rep("B",length(CycleLLIGC[16,2]:CycleLLIGC[16,3])),rep("C",length(CycleLLIGC[17,2]:CycleLLIGC[17,3])),rep("D",length(CycleLLIGC[18,2]:CycleLLIGC[18,3])),rep("E",length(CycleLLIGC[19,2]:CycleLLIGC[19,3])),rep("F",length(CycleLLIGC[20,2]:CycleLLIGC[20,3])),rep("G",length(CycleLLIGC[21,2]:CycleLLIGC[21,3])))
CutL20LIGC=c(CycleLLIGC[22,2]:CycleLLIGC[22,3],CycleLLIGC[23,2]:CycleLLIGC[23,3],CycleLLIGC[24,2]:CycleLLIGC[24,3],CycleLLIGC[25,2]:CycleLLIGC[25,3],CycleLLIGC[26,2]:CycleLLIGC[26,3],CycleLLIGC[27,2]:CycleLLIGC[27,3],CycleLLIGC[28,2]:CycleLLIGC[28,3])+t3
CycleL20LIGC=c(rep("A",length(CycleLLIGC[22,2]:CycleLLIGC[22,3])),rep("B",length(CycleLLIGC[23,2]:CycleLLIGC[23,3])),rep("C",length(CycleLLIGC[24,2]:CycleLLIGC[24,3])),rep("D",length(CycleLLIGC[25,2]:CycleLLIGC[25,3])),rep("E",length(CycleLLIGC[26,2]:CycleLLIGC[26,3])),rep("F",length(CycleLLIGC[27,2]:CycleLLIGC[27,3])),rep("G",length(CycleLLIGC[28,2]:CycleLLIGC[28,3])))
CutL25LIGC=c(CycleLLIGC[29,2]:CycleLLIGC[29,3],CycleLLIGC[30,2]:CycleLLIGC[30,3],CycleLLIGC[31,2]:CycleLLIGC[31,3],CycleLLIGC[32,2]:CycleLLIGC[32,3],CycleLLIGC[33,2]:CycleLLIGC[33,3],CycleLLIGC[34,2]:CycleLLIGC[34,3],CycleLLIGC[35,2]:CycleLLIGC[35,3])+t4
CycleL25LIGC=c(rep("A",length(CycleLLIGC[29,2]:CycleLLIGC[29,3])),rep("B",length(CycleLLIGC[30,2]:CycleLLIGC[30,3])),rep("C",length(CycleLLIGC[31,2]:CycleLLIGC[31,3])),rep("D",length(CycleLLIGC[32,2]:CycleLLIGC[32,3])),rep("E",length(CycleLLIGC[33,2]:CycleLLIGC[33,3])),rep("F",length(CycleLLIGC[34,2]:CycleLLIGC[34,3])),rep("G",length(CycleLLIGC[35,2]:CycleLLIGC[35,3])))

CutLLIGC=TestLIGC[c(CutL5LIGC,CutL10LIGC,CutL15LIGC,CutL20LIGC,CutL25LIGC),]
CutLLIGC$Cycle=c(CycleL5LIGC,CycleL10LIGC,CycleL15LIGC,CycleL20LIGC,CycleL25LIGC)

CutG5LIGC=c(CycleGLIGC[1,2]:CycleGLIGC[1,3],CycleGLIGC[2,2]:CycleGLIGC[2,3],CycleGLIGC[3,2]:CycleGLIGC[3,3],CycleGLIGC[4,2]:CycleGLIGC[4,3],CycleGLIGC[5,2]:CycleGLIGC[5,3],CycleGLIGC[6,2]:CycleGLIGC[6,3],CycleGLIGC[7,2]:CycleGLIGC[7,3])+t0
CycleG5LIGC=c(rep("A",length(CycleGLIGC[1,2]:CycleGLIGC[1,3])),rep("B",length(CycleGLIGC[2,2]:CycleGLIGC[2,3])),rep("C",length(CycleGLIGC[3,2]:CycleGLIGC[3,3])),rep("D",length(CycleGLIGC[4,2]:CycleGLIGC[4,3])),rep("E",length(CycleGLIGC[5,2]:CycleGLIGC[5,3])),rep("F",length(CycleGLIGC[6,2]:CycleGLIGC[6,3])),rep("G",length(CycleGLIGC[7,2]:CycleGLIGC[7,3])))
CutG10LIGC=c(CycleGLIGC[8,2]:CycleGLIGC[8,3],CycleGLIGC[9,2]:CycleGLIGC[9,3],CycleGLIGC[10,2]:CycleGLIGC[10,3],CycleGLIGC[11,2]:CycleGLIGC[11,3],CycleGLIGC[12,2]:CycleGLIGC[12,3],CycleGLIGC[13,2]:CycleGLIGC[13,3],CycleGLIGC[14,2]:CycleGLIGC[14,3])+t1
CycleG10LIGC=c(rep("A",length(CycleGLIGC[8,2]:CycleGLIGC[8,3])),rep("B",length(CycleGLIGC[9,2]:CycleGLIGC[9,3])),rep("C",length(CycleGLIGC[10,2]:CycleGLIGC[10,3])),rep("D",length(CycleGLIGC[11,2]:CycleGLIGC[11,3])),rep("E",length(CycleGLIGC[12,2]:CycleGLIGC[12,3])),rep("F",length(CycleGLIGC[13,2]:CycleGLIGC[13,3])),rep("G",length(CycleGLIGC[14,2]:CycleGLIGC[14,3])))
CutG15LIGC=c(CycleGLIGC[15,2]:CycleGLIGC[15,3],CycleGLIGC[16,2]:CycleGLIGC[16,3],CycleGLIGC[17,2]:CycleGLIGC[17,3],CycleGLIGC[18,2]:CycleGLIGC[18,3],CycleGLIGC[19,2]:CycleGLIGC[19,3],CycleGLIGC[20,2]:CycleGLIGC[20,3],CycleGLIGC[21,2]:CycleGLIGC[21,3])+t2
CycleG15LIGC=c(rep("A",length(CycleGLIGC[15,2]:CycleGLIGC[15,3])),rep("B",length(CycleGLIGC[16,2]:CycleGLIGC[16,3])),rep("C",length(CycleGLIGC[17,2]:CycleGLIGC[17,3])),rep("D",length(CycleGLIGC[18,2]:CycleGLIGC[18,3])),rep("E",length(CycleGLIGC[19,2]:CycleGLIGC[19,3])),rep("F",length(CycleGLIGC[20,2]:CycleGLIGC[20,3])),rep("G",length(CycleGLIGC[21,2]:CycleGLIGC[21,3])))
CutG20LIGC=c(CycleGLIGC[22,2]:CycleGLIGC[22,3],CycleGLIGC[23,2]:CycleGLIGC[23,3],CycleGLIGC[24,2]:CycleGLIGC[24,3],CycleGLIGC[25,2]:CycleGLIGC[25,3],CycleGLIGC[26,2]:CycleGLIGC[26,3],CycleGLIGC[27,2]:CycleGLIGC[27,3],CycleGLIGC[28,2]:CycleGLIGC[28,3])+t3
CycleG20LIGC=c(rep("A",length(CycleGLIGC[22,2]:CycleGLIGC[22,3])),rep("B",length(CycleGLIGC[23,2]:CycleGLIGC[23,3])),rep("C",length(CycleGLIGC[24,2]:CycleGLIGC[24,3])),rep("D",length(CycleGLIGC[25,2]:CycleGLIGC[25,3])),rep("E",length(CycleGLIGC[26,2]:CycleGLIGC[26,3])),rep("F",length(CycleGLIGC[27,2]:CycleGLIGC[27,3])),rep("G",length(CycleGLIGC[28,2]:CycleGLIGC[28,3])))
CutG25LIGC=c(CycleGLIGC[29,2]:CycleGLIGC[29,3],CycleGLIGC[30,2]:CycleGLIGC[30,3],CycleGLIGC[31,2]:CycleGLIGC[31,3],CycleGLIGC[32,2]:CycleGLIGC[32,3],CycleGLIGC[33,2]:CycleGLIGC[33,3],CycleGLIGC[34,2]:CycleGLIGC[34,3],CycleGLIGC[35,2]:CycleGLIGC[35,3])+t4
CycleG25LIGC=c(rep("A",length(CycleGLIGC[29,2]:CycleGLIGC[29,3])),rep("B",length(CycleGLIGC[30,2]:CycleGLIGC[30,3])),rep("C",length(CycleGLIGC[31,2]:CycleGLIGC[31,3])),rep("D",length(CycleGLIGC[32,2]:CycleGLIGC[32,3])),rep("E",length(CycleGLIGC[33,2]:CycleGLIGC[33,3])),rep("F",length(CycleGLIGC[34,2]:CycleGLIGC[34,3])),rep("G",length(CycleGLIGC[35,2]:CycleGLIGC[35,3])))

CutGLIGC=TestLIGC[c(CutG5LIGC,CutG10LIGC,CutG15LIGC,CutG20LIGC,CutG25LIGC),]
CutGLIGC$Cycle=c(CycleG5LIGC,CycleG10LIGC,CycleG15LIGC,CycleG20LIGC,CycleG25LIGC)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLLIGC=CutLLIGC[!CutLLIGC$Cycle=="A",]; CutGLIGC=CutGLIGC[!CutGLIGC$Cycle=="A",]
MinLLIGC=setDT(CutLLIGC)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLLIGC=setDT(CutLLIGC)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGLIGC=setDT(CutGLIGC)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGLIGC=setDT(CutGLIGC)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLLIGC=rbind(MinLLIGC[,c(3,1,2,4)],MaxLLIGC[,c(3,1,2,4)])
BiomGLIGC=rbind(MinGLIGC[,c(3,1,2,5)],MaxGLIGC[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLLIGC=as.data.frame(BiomLLIGC %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLLIGC=SlopLLIGC[order(as.numeric(SlopLLIGC$Temp)),]; SlopLLIGC$Slope=unlist(SlopLLIGC$Slope)
BiomLLIGC=cbind(SlopLLIGC,Max=BiomLLIGC[BiomLLIGC$L>100000]$L)
MaxiLLIGC=setDT(BiomLLIGC)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLLIGC=setDT(BiomLLIGC)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGLIGC=as.data.frame(BiomGLIGC %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGLIGC=SlopGLIGC[order(as.numeric(SlopGLIGC$Temp)),]; SlopGLIGC$Slope=unlist(SlopGLIGC$Slope)
BiomGLIGC=cbind(SlopGLIGC,Max=BiomGLIGC[BiomGLIGC$G>10000]$G)
MaxiGLIGC=setDT(BiomGLIGC)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGLIGC=setDT(BiomGLIGC)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL52=data.frame(MeanLLIGC[,1],MeanLLIGC[,2],MeanLLIGC[,3],TimeLLIGC[,2],TimeLLIGC[,3])
colnames(DataL52)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL")
DataL52$Scenario=rep("LIGC",5)
DataL52$Subscenario=rep("GC",5)

DataG52=data.frame(MeanGLIGC[,1],MeanGLIGC[,2],MeanGLIGC[,3],TimeGLIGC[,2],TimeGLIGC[,3])
colnames(DataG52)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG")
DataG52$Scenario=rep("LIGC",5)
DataG52$Subscenario=rep("GC",5)


################################################
### SUBSCENARIO 3: GAMMARUS BIOMASS INCREASE ###
################################################

# Leaf fall = 900 gC/m2/an = 900 000 mgC/m2/an
# Gammarus density = 90 mgDM/m2 = 45 mgC/m2
# Gammarus mean body mass = 4.26 mgDM

# Duration of the leaves fall
Days=15

# Plot the dynamics for 5 temperatures
par(mfrow=c(2,3))
par(oma=c(4,4,1,1))
par(mar=c(3,3,0,0))

Test5LIGI=GammLeafModel(Temp=5,GammMass=TSRA(5,4.26),Leaf=900000/Days,Gamm=45)
Test5LIGI=Test5LIGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test5LIGI=Test5LIGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test5LIGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("5°C", cex=1.5, side=3, line=-3, at=300)

Test10LIGI=GammLeafModel(Temp=10,GammMass=TSRA(10,4.26),Leaf=900000/Days,Gamm=45)
Test10LIGI=Test10LIGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test10LIGI=Test10LIGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test10LIGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("10°C", cex=1.5, side=3, line=-3, at=300)

Test15LIGI=GammLeafModel(Temp=15,GammMass=TSRA(15,4.26),Leaf=900000/Days,Gamm=45)
Test15LIGI=Test15LIGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test15LIGI=Test15LIGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test15LIGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("15°C", cex=1.5, side=3, line=-3, at=300)

Test20LIGI=GammLeafModel(Temp=20,GammMass=TSRA(20,4.26),Leaf=900000/Days,Gamm=45)
Test20LIGI=Test20LIGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test20LIGI=Test20LIGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test20LIGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("20°C", cex=1.5, side=3, line=-3, at=300)

Test25LIGI=GammLeafModel(Temp=25,GammMass=TSRA(25,4.26),Leaf=900000/Days,Gamm=45)
Test25LIGI=Test25LIGI %>% mutate(L=if_else(L<10^-3, 0, L))
Test25LIGI=Test25LIGI %>% mutate(G=if_else(G<10^-3, 0, G))
matplot(Test25LIGI[,-1]/10^5, type="l", ylim=c(0,10), ylab="", xlab="", cex.axis=2.0, las=1, col=c("black","tomato2"), lty=c(1,1,1), lwd=c(1.5,2))
mtext("25°C", cex=1.5, side=3, line=-3, at=300)

mtext("Time (d)", side=1, outer=T, line=1, cex=1.5)
mtext(expression('Biomass'~'('*10^5~mg~C~m^-2*')'), side=2, outer=T, line=1, cex=1.5)
par(mfrow=c(1,1))


# Calculate population metabolism and leaf ingestion
Test5LIGI$M=(MetaQuadra(5,TSRM(5,4.26))/1000)*Test5LIGI[,3]
Test5LIGI$I=(AttackQuadra(5,TSRM(5,4.26))*Test5LIGI[,2]/(1+AttackQuadra(5,TSRM(5,4.26))*1/(IngQuadra(5,TSRM(5,4.26))/1000)*Test5LIGI[,2]))*0.30*Test5LIGI[,3]
Test10LIGI$M=(MetaQuadra(10,TSRM(5,4.26))/1000)*Test10LIGI[,3]
Test10LIGI$I=(AttackQuadra(10,TSRM(5,4.26))*Test10LIGI[,2]/(1+AttackQuadra(10,TSRM(5,4.26))*1/(IngQuadra(10,TSRM(5,4.26))/1000)*Test10LIGI[,2]))*0.30*Test10LIGI[,3]
Test15LIGI$M=(MetaQuadra(15,TSRM(5,4.26))/1000)*Test15LIGI[,3]
Test15LIGI$I=(AttackQuadra(15,TSRM(5,4.26))*Test15LIGI[,2]/(1+AttackQuadra(15,TSRM(5,4.26))*1/(IngQuadra(15,TSRM(5,4.26))/1000)*Test15LIGI[,2]))*0.30*Test15LIGI[,3]
Test20LIGI$M=(MetaQuadra(20,TSRM(5,4.26))/1000)*Test20LIGI[,3]
Test20LIGI$I=(AttackQuadra(20,TSRM(5,4.26))*Test20LIGI[,2]/(1+AttackQuadra(20,TSRM(5,4.26))*1/(IngQuadra(20,TSRM(5,4.26))/1000)*Test20LIGI[,2]))*0.30*Test20LIGI[,3]
Test25LIGI$M=(MetaQuadra(25,TSRM(5,4.26))/1000)*Test25LIGI[,3]
Test25LIGI$I=(AttackQuadra(25,TSRM(5,4.26))*Test25LIGI[,2]/(1+AttackQuadra(25,TSRM(5,4.26))*1/(IngQuadra(25,TSRM(5,4.26))/1000)*Test25LIGI[,2]))*0.30*Test25LIGI[,3]

# Combine the dataframes
TestLIGI=rbind(Test5LIGI,Test10LIGI,Test15LIGI,Test20LIGI,Test25LIGI)
TestLIGI$Temp=rep(c("5","10","15","20","25"), each=2556)
TestLIGI$Temp=factor(TestLIGI$Temp, levels=unique(TestLIGI$Temp))
colnames(TestLIGI)[1]="Time"


#############################
### Mean annual biomasses ### 
#############################

# Calculate litter and Gammarus mean annual biomasses
TestLIGI=subset(TestLIGI, !Time==2555); TestLIGI$Year=rep(rep(1:(nrow(TestLIGI)/(365*5)), each=365),5)
CutLIGI=as.data.frame(setDT(TestLIGI)[TestLIGI[, tail(.I, -16), by=Year]$V1]); CutLIGI=CutLIGI[!CutLIGI$Year=="1",]
MeanL=as.data.frame(setDT(CutLIGI)[, .(MeanL=mean(L)), by=list(Temp,Year)])
MeanG=as.data.frame(setDT(CutLIGI)[, .(MeanG=mean(G)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for annual biomasses
MeanLLIGI=as.data.frame(setDT(MeanL)[, .(MeanL=mean(MeanL), SdL=sd(MeanL)), by=list(Temp)])
MeanGLIGI=as.data.frame(setDT(MeanG)[, .(MeanG=mean(MeanG), SdG=sd(MeanG)), by=list(Temp)])


##############################################################
### Mean annual persistence times above biomass thresholds ###
##############################################################

# Calculate litter and Gammarus mean annual persistence times
TestLIGI=subset(TestLIGI, !Time==2555); TestLIGI$Year=rep(rep(1:(nrow(TestLIGI)/(365*5)), each=365),5)
CutLIGI=as.data.frame(setDT(TestLIGI)[TestLIGI[, tail(.I, -16), by=Year]$V1]); CutLIGI=CutLIGI[!CutLIGI$Year=="1",]
ThLLIGI=CutLIGI[CutLIGI$L>60000,]; ThGLIGI=CutLIGI[CutLIGI$G>5000,]
TimeLLIGI=as.data.frame(setDT(ThLLIGI)[, .(TimeL=length(Time)), by=list(Temp,Year)])
TimeGLIGI=as.data.frame(setDT(ThGLIGI)[, .(TimeG=length(Time)), by=list(Temp,Year)])

# Calculate means and deviations over 6 years for the thresholds
TimeLLIGI=as.data.frame(setDT(TimeLLIGI)[, .(MeanTimeL=mean(TimeL), SdTimeL=sd(TimeL)), by=list(Temp)])
TimeGLIGI=as.data.frame(setDT(TimeGLIGI)[, .(MeanTimeG=mean(TimeG), SdTimeG=sd(TimeG)), by=list(Temp)])


###################################################
### Mean biomass maximums and biomass decreases ###
###################################################

# Find biomass maximums and minimums
CycleLLIGI=as.data.frame(setDT(TestLIGI)[, .(MaxLLIGI=findPeaks(L), MinLLIGI=findValleys(L)[seq(2,14,2)]), by=list(Temp)])
CycleGLIGI=as.data.frame(setDT(TestLIGI)[, .(MaxGLIGI=findPeaks(G), MinGLIGI=c(findValleys(G),2555)), by=list(Temp)])

# Define litter and Gammarus biomass cycles
t0=2555*0; t1=2555*1; t2=2555*2; t3=2555*3; t4=2555*4

CutL5LIGI=c(CycleLLIGI[1,2]:CycleLLIGI[1,3],CycleLLIGI[2,2]:CycleLLIGI[2,3],CycleLLIGI[3,2]:CycleLLIGI[3,3],CycleLLIGI[4,2]:CycleLLIGI[4,3],CycleLLIGI[5,2]:CycleLLIGI[5,3],CycleLLIGI[6,2]:CycleLLIGI[6,3],CycleLLIGI[7,2]:CycleLLIGI[7,3])+t0
CycleL5LIGI=c(rep("A",length(CycleLLIGI[1,2]:CycleLLIGI[1,3])),rep("B",length(CycleLLIGI[2,2]:CycleLLIGI[2,3])),rep("C",length(CycleLLIGI[3,2]:CycleLLIGI[3,3])),rep("D",length(CycleLLIGI[4,2]:CycleLLIGI[4,3])),rep("E",length(CycleLLIGI[5,2]:CycleLLIGI[5,3])),rep("F",length(CycleLLIGI[6,2]:CycleLLIGI[6,3])),rep("G",length(CycleLLIGI[7,2]:CycleLLIGI[7,3])))
CutL10LIGI=c(CycleLLIGI[8,2]:CycleLLIGI[8,3],CycleLLIGI[9,2]:CycleLLIGI[9,3],CycleLLIGI[10,2]:CycleLLIGI[10,3],CycleLLIGI[11,2]:CycleLLIGI[11,3],CycleLLIGI[12,2]:CycleLLIGI[12,3],CycleLLIGI[13,2]:CycleLLIGI[13,3],CycleLLIGI[14,2]:CycleLLIGI[14,3])+t1
CycleL10LIGI=c(rep("A",length(CycleLLIGI[8,2]:CycleLLIGI[8,3])),rep("B",length(CycleLLIGI[9,2]:CycleLLIGI[9,3])),rep("C",length(CycleLLIGI[10,2]:CycleLLIGI[10,3])),rep("D",length(CycleLLIGI[11,2]:CycleLLIGI[11,3])),rep("E",length(CycleLLIGI[12,2]:CycleLLIGI[12,3])),rep("F",length(CycleLLIGI[13,2]:CycleLLIGI[13,3])),rep("G",length(CycleLLIGI[14,2]:CycleLLIGI[14,3])))
CutL15LIGI=c(CycleLLIGI[15,2]:CycleLLIGI[15,3],CycleLLIGI[16,2]:CycleLLIGI[16,3],CycleLLIGI[17,2]:CycleLLIGI[17,3],CycleLLIGI[18,2]:CycleLLIGI[18,3],CycleLLIGI[19,2]:CycleLLIGI[19,3],CycleLLIGI[20,2]:CycleLLIGI[20,3],CycleLLIGI[21,2]:CycleLLIGI[21,3])+t2
CycleL15LIGI=c(rep("A",length(CycleLLIGI[15,2]:CycleLLIGI[15,3])),rep("B",length(CycleLLIGI[16,2]:CycleLLIGI[16,3])),rep("C",length(CycleLLIGI[17,2]:CycleLLIGI[17,3])),rep("D",length(CycleLLIGI[18,2]:CycleLLIGI[18,3])),rep("E",length(CycleLLIGI[19,2]:CycleLLIGI[19,3])),rep("F",length(CycleLLIGI[20,2]:CycleLLIGI[20,3])),rep("G",length(CycleLLIGI[21,2]:CycleLLIGI[21,3])))
CutL20LIGI=c(CycleLLIGI[22,2]:CycleLLIGI[22,3],CycleLLIGI[23,2]:CycleLLIGI[23,3],CycleLLIGI[24,2]:CycleLLIGI[24,3],CycleLLIGI[25,2]:CycleLLIGI[25,3],CycleLLIGI[26,2]:CycleLLIGI[26,3],CycleLLIGI[27,2]:CycleLLIGI[27,3],CycleLLIGI[28,2]:CycleLLIGI[28,3])+t3
CycleL20LIGI=c(rep("A",length(CycleLLIGI[22,2]:CycleLLIGI[22,3])),rep("B",length(CycleLLIGI[23,2]:CycleLLIGI[23,3])),rep("C",length(CycleLLIGI[24,2]:CycleLLIGI[24,3])),rep("D",length(CycleLLIGI[25,2]:CycleLLIGI[25,3])),rep("E",length(CycleLLIGI[26,2]:CycleLLIGI[26,3])),rep("F",length(CycleLLIGI[27,2]:CycleLLIGI[27,3])),rep("G",length(CycleLLIGI[28,2]:CycleLLIGI[28,3])))
CutL25LIGI=c(CycleLLIGI[29,2]:CycleLLIGI[29,3],CycleLLIGI[30,2]:CycleLLIGI[30,3],CycleLLIGI[31,2]:CycleLLIGI[31,3],CycleLLIGI[32,2]:CycleLLIGI[32,3],CycleLLIGI[33,2]:CycleLLIGI[33,3],CycleLLIGI[34,2]:CycleLLIGI[34,3],CycleLLIGI[35,2]:CycleLLIGI[35,3])+t4
CycleL25LIGI=c(rep("A",length(CycleLLIGI[29,2]:CycleLLIGI[29,3])),rep("B",length(CycleLLIGI[30,2]:CycleLLIGI[30,3])),rep("C",length(CycleLLIGI[31,2]:CycleLLIGI[31,3])),rep("D",length(CycleLLIGI[32,2]:CycleLLIGI[32,3])),rep("E",length(CycleLLIGI[33,2]:CycleLLIGI[33,3])),rep("F",length(CycleLLIGI[34,2]:CycleLLIGI[34,3])),rep("G",length(CycleLLIGI[35,2]:CycleLLIGI[35,3])))

CutLLIGI=TestLIGI[c(CutL5LIGI,CutL10LIGI,CutL15LIGI,CutL20LIGI,CutL25LIGI),]
CutLLIGI$Cycle=c(CycleL5LIGI,CycleL10LIGI,CycleL15LIGI,CycleL20LIGI,CycleL25LIGI)

CutG5LIGI=c(CycleGLIGI[1,2]:CycleGLIGI[1,3],CycleGLIGI[2,2]:CycleGLIGI[2,3],CycleGLIGI[3,2]:CycleGLIGI[3,3],CycleGLIGI[4,2]:CycleGLIGI[4,3],CycleGLIGI[5,2]:CycleGLIGI[5,3],CycleGLIGI[6,2]:CycleGLIGI[6,3],CycleGLIGI[7,2]:CycleGLIGI[7,3])+t0
CycleG5LIGI=c(rep("A",length(CycleGLIGI[1,2]:CycleGLIGI[1,3])),rep("B",length(CycleGLIGI[2,2]:CycleGLIGI[2,3])),rep("C",length(CycleGLIGI[3,2]:CycleGLIGI[3,3])),rep("D",length(CycleGLIGI[4,2]:CycleGLIGI[4,3])),rep("E",length(CycleGLIGI[5,2]:CycleGLIGI[5,3])),rep("F",length(CycleGLIGI[6,2]:CycleGLIGI[6,3])),rep("G",length(CycleGLIGI[7,2]:CycleGLIGI[7,3])))
CutG10LIGI=c(CycleGLIGI[8,2]:CycleGLIGI[8,3],CycleGLIGI[9,2]:CycleGLIGI[9,3],CycleGLIGI[10,2]:CycleGLIGI[10,3],CycleGLIGI[11,2]:CycleGLIGI[11,3],CycleGLIGI[12,2]:CycleGLIGI[12,3],CycleGLIGI[13,2]:CycleGLIGI[13,3],CycleGLIGI[14,2]:CycleGLIGI[14,3])+t1
CycleG10LIGI=c(rep("A",length(CycleGLIGI[8,2]:CycleGLIGI[8,3])),rep("B",length(CycleGLIGI[9,2]:CycleGLIGI[9,3])),rep("C",length(CycleGLIGI[10,2]:CycleGLIGI[10,3])),rep("D",length(CycleGLIGI[11,2]:CycleGLIGI[11,3])),rep("E",length(CycleGLIGI[12,2]:CycleGLIGI[12,3])),rep("F",length(CycleGLIGI[13,2]:CycleGLIGI[13,3])),rep("G",length(CycleGLIGI[14,2]:CycleGLIGI[14,3])))
CutG15LIGI=c(CycleGLIGI[15,2]:CycleGLIGI[15,3],CycleGLIGI[16,2]:CycleGLIGI[16,3],CycleGLIGI[17,2]:CycleGLIGI[17,3],CycleGLIGI[18,2]:CycleGLIGI[18,3],CycleGLIGI[19,2]:CycleGLIGI[19,3],CycleGLIGI[20,2]:CycleGLIGI[20,3],CycleGLIGI[21,2]:CycleGLIGI[21,3])+t2
CycleG15LIGI=c(rep("A",length(CycleGLIGI[15,2]:CycleGLIGI[15,3])),rep("B",length(CycleGLIGI[16,2]:CycleGLIGI[16,3])),rep("C",length(CycleGLIGI[17,2]:CycleGLIGI[17,3])),rep("D",length(CycleGLIGI[18,2]:CycleGLIGI[18,3])),rep("E",length(CycleGLIGI[19,2]:CycleGLIGI[19,3])),rep("F",length(CycleGLIGI[20,2]:CycleGLIGI[20,3])),rep("G",length(CycleGLIGI[21,2]:CycleGLIGI[21,3])))
CutG20LIGI=c(CycleGLIGI[22,2]:CycleGLIGI[22,3],CycleGLIGI[23,2]:CycleGLIGI[23,3],CycleGLIGI[24,2]:CycleGLIGI[24,3],CycleGLIGI[25,2]:CycleGLIGI[25,3],CycleGLIGI[26,2]:CycleGLIGI[26,3],CycleGLIGI[27,2]:CycleGLIGI[27,3],CycleGLIGI[28,2]:CycleGLIGI[28,3])+t3
CycleG20LIGI=c(rep("A",length(CycleGLIGI[22,2]:CycleGLIGI[22,3])),rep("B",length(CycleGLIGI[23,2]:CycleGLIGI[23,3])),rep("C",length(CycleGLIGI[24,2]:CycleGLIGI[24,3])),rep("D",length(CycleGLIGI[25,2]:CycleGLIGI[25,3])),rep("E",length(CycleGLIGI[26,2]:CycleGLIGI[26,3])),rep("F",length(CycleGLIGI[27,2]:CycleGLIGI[27,3])),rep("G",length(CycleGLIGI[28,2]:CycleGLIGI[28,3])))
CutG25LIGI=c(CycleGLIGI[29,2]:CycleGLIGI[29,3],CycleGLIGI[30,2]:CycleGLIGI[30,3],CycleGLIGI[31,2]:CycleGLIGI[31,3],CycleGLIGI[32,2]:CycleGLIGI[32,3],CycleGLIGI[33,2]:CycleGLIGI[33,3],CycleGLIGI[34,2]:CycleGLIGI[34,3],CycleGLIGI[35,2]:CycleGLIGI[35,3])+t4
CycleG25LIGI=c(rep("A",length(CycleGLIGI[29,2]:CycleGLIGI[29,3])),rep("B",length(CycleGLIGI[30,2]:CycleGLIGI[30,3])),rep("C",length(CycleGLIGI[31,2]:CycleGLIGI[31,3])),rep("D",length(CycleGLIGI[32,2]:CycleGLIGI[32,3])),rep("E",length(CycleGLIGI[33,2]:CycleGLIGI[33,3])),rep("F",length(CycleGLIGI[34,2]:CycleGLIGI[34,3])),rep("G",length(CycleGLIGI[35,2]:CycleGLIGI[35,3])))

CutGLIGI=TestLIGI[c(CutG5LIGI,CutG10LIGI,CutG15LIGI,CutG20LIGI,CutG25LIGI),]
CutGLIGI$Cycle=c(CycleG5LIGI,CycleG10LIGI,CycleG15LIGI,CycleG20LIGI,CycleG25LIGI)

# Calculate litter and Gammarus mean maximum biomasses and decreasing slopes
CutLLIGI=CutLLIGI[!CutLLIGI$Cycle=="A",]; CutGLIGI=CutGLIGI[!CutGLIGI$Cycle=="A",]
MinLLIGI=setDT(CutLLIGI)[, .SD[head(which(L<60000),1)], by=list(Temp,Cycle)]
MaxLLIGI=setDT(CutLLIGI)[, .SD[which.max(L)], by=list(Temp,Cycle)]
MinGLIGI=setDT(CutGLIGI)[, .SD[head(which(G<5000),1)], by=list(Temp,Cycle)]
MaxGLIGI=setDT(CutGLIGI)[, .SD[which.max(G)], by=list(Temp,Cycle)]

BiomLLIGI=rbind(MinLLIGI[,c(3,1,2,4)],MaxLLIGI[,c(3,1,2,4)])
BiomGLIGI=rbind(MinGLIGI[,c(3,1,2,5)],MaxGLIGI[,c(3,1,2,5)])

# Calculate means of cycles for the maximum biomasses and decreasing slopes
SlopLLIGI=as.data.frame(BiomLLIGI %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(L~Time, data=.))[[2]]))
SlopLLIGI=SlopLLIGI[order(as.numeric(SlopLLIGI$Temp)),]; SlopLLIGI$Slope=unlist(SlopLLIGI$Slope)
BiomLLIGI=cbind(SlopLLIGI,Max=BiomLLIGI[BiomLLIGI$L>100000]$L)
MaxiLLIGI=setDT(BiomLLIGI)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniLLIGI=setDT(BiomLLIGI)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]

SlopGLIGI=as.data.frame(BiomGLIGI %>% group_by(Cycle,Temp) %>% do(Slope=coef(lm(G~Time, data=.))[[2]]))
SlopGLIGI=SlopGLIGI[order(as.numeric(SlopGLIGI$Temp)),]; SlopGLIGI$Slope=unlist(SlopGLIGI$Slope)
BiomGLIGI=cbind(SlopGLIGI,Max=BiomGLIGI[BiomGLIGI$G>10000]$G)
MaxiGLIGI=setDT(BiomGLIGI)[, .(MeanMax=mean(Max), SdMax=sd(Max)), by=list(Temp)]
MiniGLIGI=setDT(BiomGLIGI)[, .(MeanSlope=mean(Slope), SdSlope=sd(Slope)), by=list(Temp)]


#########################
### Create dataframes ###
#########################

DataL53=data.frame(MeanLLIGI[,1],MeanLLIGI[,2],MeanLLIGI[,3],TimeLLIGI[,2],TimeLLIGI[,3])
colnames(DataL53)=c("Temperature","MeanBiomL","SdBiomL","MeanTimeL","SdTimeL")
DataL53$Scenario=rep("LIGI",5)
DataL53$Subscenario=rep("GI",5)

DataG53=data.frame(MeanGLIGI[,1],MeanGLIGI[,2],MeanGLIGI[,3],TimeGLIGI[,2],TimeGLIGI[,3])
colnames(DataG53)=c("Temperature","MeanBiomG","SdBiomG","MeanTimeG","SdTimeG")
DataG53$Scenario=rep("LIGI",5)
DataG53$Subscenario=rep("GI",5)




####################################################################################################




############################################################
### LEAF LITTER AND GAMMARUS BIOMASSES SCENARIOS OUTPUTS ###
############################################################

# Bind dataframes of all scenarios
DataLP=rbind(DataL31,DataL32,DataL33,DataL41,DataL42,DataL43,DataL51,DataL52,DataL53)
DataGP=rbind(DataG31,DataG32,DataG33,DataG41,DataG42,DataG43,DataG51,DataG52,DataG53)
DataLP$Temperature=as.numeric(as.character(DataLP$Temperature))
DataGP$Temperature=as.numeric(as.character(DataGP$Temperature))

# Create a predicted dataframe
Temperature=expand.grid(Temperature=seq(5,25,by=0.5))
Scenario=c(rep("LRGR",41),rep("LRGC",41),rep("LRGI",41),rep("LCGR",41),rep("LCGC",41),rep("LCGI",41),rep("LIGR",41),rep("LIGC",41),rep("LIGI",41))
Subscenario=c(rep("GR",41),rep("GC",41),rep("GI",41),rep("GR",41),rep("GC",41),rep("GI",41),rep("GR",41),rep("GC",41),rep("GI",41))
PredL=data.frame(Temperature,Scenario,Subscenario)
PredG=data.frame(Temperature,Scenario,Subscenario)

# Predict on short temperature intervals
PredL$PredBiomL=data.frame(predict(lmList(MeanBiomL ~ poly(Temperature,2)|Scenario, data=DataLP),PredL))[[1]]
PredL$PredTimeL=data.frame(predict(lmList(MeanTimeL ~ poly(Temperature,2)|Scenario, data=DataLP),PredL))[[1]]

PredG$PredBiomG=data.frame(predict(lmList(MeanBiomG ~ poly(Temperature,2)|Scenario, data=DataGP),PredG))[[1]]
PredG$PredTimeG=data.frame(predict(lmList(MeanTimeG ~ poly(Temperature,2)|Scenario, data=DataGP),PredG))[[1]]


##########################
### Plotting biomasses ###
##########################

# Plot leaf litter mean annual biomasses
PlotBiomL=ggplot(DataLP, aes(x=Temperature, y=MeanBiomL/10^5, group=interaction(Scenario,Subscenario))) +
  geom_point(aes(color=Scenario), position=position_dodge(0.5), size=4, pch=16) +
  geom_line(data=PredL, aes(x=Temperature, y=PredBiomL/10^5, color=Scenario, linetype=Scenario), size=1.0, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=(MeanBiomL-SdBiomL)/10^5, ymax=(MeanBiomL+SdBiomL)/10^5), color="grey50", size=1.0, width=1.0, position=position_dodge(0.5)) +
  ylab(expression('Mean leaf litter biomass'~'('*10^5~mg~C~m^-2*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,2.0,by=0.5), limits=c(0,2.0)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,25,by=5), limits=c(4,26)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  scale_color_manual(values=c(LCGC="black", LCGR="black", LCGI="black", LRGC="tomato2", LRGR="tomato2", LRGI="tomato2", LIGC="steelblue2", LIGR="steelblue2", LIGI="steelblue2")) +
  scale_linetype_manual(values=c(LCGC="dotted", LCGR="dotted", LCGI="dotted", LRGC="solid", LRGR="solid", LRGI="solid", LIGC="solid", LIGR="solid", LIGI="solid")) +
  theme(strip.text.x=element_text(face="plain", color="black", size=20)) +
  theme(strip.background=element_blank()) + theme(axis.line=element_line()) +
  facet_wrap(~Subscenario, scales="free", ncol=1, nrow=3) +
  theme(panel.spacing=unit(2,"lines")) +
  theme(legend.position="none")

# Plot Gammarus mean annual biomasses
PlotBiomG=ggplot(DataGP, aes(x=Temperature, y=MeanBiomG/10^4, group=interaction(Scenario,Subscenario))) +
  geom_point(aes(color=Scenario), position=position_dodge(0.5), size=4, pch=16) +
  geom_line(data=PredG, aes(x=Temperature, y=PredBiomG/10^4, color=Scenario, linetype=Scenario), size=1.0, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=(MeanBiomG-SdBiomG)/10^4, ymax=(MeanBiomG+SdBiomG)/10^4), color="grey50", size=1.0, width=1.0, position=position_dodge(0.5)) +
  ylab(expression('Mean'*~italic(Gammarus)*~'biomass'~'('*10^4~mg~C~m^-2*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,25,by=5), limits=c(4,26)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  scale_color_manual(values=c(LCGC="black", LCGR="black", LCGI="black", LRGC="tomato2", LRGR="tomato2", LRGI="tomato2", LIGC="steelblue2", LIGR="steelblue2", LIGI="steelblue2")) +
  scale_linetype_manual(values=c(LCGC="dotted", LCGR="dotted", LCGI="dotted", LRGC="solid", LRGR="solid", LRGI="solid", LIGC="solid", LIGR="solid", LIGI="solid")) +
  theme(strip.text.x=element_text(face="plain", color="black", size=20)) +
  theme(strip.background=element_blank()) + theme(axis.line=element_line()) +
  facet_wrap(~Subscenario, scales="free", ncol=1, nrow=3) +
  theme(panel.spacing=unit(2,"lines")) +
  theme(legend.position="none")

tiff('Mean Biomass Appendix.tiff', units="in", width=15, height=15, res=300)
Panel=plot_grid(PlotBiomL, PlotBiomG, align="h", vjust=1, nrow=1, ncol=2)
Xaxis=textGrob(expression('Temperature (°C)'), gp=gpar(fontface="bold", fontsize=20))
grid.arrange(arrangeGrob(Panel, bottom=Xaxis))
dev.off()


##################################
### Plotting persistence times ###
##################################

# Plot leaf litter mean annual persistence time
PlotTimeL=ggplot(DataLP, aes(x=Temperature, y=MeanTimeL, group=interaction(Scenario,Subscenario))) +
  geom_point(aes(color=Scenario), position=position_dodge(0.5), size=4, pch=16) +
  geom_line(data=PredL, aes(x=Temperature, y=PredTimeL, color=Scenario, linetype=Scenario), size=1.0, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=MeanTimeL-SdTimeL, ymax=MeanTimeL+SdTimeL), color="grey50", size=1.0, width=1.0, position=position_dodge(0.5)) +
  ylab(expression('Leaf litter persistence time (d)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +  
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,100,by=25), limits=c(0,100)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,25,by=5), limits=c(4,26)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  scale_color_manual(values=c(LCGC="black", LCGR="black", LCGI="black", LRGC="tomato2", LRGR="tomato2", LRGI="tomato2", LIGC="steelblue2", LIGR="steelblue2", LIGI="steelblue2")) +
  scale_linetype_manual(values=c(LCGC="dotted", LCGR="dotted", LCGI="dotted", LRGC="solid", LRGR="solid", LRGI="solid", LIGC="solid", LIGR="solid", LIGI="solid")) +
  theme(strip.text.x=element_text(face="plain", color="black", size=20)) +
  theme(strip.background=element_blank()) + theme(axis.line=element_line()) +
  facet_wrap(~Subscenario, scales="free", ncol=1, nrow=3) +
  theme(panel.spacing=unit(2,"lines")) +
  theme(legend.position="none")

# Plot Gammarus mean annual persistence time
PlotTimeG=ggplot(DataGP, aes(x=Temperature, y=MeanTimeG, group=interaction(Scenario,Subscenario))) +
  geom_point(aes(color=Scenario), position=position_dodge(0.5), size=4, pch=16) +
  geom_line(data=PredG, aes(x=Temperature, y=PredTimeG, color=Scenario, linetype=Scenario), size=1.0, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=MeanTimeG-SdTimeG, ymax=MeanTimeG+SdTimeG), color="grey50", size=1.0, width=1.0, position=position_dodge(0.5)) +
  ylab(expression(italic(Gammarus)*~'persistence time (d)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=20)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,400,by=100), limits=c(0,400)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,25,by=5), limits=c(4,26)) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +
  scale_color_manual(values=c(LCGC="black", LCGR="black", LCGI="black", LRGC="tomato2", LRGR="tomato2", LRGI="tomato2", LIGC="steelblue2", LIGR="steelblue2", LIGI="steelblue2")) +
  scale_linetype_manual(values=c(LCGC="dotted", LCGR="dotted", LCGI="dotted", LRGC="solid", LRGR="solid", LRGI="solid", LIGC="solid", LIGR="solid", LIGI="solid")) +
  theme(strip.text.x=element_text(face="plain", color="black", size=20)) +
  theme(strip.background=element_blank()) + theme(axis.line=element_line()) +
  facet_wrap(~Subscenario, scales="free", ncol=1, nrow=3) +
  theme(panel.spacing=unit(2,"lines")) +
  theme(legend.position="none")

tiff('Persistence Time Appendix.tiff', units="in", width=15, height=15, res=300)
Panel=plot_grid(PlotTimeL, PlotTimeG, align="h", vjust=1, nrow=1, ncol=2)
Xaxis=textGrob(expression('Temperature (°C)'), gp=gpar(fontface="bold", fontsize=20))
grid.arrange(arrangeGrob(Panel, bottom=Xaxis))
dev.off()
