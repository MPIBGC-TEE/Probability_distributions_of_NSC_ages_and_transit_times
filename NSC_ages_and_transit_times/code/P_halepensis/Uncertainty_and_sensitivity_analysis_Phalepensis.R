##### Uncertainty and Sensitivity analysis for the NSC mean age and mean transti time of Pinus halepensis

## load the necesary packages 
library(SoilR)
library(sensitivity)
library(boot)
library(FME) 
library(ggplot2)
setwd("/Users/_dherrera/NSC_ages_and_transit_times/code/P_halepensis/")

source("/Users/_dherrera/NSC_ages_and_transit_times/code/P_halepensis/Functions_Phalepensis.R")

## we load the files we need for the computations in this script. This parameters and their
## upper and lower limits where calculates in the script "Klain_2015_annual_parameter_estimation.R"

load("pars1.Rda")
load("LWLImeanannualrates4.Rda")
load("UPLImeanannualrates4.Rda")



#######################################################################################
############# Sensitivity Analysis#######################################################
############################################################################################

### sensibility analysis was run based on the Elementary Effects technique (Morris 1992). This technique
### is implemented in the sensitivity package. 

fluxnames= c("Rf","Rs","Rr","Gf","Gs","Gr",  
             "Lf","Ls","Lr","FtoS","StoF","Stor",
             "rtoS","Sf","Ss","Sr","Cf","Cs","Cr")

sensmeanage1 <- morris(model = systage, factors = fluxnames, r = 50, binf = LWLImeanannualrates4, 
                       bsup = UPLImeanannualrates4, design = list(type = "oat", levels = 100, grid.jump = 50))

print(sensmeanage1)


mu <- apply(sensmeanage1$ee, 3, function(M){
  apply(M, 2, mean)
})

mu.star <- apply(abs(sensmeanage1$ee), 3, function(M){
  apply(M, 2, mean)
})
sigma <- apply(sensmeanage1$ee, 3, function(M){
  apply(M, 2, sd)
})

png("phalepinsis_Mean_Age.png", width=400, height =400)
par(mar=c(5, 5,4,4))
plot(mu.star[,1], sigma[,1], main= "Mean Age",
     xlab=expression(paste(mu^"*")), ylab = expression(sigma),
     ylim=c(0, 0.3),xlim=c(0, 0.4),
     cex=2, cex.axis=2, cex.lab=2, cex.main=3, pch=18 )
text(mu.star[,1], sigma[,1], row.names(mu),pos=3, cex=1.3)
dev.off()


png("phalepensis_Mean_transit_time.png", width=400, height =400)
par(mar=c(5, 5,4,4))
plot(mu.star[,2], sigma[,2], main= "Mean Transit Time",
     xlab=expression(paste(mu^"*")), ylab = expression(sigma),
     ylim=c(0, 0.05),xlim=c(0, max(mu.star[,2])+0.01),
     cex=2, cex.axis=2, cex.lab=2, cex.main=3, pch=18 )
text(mu.star[,2], sigma[,2], row.names(mu),pos=3, cex=1.3)
dev.off()


### we plot the trend of the sensitivity of each of the most 
### influential parameters 


x=as.vector(unlist(sensmeanage1$X[,"Cs"]))

df=data.frame(x=x,
              y=sensmeanage1$y[,1])


m3=nls(y~I(a+b*x^-z), data=df,  start=list(a=0.04,b=0.01, z=1), control = list(maxiter=10000), trace = T)
summary(m3)
y=coef(m3)[1]+coef(m3)[2]*sort(x)^-coef(m3)[3]


png("Phalepensis_MeanAge_Vs_Cs", width = 350, height = 400)
par(mar=c(5,5,2,2))
plot(sensmeanage1$X[,"Cs"], sensmeanage1$y[,1], ylab="Man Age",
     xlab= "Cs", cex.axis=2, cex.lab=2, pch=18)
lines(sort(x),y, col="red", lwd=2)
dev.off()

png("Phalepensis_MeanAge_Vs_Sr", width = 350, height = 400)
par(mar=c(5,5,2,2))
plot(sensmeanage1$X[,"Sr"], sensmeanage1$y[,1], ylab = "Mean Age",
     xlab = "Sr",  cex.axis=2, cex.lab=2, pch=18)
abline(lm(sensmeanage1$y[,1] ~ sensmeanage1$X[,"Sr"]), col="red", lwd=2)
dev.off()


png("Phalepensis_Meantransittime_Vs_Cs", width = 350, height = 400)
par(mar=c(5,5,2,2))
plot(sensmeanage1$X[,"Cs"], sensmeanage1$y[,2], ylab = "Mean Transit Time",
     xlab = "Cs",  cex.axis=2, cex.lab=2, pch=18)
abline(lm(sensmeanage1$y[,2] ~ sensmeanage1$X[,"Cs"]), col="red", lwd=2)
dev.off()

png("Phalepensis_Meantransittime_Vs_Sr", width = 350, height = 400)
par(mar=c(5,5,2,2))
plot(sensmeanage1$X[,"Sr"], sensmeanage1$y[,2], ylab = "Mean Transit Time",
     xlab = "Sr",  cex.axis=2, cex.lab=2, pch=18)
abline(lm(sensmeanage1$y[,2] ~ sensmeanage1$X[,"Sr"]), col="red", lwd=2)
dev.off()


################################################################################################
############# Uncertainty Analysis ##################################
################################################################################################
################################################################################################

#We defined the fluxes names as they appear in the manuscript
fluxnames= c("Rf","Rs","Rr","Gf","Gs","Gr",  
             "Lf","Ls","Lr","FtoS","StoF","Stor",
             "rtoS","Sf","Ss","Sr","Cf","Cs","Cr")

# We calculated the variability limits for the parameters reported by Klein and Hoch 2015
pars_range=(pars1- LWLImeanannualrates4)

## the parameters should be trasnformed into a list for the MCS function from FME package 
pars1_l=list()
for(i in 1:length(pars1)){
pars1_l[i]=pars1[i]  
}
names(pars1_l)=fluxnames

pars1_mean=as.numeric(pars1)
names(pars1_mean)=fluxnames


systage_MCS(pars1_l)## test for the mean ages and transit time of the mean parameter values
lw=as.numeric(pars1-pars_range)
lw[lw<0]=0
up=as.numeric(pars1 + pars_range)
parRanges=cbind(lw, up)

rownames(parRanges) = fluxnames

## The Monte Carlo simulation was run to propagete the uncertainty of the model parametes in the 
## mean age and mean transti time of P. halepensis. This analysis is run just for the 4 parameters
## that showed the biggest influence in the model output. 

## This function calculate 1000 model realization with different combination of the model parameters
## selected. It takes a very long time to run and sometimes fail when some negative values for the 
## parameters are ramdomly chosen.

# 
# CRL4 <- modCRL(func = systage_MCS, parms = pars1_l, parMean = pars1_mean[c(3,12,16,18)],
#                parCovar =diag(variance_pars[c(3,12,16,18)]*0.7),
#                dist = "norm", num = 1000)
# 
# uncertianties_meanarebypool=CRL4
# 
# save(uncertianties_meanarebypool, file = "MCS_meanagebypool.Rda")

# the results of the function above can be load as. 

load("MCS_meanagebypool.Rda")
sd_Phalep_meanage=apply(uncertianties_meanarebypool,2,sd)
mean_Phalep_meanage=apply(uncertianties_meanarebypool,2,mean)
rbind(mean_Phalep_meanage, sd_Phalep_meanage)


CRL4=uncertianties_meanarebypool
png("phalepensis_meanSysAgebypool_uncertianty.png",
    width = 800,
    height = 700)
par(mfrow=c(2,4), mar=c(5,5,2,2))
hist(CRL4$meansysage, main = "", xlab = "Mean system age", breaks = 50, cex.lab=2, cex.axis=2)
hist(CRL4$meanTT, main = "", xlab = "Mean transit time", breaks = 50, cex.lab=2, cex.axis=2)
hist(CRL4$FANSC, main = "", xlab = "FANSC", breaks = 50, cex.lab=2, cex.axis=2)
hist(CRL4$FSNSC, main = "", xlab = "FSNSC", breaks = 50, cex.lab=2, cex.axis=2)
hist(CRL4$SANSC, main = "", xlab = "SANSC", breaks = 50, cex.lab=2, cex.axis=2)
hist(CRL4$SSNSC, main = "", xlab = "SSNSC", breaks = 50, cex.lab=2, cex.axis=2)
hist(CRL4$RANSC, main = "", xlab = "RANSC", breaks = 50, cex.lab=2, cex.axis=2)
hist(CRL4$RSNSC, main = "", xlab = "RSNSC", breaks = 50, cex.lab=2, cex.axis=2)
dev.off()

