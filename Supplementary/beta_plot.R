# supplementary materials FIGURES SM11 SM12 SM13: additional figures from APAFA
setwd("APAFA/Simulations/long/SCEN1")

load("scen1_1.RData")
par(mfrow=c(5,3))
par(mar=c(2,4,1,1))
matplot(ris_beta1[,], type="l",ylab=expression(beta[1]),xlab="iteration", lty=1:3)
matplot(ris_beta2[,], type="l",ylab=expression(beta[2]),xlab="iteration", lty=1:3)
matplot(ris_beta3[,], type="l",ylab=expression(beta[3]),xlab="iteration", lty=1:3)

load("scen1_2.RData")
 
matplot(ris_beta1[,], type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(ris_beta2[,], type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(ris_beta3[,], type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)


load("scen1_4.RData")

matplot(ris_beta1[,], type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(ris_beta2[,], type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(ris_beta3[,], type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)

par(mfrow=c(5,3))
load("scen1_5.RData")
 
par(mar=c(2,4,1,1))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)

load("scen1_6.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)

load("scen1_7.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)

load("scen1_8.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)

load("scen1_9.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)
 

 





#### 
# additional figures from APAFA
setwd("APAFA/Simulations/long/SCEN2")

load("scen2_5.RData")
par(mfrow=c(5,3))
par(mar=c(2,4,1,1))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)

load("scen2_6.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)

load("scen2_7.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)

load("scen2_8.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)

load("scen2_9.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1} *beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1} *beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1} *beta[3]),xlab="iteration", lty=1:3)

load("scen2_10.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(ris_beta1[,], type="l",ylab=expression(beta[1]),xlab="iteration", lty=1:3)
matplot(ris_beta2[,], type="l",ylab=expression(beta[2]),xlab="iteration", lty=1:3)
matplot(ris_beta3[,], type="l",ylab=expression(beta[3]),xlab="iteration", lty=1:3)



#####
# additional figures from APAFA
setwd("APAFA/Simulations/long/SCEN3")

load("scen3_5.RData")
par(mfrow=c(5,3))
par(mar=c(2,4,1,1))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen3_6.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen3_7.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen3_8.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen3_9.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen3_10.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)


# additional figures from APAFA
setwd("APAFA/Simulations/long/SCEN4a")

load("scen4a_5.RData")
par(mfrow=c(5,3))
par(mar=c(2,4,1,1))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen4a_6.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen4a_7.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen4a_8.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen4a_9.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen4a_10.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(ris_beta1[,], type="l",ylab=expression(beta[1]),xlab="iteration", lty=1:3)
matplot(ris_beta2[,], type="l",ylab=expression(beta[2]),xlab="iteration", lty=1:3)
matplot(ris_beta3[,], type="l",ylab=expression(beta[3]),xlab="iteration", lty=1:3)


# additional figures from APAFA
setwd("APAFA/Simulations/long/SCEN4B")

load("scen4b_5.RData")
par(mfrow=c(5,3))
par(mar=c(2,4,1,1))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen4b_6.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen4b_7.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen4b_8.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)

load("scen4b_9.RData")
#par(mfrow=c(3,1))
#par(mar=c(3,4,3,3))
matplot(plogis(ris_beta1[,]), type="l",ylab=expression(logit^{-1}*beta[1]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta2[,]), type="l",ylab=expression(logit^{-1}*beta[2]),xlab="iteration", lty=1:3)
matplot(plogis(ris_beta3[,]), type="l",ylab=expression(logit^{-1}*beta[3]),xlab="iteration", lty=1:3)
 
