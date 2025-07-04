rm(list=ls())
memory.limit(999999)
load("immune_res_imputation.RData")
y=readRDS(file="immune_data.RDS")

# compute MSE and predictive distributions
  
MSE=rep(0,10000)
for(iter in iter:10000) MSE[iter]=  sum((ris_y[iter,,]-y)^2)/30
plot(MSE, type="l")

#boxplot((MSE[1:8000]))
#boxplot((MSE[10000:8000]))
 
# index of out of sample units
idx_imp

# to obtain Figure 8 of the paper
pred=matrix(0, 2000,30)
for(iter in 8001:10000){
  for(ii in 1:30) pred[iter-8000,ii]=(ris_y[iter,idx_imp[ii,1],idx_imp[ii,2]]) 
}
yobs=matrix(0,1,30)
 
for(ii in 1:30) yobs[1,ii]=(y[idx_imp[ii,1],idx_imp[ii,2]]) 

# to obtain Figure 8 of the paper
par(mfrow=c(1,2))


plot(MSE, type="l", xlab="iteration")
mean(MSE)
boxplot((pred), ylab="prediction",ylim=c(-35,45))

points(1:30,yobs, col=2,pch=8,lwd=0.5,cex=01.8)

 
