rm(list=ls())
memory.limit(999999)
# for APAFA
load("immune_res_imputation.RData")
# for TETRIS
#load("immune_tetris_mse.RData")
y=readRDS(file="immune_data.RDS")
iter
  
MSE=rep(0,10000)
for(iter in iter:10000) MSE[iter]=  sum((ris_y[iter,,]-y)^2)/30
plot(MSE, type="l")

boxplot((MSE[1:8000]))
boxplot((MSE[10000:8000]))
 

diff=matrix(0, 2000,30)
for(iter in 8001:10000){
for(ii in 1:30) diff[iter-8000,ii]=(y[idx_imp[ii,1],idx_imp[ii,2]]-ris_y[iter,idx_imp[ii,1],idx_imp[ii,2]]) 
}
boxplot((diff))





pred=matrix(0, 2000,30)
for(iter in 8001:10000){
  for(ii in 1:30) pred[iter-8000,ii]=(ris_y[iter,idx_imp[ii,1],idx_imp[ii,2]]) 
}
yobs=matrix(0,1,30)
 
  for(ii in 1:30) yobs[1,ii]=(y[idx_imp[ii,1],idx_imp[ii,2]]) 

# figure 9
par(mfrow=c(1,2))

plot(MSE, type="l", xlab="iteration")
mean(MSE)
boxplot((pred), ylab="prediction")
points(1:30,yobs, col=2,pch=8,lwd=0.5,cex=01.8)
