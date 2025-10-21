#########################21.12 5528##############################################################
#LOAD DATA
memory.limit(100000)
load("immune_res.RData")


# plot covariance

colQuant90=function (m) apply(m,c(2), function(x) quantile(x,0.9))
colQuant10=function (m)apply(m,c(2),function(x) quantile(x,0.1))
colQuant50=function (m)apply(m,c(2),function(x) quantile(x,0.1))


L=colMeans(colQuant10(ris_tau_eta[5000:10000, ])*(ris_lambda[5000:10000,, ]))

G=colMeans(colQuant10(ris_tau_phi[5000:10000, ])*(ris_gamma[5000:10000,, ]))

#reorder specific factors according to activation %
GG1=tcrossprod(G[,1],G[,1])
GG2=tcrossprod(G[,2],G[,2])
GG3=tcrossprod(G[,3],G[,3])
GG4=tcrossprod(G[,4],G[,4])
GG5=tcrossprod(G[,5],G[,5])
GG6=tcrossprod(G[,6],G[,6])
GG7=tcrossprod(G[,7],G[,7])
GG8=tcrossprod(G[,8],G[,8])
GG9=tcrossprod(G[,9],G[,9])

#### obtain better colors
par(mfrow=c(3,3))
GG1[1,1]=max(abs(GG1))
GG1[p,p]=-max(abs(GG1))
GG2[1,1]=max(abs(GG2))
GG2[p,p]=-max(abs(GG2))
GG3[1,1]=max(abs(GG3))
GG3[p,p]=-max(abs(GG3))
GG4[1,1]=max(abs(GG4))
GG4[p,p]=-max(abs(GG4))
GG5[1,1]=max(abs(GG4))
GG5[p,p]=-max(abs(GG5))
GG6[1,1]=max(abs(GG6))
GG6[p,p]=-max(abs(GG6))
GG7[1,1]=max(abs(GG7))
GG7[p,p]=-max(abs(GG7))
GG8[1,1]=max(abs(GG8))
GG8[p,p]=-max(abs(GG8))
GG9[1,1]=max(abs(GG9))
GG9[p,p]=-max(abs(GG9))


par(mar=c(4,3,3,3))
# Define the matrix to be visualized 
# Define the color palette with red and blue
custom_colors <- c("red","white", "blue")

# Create a custom color function that maps values to colors
color_function <- colorRampPalette(custom_colors)
colors <- color_function(0.3*length(GG1))
# Plot the image with the custom color palette

# Figure 9: contribution of Gamma Gamma
par(mfrow=c(3,3))
image(GG3, axes=F,xlab=expression(Gamma[1]*Gamma[1]^T),col = colors[sort.list(t(GG3))])
image(GG4, axes=F,xlab=expression(Gamma[2]*Gamma[2]^T),col = colors[sort.list(GG4)])
image(GG5, axes=F,xlab=expression(Gamma[3]*Gamma[3]^T),col = colors[sort.list(GG1)])
image(GG1, axes=F,xlab=expression(Gamma[4]*Gamma[4]^T),col = colors[sort.list(GG2)])
image(GG2, axes=F,xlab=expression(Gamma[5]*Gamma[5]^T),col = colors[sort.list(GG5)])
image(GG6, axes=F,xlab=expression(Gamma[6]*Gamma[6]^T),col = colors[sort.list(GG6)])
image(GG7, axes=F,xlab=expression(Gamma[7]*Gamma[7]^T),col = colors[sort.list(GG7)])
image(GG8, axes=F,xlab=expression(Gamma[8]*Gamma[8]^T),col = colors[sort.list(GG8)])
image(GG9, axes=F,xlab=expression(Gamma[9]*Gamma[9]^T),col = colors[sort.list(GG9)])


dim(ris_gamma)

#### QUANTILE 0.1

par(mar=c(3.5,1.,1.,1.))
par(mfrow=c(2,9))

Giter1=sapply(8001:10000, function (j) 
  tcrossprod((ris_tau_phi[j, 1])*(ris_gamma[j,,1 ])))
Giter2=sapply(8001:10000, function (j) 
  tcrossprod((ris_tau_phi[j, 2])*(ris_gamma[j,,2 ])))
Giter3=sapply(8001:10000, function (j) 
  tcrossprod((ris_tau_phi[j, 3])*(ris_gamma[j,,3 ])))
Giter4=sapply(8001:10000, function (j) 
  tcrossprod((ris_tau_phi[j, 4])*(ris_gamma[j,,4 ])))
Giter5=sapply(8001:10000, function (j) 
  tcrossprod((ris_tau_phi[j, 5])*(ris_gamma[j,,5 ])))
Giter6=sapply(8001:10000, function (j) 
  tcrossprod((ris_tau_phi[j, 6])*(ris_gamma[j,,6 ])))
Giter7=sapply(8001:10000, function (j) 
  tcrossprod((ris_tau_phi[j, 7])*(ris_gamma[j,,7 ])))
Giter8=sapply(8001:10000, function (j) 
  tcrossprod((ris_tau_phi[j, 8])*(ris_gamma[j,,8 ])))
Giter9=sapply(8001:10000, function (j) 
  tcrossprod((ris_tau_phi[j, 9])*(ris_gamma[j,,9 ])))


MGG3=matrix(apply(Giter3,1, function(j) quantile(j,0.1)), ncol=63)
#MIN3=GG3[which.min(c(GG3))]=-max(c(GG3))
image(GG3, axes=F, col=colors,ylab="q=0.1",
      main=expression(Gamma[3]*Gamma[3]^T))
MGG4=matrix(apply(Giter4,1, function(j) quantile(j,0.1)), ncol=63)
#MIN4=GG4[which.min(c(GG4))]=-max(c(GG4))
image(GG4, axes=F, col=colors,
      main=expression(Gamma[4]*Gamma[4]^T))
MGG5=matrix(apply(Giter5,1, function(j) quantile(j,0.1)), ncol=63)
#MIN5=GG5[which.min(c(GG5))]=-max(c(GG5))
image(GG5, axes=F, col=colors,
      main=expression(Gamma[5]*Gamma[5]^T))

MGG1=matrix(apply(Giter1,1, function(j) quantile(j,0.1)), ncol=63)
#MIN1=GG1[which.min(c(GG1))]=-max(c(GG1))
image(GG1, axes=F, col=colors, 
      main=expression(Gamma[1]*Gamma[1]^T))
MGG2=matrix(apply(Giter2,1, function(j) quantile(j,0.1)), ncol=63)
#MIN2=GG2[which.min(c(GG2))]=-max(c(GG2))
image(GG2, axes=F, col=colors,
      main=expression(Gamma[2]*Gamma[2]^T))
MGG6=matrix(apply(Giter6,1, function(j) quantile(j,0.1)), ncol=63)
#MIN6=GG6[which.min(c(GG6))]=-max(c(GG6))
image(GG6, axes=F, col=colors,
      main=expression(Gamma[6]*Gamma[6]^T))
MGG7=matrix(apply(Giter7,1, function(j) quantile(j,0.1)), ncol=63)
#MIN7=GG7[which.min(c(GG7))]=-max(c(GG7))
image(GG7, axes=F, col=colors,
      main=expression(Gamma[7]*Gamma[7]^T))
MGG8=matrix(apply(Giter8,1, function(j) quantile(j,0.1)), ncol=63)
##MIN8=GG8[which.min(c(GG8))]=-max(c(GG8))
image(GG8, axes=F, col=colors,
      main=expression(Gamma[8]*Gamma[8]^T))
MGG9=matrix(apply(Giter9,1, function(j) quantile(j,0.1)), ncol=63)
#MIN9=GG9[which.min(c(GG9))]=-max(c(GG9))
image(GG9, axes=F, col=colors,
      main=expression(Gamma[9]*Gamma[9]^T))


# quantile 0.9
GG3=matrix(apply(Giter3,1, function(j) quantile(j,0.9)), ncol=63)

image(GG3, axes=F, col=colors,
      main=expression(Gamma[1]*Gamma[1]^T))
GG4=matrix(apply(Giter4,1, function(j) quantile(j,0.9)), ncol=63)
#GG4[which.min(c(GG4))]=MIN4


image(GG4, axes=F, col=colors,
      main=expression(Gamma[2]*Gamma[2]^T))
GG5=matrix(apply(Giter5,1, function(j) quantile(j,0.9)), ncol=63)
#GG5[which.min(c(GG5))]=-max(c(GG5))
image(GG5, axes=F, col=colors,
      main=expression(Gamma[3]*Gamma[3]^T))

GG1=matrix(apply(Giter1,1, function(j) quantile(j,0.9)), ncol=63)
#GG1[which.min(c(GG1))]=-max(c(GG1))
image(GG1, axes=F, col=colors, ylab="q=0.9",
      main=expression(Gamma[4]*Gamma[4]^T))
GG2=matrix(apply(Giter2,1, function(j) quantile(j,0.9)), ncol=63)
#GG2[which.min(c(GG2))]=-max(c(GG2))
image(GG2, axes=F, col=colors,
      main=expression(Gamma[5]*Gamma[5]^T))
GG6=matrix(apply(Giter6,1, function(j) quantile(j,0.9)), ncol=63)
#GG6[which.min(c(GG6))]=-max(c(GG6))
image(GG6, axes=F, col=colors,
      main=expression(Gamma[6]*Gamma[6]^T))
GG7=matrix(apply(Giter7,1, function(j) quantile(j,0.9)), ncol=63)
#GG7[which.min(c(GG7))]=-max(c(GG7))
image(GG7, axes=F, col=colors,
      main=expression(Gamma[7]*Gamma[7]^T))
GG8=matrix(apply(Giter8,1, function(j) quantile(j,0.9)), ncol=63)
#GG8[which.min(c(GG8))]=-max(c(GG8))
image(GG8, axes=F, col=colors,
      main=expression(Gamma[8]*Gamma[8]^T))
GG9=matrix(apply(Giter9,1, function(j) quantile(j,0.9)), ncol=63)
#GG9[which.min(c(GG9))]=-max(c(GG9))
image(GG9, axes=F, col=colors,
      main=expression(Gamma[9]*Gamma[9]^T))


library(rgl)

lower_bound <- MGG3
upper_bound <- GG3
zmean <- ifelse(lower_bound>0, "red", "white")
zmean <- ifelse(upper_bound<0, "blue", zmean ) 
open3d()
planes3d(0, 0, 1, 0, col = "gray", alpha = 0.25)
for (i in 1:63) {
  for (j in 1:63) {
    z1 <- lower_bound[i, j]
    z2 <- upper_bound[i, j]
    segments3d(rbind(
      c(i / 63, j / 63, z1 / 1),
      c(i / 63, j / 63, z2 / 1)
    ), col = zmean[i,j], lwd=2)
  }
}
legend3d("topright", legend =c("l,u<0",  "l,u>0"), lty=1,lwd=2, cex=1,
         col =c("blue",   "red"))
title3d(main = expression(Gamma[1]*Gamma[1]^T), cex=1)
snapshot3d(filename = 'GAMMA1.png', fmt = 'png')



library(rgl)
x_coords=1:63
y_coords=x_coords

lower_bound <- MGG4
upper_bound <- GG4
zmean <- ifelse(lower_bound>0, "red", "white")
zmean <- ifelse(upper_bound<0, "blue", zmean ) 
open3d()
planes3d(0, 0, 1, 0, col = "gray", alpha = 0.25)
for (i in 1:63) {
  for (j in 1:63) {
    z1 <- lower_bound[i, j]
    z2 <- upper_bound[i, j]
    segments3d(rbind(
      c(i / 63, j / 63, z1 / 1),
      c(i / 63, j / 63, z2 / 1)
    ), col = zmean[i,j], lwd=2)
  }
}
legend3d("topright", legend =c("l,u<0",  "l,u>0"), lty=1,lwd=2, cex=1,
         col =c("blue",   "red"))
title3d(main = expression(Gamma[2]*Gamma[2]^T), cex=1)
snapshot3d(filename = 'GAMMA2.png', fmt = 'png')



library(rgl)
x_coords=1:63
y_coords=x_coords

lower_bound <- MGG5
upper_bound <- GG5
zmean <- ifelse(lower_bound>0, "red", "white")
zmean <- ifelse(upper_bound<0, "blue", zmean ) 
open3d()
planes3d(0, 0, 1, 0, col = "gray", alpha = 0.25)
for (i in 1:63) {
  for (j in 1:63) {
    z1 <- lower_bound[i, j]
    z2 <- upper_bound[i, j]
    segments3d(rbind(
      c(i / 63, j / 63, z1 / 1),
      c(i / 63, j / 63, z2 / 1)
    ), col = zmean[i,j], lwd=2)
  }
}
legend3d("topright", legend =c("l,u<0",  "l,u>0"), lty=1,lwd=2, cex=1,
         col =c("blue",   "red"))
title3d(main = expression(Gamma[3]*Gamma[3]^T), cex=1)
snapshot3d(filename = 'GAMMA3.png', fmt = 'png')

library(rgl)
x_coords=1:63
y_coords=x_coords


lower_bound <- MGG1
upper_bound <- GG1
zmean <- ifelse(lower_bound>0, "red", "white")
zmean <- ifelse(upper_bound<0, "blue", zmean ) 
open3d()
planes3d(0, 0, 1, 0, col = "gray", alpha = 0.25)
for (i in 1:63) {
  for (j in 1:63) {
    z1 <- lower_bound[i, j]
    z2 <- upper_bound[i, j]
    segments3d(rbind(
      c(i / 63, j / 63, z1 / 1),
      c(i / 63, j / 63, z2 / 1)
    ), col = zmean[i,j], lwd=2)
  }
}
legend3d("topright", legend =c("l,u<0",  "l,u>0"), lty=1,lwd=2, cex=1,
         col =c("blue",   "red"))
title3d(main = expression(Gamma[4]*Gamma[4]^T), cex=1)
snapshot3d(filename = 'GAMMA4.png', fmt = 'png')



library(rgl)
x_coords=1:63
y_coords=x_coords

lower_bound <- MGG2
upper_bound <- GG2
zmean <- ifelse(lower_bound>0, "red", "white")
zmean <- ifelse(upper_bound<0, "blue", zmean ) 
open3d()
planes3d(0, 0, 1, 0, col = "gray", alpha = 0.25)
for (i in 1:63) {
  for (j in 1:63) {
    z1 <- lower_bound[i, j]
    z2 <- upper_bound[i, j]
    segments3d(rbind(
      c(i / 63, j / 63, z1 / 1),
      c(i / 63, j / 63, z2 / 1)
    ), col = zmean[i,j], lwd=2)
  }
}
legend3d("topright", legend =c("l,u<0",  "l,u>0"), lty=1,lwd=2, cex=1,
         col =c("blue",   "red"))
title3d(main = expression(Gamma[5]*Gamma[5]^T), cex=1)
snapshot3d(filename = 'GAMMA5.png', fmt = 'png')




library(rgl)
x_coords=1:63
y_coords=x_coords

lower_bound <- MGG6
upper_bound <- GG6
zmean <- ifelse(lower_bound>0, "red", "white")
zmean <- ifelse(upper_bound<0, "blue", zmean ) 
open3d()
planes3d(0, 0, 1, 0, col = "gray", alpha = 0.25)
for (i in 1:63) {
  for (j in 1:63) {
    z1 <- lower_bound[i, j]
    z2 <- upper_bound[i, j]
    segments3d(rbind(
      c(i / 63, j / 63, z1 / 1),
      c(i / 63, j / 63, z2 / 1)
    ), col = zmean[i,j], lwd=2)
  }
}
legend3d("topright", legend =c("l,u<0",  "l,u>0"), lty=1,lwd=2, cex=1,
         col =c("blue",   "red"))
title3d(main = expression(Gamma[6]*Gamma[6]^T), cex=1)
snapshot3d(filename = 'GAMMA6.png', fmt = 'png')

library(rgl)
x_coords=1:63
y_coords=x_coords

lower_bound <- MGG7
upper_bound <- GG7
zmean <- ifelse(lower_bound>0, "red", "white")
zmean <- ifelse(upper_bound<0, "blue", zmean ) 
open3d()
planes3d(0, 0, 1, 0, col = "gray", alpha = 0.25)
for (i in 1:63) {
  for (j in 1:63) {
    z1 <- lower_bound[i, j]
    z2 <- upper_bound[i, j]
    segments3d(rbind(
      c(i / 63, j / 63, z1 / 1),
      c(i / 63, j / 63, z2 / 1)
    ), col = zmean[i,j], lwd=2)
  }
}
legend3d("topright", legend =c("l,u<0",  "l,u>0"), lty=1,lwd=2, cex=1,
         col =c("blue",   "red"))
title3d(main = expression(Gamma[7]*Gamma[7]^T), cex=1)
snapshot3d(filename = 'GAMMA7.png', fmt = 'png')



library(rgl)
x_coords=1:63
y_coords=x_coords
lower_bound <- MGG8
upper_bound <- GG8
zmean <- ifelse(lower_bound>0, "red", "white")
zmean <- ifelse(upper_bound<0, "blue", zmean ) 
open3d()
planes3d(0, 0, 1, 0, col = "gray", alpha = 0.25)
for (i in 1:63) {
  for (j in 1:63) {
    z1 <- lower_bound[i, j]
    z2 <- upper_bound[i, j]
    segments3d(rbind(
      c(i / 63, j / 63, z1 / 1),
      c(i / 63, j / 63, z2 / 1)
    ), col = zmean[i,j], lwd=2)
  }
}
legend3d("topright", legend =c("l,u<0",  "l,u>0"), lty=1,lwd=2, cex=1,
         col =c("blue",   "red"))
title3d(main = expression(Gamma[8]*Gamma[8]^T), cex=1)
snapshot3d(filename = 'GAMMA8.png', fmt = 'png')
