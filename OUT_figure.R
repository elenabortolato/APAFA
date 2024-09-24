#this code reproduces Figure 4:
#boxplots of RV coefficient for estimated Lambda in the simulation study

library(ggplot2)
library(ggpubr)

rm(list=ls())
#long scenario n>p
setwd("C:/flexMSFA/simulation/long")
 
RIS1=readRDS("RIS_1.RDS") 
RIS2=readRDS("RIS_2.RDS") 
RIS3=readRDS("RIS_3.RDS") 
RIS4A=readRDS("RIS_4A.RDS") 
RIS4B=readRDS("RIS_4B.RDS")
  
RIS1TETRIS=readRDS("RIS1_TETRIS.RDS")         
RIS2TETRIS=readRDS("RIS2_TETRIS.RDS") 
RIS4ATETRIS=readRDS("RIS4A_TETRIS.RDS")

OUT_ALL = rbind(RIS1,RIS3, RIS2, RIS4A, RIS4B, RIS1TETRIS,RIS2TETRIS,RIS4ATETRIS)
 
len=c(dim(RIS1)[1], dim(RIS3)[1], dim(RIS2)[1], dim(RIS4A)[1], dim(RIS4B)[1],
 dim(RIS1TETRIS)[1], dim(RIS2TETRIS)[1], dim(RIS4ATETRIS)[1])
labels=c("A",  "A*","B", "C", "D", "A", "B","C")
setting=unlist(sapply(1:8, function (x) rep(labels[x], each=len[x])))
method=c("flexMSFA", "flexMSFA","flexMSFA","flexMSFA","flexMSFA", 
         "TETRIS", "TETRIS","TETRIS")
method=c("APAFA", "APAFA","APAFA","APAFA","APAFA", 
         "TETRIS", "TETRIS","TETRIS")
model=unlist(sapply(1:8, function (x) rep(method[x], each=len[x])))
OUT_ALL=cbind(OUT_ALL,setting,model)
head(OUT_ALL)  
 
OUT_ALL=as.data.frame(OUT_ALL)
str(OUT_ALL)
OUT_ALL$norm_eta=as.numeric(OUT_ALL$norm_eta)
plot1=ggplot2::ggplot(OUT_ALL,aes(x=setting, y=norm_eta, color=model))+
  geom_boxplot()+
  coord_cartesian(ylim = c(0, 1.0))+
  theme_minimal()+
  scale_color_grey() +
  ylab(expression("RV("*Lambda*hat(Lambda)*")" ))

plot1 
 


#large scenario n<p

setwd("flexMSFA/simulation/large")
 
RIS1LARGE=readRDS("RIS1LARGE.RDS") 
RIS2LARGE=readRDS("RIS_2_large.RDS") 
RIS3LARGE=readRDS("RIS_3_large.RDS") 
RIS4ALARGE=readRDS("RIS_4A_large.RDS") 
RIS4BLARGE=readRDS("RIS_4B_large.RDS")

setwd("C:/Users/Asus/Desktop/SCEN1/TETRIS1")
RIS1TETRIS=readRDS("RIS1TETRIS_large.RDS") 
 setwd("C:/Users/Asus/Dropbox/FactorModels/SCEN2")
RIS2TETRIS=readRDS("RIS2TETRIS_large.RDS")
setwd("C:/Users/Asus/Dropbox/FactorModels/SCEN4A")
RIS4ATETRIS=readRDS("RIS4ATETRIS_large.RDS")

OUT_ALL = rbind(RIS1LARGE,RIS3LARGE, RIS2LARGE, RIS4ALARGE, RIS4BLARGE, 
                RIS1TETRIS,RIS2TETRIS,RIS4ATETRIS)



len=c(dim(RIS1LARGE)[1], dim(RIS3LARGE)[1], dim(RIS2LARGE)[1], dim(RIS4ALARGE)[1], dim(RIS4BLARGE)[1],
      dim(RIS1TETRIS)[1], dim(RIS2TETRIS)[1], dim(RIS4ATETRIS)[1])
labels=c("A",  "A*","B", "C", "D", "A", "B","C")
setting=rep(labels ,len)
method=c("flexMSFA", "flexMSFA","flexMSFA","flexMSFA","flexMSFA", 
         "TETRIS", "TETRIS","TETRIS")
method=c("APAFA", "APAFA","APAFA","APAFA","APAFA", 
         "TETRIS", "TETRIS","TETRIS")
model=rep(method ,len)
OUT_ALL=cbind(OUT_ALL,setting,model)
head(OUT_ALL) 

OUT_ALL=as.data.frame(OUT_ALL)
str(OUT_ALL)
OUT_ALL$norm_eta=as.numeric(OUT_ALL$norm_eta)
plot2=ggplot2::ggplot(OUT_ALL,aes(x=setting, y=norm_eta, color=model))+
  geom_boxplot()+
  theme_minimal()+scale_color_grey() +
  coord_cartesian(ylim = c(0, 1))+
  ylab(expression("RV("*Lambda*hat(Lambda)*")" )) 
 
ggpubr::ggarrange(plot1,plot2, common.legend = T, ncol=2,labels = c("n>p", "n<p"),
                    vjust=0.6,hjust = -6.9, font.label = 1)
 
 