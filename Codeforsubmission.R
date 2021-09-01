###############################
#####Code for dissertation#####
###############################

#Library packages
library(readxl)
library(dplyr)
library(lubridate)
library(igraph)
library(ggplot2)
library(evir)
library(latex2exp)
library(POT)
library(fastDummies)
library(ismev)
library(evmix)
library(tea)
library(mev)
library(evd)
library(itertools)
library(Rlab)
defaultW <- getOption("warn") 

#Define the norm
sumnorm <- function(x) sum(x)
l2norm <- function(x) sqrt(sum(x^2))

#Install data
data <-  read_excel("C:/Users/Lucas/Desktop/Dissertation/Code/OpRisk_2011_2020.xlsx")
data 
summary(data)

#Group data by week
data.week <-data %>% mutate(day_number = (as.integer(Date-min(Date))+1)%/%86400,week = day_number%/%7) %>%
  arrange(Date)

data.week.sum <- data.week %>%
  group_by(week,BU,`Risk type`) %>%
  summarise(Loss = sum(Amount))

#Get the name of Risk Types and Business Units
nodesSet1 <- unique(data$`Risk type`)
nodesSet2 <- unique(data$BU)
nodesSet1
nodesSet2
d <- length(nodesSet1)
q <- length(nodesSet2)

#Compute the fraction matrix
n <- max(data.week.sum$week)+1
L <- array(0,dim = c(n,q,d),dimnames = list(1:n,nodesSet2,nodesSet1))
for (i in 1:n){
  test <- data.week.sum%>% dplyr::filter(week==i-1)
  testA <- matrix(0,nrow = q,ncol = d,dimnames = list(nodesSet2,nodesSet1))
  testA[cbind(test$BU,test$`Risk type`)]<-test$Loss
  L[i,,] <- testA
}

RT <- matrix(0,nrow=n, ncol=9 ,dimnames = list(1:n,nodesSet1))
BU <- matrix(0,nrow = n,ncol =4,dimnames = list(1:n,nodesSet2))
for (i in 1:n){
  RT[i,] <- colSums(L[i,,])
  BU[i,] <- rowSums(L[i,,])
}

A <- array(0,dim = c(n,q,d),dimnames = list(1:n,nodesSet2,nodesSet1))
for (i in 1:n){
  A[i,,] <- apply(L[i,,],2,function(x) x/sum(x))
}
A[is.nan(A)] <- 0
A[1,,]

#Compute the marginal loss data aggregated by risk types or business units
data.week.sumET <-data.week%>% group_by(week,`Risk type`)%>%
  summarise(Loss = sum(Amount))
data.week.sumBU <- data.week%>% group_by(week,BU)%>%
  summarise(Loss = sum(Amount))

#Can also compute using the following method
#Notice that RT is used before, so we define ET here which should be the same as RT
ET <- matrix(0,nrow = n,ncol = d,dimnames = list(1:n,nodesSet1))
BU <- matrix(0,nrow = n,ncol = q,dimnames = list(1:n,nodesSet2))
for (i in 1:n){
  test <- data.week.sumET%>% dplyr::filter(week==i-1)
  ET[i,][test$`Risk type`] <- test$Loss
  test <- data.week.sumBU %>%dplyr::filter(week==i-1)
  BU[i,][test$BU] <- test$Loss
}
head(ET)
head(BU)

#Calculate the percentage of severity corresponding to different BU and RT
ETsum <- apply(ET,2,sum)
BUsum <- apply(BU,2,sum)
ETsum/sum(ETsum)
BUsum/sum(BUsum)


apply(ET ==0,2,mean)
apply(BU ==0,2,mean)

#Plot of Figure 3.1 and 3.2
boxplot(log(data.week.sumET$Loss)~data.week.sumET$`Risk type`,ylab = "Logloss",xlab="Risk types")
boxplot(log(data.week.sumBU$Loss)~data.week.sumBU$BU,ylab = "Logloss",xlab = "Business units")

#Compute Marginal data
EDPM1 <- ET[,'EDPM1']
EF1 <- ET[,'EF1']
CPBP <- ET[,'CPBP']
EDPM2 <- ET[,'EDPM2']
BDSF <- ET[,'BDSF']
EPWS <- ET[,'EPWS']
DPA <- ET[,'DPA']
EF2 <- ET[,'EF2']
IF <- ET[,'IF']
Retail <- BU[,"Retail"]
Associates <- BU[,"Associates"]
IB <- BU[,"IB"]
Loan <- BU[,"Loan"]

#Use mrlplot to have a general idea of thresholds
mrlplot(EDPM1)
mrlplot(EF1)
mrlplot(CPBP)
mrlplot(EDPM2)
mrlplot(BDSF)
mrlplot(EPWS)
mrlplot(DPA)
mrlplot(EF2)
mrlplot(IF)
##################
#Independence test
##################
library(energy)
#For the first part, it gives p-value but nan for statistic
Avector = matrix(0,nrow=470,ncol=36)
for (i in 1:470){
  Avector[i,] <- as.vector(A[i,,])
}
indep.test(Avector,RT,method = "mvI",R=199)
#Can scale the data to have the statistic
RT1 <- RT/1000000
indep.test(Avector,RT1,method = "mvI",R=199)


##############################################################
#Use automatic method which considered Bias-Variance tradeoff#
##############################################################
#Define vector for auto thresholds
thresholdsauto <- rep(0,9)

#EDPM
SEDPM1 <- sort(EDPM1,decreasing = TRUE)
krange <- (min(40,floor(0.02*length(SEDPM1))):floor(0.2*length(SEDPM1)))
paramsEDPM1 <- matrix(0,nrow=2,ncol = length(krange))
#Compute estimates for different k corresponding to different thresholds
for (j in seq_along(krange)){
  a <- fit.gpd(EDPM1,threshold = SEDPM1[krange[j]])
  paramsEDPM1[,j] <- a$estimate
}
#Estimate the prediction error for different k
e <- rep(0,length(krange))
for (j in 1:length(krange)){
  k = krange[j]
  beta = paramsEDPM1[1,j]
  xi = paramsEDPM1[2,j]
  e1 = 0
  for (i in 1:k){
    e1 = e1 + 1/(beta^2)*(i/(k+1))^(2*xi)/((k+1)/i-1)*(SEDPM1[i]-SEDPM1[k+1]-(beta/xi)*((i/(k+1))^(-xi)-1))^2
    e1 = e1 + 2/k*(i/(k+1))^(2*xi)/((k+1)/i-1)*(1+xi)^2*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)^2
    e1 = e1 + 4/k*(i/(k+1))^(xi)/((k+1)/i-1)*(1+xi)*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)*log(i/(k+1))
    e1 = e1 + 2/k*((1/((k+1)/i-1))*((1+xi)/xi)^2*(log(i/(k+1)))^2)
    e1 = e1 - 1/k 
  }
  e[j] = e1
}
#Compute k and thresholds based on prediction error
which.min(e)
krange[which.min(e)]
thresholdsauto[1] <- SEDPM1[krange[which.min(e)]]
thresholdsauto

#EF1
SEF1 <- sort(EF1,decreasing = TRUE)
krange <- (min(40,floor(0.01*length(SEF1))):floor(0.3*length(SEF1)))
paramsEF1 <- matrix(0,nrow=2,ncol = length(krange))
#Compute estimates for different k corresponding to different thresholds
for (j in seq_along(krange)){
  a <- fit.gpd(EF1,threshold = EF1[krange[j]])
  paramsEF1[,j] <- a$estimate
}
paramsEF1
#Estimate the prediction error for different k
e <- rep(0,length(krange))
for (j in 1:length(krange)){
  k = krange[j]
  beta = paramsEF1[1,j]
  xi = paramsEF1[2,j]
  e1 = 0
  for (i in 1:k){
    e1 = e1 + 1/(beta^2)*(i/(k+1))^(2*xi)/((k+1)/i-1)*(SEF1[i]-SEF1[k+1]-(beta/xi)*((i/(k+1))^(-xi)-1))^2
    e1 = e1 + 2/k*(i/(k+1))^(2*xi)/((k+1)/i-1)*(1+xi)^2*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)^2
    e1 = e1 + 4/k*(i/(k+1))^(xi)/((k+1)/i-1)*(1+xi)*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)*log(i/(k+1))
    e1 = e1 + 2/k*((1/((k+1)/i-1))*((1+xi)/xi)^2*(log(i/(k+1)))^2)
    e1 = e1 - 1/k 
  }
  e[j] = e1
}
#Compute k and thresholds based on prediction error
which.min(e)
krange[which.min(e)]
thresholdsauto[2]<-SEF1[krange[which.min(e)]]
thresholdsauto

#CPBP
SCPBP <- sort(CPBP,decreasing = TRUE)
krange <- (min(40,floor(0.02*length(SCPBP))):floor(0.2*length(SCPBP)))
paramsCPBP <- matrix(0,nrow=2,ncol = length(krange))
#Compute estimates for different k corresponding to different thresholds
for (j in seq_along(krange)){
  a <- fit.gpd(CPBP,threshold = SCPBP[krange[j]])
  paramsCPBP[,j] <- a$estimate
}
paramsCPBP
#Estimate the prediction error for different k
e <- rep(0,length(krange))
for (j in 1:length(krange)){
  k = krange[j]
  beta = paramsCPBP[1,j]
  xi = paramsCPBP[2,j]
  e1 = 0
  for (i in 1:k){
    e1 = e1 + 1/(beta^2)*(i/(k+1))^(2*xi)/((k+1)/i-1)*(SCPBP[i]-SCPBP[k+1]-(beta/xi)*((i/(k+1))^(-xi)-1))^2
    e1 = e1 + 2/k*(i/(k+1))^(2*xi)/((k+1)/i-1)*(1+xi)^2*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)^2
    e1 = e1 + 4/k*(i/(k+1))^(xi)/((k+1)/i-1)*(1+xi)*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)*log(i/(k+1))
    e1 = e1 + 2/k*((1/((k+1)/i-1))*((1+xi)/xi)^2*(log(i/(k+1)))^2)
    e1 = e1 - 1/k 
  }
  e[j] = e1
}
#Compute k and thresholds based on prediction error
which.min(e)
krange[which.min(e)]
thresholdsauto[3]<-SCPBP[krange[which.min(e)]]
thresholdsauto

#EDPM2
SEDPM2 <- sort(EDPM2,decreasing = TRUE)
krange <- (min(40,floor(0.02*length(SEDPM2))):floor(0.2*length(SEDPM2)))
paramsEDPM2 <- matrix(0,nrow=2,ncol = length(krange))
#Compute estimates for different k corresponding to different thresholds
for (j in seq_along(krange)){
  a <- fit.gpd(EDPM2,threshold = SEDPM2[krange[j]])
  paramsEDPM2[,j] <- a$estimate
}
paramsEDPM2
#Estimate the prediction error for different k
e <- rep(0,length(krange))
for (j in 1:length(krange)){
  k = krange[j]
  beta = paramsEDPM2[1,j]
  xi = paramsEDPM2[2,j]
  e1 = 0
  for (i in 1:k){
    e1 = e1 + 1/(beta^2)*(i/(k+1))^(2*xi)/((k+1)/i-1)*(SEDPM2[i]-SEDPM2[k+1]-(beta/xi)*((i/(k+1))^(-xi)-1))^2
    e1 = e1 + 2/k*(i/(k+1))^(2*xi)/((k+1)/i-1)*(1+xi)^2*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)^2
    e1 = e1 + 4/k*(i/(k+1))^(xi)/((k+1)/i-1)*(1+xi)*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)*log(i/(k+1))
    e1 = e1 + 2/k*((1/((k+1)/i-1))*((1+xi)/xi)^2*(log(i/(k+1)))^2)
    e1 = e1 - 1/k 
  }
  e[j] = e1
}
#Compute k and thresholds based on prediction error
which.min(e)
krange[which.min(e)]
thresholdsauto[4]<-SEDPM2[krange[which.min(e)]]
thresholdsauto

#BDSF
SBDSF <- sort(BDSF,decreasing = TRUE)
krange <- (min(40,floor(0.02*length(SBDSF))):floor(0.2*length(SBDSF)))
paramsBDSF <- matrix(0,nrow=2,ncol = length(krange))
#Compute estimates for different k corresponding to different thresholds
for (j in seq_along(krange)){
  a <- fit.gpd(BDSF,threshold = SBDSF[krange[j]])
  paramsBDSF[,j] <- a$estimate
}
paramsBDSF
#Estimate the prediction error for different k
e <- rep(0,length(krange))
for (j in 1:length(krange)){
  k = krange[j]
  beta = paramsBDSF[1,j]
  xi = paramsBDSF[2,j]
  e1 = 0
  for (i in 1:k){
    e1 = e1 + 1/(beta^2)*(i/(k+1))^(2*xi)/((k+1)/i-1)*(SBDSF[i]-SBDSF[k+1]-(beta/xi)*((i/(k+1))^(-xi)-1))^2
    e1 = e1 + 2/k*(i/(k+1))^(2*xi)/((k+1)/i-1)*(1+xi)^2*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)^2
    e1 = e1 + 4/k*(i/(k+1))^(xi)/((k+1)/i-1)*(1+xi)*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)*log(i/(k+1))
    e1 = e1 + 2/k*((1/((k+1)/i-1))*((1+xi)/xi)^2*(log(i/(k+1)))^2)
    e1 = e1 - 1/k 
  }
  e[j] = e1
}
#Compute k and thresholds based on prediction error
which.min(e)
krange[which.min(e)]
thresholdsauto[5]<-SBDSF[krange[which.min(e)]]
thresholdsauto

#EPWS
SEPWS <- sort(EPWS,decreasing = TRUE)
krange <- (min(40,floor(0.01*length(SEPWS))):floor(0.4*length(SEPWS)))
paramsEPWS <- matrix(0,nrow=2,ncol = length(krange))
#Compute estimates for different k corresponding to different thresholds
for (j in seq_along(krange)){
  a <- fit.gpd(EPWS,threshold = SEPWS[krange[j]])
  paramsEPWS[,j] <- a$estimate
}
paramsEPWS
#Estimate the prediction error for different k
e <- rep(0,length(krange))
for (j in 1:length(krange)){
  k = krange[j]
  beta = paramsEPWS[1,j]
  xi = paramsEPWS[2,j]
  e1 = 0
  for (i in 1:k){
    e1 = e1 + 1/(beta^2)*(i/(k+1))^(2*xi)/((k+1)/i-1)*(SEPWS[i]-SEPWS[k+1]-(beta/xi)*((i/(k+1))^(-xi)-1))^2
    e1 = e1 + 2/k*(i/(k+1))^(2*xi)/((k+1)/i-1)*(1+xi)^2*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)^2
    e1 = e1 + 4/k*(i/(k+1))^(xi)/((k+1)/i-1)*(1+xi)*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)*log(i/(k+1))
    e1 = e1 + 2/k*((1/((k+1)/i-1))*((1+xi)/xi)^2*(log(i/(k+1)))^2)
    e1 = e1 - 1/k 
  }
  e[j] = e1
}
#Compute k and thresholds based on prediction error
which.min(e)
krange[which.min(e)]
thresholdsauto[6]<-SEPWS[krange[which.min(e)]]
thresholdsauto

#DPA
SDPA <- sort(DPA,decreasing = TRUE)
krange <- (min(40,floor(0.02*length(DPA))):floor(0.2*length(SDPA)))
paramsDPA <- matrix(0,nrow=2,ncol = length(krange))
#Compute estimates for different k corresponding to different thresholds
for (j in seq_along(krange)){
  a <- fit.gpd(DPA,threshold = SDPA[krange[j]])
  paramsDPA[,j] <- a$estimate
}
paramsDPA
#Estimate the prediction error for different k
e <- rep(0,length(krange))
for (j in 1:length(krange)){
  k = krange[j]
  beta = paramsDPA[1,j]
  xi = paramsDPA[2,j]
  e1 = 0
  for (i in 1:k){
    e1 = e1 + 1/(beta^2)*(i/(k+1))^(2*xi)/((k+1)/i-1)*(SDPA[i]-SDPA[k+1]-(beta/xi)*((i/(k+1))^(-xi)-1))^2
    e1 = e1 + 2/k*(i/(k+1))^(2*xi)/((k+1)/i-1)*(1+xi)^2*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)^2
    e1 = e1 + 4/k*(i/(k+1))^(xi)/((k+1)/i-1)*(1+xi)*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)*log(i/(k+1))
    e1 = e1 + 2/k*((1/((k+1)/i-1))*((1+xi)/xi)^2*(log(i/(k+1)))^2)
    e1 = e1 - 1/k 
  }
  e[j] = e1
}
#Compute k and thresholds based on prediction error
which.min(e)
krange[which.min(e)]
thresholdsauto[7]<-SDPA[krange[which.min(e)]]
thresholdsauto

#EF2
SEF2 <- sort(EF2,decreasing = TRUE)
krange <- (min(40,floor(0.02*length(EF2))):floor(0.2*length(SEF2)))
paramsEF2 <- matrix(0,nrow=2,ncol = length(krange))
#Compute estimates for different k corresponding to different thresholds
for (j in seq_along(krange)){
  a <- fit.gpd(EF2,threshold = SEF2[krange[j]])
  paramsEF2[,j] <- a$estimate
}
paramsEF2
#Estimate the prediction error for different k
e <- rep(0,length(krange))
for (j in 1:length(krange)){
  k = krange[j]
  beta = paramsEF2[1,j]
  xi = paramsEF2[2,j]
  e1 = 0
  for (i in 1:k){
    e1 = e1 + 1/(beta^2)*(i/(k+1))^(2*xi)/((k+1)/i-1)*(SEF2[i]-SEF2[k+1]-(beta/xi)*((i/(k+1))^(-xi)-1))^2
    e1 = e1 + 2/k*(i/(k+1))^(2*xi)/((k+1)/i-1)*(1+xi)^2*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)^2
    e1 = e1 + 4/k*(i/(k+1))^(xi)/((k+1)/i-1)*(1+xi)*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)*log(i/(k+1))
    e1 = e1 + 2/k*((1/((k+1)/i-1))*((1+xi)/xi)^2*(log(i/(k+1)))^2)
    e1 = e1 - 1/k 
  }
  e[j] = e1
}
#Compute k and thresholds based on prediction error
which.min(e)
krange[which.min(e)]
thresholdsauto[8]<-SEF2[krange[which.min(e)]]
thresholdsauto

#IF
SIF <- sort(IF,decreasing = TRUE)
krange <- (min(40,floor(0.02*length(IF))):floor(0.2*length(IF)))
paramsIF <- matrix(0,nrow=2,ncol = length(krange))
#Compute estimates for different k corresponding to different thresholds
for (j in seq_along(krange)){
  a <- fit.gpd(IF,threshold = IF[krange[j]])
  paramsIF[,j] <- a$estimate
}
paramsIF
#Estimate the prediction error for different k
e <- rep(0,length(krange))
for (j in 1:length(krange)){
  k = krange[j]
  beta = paramsIF[1,j]
  xi = paramsIF[2,j]
  e1 = 0
  for (i in 1:k){
    e1 = e1 + 1/(beta^2)*(i/(k+1))^(2*xi)/((k+1)/i-1)*(SIF[i]-SIF[k+1]-(beta/xi)*((i/(k+1))^(-xi)-1))^2
    e1 = e1 + 2/k*(i/(k+1))^(2*xi)/((k+1)/i-1)*(1+xi)^2*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)^2
    e1 = e1 + 4/k*(i/(k+1))^(xi)/((k+1)/i-1)*(1+xi)*(1+2*xi)/(xi^2)*(((i/(k+1))^(-xi)-1)/xi)*log(i/(k+1))
    e1 = e1 + 2/k*((1/((k+1)/i-1))*((1+xi)/xi)^2*(log(i/(k+1)))^2)
    e1 = e1 - 1/k 
  }
  e[j] = e1
}
#Compute k and thresholds based on prediction error
which.min(e)
krange[which.min(e)]
thresholdsauto[9]<-SIF[krange[which.min(e)]]
thresholdsauto

#Compute the results
thresholdsauto <-matrix(thresholdsauto,nrow = 1,ncol = 9,dimnames = list(c("thresholds"),nodesSet1))
thresholdsauto

#Compute the parameter estimation which gives 2 negative shape estiamates
params <- matrix(0,nrow = 3,ncol = 9,dimnames = (list(c('u','scale','shape'),nodesSet1)))
params[1,] <- thresholdsauto
params[2:3,1] <- fitgpd(EDPM1,thresholdsauto[1],est = "mple")$param
params[2:3,2] <- fitgpd(EF1,thresholdsauto[2],est = "mple")$param
params[2:3,3] <- fitgpd(CPBP,thresholdsauto[3],est = "mle")$param
params[2:3,4] <- fitgpd(EDPM2,thresholdsauto[4],est = "mple")$param
params[2:3,5] <- fitgpd(BDSF,thresholdsauto[5],est = "mple")$param
params[2:3,6] <- fitgpd(EPWS,thresholdsauto[6],est = "mple")$param
params[2:3,7] <- fitgpd(DPA,thresholdsauto[7],est = "mple")$param
params[2:3,8] <- fitgpd(EF2,thresholdsauto[8],est = "mple")$param
params[2:3,9] <- fitgpd(IF,thresholdsauto[9],est = "mple")$param
params

####################################
#Threshold Selection based on dAMSE#
####################################
thresholdsmse <- matrix(0,nrow = 1,ncol = 9,dimnames = list(c('thresholdsmse'),nodesSet1))
hillshape <- matrix(0,nrow = 2,ncol = 9,dimnames = list(c('hill shape',"hill scale"),nodesSet1))
thresholdsmse[1,'EDPM1'] <- dAMSE(EDPM1)$threshold
thresholdsmse[1,'EF1'] <- dAMSE(EF1)$threshold
thresholdsmse[1,'CPBP'] <- dAMSE(CPBP)$threshold
thresholdsmse[1,'EDPM2'] <- dAMSE(EDPM2)$threshold
thresholdsmse[1,'BDSF'] <- dAMSE(BDSF)$threshold
thresholdsmse[1,'EPWS'] <- dAMSE(EPWS[EPWS>0])$threshold
thresholdsmse[1,'DPA'] <- dAMSE(DPA[DPA>0])$threshold
thresholdsmse[1,'EF2'] <- dAMSE(EF2[EF2>0])$threshold
thresholdsmse[1,'IF'] <- dAMSE(IF[IF>0])$threshold

#It can also prvides hill estimators
hillshape[1,'EDPM1'] <- 1/dAMSE(EDPM1)$tail.index
hillshape[1,'EF1'] <- 1/dAMSE(EF1)$tail.index
hillshape[1,'CPBP'] <- 1/dAMSE(CPBP)$tail.index
hillshape[1,'EDPM2'] <- 1/dAMSE(EDPM2)$tail.index
hillshape[1,'BDSF'] <- 1/dAMSE(BDSF)$tail.index
hillshape[1,'EPWS'] <- 1/dAMSE(EPWS[EPWS>0])$tail.index
hillshape[1,'DPA'] <- 1/dAMSE(DPA[DPA>0])$tail.index
hillshape[1,'EF2'] <- 1/dAMSE(EF2[EF2>0])$tail.index
hillshape[1,'IF'] <- 1/dAMSE(IF[IF>0])$tail.index
thresholdsmse
hillshape

#Fit the gpd distribution 
paramsmse <- matrix(0,nrow = 3,ncol = 9,dimnames = (list(c('u','scale','shape'),nodesSet1)))
paramsmse[1,] <- thresholdsmse
for (i in 1:9){
  paramsmse[2:3,i] <- fitgpd(eval(parse(text = nodesSet1[i])),thresholdsmse[i],'mple')$param
}

#Compute the confidence interval
paramsse <- matrix(0,nrow = 1,ncol = 9,dimnames = list(c("se"),nodesSet1))
paramsse[1] <- fitgpd(EDPM1,thresholdsmse[1],est = "mple")$std.err[2]
paramsse[2] <- fitgpd(EF1,thresholdsmse[2],est = "mple")$std.err[2]
paramsse[3] <- fitgpd(CPBP,thresholdsmse[3],est = "mple")$std.err[2]
paramsse[4] <- fitgpd(EDPM2,thresholdsmse[4],est = "mple")$std.err[2]
paramsse[5] <- fitgpd(BDSF,thresholdsmse[5],est = "mple")$std.err[2]
paramsse[6] <- fitgpd(EPWS,thresholdsmse[6],est = "mple")$std.err[2]
paramsse[7] <- fitgpd(DPA,thresholdsmse[7],est = "mple")$std.err[2]
paramsse[8] <- fitgpd(EF2,thresholdsmse[8],est = "mple")$std.err[2]
paramsse[9] <- fitgpd(IF,thresholdsmse[9],est = "mple")$std.err[2]
paramsse
paramsmse[3,] - 1.96*paramsse
paramsmse[3,] + 1.96*paramsse

############
#Mean Alpha#
############
#The parameter estimation for gpd method need to mannually change the threshold to make fits available
paramsmse1 <- matrix(0,nrow = 3,ncol = 9,dimnames = (list(c('u','scale','shape'),nodesSet1)))
paramsmse1[1,] <- paramsmse[1,]
shapemeanmse <- mean(paramsmse[3,]) 
paramsmse1[2:3,1] <-fitgpd(EDPM1,paramsmse[1,1]+100,shape = shapemeanmse)$param
paramsmse1[1,1] <- paramsmse1[1,1]
paramsmse1[2:3,2] <-fitgpd(EF1,paramsmse1[1,2],shape = shapemeanmse)$param
paramsmse1[2:3,3] <-fitgpd(CPBP,paramsmse[1,3],shape = shapemeanmse)$param
paramsmse1[2:3,4] <-fitgpd(EDPM2,paramsmse[1,4]+100,shape=shapemeanmse)$param
paramsmse1[1,4] <- paramsmse1[1,4] +100
paramsmse1[2:3,5] <-fitgpd(BDSF,paramsmse1[1,5]+100,shape=shapemeanmse)$param
paramsmse1[1,5] <- paramsmse1[1,5]+100
paramsmse1[2:3,6] <-fitgpd(EPWS,paramsmse1[1,6]+300,shape=shapemeanmse)$param
paramsmse1[1,6] <- paramsmse1[1,6]+300
paramsmse1[2:3,7] <- fitgpd(DPA,paramsmse1[1,7]+200,shape = shapemeanmse)$param
paramsmse1[1,7] <- paramsmse1[1,7]+200
paramsmse1[2:3,8] <-fitgpd(EF2, paramsmse[1,8]+500,shape = shapemeanmse)$param
paramsmse1[1,8] <- paramsmse1[1,8]+500
paramsmse1[2:3,9] <-fitgpd(IF, paramsmse1[1,9],shape=shapemeanmse)$param

#KS TEST
library(fExtremes)
ks.test(SEDPM1[1:70],rgpd(n=50,mu=paramsmse1[1,1], beta= paramsmse1[2,1],xi = paramsmse1[3,1]))
ks.test(SEF1[1:200],rgpd(n=50,mu=paramsmse1[1,2],beta = paramsmse1[2,2],xi = paramsmse1[3,2]))
ks.test(SCPBP[1:50],rgpd(n=50,mu=paramsmse1[1,3],beta = paramsmse1[2,3],xi = paramsmse1[3,3]))
ks.test(SEDPM2[1:40],rgpd(n=50,mu=paramsmse1[1,4],beta = paramsmse1[2,4],xi = paramsmse1[3,4]))
ks.test(SBDSF[1:40],rgpd(n=50,mu=paramsmse1[1,5],beta = paramsmse1[2,5],xi = paramsmse1[3,5]))
ks.test(SEPWS[1:40],rgpd(n=50,mu=paramsmse1[1,6],beta = paramsmse1[2,6],xi = paramsmse1[3,6]))
ks.test(SDPA[1:50],rgpd(n=50,mu=paramsmse1[1,7],beta = paramsmse1[2,7],xi = paramsmse1[3,7]))
ks.test(SEF2[1:50],rgpd(n=50,mu=paramsmse1[1,8],beta = paramsmse1[2,8],xi = paramsmse1[3,8]))
ks.test(SIF[1:50],rgpd(n=50,mu=paramsmse1[1,9],beta = paramsmse1[2,9],xi = paramsmse1[3,9]))

#QQplot Figure 3.3
outer = FALSE
line = -2
cex = 2
adj  = 0.025
par(mfrow = c(3,3))
qqgpd(EDPM1,30,scale = paramsmse1[2,1],shape = paramsmse1[3,1])
title(outer=outer,adj=adj,main="RT1",cex.main=cex,col="black",font=2,line=line)
qqgpd(EF1,30,scale = paramsmse1[2,2],shape = paramsmse1[3,2])
title(outer=outer,adj=adj,main="RT2",cex.main=cex,col="black",font=2,line=line)
qqgpd(CPBP,30,scale = paramsmse1[2,3],shape = paramsmse1[3,3])
title(outer=outer,adj=adj,main="RT3",cex.main=cex,col="black",font=2,line=line)
qqgpd(BDSF,30,scale = paramsmse1[2,5],shape = paramsmse1[3,5])
title(outer=outer,adj=adj,main="RT4",cex.main=cex,col="black",font=2,line=line)
qqgpd(EPWS,30,scale = paramsmse1[2,6],shape = paramsmse1[3,6])
title(outer=outer,adj=adj,main="RT5",cex.main=cex,col="black",font=2,line=line)
qqgpd(DPA,30,scale = paramsmse1[2,7],shape = paramsmse1[3,7])
title(outer=outer,adj=adj,main="RT6",cex.main=cex,col="black",font=2,line=line)
qqgpd(EF2,30,scale = paramsmse1[2,8],shape = paramsmse1[3,8])
title(outer=outer,adj=adj,main="RT7",cex.main=cex,col="black",font=2,line=line)
qqgpd(IF,30,scale = paramsmse1[2,9],shape = paramsmse1[3,9])
title(outer=outer,adj=adj,main="RT8",cex.main=cex,col="black",font=2,line=line)
qqgpd(EDPM2,30,scale = paramsmse1[2,4],shape = paramsmse1[3,4])
title(outer=outer,adj=adj,main="RT9",cex.main=cex,col="black",font=2,line=line)

#Compute alpha and K_i
mparam1 <- matrix(0,nrow = 3,ncol = 9,dimnames = (list(c('tail.index','scaling1','scaling2'),nodesSet1)))
mparam1[1,] <-1/paramsmse1[3,]
for (j in 1:9){
  mparam1[2,j] <- (mparam1[1,j]*paramsmse1[2,j])^mparam1[1,j]*mean(eval(parse(text = nodesSet1[j]))>paramsmse1[1,j])
}
mparam1
#Compute the risk constants and risk measures
alpha <- 1/paramsmse1[3,1]
K <-mparam1[2,]
K
Kmat <-diag(mparam1[2,]^(1/mparam1[1,]))
Kmati <- solve(Kmat)
k <-50

#Scale y and compute order statistics
y <- matrix(0, nrow = n, ncol = 9)
for (i in 1:n){
  y[i,] <-Kmati%*% RT[i,]    
}

normy <- matrix(0,nrow=1,ncol=n,dimnames = list(c('norm'),0:(n-1)))
for (i in 1:n){
  normy[i] <- sumnorm(RT[i,])
}
onormy <-order(normy,decreasing = TRUE)

#risk constants for asymptotically independent case
#Business units
k=50
Amean = matrix(0,nrow = 4,ncol=9)
for (i in 1:k){
  Amean = Amean+A[onormy[i],,]
}
Amean = Amean/k
Ciind = matrix(0,nrow=4,ncol=1,dimnames = list(nodesSet2,c("Ci")))
for (i in 1:4){
  for(j in 1:9){
    Ciind[i] = Ciind[i] + K[j]*(Amean[i,j]^(alpha))
  }
}
Ciind
#Systematic
Aemean = rep(0,9)
for (i in 1:9){
  for (j in 1:k){
    Aemean[i] = Aemean[i] + sumnorm(A[onormy[j],,][,i])^(alpha)
  }
}
Aemean = Aemean/k
Csind <- 0
for (i in 1:9){
  Csind = Csind+K[i]*Aemean[i]
}
Csind

#risk constants for asymptotically fully dependent case
#Business units
Cidep = matrix(0,nrow=4,ncol=1,dimnames = list(nodesSet2,c("Ci")))
ones <-matrix(1,nrow=9,ncol=1)
for (i in 1:4){
  for (j in 1:k){
    Cidep[i] = Cidep[i] + ((A[onormy[j],,]%*%Kmat%*%ones)[i])^(alpha)
  }
}
Cidep = Cidep/k
Cidep
#Systematic
Csdep = 0
for (i in 1:k){
  Csdep  = Csdep + sumnorm(A[onormy[i],,]%*%Kmat%*%ones)^(alpha)
}
Csdep = Csdep/k
Csdep

#Use risk constants to compute VaR and ES for independent case
variind <- Ciind^(1/alpha)* 0.001^(-1/alpha)
variind
varsind <- Csind^(1/alpha)* 0.001^(-1/alpha)
varsind
cotiind <- alpha/(alpha-1)*variind
cotsind <- alpha/(alpha-1)*varsind
cotiind
cotsind
#Use risk constants to compute VaR and ES for dependent case
varidep <- Cidep^(1/alpha)* 0.001^(-1/alpha)
varidep
varsdep <- Csdep^(1/alpha)* 0.001^(-1/alpha)
varsdep
cotidep <- alpha/(alpha-1)*varidep
cotidep
cotsdep <- alpha/(alpha-1)*varsdep
cotsdep

#Compute estiamtes based on spectral measure
ynew <- matrix(0,nrow=470,ncol = 9)
for (i in 1:470){
  ynew[i,] <- y[i,]/sumnorm(y[i,])
}

#Compute denominator
denomC <- 0
for (i in 1:k){
  denomC <- denomC + ynew[onormy[i],1]^(alpha) 
}
denomC

#Systemic constant
CsA <-0
for (i in 1:470){
  temp = 0
  for (j in 1:k){
    temp = temp + (sumnorm(A[i,,]%*%Kmat %*% ynew[onormy[j],]))^(alpha)/denomC
  }
  CsA <- CsA +temp
}
CsA <- CsA/470
CsA

#ind constant
CiA <-matrix(0,nrow=1,ncol=4,dimnames = list(c('ci'),nodesSet2))

for (i in 1:470){
  temp = rep(0,4)
  for (j in 1:k){
    temp = temp + ((A[i,,]%*%Kmat %*% ynew[onormy[j],]))^(alpha)/denomC
  }
  CiA <- CiA +t(temp)
}


CiA <- CiA/470
CiA

#Compute VaR and ES based on spectral measure
varcsa <- CsA^(1/alpha)* 0.001^(-1/alpha)
varcsa
varcia <- CiA^(1/alpha)* 0.001^(-1/alpha)
varcia
coti <- alpha/(alpha-1)*varcia
cots <- alpha/(alpha-1)*varcsa
coti
cots

#Capital alllocation
CAi <-matrix(0,nrow=1,ncol=4,dimnames = list(c('ci'),nodesSet2))
for (i in 1:470){
  temp = rep(0,4)
  for (j in 1:k){
    temp = temp + ((A[i,,]%*%Kmat %*% ynew[onormy[j],]))*(l2norm(A[i,,]%*%Kmat %*% ynew[onormy[j],]))^(alpha-1)/denomC
  }
  CAi <- CAi +t(temp)
}


CAi <- CAi/470
CAi

varca <- CsA^(1/alpha-1)*CAi* 0.001^(-1/alpha)
coteca <- CsA^(1/alpha-1)*CAi* 0.001^(-1/alpha)*alpha/(alpha-1)
varca
coteca

#################
#Minimum alpha###
#################
paramsmse2 <- matrix(0,nrow = 3,ncol = 9,dimnames = (list(c('u','scale','shape'),nodesSet1)))
paramsmse2[1,] <- paramsmse[1,]
shapemaxmse <- max(paramsmse[3,]) 
paramsmse2[2:3,1] <-fitgpd(EDPM1,paramsmse2[1,1],shape = shapemaxmse)$param
paramsmse2[1,1] <- paramsmse2[1,1]
paramsmse2[2:3,2] <-fitgpd(EF1,paramsmse2[1,2],shape = shapemaxmse)$param
paramsmse2[2:3,3] <-fitgpd(CPBP,paramsmse2[1,3]+200,shape = shapemaxmse)$param
paramsmse2[1,1] <- paramsmse2[1,3] +200
paramsmse2[2:3,4] <-fitgpd(EDPM2,paramsmse2[1,4],shape=shapemaxmse)$param
paramsmse2[2:3,5] <-fitgpd(BDSF,paramsmse2[1,5],shape=shapemaxmse)$param
paramsmse2[2:3,6] <-fitgpd(EPWS,paramsmse2[1,6],shape=shapemaxmse)$param
paramsmse2[2:3,7] <- fitgpd(DPA,paramsmse2[1,7],shape = shapemaxmse)$param
paramsmse2[2:3,8] <-fitgpd(EF2, paramsmse2[1,8]+1500,shape = shapemaxmse)$param
paramsmse2[1,8] <- paramsmse2[1,8]+1500
paramsmse2[2:3,9] <-fitgpd(IF, paramsmse2[1,9]+400,shape=shapemaxmse)$param
paramsmse2[1,9] <- paramsmse2[1,8]+400

#For the estimation, it follows exactly the same procedures as above
#The only difference is alph2,K2,Kmat2 and Kmati2 defined below
alpha2 <- 1/paramsmse2[3,1]
K2 <-mparam2[2,]
K2
Kmat2 <-diag(mparam2[2,]^(1/mparam2[1,]))
Kmati2 <- solve(Kmat2)
k <-50

######################################################
#Use LDA to compute VaR and ES estiamtes as benchmark#
######################################################
thresholdsBU <- matrix(0,nrow = 1,ncol = 4,dimnames = list(c('thresholdsBU'),nodesSet2))
thresholdsBU[1,'Retail'] <- dAMSE(Retail)$threshold
thresholdsBU[1,'Associates'] <- dAMSE(Associates)$threshold
thresholdsBU[1,'IB'] <- dAMSE(IB)$threshold
thresholdsBU[1,'Loan'] <- dAMSE(Loan)$threshold

#Estimate the parameters
paramsBU <- matrix(0,nrow = 3,ncol = 4,dimnames = (list(c('u','scale','shape'),nodesSet2)))
paramsBU[1,] <- thresholdsBU
paramsBU[2:3,1] <- fitgpd(Retail,11698797,'mle')$param
paramsBU[2:3,2] <- fitgpd(Associates,11773209,'mle')$param
paramsBU[2:3,3] <- fitgpd(IB,thresholdsBU[3],'mle')$param
paramsBU[2:3,4] <- fitgpd(Loan,thresholdsBU[4]+100,'mle')$param
paramsBU[1,1] <-11698797
paramsBU[1,2] <-11773209
paramsBU[1,4] <-thresholdsBU[4]+100

#Functions for LDA

#####################################################
## Fit a distribution

# This is the Maximum Likelihood function
MLE <- function(Params, InData, Type, Threshold)
{
  
  if (Type=="LN")
  {
    Mu1          <- Params[1]
    Sigma1       <- Params[2]
    Mu1          <- max(0.001,Mu1)
    Sigma1       <- max(0.001,Sigma1)
    Mu1          <- min(11, Mu1)
    Sigma1       <- min(3,Sigma1)
    
    TmpLikeFunc <- -sum(log(dlnorm(InData, Mu1, Sigma1) / (1-plnorm(Threshold, Mu1, Sigma1))))    
  }
  
  return(TmpLikeFunc)
}

FitDistribution <- function(data, Initialparameters, Threshold, Type)
{
  data <- data[data > Threshold]
  
  if (Type =="LN") 
  {
    shape       <- Initialparameters[1]
    scale       <- Initialparameters[2]
    
    model       <- optim(c(shape, scale), MLE, method="Nelder-Mead", control=list(abstol=1e-5), InData=data, Threshold=Threshold, Type="LN")
    if(model$par[1]<=0||model$par[2]<=0)
    {
      modelpar <- Initialparameters
    }
    else
    {
      modelpar <- model$par
    }
  }
  
  return(modelpar)
}

############################################################
## Graph the CDF

DisplayCDF1 <- function(data, params, dist, graphTitle,loc=10)
{
  #data <- data[data > Threshold]  # if you want to only display data above the threshold
  #data <- data$Loss
  data <- sort(data)  # important!
  len <- length(data)
  p <- ((1:len)-0.5)/len
  points <- matrix(c(data,p), ncol=2)
  maxdata <- max(data)
  
  if (dist=="LN")
  {
    probCurve <- plnorm(data, params[1], params[2])
    curve <- matrix(c(data, probCurve), ncol=2) 
    problim <- range(c(probCurve, p))
  }
  if (dist=="GPD")
  {
    probCurve <- pgpd(data,mu=loc,beta =  params[1],xi =  params[2])
    curve <- matrix(c(data, probCurve), ncol=2) 
    problim <- range(c(probCurve, p))
  }
  
  plot(points, pch=16, col="red", cex=0.5, xlab="Loss", ylab="Cumulative Probability", main=graphTitle, ylim=c(0,1), xlim=c(0,maxdata))
  par(new=T)
  plot(curve, type="l", col="blue", xlab="", ylab="", main="", ylim=c(0,1), xlim=c(0,maxdata))
  
}

####################################################################
## Loss Distribution
LDA2 <- function(params, mcTrials, data, dist, years, threshold)
{
  # generate the annual freq.
  #data <- data$Loss
  freq <- length(data)/years
  lossesSum <- numeric(mcTrials)
  for (t in 1:mcTrials)
  {
    # get a random sample for the annual size
    #n <- rpois(1,freq)
    
    if (dist=="LN")
    {
      # Apply a frequency correction: freq1 = freq/Prob(loss>threshold) 
      cdf <- plnorm(threshold, params[1], params[2])
      freq1 <- freq/(1 - cdf)
      n <- rpois(1,freq1)
      losses <- rlnorm(n, params[1], params[2])
    }
    if (dist=="GPD")
    {
      # Apply a frequency correction: freq1 = freq/Prob(loss>threshold) 
      cdf <- pgpd(threshold,mu=threshold,beta =  params[1],xi =  params[2])
      freq1 <- freq/(1 - cdf)
      n <- rpois(1,freq1)
      losses <- rgpd(n,mu=threshold,beta =  params[1],xi =  params[2])
    }
    
    lossesSum[t] <- sum(losses)
    
  }
  
  # 99.9% VaR and expected value
  lossesSum <- sort(lossesSum)
  var <- lossesSum[ceiling(0.999*length(lossesSum))]
  es <- mean(lossesSum[(ceiling(0.999*length(lossesSum))+1):mcTrials])
  #print(paste("VaR = ", var, sep=""))
  ev <- mean(lossesSum)
  ret <- c(var, ev,es)/1e6   # in millions
  return(ret)
}

#Compute results using weekly data aggregated by BU
paramstest123 <-fitdistr(Retail, densfun="lognormal")$estimate
DisplayCDF1(Retail, paramstest123, "LN", "Lognormal fit")
LDA2(paramstest123, mcTrials=1000000, data=Retail, dist="LN", years=10, threshold=10)

paramstest123 <-fitdistr(Associates, densfun="lognormal")$estimate
DisplayCDF1(Associates, paramstest123, "LN", "Lognormal fit")
LDA2(paramstest123, mcTrials=1000000, data=Associates, dist="LN", years=10, threshold=10)

paramstest123 <-fitdistr(IB, densfun="lognormal")$estimate
DisplayCDF1(IB, paramstest123, "LN", "Lognormal fit")
LDA2(paramstest123, mcTrials=1000000, data=IB, dist="LN", years=10, threshold=10)

paramstest123 <-fitdistr(Loan, densfun="lognormal")$estimate
DisplayCDF1(Loan, paramstest123, "LN", "Lognormal fit")
LDA2(paramstest123, mcTrials=1000000, data=Loan, dist="LN", years=10, threshold=10)

#Compute results using weekly data NOT AGGREGATED
paramstest123 <-fitdistr(Associates, densfun="lognormal")$estimate
DisplayCDF1(Associates, paramstest123, "LN", "Lognormal fit")
a <-LDA2(paramstest123, mcTrials=200, data=Associates, dist="LN", years=10, threshold=10)[1]
varaggregate <- matrix(0,nrow=4,ncol=9,dimnames = list(nodesSet2,nodesSet1))
coteaggregate <- matrix(0,nrow=4,ncol=9,dimnames = list(nodesSet2,nodesSet1))
for(i in 1:4){
  for(j in 1:9){
    if (length(L[,i,j][L[,i,j]!=0])!=0){
      paramstest123 <-fitdistr(L[,i,j][L[,i,j]!=0], densfun="lognormal")$estimate
      temp <- LDA2(paramstest123, mcTrials=200000, data=L[,i,j], dist="LN", years=10, threshold=10)
      varaggregate[i,j] <- temp[1] 
      coteaggregate[i,j] <- temp[3] 
    }
  }
}

apply(varaggregate,MARGIN = 1,sum)
apply(coteaggregate,MARGIN = 1,sum)

#Compute results using daily data aggregated by BU
#Same method but different data structure
data.day.sum <- data %>% mutate(startdate = min(Date))%>%
  mutate(numdays = difftime(Date,startdate,units="days"))%>% 
  group_by(numdays,BU,`Risk type`)%>%
  summarise(Loss =sum(Amount))


#Non aggregation
nday <- length(unique(data.day.sum$numdays))
Lday <- array(0,dim = c(nday,q,d),dimnames = list(1:nday,nodesSet2,nodesSet1))
for (i in 1:nday){
  test <- data.day.sum%>% dplyr::filter(numdays==unique(data.day.sum$numdays)[i])
  testA <- matrix(0,nrow = q,ncol = d,dimnames = list(nodesSet2,nodesSet1))
  testA[cbind(test$BU,test$`Risk type`)]<-test$Loss
  Lday[i,,] <- testA
}
Lday[1,,]

varaggregateday <- matrix(0,nrow=4,ncol=9,dimnames = list(nodesSet2,nodesSet1))
coteaggregateday <- matrix(0,nrow=4,ncol=9,dimnames = list(nodesSet2,nodesSet1))
for(i in 1:4){
  for(j in 1:9){
    if (length(Lday[,i,j][Lday[,i,j]!=0])>5){
      paramstest123 <-fitdistr(Lday[,i,j][Lday[,i,j]!=0], densfun="lognormal")$estimate
      paramstest123
      temp <- LDA2(params = paramstest123, mcTrials=1000000, data= Lday[,i,j], dist="LN", years=10, threshold=10)
      varaggregateday[i,j] <- temp[1] 
      coteaggregateday[i,j] <- temp[3] 
    }
  }
}
apply(varaggregateday,MARGIN = 1,sum)
apply(coteaggregateday,MARGIN = 1,sum)


#Using aggregation method
data.day.sumBU <- data.day.sum%>% group_by(numdays,BU)%>%
  summarise(LossBU = sum(Loss))

Retailday <- data.day.sumBU %>%filter(BU=="Retail")%>%pull(LossBU)
Associatesday <- data.day.sumBU %>%filter(BU=="Associates")%>%pull(LossBU)
IBday <- data.day.sumBU %>%filter(BU=="IB")%>%pull(LossBU)
Loanday <- data.day.sumBU %>%filter(BU=="Loan")%>%pull(LossBU)

paramstest123 <-fitdistr(Retailday, densfun="lognormal")$estimate
DisplayCDF1(Retail, paramstest123, "LN", "Lognormal fit")
LDA2(paramstest123, mcTrials=500000, data=Retailday, dist="LN", years=10, threshold=10)

paramstest123 <-fitdistr(Associatesday, densfun="lognormal")$estimate
DisplayCDF1(Associatesday, paramstest123, "LN", "Lognormal fit")
LDA2(paramstest123, mcTrials=500000, data=Associatesday, dist="LN", years=10, threshold=10)

paramstest123 <-fitdistr(IBday, densfun="lognormal")$estimate
DisplayCDF1(IBday, paramstest123, "LN", "Lognormal fit")
LDA2(paramstest123, mcTrials=500000, data=IBday, dist="LN", years=10, threshold=10)

paramstest123 <-fitdistr(Loanday, densfun="lognormal")$estimate
DisplayCDF1(Loanday, paramstest123, "LN", "Lognormal fit")
LDA2(paramstest123, mcTrials=500000, data=Loanday[Loanday>10], dist="LN", years=10, threshold=10)

##################
#Hill Estimator###
##################

#Compute the order index of thresholds
hillmse <- matrix(0,nrow=3,ncol=9,dimnames = list(c("kmse","alpha","K"),nodesSet1))
hillmse[1,1] <- which(SEDPM1==thresholdsmse[1])
hillmse[1,2] <- which(SEF1==thresholdsmse[2])
hillmse[1,3] <- which(SCPBP==thresholdsmse[3])
hillmse[1,4] <- which(SEDPM2==thresholdsmse[4])
hillmse[1,5] <- which(SBDSF==thresholdsmse[5])
hillmse[1,6] <- which(SEPWS==thresholdsmse[6])
hillmse[1,7] <- which(SDPA==thresholdsmse[7])
hillmse[1,8] <- which(SEF2==thresholdsmse[8])
hillmse[1,9] <- which(SIF==thresholdsmse[9])
hillmse

#Compute the tail index
hillmse[2,1] <- HW(SEDPM1)$tail.index
hillmse[2,2] <- HW(SEF1)$tail.index
hillmse[2,3] <- HW(SCPBP)$tail.index
hillmse[2,4] <- HW(SEDPM2)$tail.index
hillmse[2,5] <- HW(SBDSF)$tail.index
hillmse[2,6] <- HW(SEPWS)$tail.index
hillmse[2,7] <- HW(SDPA)$tail.index
hillmse[2,8] <- HW(SEF2[SEF2!=0])$tail.index
hillmse[2,9] <- HW(SIF[SIF!=0])$tail.index
hillmse
#Take mean value
hillmse[2,] <- rep(mean(hillmse[2,]),9)
hillmse

#Estimate the scaling parameter
hillmse[3,1] <- SEDPM1[hillmse[1,1]]^(hillmse[2,1])*hillmse[1,1]/n
hillmse[3,2] <- SEF1[hillmse[1,2]]^(hillmse[2,2])*hillmse[1,2]/n
hillmse[3,3] <- SCPBP[hillmse[1,3]]^(hillmse[2,3])*hillmse[1,3]/n
hillmse[3,4] <- SEDPM2[hillmse[1,4]]^(hillmse[2,4])*hillmse[1,4]/n
hillmse[3,5] <- SBDSF[hillmse[1,5]]^(hillmse[2,5])*hillmse[1,5]/n
hillmse[3,6] <- SEPWS[hillmse[1,6]]^(hillmse[2,6])*hillmse[1,6]/n
hillmse[3,7] <- SDPA[hillmse[1,7]]^(hillmse[2,7])*hillmse[1,7]/n
hillmse[3,8] <- SEF2[hillmse[1,8]]^(hillmse[2,8])*hillmse[1,8]/n
hillmse[3,9] <- SIF[hillmse[1,9]]^(hillmse[2,9])*hillmse[1,9]/n
hillmse

#The parameters are given as follows and the procedure of risk
#measure estimation are the same as aboce
alphahill <- hillmse[2,1]
Khill <-hillmse[3,]
Khill
Khillmat <-diag(hillmse[3,]^(1/alphahill))
Khillmati <- solve(Khillmat)

#Mean alpha but remove alphas<1
hillmse2 <- matrix(0,nrow=3,ncol=9,dimnames = list(c("kmse","alpha","K"),nodesSet1))
hillmse2[1,1] <- which(SEDPM1==thresholdsmse[1])
hillmse2[1,2] <- which(SEF1==thresholdsmse[2])
hillmse2[1,3] <- which(SCPBP==thresholdsmse[3])
hillmse2[1,4] <- which(SEDPM2==thresholdsmse[4])
hillmse2[1,5] <- which(SBDSF==thresholdsmse[5])
hillmse2[1,6] <- which(SEPWS==thresholdsmse[6])
hillmse2[1,7] <- which(SDPA==thresholdsmse[7])
hillmse2[1,8] <- which(SEF2==thresholdsmse[8])
hillmse2[1,9] <- which(SIF==thresholdsmse[9])

hillmse2[2,1] <- HW(SEDPM1)$tail.index
hillmse2[2,2] <- HW(SEF1)$tail.index
hillmse2[2,3] <- HW(SCPBP)$tail.index
hillmse2[2,4] <- HW(SEDPM2)$tail.index
hillmse2[2,5] <- HW(SBDSF)$tail.index
hillmse2[2,6] <- HW(SEPWS)$tail.index
hillmse2[2,7] <- HW(SDPA)$tail.index
hillmse2[2,8] <- HW(SEF2[SEF2!=0])$tail.index
hillmse2[2,9] <- HW(SIF[SIF!=0])$tail.index
hillmse2[2,] <- rep(mean(hillmse2[2,][hillmse2[2,]>1]),9)


hillmse2[3,1] <- SEDPM1[hillmse2[1,1]]^(hillmse2[2,1])*hillmse2[1,1]/n
hillmse2[3,2] <- SEF1[hillmse2[1,2]]^(hillmse2[2,2])*hillmse2[1,2]/n
hillmse2[3,3] <- SCPBP[hillmse2[1,3]]^(hillmse2[2,3])*hillmse2[1,3]/n
hillmse2[3,4] <- SEDPM2[hillmse2[1,4]]^(hillmse2[2,4])*hillmse2[1,4]/n
hillmse2[3,5] <- SBDSF[hillmse2[1,5]]^(hillmse2[2,5])*hillmse2[1,5]/n
hillmse2[3,6] <- SEPWS[hillmse2[1,6]]^(hillmse2[2,6])*hillmse2[1,6]/n
hillmse2[3,7] <- SDPA[hillmse2[1,7]]^(hillmse2[2,7])*hillmse2[1,7]/n
hillmse2[3,8] <- SEF2[hillmse2[1,8]]^(hillmse2[2,8])*hillmse2[1,8]/n
hillmse2[3,9] <- SIF[hillmse2[1,9]]^(hillmse2[2,9])*hillmse2[1,9]/n
hillmse2

alphahill2 <- hillmse2[2,1]
Khill2 <-hillmse2[3,]
Khill2
Khillmat2 <-diag(hillmse2[3,]^(1/alphahill2))
Khillmati2 <- solve(Khillmat2)
k <-50
alphahill2


#Functions that can be used to estimate VaR and ES
#varcoteca uses l2norm
varcoteca <- function(alpha,K,Kmat,Kmati,k,A){
  #Scale the data and compute order statistics.
  y <- matrix(0, nrow = n, ncol = 9)
  for (i in 1:n){
    y[i,] <-Kmati%*% RT[i,]    
  }
  
  normy <- matrix(0,nrow=1,ncol=n,dimnames = list(c('norm'),0:(n-1)))
  for (i in 1:n){
    normy[i] <- l2norm(RT[i,])
  }
  onormy <-order(normy,decreasing = TRUE)
  
  #Compute risk constants
  Amean = matrix(0,nrow = 4,ncol=9)
  for (i in 1:k){
    Amean = Amean+A[onormy[i],,]
  }
  Amean = Amean/k
  Ciind = matrix(0,nrow=4,ncol=1,dimnames = list(nodesSet2,c("Ci")))
  for (i in 1:4){
    for(j in 1:9){
      Ciind[i] = Ciind[i] + K[j]*(Amean[i,j]^(alpha))
    }
  }
  
  Aemean = rep(0,9)
  for (i in 1:9){
    for (j in 1:k){
      Aemean[i] = Aemean[i] + l2norm(A[onormy[j],,][,i])^(alpha)
    }
  }
  Aemean = Aemean/k
  Csind <- 0
  for (i in 1:9){
    Csind = Csind+K[i]*Aemean[i]
  }
  
  Cidep = matrix(0,nrow=4,ncol=1,dimnames = list(nodesSet2,c("Ci")))
  ones <-matrix(1,nrow=9,ncol=1)
  for (i in 1:4){
    for (j in 1:k){
      Cidep[i] = Cidep[i] + ((A[onormy[j],,]%*%Kmat%*%ones)[i])^(alpha)
    }
  }
  Cidep = Cidep/k
  
  Csdep = 0
  for (i in 1:k){
    Csdep  = Csdep + l2norm(A[onormy[i],,]%*%Kmat%*%ones)^(alpha)
  }
  Csdep = Csdep/k
  #use risk constants to compute risk measures
  variind <- Ciind^(1/alpha)* 0.001^(-1/alpha)
  varsind <- Csind^(1/alpha)* 0.001^(-1/alpha)
  
  cotiind <- alpha/(alpha-1)*variind
  cotsind <- alpha/(alpha-1)*varsind
  
  varidep <- Cidep^(1/alpha)* 0.001^(-1/alpha)
  varsdep <- Csdep^(1/alpha)* 0.001^(-1/alpha)
  cotidep <- alpha/(alpha-1)*varidep
  cotsdep <- alpha/(alpha-1)*varsdep
  
  ynew <- matrix(0,nrow=470,ncol = 9)
  for (i in 1:470){
    ynew[i,] <- y[i,]/l2norm(y[i,])
  }
  
  #denominator
  denomC <- 0
  for (i in 1:k){
    denomC <- denomC + ynew[onormy[i],1]^(alpha) 
  }
  denomC
  
  #Systemic constant
  CsA <-0
  for (i in 1:470){
    temp = 0
    for (j in 1:k){
      temp = temp + (l2norm(A[i,,]%*%Kmat %*% ynew[onormy[j],]))^(alpha)/denomC
    }
    CsA <- CsA +temp
  }
  CsA <- CsA/470
  
  #ind constant
  CiA <-matrix(0,nrow=1,ncol=4,dimnames = list(c('ci'),nodesSet2))
  
  for (i in 1:470){
    temp = rep(0,4)
    for (j in 1:k){
      temp = temp + ((A[i,,]%*%Kmat %*% ynew[onormy[j],]))^(alpha)/denomC
    }
    CiA <- CiA +t(temp)
  }
  CiA <- CiA/470
  
  #Calculate var
  varcsa <- CsA^(1/alpha)* 0.001^(-1/alpha)
  varcia <- CiA^(1/alpha)* 0.001^(-1/alpha)
  coti <- alpha/(alpha-1)*varcia
  cots <- alpha/(alpha-1)*varcsa
  
  #Capital allocation
  CAi <-matrix(0,nrow=1,ncol=4,dimnames = list(c('ci'),nodesSet2))
  for (i in 1:470){
    temp = rep(0,4)
    for (j in 1:k){
      temp = temp + ((A[i,,]%*%Kmat %*% ynew[onormy[j],]))*(l2norm(A[i,,]%*%Kmat %*% ynew[onormy[j],]))^(alpha-1)/denomC
    }
    CAi <- CAi +t(temp)
  }
  
  
  CAi <- CAi/470
  
  varca <- CsA^(1/alpha-1)*CAi* 0.001^(-1/alpha)
  coteca <- CsA^(1/alpha-1)*CAi* 0.001^(-1/alpha)*alpha/(alpha-1)
  
  return(list(variind,varsind,cotiind,cotsind,varidep,varsdep,cotidep,cotsdep,varcia,varcsa,coti,cots,varca,coteca))
}

#varcotecasum uses sum norm
varcotecasum <- function(alpha,K,Kmat,Kmati,k,A){
  #Scale the data and compute order statistics.
  y <- matrix(0, nrow = n, ncol = 9)
  for (i in 1:n){
    y[i,] <-Kmati%*% RT[i,]    
  }
  
  normy <- matrix(0,nrow=1,ncol=n,dimnames = list(c('norm'),0:(n-1)))
  for (i in 1:n){
    normy[i] <- sumnorm(RT[i,])
  }
  onormy <-order(normy,decreasing = TRUE)
  
  #Compute risk constants
  Amean = matrix(0,nrow = 4,ncol=9)
  for (i in 1:k){
    Amean = Amean+A[onormy[i],,]
  }
  Amean = Amean/k
  Ciind = matrix(0,nrow=4,ncol=1,dimnames = list(nodesSet2,c("Ci")))
  for (i in 1:4){
    for(j in 1:9){
      Ciind[i] = Ciind[i] + K[j]*(Amean[i,j]^(alpha))
    }
  }
  
  Aemean = rep(0,9)
  for (i in 1:9){
    for (j in 1:k){
      Aemean[i] = Aemean[i] + sumnorm(A[onormy[j],,][,i])^(alpha)
    }
  }
  Aemean = Aemean/k
  Csind <- 0
  for (i in 1:9){
    Csind = Csind+K[i]*Aemean[i]
  }
  
  Cidep = matrix(0,nrow=4,ncol=1,dimnames = list(nodesSet2,c("Ci")))
  ones <-matrix(1,nrow=9,ncol=1)
  for (i in 1:4){
    for (j in 1:k){
      Cidep[i] = Cidep[i] + ((A[onormy[j],,]%*%Kmat%*%ones)[i])^(alpha)
    }
  }
  Cidep = Cidep/k
  
  Csdep = 0
  for (i in 1:k){
    Csdep  = Csdep + sumnorm(A[onormy[i],,]%*%Kmat%*%ones)^(alpha)
  }
  Csdep = Csdep/k
  
  variind <- Ciind^(1/alpha)* 0.001^(-1/alpha)
  varsind <- Csind^(1/alpha)* 0.001^(-1/alpha)
  
  cotiind <- alpha/(alpha-1)*variind
  cotsind <- alpha/(alpha-1)*varsind
  
  varidep <- Cidep^(1/alpha)* 0.001^(-1/alpha)
  varsdep <- Csdep^(1/alpha)* 0.001^(-1/alpha)
  cotidep <- alpha/(alpha-1)*varidep
  cotsdep <- alpha/(alpha-1)*varsdep
  
  ynew <- matrix(0,nrow=470,ncol = 9)
  for (i in 1:470){
    ynew[i,] <- y[i,]/sumnorm(y[i,])
  }
  
  #denominator
  denomC <- 0
  for (i in 1:k){
    denomC <- denomC + ynew[onormy[i],1]^(alpha) 
  }
  denomC
  
  #Systemic constant
  CsA <-0
  for (i in 1:470){
    temp = 0
    for (j in 1:k){
      temp = temp + (sumnorm(A[i,,]%*%Kmat %*% ynew[onormy[j],]))^(alpha)/denomC
    }
    CsA <- CsA +temp
  }
  CsA <- CsA/470
  
  #ind constant
  CiA <-matrix(0,nrow=1,ncol=4,dimnames = list(c('ci'),nodesSet2))
  
  for (i in 1:470){
    temp = rep(0,4)
    for (j in 1:k){
      temp = temp + ((A[i,,]%*%Kmat %*% ynew[onormy[j],]))^(alpha)/denomC
    }
    CiA <- CiA +t(temp)
  }
  CiA <- CiA/470
  
  #Calculate var
  varcsa <- CsA^(1/alpha)* 0.001^(-1/alpha)
  varcia <- CiA^(1/alpha)* 0.001^(-1/alpha)
  coti <- alpha/(alpha-1)*varcia
  cots <- alpha/(alpha-1)*varcsa
  
  #Capital allocation
  CAi <-matrix(0,nrow=1,ncol=4,dimnames = list(c('ci'),nodesSet2))
  for (i in 1:470){
    temp = rep(0,4)
    for (j in 1:k){
      temp = temp + ((A[i,,]%*%Kmat %*% ynew[onormy[j],]))*(sumnorm(A[i,,]%*%Kmat %*% ynew[onormy[j],]))^(alpha-1)/denomC
    }
    CAi <- CAi +t(temp)
  }
  
  
  CAi <- CAi/470
  
  varca <- CsA^(1/alpha-1)*CAi* 0.001^(-1/alpha)
  coteca <- CsA^(1/alpha-1)*CAi* 0.001^(-1/alpha)*alpha/(alpha-1)
  
  return(list(variind,varsind,cotiind,cotsind,varidep,varsdep,cotidep,cotsdep,varcia,varcsa,coti,cots,varca,coteca))
}


#The risk constants estimates under different value of k, the number of exceedances
#Define a matrix of sorted data
SRT <- matrix(0,nrow=470,ncol=9,dimnames = list(1:470,nodesSet1))
SRT[,1] <- SEDPM1
SRT[,2] <- SEF1
SRT[,3] <- SCPBP
SRT[,4] <- SEDPM2
SRT[,5] <- SBDSF
SRT[,6] <- SEPWS
SRT[,7] <- SDPA
SRT[,8] <- SEF2
SRT[,9] <- SIF

kexceed <- 1:470
Khillexceed <- matrix(0,nrow=470,ncol=9,dimnames=list(1:470,nodesSet1))
for (i in 1:9){
  Khillexceed[,i] <- kexceed/470*(SRT[,i]^alphahill)
}

#Systemic constant
CsAhillexceed <- rep(0,470)
for (x in 1:470){
  for (i in 1:470){
    temphill = 0
    for (j in 1:k){
      temphill = temphill + (sumnorm(A[i,,]%*%diag(Khillexceed[x,]) %*% ynewhill[onormyhill[j],]))^(alphahill)/denomChill
    }
    CsAhillexceed[x] <- CsAhillexceed[x] +temphill
  }
}

CsAhillexceed <- CsAhillexceed/470
plot(c(0,CsAhillexceed^(1/alphahill)),type = "l")
abline(h= Csind^(1/alpha),col = "red")
#Plot to compare the estiamtes and lower bound
plot(log(c(0,CsAhillexceed^(1/alphahill))),type = "l",xlab="number of exceedances",ylab="",main = "Systematic Risk")
abline(h=log(Csind^(1/alpha)),col = "red")


#ind constant
CiAhillexceed <-matrix(0,nrow=470,ncol=4,dimnames = list(1:470,nodesSet2))
for (x in 1:470){
  for (i in 1:470){
    temphill = rep(0,4)
    for (j in 1:k){
      temphill = temphill + (A[i,,]%*%diag(Khillexceed[x,]) %*% ynewhill[onormyhill[j],])^(alphahill)/denomChill
    }
    CiAhillexceed[x,] <- CiAhillexceed[x,] +t(temphill)
  }
}


CiAhillexceed <- CiAhillexceed/47


#Plot to compare the estiamtes and lower bound, Figure 3.6
par(mfrow=c(2,2))
plot(log(c(0,CiAhillexceed[,1]^(1/alphahill))),type = "l",xlab= "number of exceedances",main="Retail",ylab = "")
abline(h= log(Ciind[1]^(1/alpha)),col = "red")
plot(log(c(0,CiAhillexceed[,2]^(1/alphahill))),type = "l",xlab= "number of exceedances",main = "Associates",ylab = "")
abline(h= log(Ciind[2]^(1/alpha)),col = "red")
plot(log(c(0, CiAhillexceed[,3]^(1/alphahill))),type = "l",xlab= "number of exceedances",main = "IB",ylab = "")
abline(h= log(Ciind[3]^(1/alpha)),col = "red")
plot(log(c(0,CiAhillexceed[,4]^(1/alphahill))),type = "l",xlab= "number of exceedances",main = "Retail",ylab = "")
abline(h= log(Ciind[4]^(1/alpha)),col = "red")

###################################
#Fraction Matrix Scenario analysis#
###################################
#Homogenerous model
p = mean(A!=0)#0.625591
A1 <- array(0,dim = c(n,q,d),dimnames = list(1:n,nodesSet2,nodesSet1))
L1 <- array(0,dim = c(n,q,d),dimnames = list(1:n,nodesSet2,nodesSet1))

#Compute risk constants with homogeneous fraction matrix
for (i in 1:n){
  L1[i,,] <- matrix(rbern(q*d,p),nrow=q,ncol=d,dimnames = list(nodesSet2,nodesSet1))
}
for (i in 1:n){
  A1[i,,] <- apply(L1[i,,],2,function(x) x/sum(x))
}
A1[is.nan(A1)] <- 0

CiAhillA1 <-matrix(0,nrow=470,ncol=4,dimnames = list(1:470,nodesSet2))
for (l in 1:4){
  for (i in 1:470){
    temphill = rep(0,4)
    for (j in 1:k){
      temphill = temphill + ((A1[i,,]%*%Khillmat %*% ynewhill[onormyhill[j],]))^(alphahill)/denomChill
    }
    CiAhillA1[i,] <- t(temphill)
  }
}

#Figure 4.2
boxplot(log(CiAhillA1^(1/alphahill)),ylab = TeX("$log((C^{i})^{1/\\alpha})$"))

#Compute risk constants with empricial fraction matrix
CiAhillA <-matrix(0,nrow=470,ncol=4,dimnames = list(1:470,nodesSet2))
for (l in 1:4){
  for (i in 1:470){
    temphill = rep(0,4)
    for (j in 1:k){
      temphill = temphill + ((A[i,,]%*%Khillmat %*% ynewhill[onormyhill[j],]))^(alphahill)/denomChill
    }
    CiAhillA[i,] <- t(temphill)
  }
}
#Figure 4.1
boxplot(log(CiAhillA^(1/alphahill)),ylab = TeX("$log((C^{i})^{1/\\alpha})$"))


#############
#Rasch Model#
#############
#Logistic Regression method
#construct data
rasch_data1$binary <- as.integer(rasch_data1$Freq>0)
rasch_data1$binary
#Fit the model
rasch_data.glm <- glm(binary ~ .-Var1-Var2-Freq-1-log_freq,data = rasch_data1,family = binomial)
summary(rasch_data.glm)
#Results are not reasonable

##################################
#Proposed method to fit the model#
##################################
#Define the logit function and transform the response
logit <- function(x) log(x/(1-x))
logiti <- function(x) 1-1/(1+exp(x))
rasch_data <-  as.data.frame(as.table(A[1,,]))
for (i in 2:470){
  rasch_data <- rbind(rasch_data,as.data.frame(as.table(A[i,,])))
}
#set the constants a and b 
rasch_data$log_freq <- logit(rasch_data$Freq)
rasch_data$log_freq[rasch_data$log_freq == -Inf] = -2.5
rasch_data$log_freq[rasch_data$log_freq == Inf] = 1.5
rasch_data1 <-dummy_cols(rasch_data, select_columns = 'Var1')
rasch_data1 <-dummy_cols(rasch_data1, select_columns = 'Var2')

rasch_data1$Var1_Loan[839] <- 1
#Fit the model
rasch_data1.fit <- lm(log_freq ~ .-Var1-Var2-Freq-1,data = rasch_data1)
summary(rasch_data1.fit)

#Use the estimates to compute p_ij matrix
p_rasch <- matrix(0,nrow = q,ncol = d,dimnames = list(nodesSet2,nodesSet1))
logbeta_rasch <- rasch_data1.fit$coefficients[1:4]
logdelta_rasch <- rasch_data1.fit$coefficients[5:13]
for (i in 1:4){
  for (j in 1:9){
    p_rasch[i,j] <- logiti(logbeta_rasch[i]+logdelta_rasch[j])
  }
}
p_rasch

#Calculate and plot the results
Afit <- array(0,dim = c(n,q,d),dimnames = list(1:n,nodesSet2,nodesSet1))
Lfit <- array(0,dim = c(n,q,d),dimnames = list(1:n,nodesSet2,nodesSet1))
 
for (i in 1:q){
  for(j in 1:d){
    Lfit[,i,j] <- rbern(470,p_rasch[i,j])
  }
}
 
for (i in 1:n){
  Afit[i,,] <- apply(Lfit[i,,],2,function(x) x/sum(x))
}
Afit[is.nan(Afit)] <- 0
Afit[1,,]
 
CiAhillAfit <-matrix(0,nrow=470,ncol=4,dimnames = list(1:470,nodesSet2))
for (l in 1:4){
  for (i in 1:470){
    temphill = rep(0,4)
    for (j in 1:k){
      temphill = temphill + ((Afit[i,,]%*%Khillmat %*% ynewhill[onormyhill[j],]))^(alphahill)/denomChill
    }
    CiAhillAfit[i,] <- t(temphill)
  }
}
box_rasch <-log(CiAhillAfit^(1/alphahill))
box_rasch[box_rasch==-Inf] = NaN
boxplot(box_rasch,ylab = TeX("$log((C^{i})^{1/\\alpha})$"))



#Compute the estimates corresponding to different matrix
resultsA1 <- varcoteca(alpha,K,Kmat,Kmati,k,A1)
resultsA1[9:14]

resultsA1sum <- varcotecasum(alpha,K,Kmat,Kmati,k,A1)
resultsA1sum[9:14]

resultsAfit <- varcoteca(alpha,K,Kmat,Kmati,k,Afit)
resultsAfit[9:14]

resultsAfitsum <- varcotecasum(alpha,K,Kmat,Kmati,k,Afit)
resultsAfitsum[9:14]
