
##########################################################################
#                                                                        #
# This program calculate the confidence interval of the implementation   #
# of the model and the posterior predictive distribution of the          #
# manuscript "Effect of daily periodic human movement on dengue          #
# dynamics: the case of the 2010 outbreak in Hermosillo, Mexico."  We    #
# employ txt files  which contain 50000 solutions of new weekly infected #
# individuals obtained solving the mathematical model with 50000 samples # 
# of the model parameters for each Tl=1, Tl=0.75 and Tl=0.5.             #                                            #
#                                                                        #
##########################################################################

###########################################################
###########################################################
#                     MODEL                               #
###########################################################
###########################################################

#############
# PACKAGES  #
#############
library(MASS)
library(deSolve)

#############
# FUNCTIONS #
#############
source("Model_la.R",local = FALSE)
source("Model_ha.R",local = FALSE)
source("SolModel.R",local = FALSE)

##############
#  SCENARIO  #
##############

Population_n<-374102
Population_s<-340959

rho<-2
NumInfectedTime0_n<-1
NumInfectedTime0_s<-1
NumSuceptibleTime0_n<-Population_n-NumInfectedTime0_n
NumSuceptibleTime0_s<-Population_s-NumInfectedTime0_s
NumRecoveredTime0_n<-0
NumRecoveredTime0_s<-0
NumSuceptibleVecTime0_n<-rho*Population_n
NumSuceptibleVecTime0_s<-rho*Population_s
NumInfectedVecTime0_n<-0
NumInfectedVecTime0_s<-0
NumCumInfectedTime0_n<-NumInfectedTime0_n
NumCumInfectedTime0_s<-NumInfectedTime0_s
t0<-32
T0<-40
TimeObs<-seq(t0+1,T0,1)

yini_n<-c(S1=NumSuceptibleTime0_n,I1=NumInfectedTime0_n,R1=NumRecoveredTime0_n,
          Vs1=NumSuceptibleVecTime0_n,Vi1=NumInfectedVecTime0_n)
yini_s<-c(S2=NumSuceptibleTime0_s,I2=NumInfectedTime0_s,R2=NumRecoveredTime0_s,
          Vs2=NumSuceptibleVecTime0_s,Vi2=NumInfectedVecTime0_s)
yiniC<-c(C1=NumCumInfectedTime0_n,C2=NumCumInfectedTime0_s)
Yinitial<-c(yini_n,yini_s,yiniC)

Tl<-0.5
Period<-c()
Period[1]<-0
dt<-0.01
sizePeriod<-(length(TimeObs))*14+1
CumTimeWeekly<-seq(0,(length(TimeObs))*7,7)
for (i in 2:sizePeriod)
{
  if((i%%2)==0) 
  {
    Period[i] <- Period[i-1]+Tl
  }
  else
  {
    Period[i] <- Period[i-1]+(1.0-Tl)
  }
}  

###################
#    SIMULATION   #
###################

VecPar<-c(LogAlphah1=-1.9,LogGamma1=-1.3,LogAlphav1=-1.3,
          LogAlphah2=-2.2,LogGamma2=-1.4,LogAlphav2=-1.4,
          Delta=0.4,theta1=0.45,theta2=0.1)

OutModel<-SolModel(VecPar,Yinitial,Period,dt)
CumInfectedPatch1<-OutModel[,12][OutModel[,1] %in% CumTimeWeekly]
CumInfectedPatch2<-OutModel[,13][OutModel[,1] %in% CumTimeWeekly]
InfectedPatch1<-CumInfectedPatch1[2:(length(TimeObs)+1)]-CumInfectedPatch1[1:length(TimeObs)]
InfectedPatch2<-CumInfectedPatch2[2:(length(TimeObs)+1)]-CumInfectedPatch2[1:length(TimeObs)]

par(mfrow=c(1,2))
plot(TimeObs,InfectedPatch1,type="l",col=2,xlab="Time (Weeks)",ylab="Infected of Patch 1",
     lwd=2)
plot(TimeObs,InfectedPatch2,type="l",col=3,xlab="Time (Weeks)",ylab="Infected of Patch 2",
     lwd=2)


###########################################################
###########################################################
#             Posterior predictive distribution           #
###########################################################
###########################################################

#############
# FUNCTIONS #
#############

source("PosteriorPredictiveDensity.R",local = FALSE)

#################
# READING FILES #
#################

MatrixSolInfected1<-list()
MatrixSolInfected2<-list()

MatrixSolInfected1[[1]] <- read.table("SolInfectedPatch1_Tl1.txt",header=FALSE,sep=" ")
MatrixSolInfected2[[1]] <- read.table("SolInfectedPatch2_Tl1.txt",header=FALSE,sep=" ")

MatrixSolInfected1[[2]] <- read.table("SolInfectedPatch1_Tl075.txt",header=FALSE,sep=" ")
MatrixSolInfected2[[2]] <- read.table("SolInfectedPatch2_Tl075.txt",header=FALSE,sep=" ")

MatrixSolInfected1[[3]] <- read.table("SolInfectedPatch1_Tl05.txt",header=FALSE,sep=" ")
MatrixSolInfected2[[3]] <- read.table("SolInfectedPatch2_Tl05.txt",header=FALSE,sep=" ")

########################################################
# Calculation of the posterior predictive distribution #  
########################################################

t0<-32
T0<-40
TimeObs<-seq(t0+1,T0,1)
vecvalInfected1<-seq(0,420,1)
vecvalInfected2<-seq(0,85,1)

DensityPoissonPred1<-list()
DensityPoissonPred2<-list()

DensityPoissonPred1[[1]]<-PosteriorPredictiveDensity(MatrixSolInfected1[[1]],
                                                             vecvalInfected1,
                                                             length(TimeObs))
DensityPoissonPred2[[1]]<-PosteriorPredictiveDensity(MatrixSolInfected2[[1]],
                                                             vecvalInfected2,
                                                             length(TimeObs))
DensityPoissonPred1[[2]]<-PosteriorPredictiveDensity(MatrixSolInfected1[[2]],
                                                             vecvalInfected1,
                                                             length(TimeObs))
DensityPoissonPred2[[2]]<-PosteriorPredictiveDensity(MatrixSolInfected2[[2]],
                                                             vecvalInfected2,
                                                             length(TimeObs))
DensityPoissonPred1[[3]]<-PosteriorPredictiveDensity(MatrixSolInfected1[[3]],
                                                             vecvalInfected1,
                                                             length(TimeObs))
DensityPoissonPred2[[3]]<-PosteriorPredictiveDensity(MatrixSolInfected2[[3]],
                                                             vecvalInfected2,
                                                             length(TimeObs))
DistributionPoissonPred1<-list()
DistributionPoissonPred2<-list()

DistributionPoissonPred1[[1]]<-apply(DensityPoissonPred1[[1]],2,cumsum)
DistributionPoissonPred2[[1]]<-apply(DensityPoissonPred2[[1]],2,cumsum)
DistributionPoissonPred1[[2]]<-apply(DensityPoissonPred1[[2]],2,cumsum)
DistributionPoissonPred2[[2]]<-apply(DensityPoissonPred2[[2]],2,cumsum)
DistributionPoissonPred1[[3]]<-apply(DensityPoissonPred1[[3]],2,cumsum)
DistributionPoissonPred2[[3]]<-apply(DensityPoissonPred2[[3]],2,cumsum)

#######################################
# 95% PREDICTIVE CONFIDENCE INTERVAL  #
#######################################

CI95_1<-list()
CI95_2<-list()
CI95_1[[1]]<-matrix(NA,2,length(TimeObs)) 
CI95_1[[2]]<-matrix(NA,2,length(TimeObs)) 
CI95_1[[3]]<-matrix(NA,2,length(TimeObs)) 
CI95_2[[1]]<-matrix(NA,2,length(TimeObs)) 
CI95_2[[2]]<-matrix(NA,2,length(TimeObs)) 
CI95_2[[3]]<-matrix(NA,2,length(TimeObs)) 

for(i in 1:3)
{
  for(j in 1:(dim(DistributionPoissonPred1[[1]])[2]))
  {
    CI95_1[[i]][,j]<-c(vecvalInfected1[which(DistributionPoissonPred1[[i]][,j]>=0.025)[1]],
                       vecvalInfected1[which(DistributionPoissonPred1[[i]][,j]>=0.975)[1]])
  }
}
for(i in 1:3)
{
  for(j in 1:(dim(DistributionPoissonPred2[[1]])[2]))
  {
    CI95_2[[i]][,j]<-c(vecvalInfected2[which(DistributionPoissonPred2[[i]][,j]>=0.025)[1]],
                       vecvalInfected2[which(DistributionPoissonPred2[[i]][,j]>=0.975)[1]])
  }
}

#PLOT
col=c("black","blue","red")

par(mfrow=c(1,2))
par(mar = c(3.2, 2.3, 0.3, 0.3),mgp = c(2.2, 1, 0))
plot(x=c(32.5,40.5),y=c(0,425), type="n",xlab="Time (Weeks)",ylab="",
     xaxt="n",yaxt="n")
axis(1, at=seq(33,40,2))
axis(2, at=seq(0,400,100))
for(i in 1:3)
{
  for(j in 1:length(TimeObs))
  {
    arrows(x0=TimeObs[j],y0=CI95_1[[i]][1,j],x1=TimeObs[j],
           y1=CI95_1[[i]][2,j],col=col[i],lwd=1,code=3,angle=90,
           length=0.05,lty=1)
  }
}
text(32.7,415,"A",cex=1.5)
legend(33,425,c(expression(paste("95%CI, ",T[l]==1)),
                expression(paste("95%CI, ",T[l]==0.75)),
                expression(paste("95%CI, ",T[l]==0.5))),
       lty=c(1,1,1),col=col,lwd=c(1,1,1),
       pch=c(NA,NA,NA))

par(mar = c(3.2, 2.3, 0.3, 0.3),mgp = c(2.2, 1, 0))
plot(x=c(32.5,40.5),y=c(0,65), type="n",xlab="Time (Weeks)",ylab="",
     xaxt="n",yaxt="n")
axis(1, at=seq(33,40,2))
axis(2, at=seq(0,60,15))
for(i in 1:3)
{
  for(j in 1:length(TimeObs))
  {
    arrows(x0=TimeObs[j],y0=CI95_2[[i]][1,j],x1=TimeObs[j],
           y1=CI95_2[[i]][2,j],col=col[i],lwd=1,code=3,angle=90,
           length=0.05,lty=1)
  }
}
text(32.7,63,"B",cex=1.5)
