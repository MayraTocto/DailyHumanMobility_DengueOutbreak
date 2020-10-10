##########################################################################
#                                                                        #
# This program calculate the confidence interval of the posterior        #
# predictive distribution of the manuscript "Effect of daily periodic    #
# human movement on dengue dynamics: the case of the 2010 outbreak in    #
# Hermosillo, Mexico."  We employ txt files  which contain 50000         #
# solutions of new weekly infected individuals obtained solving the      #
# mathematical model with 50000 samples of the model parameters for each #
# Tl=1, Tl=0.75 and Tl=0.5.                                              #
#                                                                        #
##########################################################################

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

