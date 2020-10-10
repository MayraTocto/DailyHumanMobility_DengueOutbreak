PosteriorPredictiveDensity<-function(MatrixSolInfected,vecvalInfected,
                                     LengthTimeObs)
{
  n<-dim(MatrixSolInfected)[1]
  m<-length(vecvalInfected)
  DensityPoissonPred_temp<-c()
  DensityPoissonPred<-matrix(NA,nrow=m,ncol=LengthTimeObs)
  MatvecvalInfected<-matrix(rep(vecvalInfected,n), ncol = n)
  
  for(j in 1:LengthTimeObs)
  {
    for(i in 1:m)
    {
      DensityPoissonPred_temp[i]<-mean(dpois(MatvecvalInfected[i,],
                                      lambda=MatrixSolInfected[,j]))
    }
    DensityPoissonPred[,j]<-DensityPoissonPred_temp/sum(DensityPoissonPred_temp)
  }
  return(DensityPoissonPred)
}