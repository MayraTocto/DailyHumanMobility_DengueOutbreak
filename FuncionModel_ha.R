FuncionModel_ha<-function(Time,State,Pars)
{
  with(as.list(c(State, Pars)), {
    mu<-1/14
    rho<-2
    N1<-374102
    N2<-340959
    dS11 <- -exp(LogAlphah1) * S11 * Vi1 / ((1-theta1)*N1+theta2*N2)
    dI11 <- exp(LogAlphah1) * S11 * Vi1 / ((1-theta1)*N1+theta2*N2) - exp(LogGamma1) * I11
    dR11 <- exp(LogGamma1) * I11
    dS21 <- -exp(LogAlphah1) * S21 * Vi1 / ((1-theta1)*N1+theta2*N2)
    dI21 <- exp(LogAlphah1) * S21 * Vi1 / ((1-theta1)*N1+theta2*N2) - exp(LogGamma1) * I21
    dR21 <- exp(LogGamma1) * I21
    dVs1 <- mu * (rho*N1) - exp(LogAlphav1) * Vs1 * (I11+I21) / ((1-theta1)*N1+theta2*N2) - mu * Vs1
    dVi1 <- exp(LogAlphav1) * Vs1 * (I11+I21) / ((1-theta1)*N1+theta2*N2) - mu * Vi1
    dS22 <- -exp(LogAlphah2) * S22 * Vi2 / ((1-theta2)*N2+theta1*N1)
    dI22 <- exp(LogAlphah2) * S22 * Vi2 / ((1-theta2)*N2+theta1*N1) - exp(LogGamma2) * I22
    dR22 <- exp(LogGamma2) * I22
    dS12 <- -exp(LogAlphah2) * S12 * Vi2 / ((1-theta2)*N2+theta1*N1)
    dI12 <- exp(LogAlphah2) * S12 * Vi2 / ((1-theta2)*N2+theta1*N1) - exp(LogGamma2) * I12
    dR12 <- exp(LogGamma2) * I12
    dVs2 <- mu * (rho*N2) - exp(LogAlphav2) * Vs2 * (I22+I12) / ((1-theta2)*N2+theta1*N1) - mu * Vs2
    dVi2 <- exp(LogAlphav2) * Vs2 * (I22+I12) / ((1-theta2)*N2+theta1*N1) - mu * Vi2
    dC1  <- Delta * ( exp(LogAlphah1) * S11 * Vi1 / ((1-theta1)*N1+theta2*N2) + exp(LogAlphah2) * S12 * Vi2 / ((1-theta2)*N2+theta1*N1))
    dC2  <- Delta * ( exp(LogAlphah2) * S22 * Vi2 / ((1-theta2)*N2+theta1*N1) + exp(LogAlphah1) * S21 * Vi1 / ((1-theta1)*N1+theta2*N2))
    list(c(dS11, dI11, dR11, dS21, dI21, dR21, dVs1, dVi1, dS22, dI22, dR22, dS12, dI12, dR12, dVs2, dVi2, dC1, dC2)) })
}
