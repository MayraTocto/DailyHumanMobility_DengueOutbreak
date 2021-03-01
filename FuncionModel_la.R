FuncionModel_la<-function(Time,State,Pars)
{
  with(as.list(c(State, Pars)), {
    mu<-1/14
    rho<-2
    N1<-374102
    N2<-340959
    dS1 <- -exp(LogAlphah1) * S1 * Vi1 / N1
    dI1 <- exp(LogAlphah1) * S1 * Vi1 / N1 - exp(LogGamma1) * I1
    dR1 <- exp(LogGamma1) * I1
    dVs1 <- mu * (rho*N1) - exp(LogAlphav1) * Vs1 * I1 / N1 - mu * Vs1
    dVi1 <- exp(LogAlphav1) * Vs1 * I1 / N1 - mu * Vi1
    dS2 <- -exp(LogAlphah2) * S2 * Vi2 / N2
    dI2 <- exp(LogAlphah2) * S2 * Vi2 / N2 - exp(LogGamma2) * I2
    dR2 <- exp(LogGamma2) * I2
    dVs2 <- mu * (rho*N2) - exp(LogAlphav2) * Vs2 * I2 / N2 - mu * Vs2
    dVi2 <- exp(LogAlphav2) * Vs2 * I2 / N2 - mu * Vi2
    dC1 <- Delta*exp(LogAlphah1) * S1 * Vi1 / N1
    dC2 <- Delta*exp(LogAlphah2) * S2 * Vi2 / N2
    list(c(dS1, dI1, dR1, dVs1, dVi1, dS2, dI2, dR2, dVs2, dVi2, dC1, dC2)) })
}
