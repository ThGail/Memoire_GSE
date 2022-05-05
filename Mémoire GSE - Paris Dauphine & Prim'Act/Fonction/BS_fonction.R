# définition des fonctions utiles du modèle Black&Scholes

# simulation prix actif en version Wasicek
PA_BS_Vas_FF_sim <- function(N, t, S0, param_Vas, param_actif, r0=TauxZC[1]){
  nb_mat <- length(t)
  mat <- t(matrix(rep(t,N),nrow=nb_mat))
  r <- taux_Vas_sim(N, t, param_Vas, r0)
  return(S0 * exp((r - 0.5 * param_actif^2) * t + param_actif * rnorm(N, 0, sqrt(t))))
}
# simulation prix actif sur plusieurs temps
PA_BS_Vas_FF_sim.t <- Vectorize(PA_BS_Vas_FF_sim,"t")

# simulation prix actif en version HW
PA_BS_HW_FF_sim <- function(N, t, S0, param_HW, param_actif, r0=TauxZC[1]){
  nb_mat <- length(t)
  mat <- t(matrix(rep(t,N),nrow=nb_mat))
  r <- taux_HW_sim(N, t, param_HW, r0)
  return(S0 * exp((r - 0.5 * param_actif^2) * t + param_actif * rnorm(N, 0, sqrt(t))))
}
# simulation prix actif sur plusieurs temps
PA_BS_HW_FF_sim.t <- Vectorize(PA_BS_HW_FF_sim,"t")

# simulation prix actif pour le test de martingalité
PAact_BS_FF_sim <- function(N,t,S0,param_Vas, param_actif, r0=TauxZC[1]){
  r <- TZC_Vas_FF_calibrage(t, param_Vas, r0)
  return(exp(-r*t)*S0*exp((r-0.5*param_actif^2)*t+param_actif*rnorm(N,0,sqrt(t))))
}
PAact_BS_FF_sim.t <- Vectorize(PAact_BS_FF_sim,"t")

# simulation prix de l'option Call actualisé pour une maturité donnée
call_BS_sim <- function(N,t,TT,S0,param_Vas,param_actif,r0=TauxZC[1]){
  tau <- TT - t
  nb_mat <- length(tau)
  theta <- t(matrix(rep(tau,N),nrow=nb_mat))
  r <- taux_Vas_sim.t(N, t, param_Vas, r0)
  return(exp(-r*theta)*pmax(PA_BS_Vas_FF_sim(N, TT, S0, param_Vas, param_actif, r0)-K,0))
}
call_BS_sim.T <- Vectorize(call_BS_sim,"TT")

call_BS_FF <- function(TT,S0,K,param_Vas,param_actif,r0=TauxZC[1]) {
  r <- TZC_Vas_FF_calibrage(TT,param_Vas,r0)
  d1 <- (log(S0/K) + (r + 0.5*param_actif^2)*TT)/(param_actif*sqrt(TT))
  d2 <- d1 - param_actif*sqrt(TT)
  return(S0*pnorm(d1) - K*exp(-r*TT)*pnorm(d2))
}
