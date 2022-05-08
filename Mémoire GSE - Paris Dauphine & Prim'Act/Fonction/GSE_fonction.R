### Définition fonctions utiles ###################################
# la fonction qui permet de simuler les taux instantanées de Vasicek
Vas_W <- function(t, param, W, r0=TauxZC[1]) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(r0 * exp(-a * t) + b * (1 - exp(-a * t)) + sigma * sqrt((1 - exp(-2 * a * t)) / (2 * a)) * W)
}

# transformation du taux instantanée en TZC de Vasicek
tauxToTZC_Vas <- function(t,TT,taux_t,param){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  Ri <- b-sigma^2/(2*a^2)
  theta <- TT-t
  return(Ri - ((Ri-taux_t)*(1-exp(-a*theta))-sigma^2/(4*a^2)*(1-exp(-a*theta))^2)/(a*theta))
}

tauxToPZC_Vas <- function(t, TT, taux_t, param){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  theta <- TT-t
  return(exp(-b * theta) * exp(-(taux_t - b) * (1 - exp(-a * theta)) / a + 0.5 * (sigma^2 * theta / a^2 - sigma^2 / a^3 * (1 - exp(-a * theta)) - sigma^2 / (2 * a^3) * (1 - exp(-a * theta))^2)))
}

# la fonction qui permet de simuler les taux instantanées de HW
HW_W <- function(t, param, W, r0=TauxZC[1]){
  a <- param[1]
  sigma <- param[2]
  return(r0*exp(-a*t)+beta_HW(t,param)-beta_HW(0,param)*exp(-a*t) + sqrt((sigma^2)/(2*a)*(1-exp(-2*a*t)))*W)
}

# transformation du taux instantanée en TZC de HW
tauxToTZC_HW <- function(t,TT,taux_t,param){
  PZC <- A_HW(t,TT,param)*exp(-B_HW(t,TT,param)*taux_t)
  return(-log(PZC)/(TT-t))
}

tauxToPZC_HW <- function(t,TT,taux_t,param){
  return(A_HW(t,TT,param)*exp(-B_HW(t,TT,param)*taux_t))
}

# la fonction qui simule le prix des actifs du BS
BS_Vas_W <- function(t, param_Vas, param_actif, W_vas, W_act, S0=1, r0=TauxZC[1]) {
  return(S0 * exp((Vas_W(t, param_Vas, W_vas, r0) - 0.5 * param_actif^2) * t + param_actif * sqrt(t) * W_act))
}

BS_HW_W <- function(t, param_HW, param_actif, W_HW, W_act, S0=1, r0=TauxZC[1]) {
  return(S0 * exp((HW_W(t, param_HW, W_HW, r0) - 0.5 * param_actif^2) * t + param_actif * sqrt(t) * W_act))
}





### GSE complet avec risque de crédit CIR en utilisant Vasicek pour le generateur de taux sans risque ####
GSE_Vas_CIR <- function(N, t, TT, param_Vas, param_action, param_immo, param_CIR, L, LGD, S0=1,r0=TauxZC[1]){
  
  # simulation de gaussiennes corrélées
  # on simule 3 car : 1 taux sans risque + 2 actifs + 6 crédits
  Z <- matrix(rnorm(N*3), nrow = N)
  W <- L %*% t(Z)
  W.vas <- W[1,]
  W.action <- W[2,]
  W.immo <- W[3,]
  
  # calcul des grandeurs de base
  vas <- Vas_W(t, param_Vas, W.vas, r0)
  action <- BS_Vas_W(t, param_Vas, param_action, W.vas, W.action,S0,r0)
  immo <- BS_Vas_W(t, param_Vas, param_immo, W.vas, W.immo,S0,r0)
  
  TZCvas <- tauxToTZC_Vas(t,TT,vas,param_Vas)
  PZCvas <- tauxToPZC_Vas(t,TT,vas,param_Vas)
  
  # déflateur
  defl <- exp(-(TT-t)*TZCvas)
  
  # calcul des spreads de crédit / proba de survie / proba de défaut / PZCr
  # on simule 6 autres pour les chacunes de rating (transformé en proba de survie)
  # les crédits sont indépendants les uns des autres (hypothèses)
  ps.AAA <- survie_CIR_sim(N, t, TT, param_CIR$AAA)
  ps.AA <- survie_CIR_sim(N, t, TT, param_CIR$AA)
  ps.A <- survie_CIR_sim(N, t, TT, param_CIR$A)
  ps.BBB <- survie_CIR_sim(N, t, TT, param_CIR$BBB)
  ps.BB <- survie_CIR_sim(N, t, TT, param_CIR$BB)
  ps.B <- survie_CIR_sim(N, t, TT, param_CIR$B)
  
  # calcul du PZR risqué
  pzcr.AAA <- PZCvas * (1-LGD+LGD*ps.AAA)
  pzcr.AA <- PZCvas * (1-LGD+LGD*ps.AA)
  pzcr.A <- PZCvas * (1-LGD+LGD*ps.A)
  pzcr.BBB <- PZCvas * (1-LGD+LGD*ps.BBB)
  pzcr.BB <- PZCvas * (1-LGD+LGD*ps.BB)
  pzcr.B <- PZCvas * (1-LGD+LGD*ps.B)
  
  return(list(deflateur_c=defl,PZC_Vas=PZCvas,TZC_Vas=TZCvas,ind_action=action,ind_immo=immo,
              ps_AAA=ps.AAA,ps_AA=ps.AA,ps_A=ps.A,ps_BBB=ps.BBB,ps_BB=ps.BB,ps_B=ps.B,
              PZCr_AAA=pzcr.AAA,PZCr_AA=pzcr.AA,PZCr_A=pzcr.A,PZCr_BBB=pzcr.BBB,PZCr_BB=pzcr.BB,PZCr_B=pzcr.B))
}





### GSE complet avec risque de crédit JLT en utilisant Vasicek pour le generateur de taux sans risque ####
GSE_Vas_JLT <- function(N, t, TT, param_Vas, param_action, param_immo, param_JLT, L, LGD, S0=1,r0=TauxZC[1]){
  
  Z <- matrix(rnorm(N*3), nrow = N)
  W <- L %*% t(Z)
  W.vas <- W[1,]
  W.action <- W[2,]
  W.immo <- W[3,]
  
  vas <- Vas_W(t, param_Vas, W.vas, r0)
  action <- BS_Vas_W(t, param_Vas, param_action, W.vas, W.action,S0,r0)
  immo <- BS_Vas_W(t, param_Vas, param_immo, W.vas, W.immo,S0,r0)
  
  TZCvas <- tauxToTZC_Vas(t,TT,vas,param_Vas)
  PZCvas <- tauxToPZC_Vas(t,TT,vas,param_Vas)
  
  defl <- exp(-(TT-t)*TZCvas)
  
  ps.AAA <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 1)
  ps.AA <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 2)
  ps.A <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 3)
  ps.BBB <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 4)
  ps.BB <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 5)
  ps.B <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 6)
  
  pzcr.AAA <- PZCvas * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 1))
  pzcr.AA <- PZCvas * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 2))
  pzcr.A <- PZCvas * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 3))
  pzcr.BBB <- PZCvas * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 4))
  pzcr.BB <- PZCvas * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 5))
  pzcr.B <- PZCvas * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 6))
  
  return(list(deflateur_c=defl,PZC_Vas=PZCvas,TZC_Vas=TZCvas,ind_action=action,ind_immo=immo,
              ps_AAA=ps.AAA,ps_AA=ps.AA,ps_A=ps.A,ps_BBB=ps.BBB,ps_BB=ps.BB,ps_B=ps.B,
              PZCr_AAA=pzcr.AAA,PZCr_AA=pzcr.AA,PZCr_A=pzcr.A,PZCr_BBB=pzcr.BBB,PZCr_BB=pzcr.BB,PZCr_B=pzcr.B))
}





### GSE complet avec risque de crédit CIR en utilisant Hull&White pour le generateur de taux sans risque ####
GSE_HW_CIR <- function(N, t, TT, param_HW, param_action, param_immo, param_CIR, L, LGD, S0=1,r0=TauxZC[1]){
  
  Z <- matrix(rnorm(N*3), nrow = N)
  W <- L %*% t(Z)
  W.hw <- W[1,]
  W.action <- W[2,]
  W.immo <- W[3,]
  
  hw <- HW_W(t, param_HW, W.hw, r0)
  action <- BS_HW_W(t, param_HW, param_action, W.hw, W.action,S0,r0)
  immo <- BS_HW_W(t, param_HW, param_immo, W.hw, W.immo,S0,r0)
  
  PZChw <- tauxToPZC_HW(t,TT,hw,param_HW)
  TZChw <- -log(PZChw)/(TT-t)
  
  defl <- exp(-(TT-t)*TZChw)
  
  ps.AAA <- survie_CIR_sim(N, t, TT, param_CIR$AAA)
  ps.AA <- survie_CIR_sim(N, t, TT, param_CIR$AA)
  ps.A <- survie_CIR_sim(N, t, TT, param_CIR$A)
  ps.BBB <- survie_CIR_sim(N, t, TT, param_CIR$BBB)
  ps.BB <- survie_CIR_sim(N, t, TT, param_CIR$BB)
  ps.B <- survie_CIR_sim(N, t, TT, param_CIR$B)
  
  pzcr.AAA <- PZChw * (1-LGD+LGD*ps.AAA)
  pzcr.AA <- PZChw * (1-LGD+LGD*ps.AA)
  pzcr.A <- PZChw * (1-LGD+LGD*ps.A)
  pzcr.BBB <- PZChw * (1-LGD+LGD*ps.BBB)
  pzcr.BB <- PZChw * (1-LGD+LGD*ps.BB)
  pzcr.B <- PZChw * (1-LGD+LGD*ps.B)
  
  return(list(deflateur_c=defl,PZC_HW=PZChw,TZC_HW=TZChw,ind_action=action,ind_immo=immo,
              ps_AAA=ps.AAA,ps_AA=ps.AA,ps_A=ps.A,ps_BBB=ps.BBB,ps_BB=ps.BB,ps_B=ps.B,
              PZCr_AAA=pzcr.AAA,PZCr_AA=pzcr.AA,PZCr_A=pzcr.A,PZCr_BBB=pzcr.BBB,PZCr_BB=pzcr.BB,PZCr_B=pzcr.B))
}





### GSE complet avec risque de crédit JLT en utilisant Hull&White pour le generateur de taux sans risque ####
GSE_HW_JLT <- function(N, t, TT, param_HW, param_action, param_immo, param_JLT, L, LGD, S0=1,r0=TauxZC[1]){
  
  Z <- matrix(rnorm(N*3), nrow = N)
  W <- L %*% t(Z)
  W.hw <- W[1,]
  W.action <- W[2,]
  W.immo <- W[3,]
  
  hw <- HW_W(t, param_HW, W.hw, r0)
  action <- BS_HW_W(t, param_HW, param_action, W.hw, W.action,S0,r0)
  immo <- BS_HW_W(t, param_HW, param_immo, W.hw, W.immo,S0,r0)
  
  PZChw <- tauxToPZC_HW(t,TT,hw,param_HW)
  TZChw <- -log(PZChw)/(TT-t)
  
  defl <- exp(-(TT-t)*TZChw)
  
  ps.AAA <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 1)
  ps.AA <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 2)
  ps.A <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 3)
  ps.BBB <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 4)
  ps.BB <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 5)
  ps.B <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 6)
  
  pzcr.AAA <- PZChw * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 1))
  pzcr.AA <- PZChw * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 2))
  pzcr.A <- PZChw * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 3))
  pzcr.BBB <- PZChw * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 4))
  pzcr.BB <- PZChw * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 5))
  pzcr.B <- PZChw * (1-LGD*proba_defaut_i_JLT(N, t, TT, param_JLT, M, D, 6))
  
  return(list(deflateur_c=defl,PZC_HW=PZChw,TZC_HW=TZChw,ind_action=action,ind_immo=immo,
              ps_AAA=ps.AAA,ps_AA=ps.AA,ps_A=ps.A,ps_BBB=ps.BBB,ps_BB=ps.BB,ps_B=ps.B,
              PZCr_AAA=pzcr.AAA,PZCr_AA=pzcr.AA,PZCr_A=pzcr.A,PZCr_BBB=pzcr.BBB,PZCr_BB=pzcr.BB,PZCr_B=pzcr.B))
}
GSE_HW_JLT.t <- Vectorize(GSE_HW_JLT,"t")
GSE_HW_JLT.T <- Vectorize(GSE_HW_JLT,"TT")
