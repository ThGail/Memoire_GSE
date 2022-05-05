rm(list = ls())
library(readxl)
library(dfoptim)
library(writexl)
library(matrixStats)



chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = paste(chemin,"/Fonction", sep="")
chemin_cal = paste(chemin,"/Calibrage", sep="")
setwd(chemin)

chemin_export = "C:/Users/mengr/Desktop/Export"


### Importation des donnees ###################################################
source(paste(chemin_cal,"/Vasicek_calibrage.R",sep=""))
source(paste(chemin_cal,"/HW_calibrage.R",sep=""))
source(paste(chemin_cal,"/BS_calibrage.R",sep=""))
source(paste(chemin_cal,"/CIR_calibrage.R",sep=""))
source(paste(chemin_cal,"/JLT_calibrage.R",sep=""))
source(paste(chemin_fct,"/GSE_fonction.R",sep=""))



sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)



Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
PrixZC <- exp(-Maturite*TauxZC)
TauxIflZC <- as.numeric(sheetDGlo$...5[7:106])
PrixCaps <- as.numeric(sheetCali$...5[3:22])
SpreadMarket <- as.numeric(sheetDGlo$...15[3:9])
Maturite10 <- Maturite[1:10]
Maturite15 <- Maturite[1:15]
Maturite30 <- Maturite[1:30]
Maturite50 <- Maturite[1:50]


### Parametrisation #######################################
# matrice de correlaion entre taux_court, action_glo et immo
Corr <- matrix(c(1,0.1,0.05,0.1,1,-0.3,0.05,-0.3,1),nrow=3)
# d?compostion de Cholesky :
L <- t(chol(Corr, pivot = FALSE))



N = 1000; TT = 50; LGD = 0.3; S0=1;



### HLG-HW_20220224_11h50m19s_QDeflateur ##############
defla_HW <- c(rep(1,N)) # pour T=0
defla_HW <- cbind(defla_HW,deflateur_HW.T(N,1/360,Maturite50,paramHW))
defla_HW.df <- as.data.frame(defla_HW)
names(defla_HW.df) <- c(0,Maturite50)
write_xlsx(defla_HW.df,paste(chemin_export,"/HLG-HW_20220224_11h50m19s_QDeflateur", ".xlsx", sep=""))

defla_Vas <- c(rep(1,N)) # pour T=0
defla_Vas <- cbind(defla_Vas,deflateur_Vas.T(N,1/360,Maturite50,paramVas))
defla_Vas.df <- as.data.frame(defla_Vas)
names(defla_Vas.df) <- c(0,Maturite50)
write_xlsx(defla_Vas.df,paste(chemin_export,"/HLG-HW_20220224_11h50m19s_QDeflateur", ".xlsx", sep=""))

### HLG-HW_20220224_11h50m19sQActionsGlobales ####
### HLG-HW_20220224_11h50m19sQImmobilier #################
action_HW <- c(rep(S0,N)) # pour T=0
immo_HW <- c(rep(S0,N))
for (t in Maturite50){
  GSE <- GSE_HW_JLT(N, t, NA, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  action_HW <- cbind(action_HW,GSE$ind_action)
  immo_HW <- cbind(immo_HW,GSE$ind_immo)
}
action_HW.df <- as.data.frame(action_HW)
immo_HW.df <- as.data.frame(immo_HW)
names(action_HW.df) <- c(0,Maturite50)
names(immo_HW.df) <- c(0,Maturite50)
write.csv2(action_HW.df,paste(chemin_export,"/HLG-HW_20220224_11h50m19sQActionsGlobales", ".csv", sep=""),row.names = FALSE)
write.csv2(immo_HW.df,paste(chemin_export,"/HLG-HW_20220224_11h50m19sQImmobilier", ".csv", sep=""),row.names = FALSE)

action_Vas <- c(rep(S0,N)) # pout T=0
immo_Vas <- c(rep(S0,N))
for (t in Maturite50){
  GSE <- GSE_Vas_JLT(N, t, NA, paramVas, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  action_Vas <- cbind(action_Vas,GSE$ind_action)
  immo_Vas <- cbind(immo_Vas,GSE$ind_immo)
}
action_Vas.df <- as.data.frame(action_Vas)
immo_Vas.df <- as.data.frame(immo_Vas)
names(action_Vas.df) <- c(0,Maturite50)
names(immo_Vas.df) <- c(0,Maturite50)
write_xlsx(action_Vas.df,paste(chemin_export,"/HLG-HW_20220224_11h50m19sQActionsGlobales", ".xlsx", sep=""))
write_xlsx(immo_Vas.df,paste(chemin_export,"/HLG-HW_20220224_11h50m19sQImmobilier", ".xlsx", sep=""))


### Courbe_Taux_dans_XX_an_HLG-HW_20220224_11h50m19s_QZCsto _ #####################################

for(t in 0:50){
  print(paste("Année de projection : ", t))
  CourbeZC <- c()
  
  #Maturité 0.083333
  GSET50.test <- GSE_HW_JLT(N, t, t+0.0833333, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  #Maturité 0.25
  GSET50.test <- GSE_HW_JLT(N, t, t+0.25, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  #Maturité 0.5
  GSET50.test <- GSE_HW_JLT(N, t, t+0.5, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  #Maturité 0.75
  GSET50.test <- GSE_HW_JLT(N, t, t+0.75, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  
  for(Tt in 1:30)
  {
    GSET50.test <- GSE_HW_JLT(N, t, t+Tt, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
    CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  }
  
  #Maturité 40
  GSET50.test <- GSE_HW_JLT(N, t, t+40, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  #Maturité 50
  GSET50.test <- GSE_HW_JLT(N, t, t+50, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  
  CourbeZC.dataframe <- as.data.frame(CourbeZC)
  
  #On renomme les colonnes
  names(CourbeZC.dataframe)[1] <- 0.0833333
  names(CourbeZC.dataframe)[2] <- 0.25
  names(CourbeZC.dataframe)[3] <- 0.5
  names(CourbeZC.dataframe)[4] <- 0.75
  names(CourbeZC.dataframe)[5:34] <- Maturite30
  names(CourbeZC.dataframe)[35] <- 40
  names(CourbeZC.dataframe)[36] <- 50
  write_xlsx(CourbeZC.dataframe,paste(chemin_export,"/Courbe_Taux_dans_",t,"_an_HLG-HW_20220224_11h50m19s_QZCsto _.xlsx", sep=""))
}



### HLG-HW_20220224_11h50m19sQXXX #####################################
projection = 3

GSET50.test <- GSE_HW_JLT(N, TT, TT, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
Proba_survie_AAA <- c(GSET50.test$PZCr_AAA)
Proba_survie_AA <- c(GSET50.test$PZCr_AA)
Proba_survie_A <- c(GSET50.test$PZCr_A)
Proba_survie_BBB <- c(GSET50.test$PZCr_BBB)
Proba_survie_BB <- c(GSET50.test$PZCr_BB)
Proba_survie_B <- c(GSET50.test$PZCr_B)


for(i in 1:TT){
  print(paste("Année de projection : ", i))
  GSET50.test <- GSE_HW_JLT(N, TT-i, TT, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  GSET50.test <- GSE_HW_JLT(N, TT-i, TT, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  Proba_survie_AAA <- cbind(Proba_survie_AAA, GSET50.test$PZCr_AAA)
  Proba_survie_AA <- cbind(Proba_survie_AA, GSET50.test$PZCr_AA)
  Proba_survie_A <- cbind(Proba_survie_A, GSET50.test$PZCr_A)
  Proba_survie_BBB <- cbind(Proba_survie_BBB, GSET50.test$PZCr_BBB)
  Proba_survie_BB <- cbind(Proba_survie_BB, GSET50.test$PZCr_BB)
  Proba_survie_B <- cbind(Proba_survie_B, GSET50.test$PZCr_B)
}

plot(colMeans(Proba_survie_))
plot(PrixZC)
Proba_survie_A[3,]
lines(t(apply(t(Proba_survie_AAA), 2, cummin))[75,])
lines(colMeans(t(apply(t(Proba_survie_AAA), 2, cummin))))
lines(t(apply(t(Proba_survie_AA), 2, cummin))[3,])
plot(1-LGD+LGD*t(apply(t(Proba_survie_BB), 2, cummin))[53,])
Proba_survie_AA[6,]

# PZCr = PZC * (1-LGD+LGD*survie)

Proba_survie_AAA.dataframe <- as.data.frame(Proba_survie_AAA)
Proba_survie_AA.dataframe <- as.data.frame(Proba_survie_AA)
Proba_survie_A.dataframe <- as.data.frame(Proba_survie_A)
Proba_survie_BBB.dataframe <- as.data.frame(Proba_survie_BBB)
Proba_survie_BB.dataframe <- as.data.frame(Proba_survie_BB)
Proba_survie_B.dataframe <- as.data.frame(Proba_survie_B)



#On renomme les colonnes
names(Proba_survie_AAA.dataframe)[1] <- 0
names(Proba_survie_AA.dataframe)[1] <- 0
names(Proba_survie_A.dataframe)[1] <- 0
names(Proba_survie_BBB.dataframe)[1] <- 0
names(Proba_survie_BB.dataframe)[1] <- 0
names(Proba_survie_B.dataframe)[1] <- 0



for(i in 2:TT)
{
  names(Proba_survie_AAA.dataframe)[i] <- i-1
  names(Proba_survie_AA.dataframe)[i] <- i-1
  names(Proba_survie_A.dataframe)[i] <- i-1
  names(Proba_survie_BBB.dataframe)[i] <- i-1
  names(Proba_survie_BB.dataframe)[i] <- i-1
  names(Proba_survie_B.dataframe)[i] <- i-1
}
write_xlsx(Proba_survie_AAA.dataframe,paste(chemin_export,"/Proba_survie_AAA", ".xlsx", sep=""))
write_xlsx(Proba_survie_AA.dataframe,paste(chemin_export,"/Proba_survie_AA", ".xlsx", sep=""))
write_xlsx(Proba_survie_A.dataframe,paste(chemin_export,"/Proba_survie_A", ".xlsx", sep=""))
write_xlsx(Proba_survie_BBB.dataframe,paste(chemin_export,"/Proba_survie_BBB", ".xlsx", sep=""))
write_xlsx(Proba_survie_BB.dataframe,paste(chemin_export,"/Proba_survie_BB", ".xlsx", sep=""))
write_xlsx(Proba_survie_B.dataframe,paste(chemin_export,"/Proba_survie_B", ".xlsx", sep=""))

