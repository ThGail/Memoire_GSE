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

N = 1000; TT = 50; LGD = 0.3; S0=1

# fonction permettant de passer du PZC au TZC 
FP <- function(TT, t, x){return(-log(x)/(TT-t))}



# code à changer pour PZCr
# code à changer pour différente combinaison de modèle pour le GSE

### Exportation déflateur - action - immo ########
defla_HW <- c(rep(1,N)) # pour T=0
action_HW <- c(rep(S0,N)) # pour T=0
immo_HW <- c(rep(S0,N))

for (t in Maturite50){
  defla_HW <- cbind(defla_HW,GSE_HW_JLT(N, 1/360, t, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)$deflateur_c)
  GSE <- GSE_HW_JLT(N, t, NA, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  action_HW <- cbind(action_HW,GSE$ind_action)
  immo_HW <- cbind(immo_HW,GSE$ind_immo)
}
defla_HW.df <- as.data.frame(defla_HW)
action_HW.df <- as.data.frame(action_HW)
immo_HW.df <- as.data.frame(immo_HW)

names(defla_HW.df) <- c(0,Maturite50)
names(action_HW.df) <- c(0,Maturite50)
names(immo_HW.df) <- c(0,Maturite50)

write.csv2(defla_HW.df,paste(chemin_export,"/HLG-HW_20220224_11h50m19s_QDeflateur", ".csv", sep=""), row.names = FALSE)
write.csv2(action_HW.df,paste(chemin_export,"/HLG-HW_20220224_11h50m19sQActionsGlobales", ".csv", sep=""),row.names = FALSE)
write.csv2(action_HW.df,paste(chemin_export,"/HLG-HW_20220224_11h50m19sQActionsAutres", ".csv", sep=""),row.names = FALSE)
write.csv2(immo_HW.df,paste(chemin_export,"/HLG-HW_20220224_11h50m19sQImmobilier", ".csv", sep=""),row.names = FALSE)





### Exportation TZC sans risque #####################################
for(t in 0:50){
  print(paste("Année de projection : ", t))
  CourbeZC <- c()
  
  #Maturité 0.083333
  GSET50.test <- GSE_HW_JLT(N, t, t+0.0833333, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  CourbeZC_AAA <- cbind(CourbeZC, f( GSET50.test$PZC_AAA))
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
  write.csv2(CourbeZC.dataframe,paste(chemin_export,"/Courbe_Taux_dans_",t,"_an_HLG-HW_20220224_11h50m19s_QZCsto _.csv", sep=""),row.names = FALSE)
}





### Exportation TZC risqué #####################################
for(t in 0:50){
  print(paste("Année de projection : ", t))
  CourbeZC_B <- c()
  
  #Maturité 0.083333
  GSET50.test <- GSE_HW_JLT(N, t, t+0.0833333, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC_B <- cbind(CourbeZC_B, FP(t, t+0.0833333, GSET50.test$PZCr_B))
  #Maturité 0.25
  GSET50.test <- GSE_HW_JLT(N, t, t+0.25, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC_B <- cbind(CourbeZC_B, FP(t, t+0.25, GSET50.test$PZCr_B))
  #Maturité 0.5
  GSET50.test <- GSE_HW_JLT(N, t, t+0.5, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC_B <- cbind(CourbeZC_B, FP(t, t+0.5, GSET50.test$PZCr_B))
  #Maturité 0.75
  GSET50.test <- GSE_HW_JLT(N, t, t+0.75, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC_B <- cbind(CourbeZC_B, FP(t, t+0.75, GSET50.test$PZCr_B))
  
  for(Tt in 1:30)
  {
    GSET50.test <- GSE_HW_JLT(N, t, t+Tt, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
    CourbeZC_B <- cbind(CourbeZC_B, FP(t, t+Tt, GSET50.test$PZCr_B))
  }
  
  #Maturité 40
  GSET50.test <- GSE_HW_JLT(N, t, t+40, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC_B <- cbind(CourbeZC_B, FP(t, t+40, GSET50.test$PZCr_B))
  #Maturité 50
  GSET50.test <- GSE_HW_JLT(N, t, t+50, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC_B <- cbind(CourbeZC_B, FP(t, t+50, GSET50.test$PZCr_B))
  
  CourbeZC_B.dataframe <- as.data.frame(CourbeZC_B)
  
  #On renomme les colonnes
  names(CourbeZC_B.dataframe)[1] <- 0.0833333
  names(CourbeZC_B.dataframe)[2] <- 0.25
  names(CourbeZC_B.dataframe)[3] <- 0.5
  names(CourbeZC_B.dataframe)[4] <- 0.75
  names(CourbeZC_B.dataframe)[5:34] <- Maturite30
  names(CourbeZC_B.dataframe)[35] <- 40
  names(CourbeZC_B.dataframe)[36] <- 50
  write.csv2(CourbeZC_B.dataframe,paste(chemin_export,"/Courbe_Taux_dans_",t,"_an_HLG-HW_20220224_11h50m19s_QZCsto _.csv", sep=""),row.names = FALSE)
}
