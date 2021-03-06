# Mémoire GSE - Best Estimate

__Mémoire de Master 1, Université Paris-Dauphine & Prim'Act__

Dossier rendu pour le mémoire de M1 sur le sujet __"Générateur de scénarios économiques pour le risque de crédit en assurance"__, 
réalisé par _Mengru CHEN_, _Thibault GAILLARD_ et _Sina MBAYE_.

Mise en place d'un générateur de scénarios économiques. Modèles mis en place :
 - Taux sans risque : Modèle Vasicek et Modèle de Hull&White
 - Modèle action et immibilier : Modèle de Black&Scholes
 - Modèle de risque de crédit : Modèle CIR et Modèle JLT

Trois dossiers contiennent différents codes
- Fonction : contient l'ensemble des fonctions que nous avons définies
- Calibrage : contient le calibrage des modèles
- Resultat : contient les graphiques des simulations et tests de martingalité

Les données sont dans le fichier excel "Input_20210118_18h41m33s" contenant les données fournies par Prim'Act.

L'ensemble des fichiers ont servi à la construction de GSE où les résultats se trouvent dans "Resultat > GSE_resultat".
Etant donné que nous avons plusieurs modèles de taux et de crédit, nous avons élaboré 4 versions :
- Vasicek / CIR
- Vasicek / JLT
- Hull&White / CIR
- Hull&White / JLT

Attention : Lors de l'exécution, pensez à modifier les chemins d'accès dans TOUS les fichiers.

Mode d'emploi plus détaillé "Mode d'emploi.txt"
