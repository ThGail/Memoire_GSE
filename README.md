# Mémoire GSE - Best Estimate

Ces dernières années sont marquées par une forte augmentation d’épargne avec un record atteint en
2019 au début de la crise sanitaire. La plupart des individus choisissent d’investir en assurance-vie
dans l’optique de financer leur retraite voire un projet futur ou bien faire profiter leurs héritiers en cas
décès. Face à cette forte demande, les assureurs doivent être en mesure à la fois de faire fructifier les
capitaux à travers des placements dans des produits financiers mais également de contrôler les risques
engendrés.
C’est dans le cadre de ces contrats d’épargne en fonds euros que l’évaluation des provisions techniques
«Best Estimate» (BE) intervient. Il s’agit d’une réserve permettant de garantir les engagements de
versement des organismes d’assurance envers leurs assurés. Sous la directive de la Solvabilité II, la
notion du BE est définit par « la meilleure estimation correspond à la moyenne pondérée par leur
probabilité des flux de trésorerie futurs, compte tenu de la valeur temporelle de l’argent (valeur actuelle
attendue des flux de trésorerie futurs) ».
Remarquons que le portefeuille d’un assureur est constitué majoritairement d’obligations souveraines
et de sociétés (58% en 2020 d’après ACPR)1. Compte tenu du poids important des actifs de l’assureur
dans les titres d’entreprise, la construction du BE doit permettre d’anticiper le « risque de crédit »
ainsi que les possibles états futurs de l’économie au moment des versements.
Pour un émetteur de titre donné, le risque de crédit se concrétise par :
• un risque de défaut
• un risque d’écartement des spreads de crédit
• un risque de dégradation de notation des entreprises
L’objectif de ce mémoire est donc de s’intéresser sur la modélisation du risque de crédit et d’évaluer
son influence sur la valeur du BE.
Pour se faire nous allons construire un générateur de scénarios économiques (GSE) permettant de
projeter et simuler la variation des taux financiers ainsi que les grandeurs économiques à différents
horizons temporels à travers plusieurs modèles de base (taux sans risque, crédit, immobilier, action).
Après avoir étudié différents modèles (de taux, action et immobilier), nous construirons un modèle
témoin qui n’est autre qu’un simple GSE sans tenir compte du risque de crédit. Ensuite nous inclurons
ce dernier dans le GSE grâce à plusieurs modélisations. Enfin nous allons étudier l’impact de ce dernier
sur la valorisation du BE, et des différents paramètres pouvant influencer sur cette provisi
