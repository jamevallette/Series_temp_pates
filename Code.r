#Projet de séries temporelles
#Etude de la fabrication de boissons

#Installation des packages necessaires
install.packages("zoo")
install.packages("tseries")
install.packages("urca")

require(zoo)
require(tseries)
library(urca)

## PARTIE 1 ##

# Importation et nettoyage rapide de la serie
datafile <- "fabricationdeboissons.csv"
data <- read.csv(datafile,sep=";", skip=3, header = TRUE)

colnames(data) <- c("Periode", "Indice","Code")
data$Periode <- as.Date(paste0(data$Periode, "-01"))

xm_source <- zoo(data$Indice, order.by = data$Periode)
T <-length(xm_source)
xm <- xm_source[1:(T-4)]

plot(xm, type = "l", col = "blue", xlab = "Date", ylab = "Indice")

#### QUESTION 2 : 

# On effectue une différentiation à l'ordre 1 car on remarque une tendance dans la serie

xm_diff <- diff(xm)

#### QUESTION 3 : 

## Ouverture d'une fenêtre à 2 graphiques
par(mfrow = c(2, 1))  # 2 lignes, 1 colonne

# Premier graphique : série log-transformée
plot(xm, type = "l", col = "blue", xlab = "Date", ylab = "Indice")

# Deuxième graphique : série différenciée
plot(xm_diff, type = "l", col = "red", xlab = "Date", ylab = "Indice différencié")

#Elle semble maintenant stationnaire avec deux valeurs abérrantes liées au COVID

## PARTIE 2 ##
#On affiche les acf et pacf pour estimer un ARMA(p,q) sur la serie differenciée qui semble stationnaire
acf(coredata(xm_diff))
pacf(coredata(xm_diff))

#3 Tests pour savoir si la série admet des racines unitaires : 

# ADF sans constante ni tendance
adf_none <- ur.df(coredata(xm_diff), type = "none", selectlags = "AIC")
summary(adf_none)
#On rejette bien H0 qui est la présence de racine unitaire

# PP sans constante ni tendance
pp_test <- ur.pp(coredata(xm_diff), type = "Z-tau", model = "constant")
summary(pp_test)
#On rejette H0 la présence de racine unitaire également (avec le parametre constante qui n'est pas significatif )

# Test KPSS
test_kpss <- kpss.test(coredata(xm_diff), null = "Level")  # Série stationnaire autour de la moyenne nulle
print(test_kpss)
# On ne rejette pas H0 ici qui est la non présence de racine unitaire. (H0 : modèle stationnaire (TS))

# On peut ainsi estimé un ARMA sur la série différenciée :
#On affiche les acf et pacf pour estimer un ARMA(p,q) sur la serie differenciée qui semble stationnaire
acf(coredata(xm_diff))
pacf(coredata(xm_diff))

#On estime q avec le graphe des acf => on prend q=2 (les lag après sont nuls hormis 1 ou 2 qui dépasse très peu la borne)
#On estime p avec le graphe des pacf => on peut commencer avec p = 7

#On estime le modèle ARMA sur la serie corrigée : 
arima(coredata(xm_diff),c(7,0,2))

