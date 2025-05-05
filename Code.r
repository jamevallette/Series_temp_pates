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
arima702 <- arima(coredata(xm_diff),c(7,0,2))
Box.test(arima702$residuals, lag=10, type="Ljung-Box", fitdf=5) #test de Ljung-Box

#test pour plus de lag : 
Qtests <- function(series, k, fitdf=0) {
pvals <- apply(matrix(1:k), 1, FUN=function(l) {
pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
return(c("lag"=l,"pval"=pval))
})
return(t(pvals))
}
Qtests(arima702$residuals, 24, 9) #tests de LB pour les ordres 1 `a 24

## Test différents ARMA : 

#fonction de test des significativit´es individuelles des coefficients
signif <- function(estim){

coef <- estim$coef
se <- sqrt(diag(estim$var.coef))
t <- coef/se
pval <- (1-pnorm(abs(t)))*2
return(rbind(coef,se,pval))
}
signif(arima702) #tests de siginificativit´e de l’ARIMA(3,0,2)

##fonction d’affichage des tests pour la sélection du modèle ARIMA
arimafit <- function(estim){
adjust <- round(signif(estim),3)
pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
pvals <- matrix(apply(matrix(1:24,nrow=12),2,function(c) round(pvals[c,],3)),nrow=6)
colnames(pvals) <- rep(c("lag", "pval"),4)
cat("tests de nullit´e des coefficients :\n")
print(adjust)
cat("\n tests d’absence d’autocorr´elation des r´esidus : \n")
print(pvals)
}

estim <- arima(coredata(xm_diff),c(7,0,2)); arimafit(estim)
