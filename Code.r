#Projet de séries temporelles
#Etude de la fabrication de boissons

#Installation des packages necessaires
install.packages("zoo")
install.packages("tseries")
install.packages("urca")
install.packages("fUnitRoots")

require(zoo)
require(tseries)
library(urca)
require(fUnitRoots)

## PARTIE 1 ##

# Importation et nettoyage rapide de la serie
datafile <- "/home/onyxia/work/Series_temp_vin/fabricationdeboissons.csv"
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

arima702 <- arima(coredata(xm_diff),c(7,0,2))

Qtests(arima702$residuals, 24, 3)
signif(arima702) #tests de siginificativit´e de l’ARIMA(3,0,2)


pmax=7
qmax=2

mat <- matrix(NA,nrow=pmax+1,ncol=qmax+1) #matrice vide `a remplir
rownames(mat) <- paste0("p=",0:pmax) #renomme les lignes
colnames(mat) <- paste0("q=",0:qmax) #renomme les colonnes
AICs <- mat #matrice des AIC non remplie
BICs <- mat #matrice des BIC non remplie
pqs <- expand.grid(0:pmax,0:qmax) #toutes les combinaisons possibles de p et q
for (row in 1:dim(pqs)[1]){ #boucle pour chaque (p,q)
p <- pqs[row,1] #r´ecup`ere p
q <- pqs[row,2] #r´ecup`ere q
estim <- try(arima(coredata(xm_diff),c(p,0,q),include.mean = F)) #tente d’estimer l’ARIMA
AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic #assigne l’AIC
BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) #assigne le BIC
}
AICs
AICs==min(AICs)
BICs
BICs == min(BICs)

arima101 <- arima(coredata(xm_diff),c(1,0,1))


adj_r2 <- function(model, data){
  # Somme des résidus au carré
  ss_res <- sum(model$residuals^2)

  # Ordres AR et MA
  p <- model$arma[1]
  q <- model$arma[2]

  # Série utilisée : data sans les premières valeurs perdues à cause des retards
  ss_tot <- sum(data[-c(1:max(p, q))]^2)

  # Nombre d'observations utiles (nobs déjà corrigé par arima)
  n <- model$nobs - max(p, q)

  # Calcul du R² ajusté
  adj_r2 <- 1 - (ss_res / (n - p - q - 1)) / (ss_tot / (n - 1))

  return(adj_r2)
}
adj_r2(arima101, coredata(xm_diff))


dev.off() #r´einitialise les param`etres graphiques
plot(arima101$residuals)




# Rechargement de la série source
T <- length(xm_source)
xm <- xm_source[1:(T-4)]  # série réelle (non différenciée)

# Ajustement du modèle ARIMA(1,1,1) sur la série réelle
model <- arima(coredata(xm), order = c(1, 1, 1))

# Prévision à 12 mois
h <- 12
forecast <- predict(model, n.ahead = h)

# Prévisions et erreur standard
pred <- forecast$pred
se <- forecast$se

# Construction de l’intervalle de confiance à 95 %
upper <- pred + 1.96 * se
lower <- pred - 1.96 * se

# Génération des dates futures
last_date <- index(xm)[length(xm)]
future_dates <- seq(from = last_date + 30, by = "month", length.out = h)

y_min <- min(xm, lower)
y_max <- max(xm, upper)

# Tracé de la série observée + prévision + IC
plot(xm, type = "l", 
     xlim = c(start(xm), future_dates[h] + 30),
     ylim = c(y_min, 110),
     main = "Prévision ARIMA(1,1,1) avec intervalle de confiance à 95%",
     ylab = "Indice de production", xlab = "Date")

lines(future_dates, pred, col = "blue", lwd = 2)
lines(future_dates, upper, col = "red", lty = 2)
lines(future_dates, lower, col = "red", lty = 2)
legend("topleft", legend = c("Prévision", "Intervalle 95%"),
       col = c("blue", "red"), lty = c(1,2), bty = "n")
