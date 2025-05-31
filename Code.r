"""
Projet de séries temporelles : Etude de la fabrication de boissons

Trois parties :
1. Préparation des données, représentations et tests de stationnarité
2. Estimation du modèle le plus adapté à la série
3. Prédiction pour deux valeurs futures
"""



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
datafile <- "fabricationdeboissons.csv"
data <- read.csv(datafile,sep=";", skip=3, header = TRUE)

colnames(data) <- c("Periode", "Indice","Code")
data$Periode <- as.Date(paste0(data$Periode, "-01"))

xm_source <- zoo(data$Indice, order.by = data$Periode)

# Nous tronquons la série sur ses deux dernières observations pour pouvoir comparer nos observations avec les valeurs réelles
T <-length(xm_source)
xm <- xm_source[1:(T-2)]

pdf("Images/serie.pdf", width = 7, height = 3.5)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
plot(xm, type = "l", col = "blue", xlab = "Date", ylab = "Indice")
dev.off()


#### QUESTION 2 : 

# On effectue une différentiation à l'ordre 1 car on remarque une tendance dans la serie

xm_diff <- diff(xm)

#### QUESTION 3 : 

## Test LJUNGBOX et KPSS sur la série iniale :
# WARNING : Test ADF valide que si les résidus sont décorrelés(LJUNGBOX valide) 

##Test LJUNGBOX pour vérifier la validité du Test ADF : 

#Vérification présence constante ou tendance linéaire : 
summary(lm(coredata(xm) ~ as.numeric(index(xm))))

# Constante positive et le coefficient associé à la tendance linéaire (Période) est bien positif. 
# On ne peut pas interpréter de façon sûr la significativité car les résidus sont peut-être autocorrélés.

# Test ADF avec constante et tendance
adf <- adfTest(coredata(xm), lag = 0, type = "ct")

# Récupération des résidus
res <- adf@test$lm$residuals

# Test de Ljung-Box sur les résidus (lags 1 à 24)
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l <= fitdf) NA else Box.test(series, lag = l, type = "Ljung-Box", fitdf = fitdf)$p.value
    return(c("lag" = l, "pval" = pval))
  })
  return(t(pvals))
}

Qtests(res, k = 24, fitdf = 3)
#Absence d'autocorrelation entre les résidus est rejetée de Q(8) à Q(24) donc on ajoute des retards pour avoir des résidus non corrélés

# Fonction testant la validité de l'ADF en vérifiant la non-corrélation des résidus
adfTest_valid <- function(series, kmax, type) { 
  k <- 0
  noautocorr <- 0
  while (noautocorr == 0 && k <= kmax) {
    cat(paste0("ADF with ", k, " lags: residuals OK? "))
    adf <- adfTest(series, lags = k, type = type)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[, 2]
    
    if (sum(pvals < 0.05, na.rm = TRUE) == 0) {
      noautocorr <- 1
      cat("OK \n")
    } else {
      cat("nope \n")
      k <- k + 1
    }
  }
  
  if (k > kmax) {
    warning("Aucun ADF avec résidus non autocorrélés trouvé jusqu'à ", kmax, " lags.")
    return(NULL)
  } else {
    return(adf)
  }
}

#Utilisation de KPSS : 
kpss_test <- kpss.test(coredata(xm), null = "Level", lshort = TRUE)
print(kpss_test)

#Graphique série différenciée et non différenciée
pdf("Images/serie_init_diff.pdf", width = 7, height = 6)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.2)
plot(xm, type = "l", col = "blue", xlab = "Date", ylab = "Indice", ylim = c(40, 110),
     main = "Série initiale")
plot(xm_diff, type = "l", col = "red", xlab = "Date", ylab = "Indice différencié", ylim = c(-30, 30),
     main = "Série différenciée")
dev.off()

#Elle semble maintenant stationnaire avec deux valeurs abérrantes liées au COVID




### PARTIE 2 ###

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

## On peut ainsi estimé un ARMA sur la série différenciée :

# On affiche les acf et pacf pour estimer un ARMA(p,q) sur la serie differenciée qui semble stationnaire

pdf("Images/acf_pacf_diff.pdf", width = 7, height = 3.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
acf(coredata(xm_diff), main = "ACF de la série différenciée")
pacf(coredata(xm_diff), main = "PACF de la série différenciée")
dev.off()


# On estime q avec le graphe des acf => on prend q=2 (les lag après sont nuls hormis 1 ou 2 qui dépasse très peu la borne)
# On estime p avec le graphe des pacf => on peut commencer avec p = 7

# On test différent modèle ARMA : 

#1 - On estime le modèle ARMA(7,2) sur la serie corrigée : 
arima702 <- arima(coredata(xm_diff),c(7,0,2))

Box.test(arima702$residuals, lag=12, type="Ljung-Box", fitdf=9) #test de Ljung-Box : fitdf = p+q et lag=12 car data mensuelle
# P-value très forte donc les residus ne sont pas autocorrélés

# Fonction de test des significativités individuelles des coefficients
signif <- function(estim){coef <- estim$coef
se <- sqrt(diag(estim$var.coef))
t <- coef/se
pval <- (1-pnorm(abs(t)))*2
return(rbind(coef,se,pval))}

# Rassemblement des commandes pour afficher le modèle ARMA(7,2)
arima702 <- arima(coredata(xm_diff),c(7,0,2))

Qtests(arima702$residuals, 24, 9)
# On accepte l'hypothèse de non corrélation des résidus

signif(arima702) 
# On remarque que le coeff MA(2) ne rejette pas l'hypothèse de nullité, le modèle est donc mal ajusté

# Estimation des modèles ARMA(p,q) valides avec pmax=7 et qmax=2
for (p in 0:7) {
  for (q in 0:2) {
    cat("\n--- ARMA(", p, ",", q, ") ---\n")
    model <- try(arima(coredata(xm_diff), order = c(p, 0, q)), silent = TRUE)
    if (!inherits(model, "try-error")) {
      assign(paste0("arma", p, q), model)
      print(signif(model))
    } else {
      cat("Échec de l'estimation\n")
    }
  }
}
# Modèles bien ajustés : MA(1), MA(2), AR(1), ARMA(1,1), ARMA(3,1), AR(2), AR(3), AR(4), AR(5), ARMA(6,2), AR(7)

# On teste l'absence d'autocorrelation des résidus pour ces modèles bien ajustés : 

## MA(1)
ma1 <- arima(coredata(xm_diff),c(0,0,1))
Qtests(ma1$residuals, 12, 1)
# Autocorrelations des résidus donc MA(1) non valide
# p-value <5% donc on rejette H0 : non autocorrelation des résidus

## MA(2)
ma2 <- arima(coredata(xm_diff),c(0,0,2))
Qtests(ma2$residuals, 12, 2)
# Non autocorrelation des résidus donc MA(2) valide
# pvalue > 5% donc on accepte H0 : non autocorrelation des résidus

## AR(1)
ar1 <- arima(coredata(xm_diff),c(1,0,0))
Qtests(ar1$residuals, 12, 1)
# Autocorrelations des résidus donc AR(1) non valide

## AR(2)
ar2 <- arima(coredata(xm_diff),c(2,0,0))
Qtests(ar2$residuals, 12, 2)
# Autocorrelations des résidus donc AR(2) non valide

## AR(3)
ar3 <- arima(coredata(xm_diff),c(3,0,0))
Qtests(ar3$residuals, 12, 3)
# Autocorrelations des résidus donc AR(3) non valide

## AR(4)
ar4 <- arima(coredata(xm_diff),c(4,0,0))
Qtests(ar4$residuals, 12, 4)
# Autocorrelations des résidus donc AR(4) non valide

## AR(5)
ar5 <- arima(coredata(xm_diff),c(5,0,0))
Qtests(ar5$residuals, 12, 5)
# Autocorrelations des résidus donc AR(5) non valide

## AR(7)
ar7 <- arima(coredata(xm_diff),c(7,0,0))
Qtests(ar7$residuals, 12, 7)
# Non Autocorrelation des résidus donc AR(7) valide

## ARMA(1,1)
arima101 <- arima(coredata(xm_diff),c(1,0,1))
Qtests(arima101$residuals, 12, 2)
# Non Autocorrelation des résidus donc ARMA(1,1) valide
# pvalue > 5% donc on accepte H0 : non autocorrelation des résidus

## ARMA(3,1)
arima301 <- arima(coredata(xm_diff),c(3,0,1))
Qtests(arima301$residuals, 12, 4)
# Autocorrelation des résidus à partir de Q(4) donc ARMA(3,1) non valide
# pvalue < 5% à partir de Q(4) donc on rejette l'hypothèse d'abscence d'autocorrelation des résidus 

## ARMA(6,2)
arima602 <- arima(coredata(xm_diff),c(6,0,2))
Qtests(arima602$residuals, 12, 8)
# Autocorrelation des résidus donc ARMA(6,2) non valide


### On recherche donc le meilleur modèle ARMA en minimisant le critère AIC et BIC: 
pmax=7
qmax=2

mat <- matrix(NA,nrow=pmax+1,ncol=qmax+1) 
rownames(mat) <- paste0("p=",0:pmax) 
colnames(mat) <- paste0("q=",0:qmax) 
AICs <- mat 
BICs <- mat 
pqs <- expand.grid(0:pmax,0:qmax) 
for (row in 1:dim(pqs)[1]){
  p <- pqs[row,1] 
  q <- pqs[row,2] 
  estim <- try(arima(coredata(xm_diff),c(p,0,q),include.mean = F)) 
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic 
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) 
  }

## Résultats :
AICs
AICs==min(AICs)
BICs
BICs == min(BICs)
### L'ARMA(1,1) sur la série corrigée est donc le modèle valide qui minimise les critères BIC et AIC



# Calcul du R² ajusté de ce modèle : 
adj_r2 <- function(model, data){
  ss_res <- sum(model$residuals^2)
  p <- model$arma[1]
  q <- model$arma[2]
  ss_tot <- sum(data[-c(1:max(p, q))]^2)
  n <- model$nobs - max(p, q)
  adj_r2 <- 1 - (ss_res / (n - p - q - 1)) / (ss_tot / (n - 1))

  return(adj_r2)
  }

adj_r2(arima101, coredata(xm_diff))

## Grape des résidus du modèle ARMA(1,1) de la série corrigée: 
plot(arima101$residuals)

### Cela revient à faire un modèle ARIMA(1,1,1) sur la série initiale

### Partie 3 ###

# On cherche dans cette partie à réaliser des prédiction pour deux valeurs de la série
# Nous réalisons les prédictions sur les deux dernières valeurs de la série afin de pourvoir les comparer avec les valeurs réelles

# Construction de la série tronquée pour estimation (jusqu'à T-2)
T <- length(xm_source)
xm_full <- xm_source
xm_train <- xm_source[1:(T - 2)]

# Modèle ARIMA(1,1,1) sur la série tronquée
model <- arima(coredata(xm_train), order = c(1, 1, 1))

# Prédiction à 2 pas
h <- 2
forecast <- predict(model, n.ahead = h)

# Prévisions et intervalles de confiance à 95 %
pred <- forecast$pred
se <- forecast$se
upper <- pred + 1.96 * se
lower <- pred - 1.96 * se

# Valeurs observées réelles pour comparaison
obs <- coredata(xm_full)[(T - 1):T]
last_train_date <- index(xm_train)[length(xm_train)]
future_dates <- seq(from = as.Date(last_train_date) + 30, by = "month", length.out = h)

# Tableau de synthèse
resultats <- data.frame(
  Date = future_dates,
  Prediction = round(pred, 2),
  IC_bas = round(lower, 2),
  IC_haut = round(upper, 2),
  Observé = round(obs, 2)
)

print(resultats)

## On peut enfin tracer la série réelle et y ajouter nos estimations pour les mois de janvier et février 2025

# Tracé de la série observée jusqu’à T-2
plot(xm_train, type = "l",
     main = "Série observée jusqu’à T-2 + Prédictions ARIMA(1,1,1)",
     xlab = "Date", ylab = "Indice de production",
     ylim = range(c(xm_train, pred)))  # pour inclure les IC

# Ajout des prédictions
points(future_dates, pred, col = "blue", pch = 19)
lines(future_dates, pred, col = "blue", lty = 1)

# Ajout des intervalles de confiance à 95 %
lines(future_dates, upper, col = "red", lty = 2)
lines(future_dates, lower, col = "red", lty = 2)


 
## Cependant cette représentation est trop large et ne permet pas de bien visualiser nos prédictions avec leur intervalle de confiance
## On affiche donc les mêmes résultats mais seulement à l'échelle des 12 derniers mois.

n_recent <- 12
xm_recent <- window(xm_train, start = index(xm_train)[length(xm_train) - n_recent + 1])

# Construire un objet zoo pour les prédictions
pred_zoo <- zoo(pred, order.by = future_dates)
upper_zoo <- zoo(upper, order.by = future_dates)
lower_zoo <- zoo(lower, order.by = future_dates)
obs_zoo <- zoo(obs, order.by = future_dates)

# Fusionner série observée récente et prédictions
xm_plot <- c(xm_recent, pred_zoo)

# Export en PDF
pdf("Images/zoom_predictions.pdf", width = 7, height = 4)

# Tracé dans le fichier PDF
plot(xm_plot, type = "l",
     xlab = "Date", ylab = "Indice de production",
     ylim = range(c(upper, lower)))

lines(upper_zoo, col = "red", lty = 2)
lines(lower_zoo, col = "red", lty = 2)
points(pred_zoo, col = "blue", pch = 19)
points(upper_zoo, col = "red", pch = 4)
points(lower_zoo, col = "red", pch = 4)
points(obs_zoo, col = "black", pch = 17)

legend("bottomleft", 
       legend = c("Observé", "Prédiction", "IC 95%", "Observé T+1, T+2"),
       col = c("black", "blue", "red", "black"),
       lty = c(1, 1, 2, NA),
       pch = c(NA, 19, NA, 17),
       bty = "n")

# Fermeture du fichier
dev.off()
