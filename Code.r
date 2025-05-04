#Projet de séries temporelles
#Etude de la fabrication de boissons

#Installation des packages necessaires
install.packages("zoo")
install.packages("tseries")

require(zoo)
require(tseries)

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

acf(coredata(xm_diff))
pacf(coredata(xm_diff))

adf.test(coredata(xm_diff), alternative = "stationary")
pp.test(coredata(xm_diff), alternative = "stationary")
