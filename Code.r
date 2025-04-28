#Projet de séries temporelles
#Etude de la production de pâtes alimentaires

install.packages("zoo")
install.packages("tseries")

require(zoo)
require(tseries)

#   QUESTION 2

##Importation et nettoyage rapide de la série
datafile <- "pates_alim.csv"
data <- read.csv(datafile,sep=";", skip=3, header = TRUE)

colnames(data) <- c("Periode", "Indice","Code")
data$Periode <- as.Date(paste0(data$Periode, "-01"))
xm_source <- zoo(data$Indice, order.by = data$Periode)
T <-length(xm_source)
xm <- xm_source[1:(T-4)]

#Tracé de la série 
plot(xm, type = "l", xlab = "Date", ylab = "Indice")

acf(coredata(xm))
pacf(coredata(xm))

#Première différenciation 
xm_diff <- diff(xm)

plot(xm_diff, type = "l", xlab = "Date", ylab = "Indice")
acf(coredata(xm_diff))
pacf(coredata(xm_diff))


#   QUESTION 3

## Ouverture d'une fenêtre à 2 graphiques
par(mfrow = c(2, 1))  # 2 lignes, 1 colonne

# Premier graphique : série log-transformée
plot(xm, type = "l", col = "blue", xlab = "Date", ylab = "Indice")

# Deuxième graphique : série différenciée
plot(xm_diff, type = "l", col = "red", xlab = "Date", ylab = "Indice différencié")

adf.test(coredata(xm_diff), alternative = "stationary")
