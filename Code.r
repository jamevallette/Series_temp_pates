#Projet de séries temporelles
#Etude de la production de pâtes alimentaires

datafile <- "pates_alim.csv"
data <- read.csv(datafile,sep=";", skip=3, header = TRUE)


install.packages("zoo")
install.packages("tseries")

require(zoo)
require(tseries)


colnames(data) <- c("Periode", "Indice", "Code")

# Conversion de la période en Date (on ajoute "-01" pour avoir un jour du mois)
data$Date <- as.Date(paste0(data$Periode, "-01"))

# Création du zoo
xm_source <- zoo(data$Indice, order.by = data$Date)

# Tracé
plot(xm_source, type = "l", xlab = "Date", ylab = "Indice", main = "Indice de Production Industrielle - Pâtes alimentaires")

head(data)