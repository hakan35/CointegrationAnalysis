#generated time series
x <- diffinv(rnorm(999))

#Generate dummy variable, which incorporates two structural breaks at position 400 and 800
myDt <- cbind(
  "Dt1" =  c(rep(0,399),rep(1,401),rep(0,200)),
  "Dt2" =  c(rep(0,801),rep(1,199)))


#Ijt - the indicator dummy which is only 1 on the break date

myIndicationD <- cbind( "IndiD1" = c(rep(0,399),rep(1,1),rep(0,600)),
                        "IndiD2" = c(rep(0,799),rep(1,1),rep(0,200)))
#Create trend variable
myTrend <-as.matrix(cbind( "Trend" = as.numeric(1:nrow(myDt))))

#Number of lags in the VAR specification
mylags <- 5

#Number of sub-sample periods = number of breaks +1
mySubSample <- 3

#Create lagged variables, which are needed for the estimation of the VAR model.
# Variables to create
# - (Linear) trend is named myTrend
# - number of lags used in the estimation k - mylags 
# - D2,t-k (where k, the maximum lag length)  
# - Trend*D2,t-k
# - I2,t ; I2,t-1; ......; I2,(t-(k-1); where again, k = 5 in our case.

#Create lagmatrix function including automatic naming of the columns
lagmatrix <- function(x,max.lag){
  x <- as.matrix(x)
  if(is.null(colnames(x))== TRUE){
    colnames(x) <- "VarCol0"
  }
  return.matrix <-  embed(c(rep(NA,max.lag),x),max.lag+1)
  dimnames(return.matrix)[[2]] <- c(colnames(x)[1, drop = FALSE], paste(colnames(x)[1,drop = FALSE],".l",1:max.lag, sep = ""))
  return(return.matrix)
}

# InterventionDummy multiplied with trend matrix

myDtTrend <- myDt * as.vector(myTrend)
colnames(myDtTrend) <- paste(colnames(myDt),colnames(myTrend), sep = "*")

# Create matrix, which includes all the lagged dummy variable
# Combination of the Intervention and Indication dummies do cause issues
# Dummy matrix according to Joyeux(2007: 11) or equation (13)
dummat <- matrix(NA,dt.myVariables[,.N], mylags + 4 + (mySubSample -1))
dummat <- cbind(#myTrend
  lagmatrix(myDtTrend[,c("Dt1*Trend"), drop = FALSE], max.lag = mylags)[,mylags +1, drop = FALSE]
  ,lagmatrix(myDtTrend[,c("Dt2*Trend"), drop = FALSE], max.lag = mylags)[,mylags +1, drop = FALSE]
  ,lagmatrix(myDt[,c("Dt1"), drop = FALSE], max.lag = mylags)[,mylags + 1, drop = FALSE]
  ,lagmatrix(myDt[,c("Dt2"), drop = FALSE], max.lag = mylags)[,mylags + 1, drop = FALSE]
  ,lagmatrix(myIndicationD[,c("IndiD1"), drop = FALSE], max.lag = mylags-1)#[,1, drop = FALSE]
  ,lagmatrix(myIndicationD[,c("IndiD2"), drop = FALSE], max.lag = mylags-1)#[,1, drop = FALSE]
)

#Replace NA with 0
dummat[is.na(dummat)] <- 0
