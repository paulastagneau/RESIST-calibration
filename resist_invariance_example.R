#------------------------------------------------------------------------------
## Example script to calculate the invariance metric of the RESIST calibration
## The GR4J model example comes from the airGR package documentation
## Please see the manuscript for details about the invariance metric
##
## R version 4.3.1
## Creation: 23-06-2025
## Author: Paul C. Astagneau (paul.astagneau@slf.ch)
#------------------------------------------------------------------------------

# Libraries
library(transport) #For the Wasserstein distance
library(airGR)     #For the GR4J model


#------------------------------------------------------
#- Function to calculate the invariance metric
#------------------------------------------------------
invarCrit <- function(theResiduals, theY, N=100) { 
  #theResiduals DOUBLE   residual time series 
  #theY         INTEGER  Years corresponding to the residual time series 
  #N            INTEGER  Number of repetition (see equations in the manuscript)          
  require(transport)
  
  critRes <- matrix(nrow = N, ncol = 2)
  iPer <- 1
  compteur <- 0
  # Loop over N repetitions
  while (iPer<=N & compteur<(2*N)) {
    compteur <- compteur+1
    
    # Years sampling (eq. 2 to 5)
    sampleSize <- as.integer(length(unique(theY))/2)
    invSampleSize <- length(unique(theY))-sampleSize
    rdVec <- sample(size=sampleSize, x=unique(theY)) 
    a <- theResiduals[theY%in%rdVec]
    b <- theResiduals[!theY%in%rdVec]
    
    # Wasserstein distance (eq. 6 and 7)
    critRes[iPer,1] <- transport::wasserstein1d(a=a[!is.na(a)], b=b[!is.na(b)], p=1) 

    iPer <- iPer+1
  }
  
  # Mean distance over N (eq. 8)
  critResFinal <- mean(critRes[,1])
  
  #Out
  return(critResFinal)
}


#------------------------------------------------------
#- GR4J model run (see airGR documentation)
#------------------------------------------------------

# Loading catchment data example
data(L0123002)

# Preparation of the InputsModel object
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = BasinObs$DatesR,
                                 Precip = BasinObs$P, PotEvap = BasinObs$E)

# run period selection
Ind_Run <- seq(which(format(BasinObs$DatesR, format = "%Y-%m-%d")=="1991-01-01"),
               which(format(BasinObs$DatesR, format = "%Y-%m-%d")=="2010-12-31"))
theYears <- as.integer(format(BasinObs$DatesR[Ind_Run],"%Y"))

# Preparation of the RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                               InputsModel = InputsModel, IndPeriod_Run = Ind_Run)

# Simulation
Param <- c(X1 = 257.238, X2 = 1.012, X3 = 88.235, X4 = 2.208)
OutputsModel <- RunModel_GR4J(InputsModel = InputsModel,
                              RunOptions = RunOptions, Param = Param)


#------------------------------------------------------
#- Invariance metric
#------------------------------------------------------

# Boxcox transformation
lambda <- 0.2 #See Berthet (2020)
obsBox <- ((BasinObs$Qmm[Ind_Run] ^ lambda - (0.01*mean(BasinObs$Qmm[Ind_Run], na.rm=TRUE))^lambda) / lambda) #Santos (2018) eq10
simBox <- ((OutputsModel$Qsim ^ lambda - (0.01*mean(OutputsModel$Qsim, na.rm=TRUE))^lambda) / lambda) #Santos (2018) eq10

# Residuals (eq. 1)
resBox <- obsBox-simBox 

# Invariance metric
invarCrit(theResiduals=resBox, theY=theYears, N=100)


