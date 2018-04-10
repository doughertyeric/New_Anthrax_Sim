library(raster)
library(survival)
library(tidyverse)
library(moveHMM)

####################################################################
############                Functions                 ##############
####################################################################

weighted_sampling <- function(input.raster, N) {
  x <- getValues(input.raster)
  x[is.na(x)] <- 0
  samp <- sample(nrow(input.raster)*ncol(input.raster), size=N, prob=x)
  samp.raster <- raster(input.raster)
  samp.raster[samp] <- 1

  points <- data.frame(rasterToPoints(samp.raster, fun=function(x){x == 1}))[,1:2]
  for (i in 1:length(points)) {
    for (j in 1:2) {
      points[i,j] <- points[i,j] + runif(1,-14.99,14.99)
    }
  }
  points <- SpatialPoints(points)
  return(points)
}

selection_functions <- function(name.list) {
  forage <- read.csv(paste0('Zebra_Data/', name.list[1], "_Foraging_Final.csv"))
  for (i in 2:length(name.list)) {
    temp <- read.csv(paste0('Zebra_Data/', name.list[i], "_Foraging_Final.csv"))
    forage <- data.frame(rbind(forage, temp))
  }
  
  directed <- read.csv(paste0('Zebra_Data/', name.list[1], "_Directed_Final.csv"))
  for (i in 2:length(name.list)) {
    temp <- read.csv(paste0('Zebra_Data/', name.list[i], "_Directed_Final.csv"))
    directed <- data.frame(rbind(directed, temp))
  }
  
  forage_clogit <- clogit(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm + Risk_Norm + 
                                   strata(fix) + strata(ID), data=forage)
  directed_clogit <- clogit(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm + Risk_Norm + 
                              strata(fix) + strata(ID), data=directed)
  
  out_list <- list(forage_clogit, directed_clogit)
  return(out_list)
}

HMM <- function(data, states = 3) {

  all_data <- prepData(trackData = data.frame(data), type = 'UTM',
                        coordNames = c('x','y'))
  if (states == 3) {
    mu0 <- c(5,50,500) # step mean (three parameters: one for each state)
    sigma0 <- c(2,20,200) # step SD
    zeromass0 <- c(0,0,0) # step zero-mass
    stepPar0 <- c(mu0,sigma0,zeromass0)
  
    angleMean0 <- c(0,0,0) # angle mean
    kappa0 <- c(1,1,0.1) # angle concentration
    anglePar0 <- c(angleMean0,kappa0)
  
    mod <- fitHMM(data=all_data, nbStates=3, 
                 stepPar0=stepPar0, anglePar0=anglePar0)
  } else {
    mu0 <- c(5,500) # step mean (two parameters: one for each state)
    sigma0 <- c(2,250) # step SD
    zeromass0 <- c(0,0) # step zero-mass
    stepPar0 <- c(mu0,sigma0,zeromass0)
    
    angleMean0 <- c(0,0) # angle mean
    kappa0 <- c(0.1,1) # angle concentration
    anglePar0 <- c(angleMean0,kappa0)
    
    mod <- fitHMM(data=all_data, nbStates=2, 
                  stepPar0=stepPar0, anglePar0=anglePar0)
  }
  return(mod)
}

####################################################################

name_list = c('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009',
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010',
              'AG253_2010', 'AG255_2010', 'AG256_2010')

risk09 <- raster('Layers/Anthrax_Risk_noNA_2009.tif')
risk10 <- raster('Layers/Anthrax_Risk_noNA_2010.tif')
risk_mean <- mean(risk09, risk10)

LIZs <- weighted_sampling(input.raster = risk_mean, N = 200)
SSFs <- selection_functions(name.list = name_list)
foraging.coeff <- SSFs[[1]]$coefficients
directed.coeff <- SSFs[[2]]$coefficients

zebra09 <- read_csv("Zebra_Data/Zebra_Anthrax_2009_Cleaned.csv") %>%
  dplyr::select(x,y,date,ID) 
zebra10 <- read_csv("Zebra_Data/Zebra_Anthrax_2010_Cleaned.csv") %>%
  dplyr::select(x,y,date,ID)
all_data <- rbind(zebra09, zebra10)

HMMs <- HMM(data = all_data, states = 3)
trans.mat <- HMMs$mle$gamma


