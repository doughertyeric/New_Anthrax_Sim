library(raster)
library(survival)
library(tidyverse)
library(moveHMM)
library(rgeos)
library(foreach)
library(doParallel)
library(sf)
library(velox)
library(truncnorm)
library(fasterize)
library(adehabitatLT)
library(lubridate)
library(fitdistrplus)

####################################################################

HMM <- function(data, states = 3) {
  
  all_data <- prepData(trackData = data.frame(data), type = 'UTM',
                       coordNames = c('x','y'))
  if (states == 3) {
    mu0 <- c(5,50,500) # step mean (three parameters: one for each state)
    sigma0 <- c(2,20,200) # step SD
    zeromass0 <- c(0,0,0) # step zero-mass
    stepPar0 <- c(mu0,sigma0)#,zeromass0)
    
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

sim.paths <- function(pt1, pt2, vmax, time.steps, crs) {
  sub.v <- vmax/time.steps
  new.pos <- st_as_sf(pt1, coords=1:2, crs=crs)
  pt2 <- st_as_sf(pt2, coords=1:2, crs=crs)
  intermediate.pts <- data.frame(matrix(0,time.steps,2))
  intermediate.pts[1,] <- st_coordinates(new.pos$geometry)
  for (n in 2:(time.steps-1)) {
    origin.circle.rad <- sub.v*n
    origin.circle <- st_buffer(new.pos, dist=origin.circle.rad)
    destination.circle.rad <- sub.v*(time.steps-n)
    destination.circle <- st_buffer(pt2, dist=destination.circle.rad)
    new.pos.region <- st_intersection(origin.circle, destination.circle)
    new.pt = st_sample(new.pos.region, size=10, type='random')
    intermediate.pts[n,] <- st_coordinates(new.pt[[1]])
  }
  intermediate.pts[time.steps,] <- st_coordinates(pt2$geometry)
  colnames(intermediate.pts) <- c('x', 'y')
  return(intermediate.pts)
}

intermediate.steps <- function(xy, state.vector) {
  path.dens <- data.frame(matrix(0,0,3))
  for (i in 1:(nrow(xy) - 1)) {
    dist.x <- abs(xy[i,1] - xy[(i+1),1])
    dist.y <- abs(xy[i,2] - xy[(i+1),2])
    current.dist <- sqrt(dist.x^2 + dist.y^2)
    
    v_max <- rtruncnorm(n=1,mean=1.1,sd=0.05,a=1.01)*current.dist
    path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], v_max, 21, crs)
    path.dens.temp <- cbind(path.dens.temp, rep(state.vector[i],21))
    colnames(path.dens.temp) <- c('x', 'y', 'state')
    path.dens <- data.frame(rbind(path.dens, path.dens.temp))
    #print(i)
  }
  return(path.dens)
}

LIZ.layer <- function(LIZs, mask.sp, crs) {
  LIZ.df <- data.frame(matrix(0,nrow(LIZs@coords),5))
  
  for (m in 1:nrow(LIZs@coords)) {
    LIZ.df[m,1] <- LIZs@coords[m,1]
    LIZ.df[m,2] <- LIZs@coords[m,2]
    LIZ.df[m,3] <- floor(runif(1,1,3.9999)) #years since death
    rand <- runif(1,0,1)
    if (rand < 0.826590) { #Probability of springbok
      LIZ.df[m,4] <- 1
      LIZ.df[m,5] <- rnorm(1, 350, (350/8)) # body size
    } else if (rand > 0.826589 && rand < 0.974950) {
      LIZ.df[m,4] <- 2
      LIZ.df[m,5] <- rnorm(1, 35, (35/8))
    } else {
      LIZ.df[m,4] <- 3
      LIZ.df[m,5] <- rnorm(1, 3500, (3500/8))
    }
  }
  
  LIZ_sf <- st_as_sf(LIZ.df, coords=1:2, crs=crs)
  
  buffers <- c()
  for (i in 1:nrow(LIZ.df)) {
    age <- LIZ.df[i,3]
    size.factor <- log(LIZ.df[i,5])/log(350)
    buffers[i] <- size.factor * (rnorm(1, (3 * (4 - age)), ((3 * (4 - age))/8)))
  }
  buffers.sp <- st_buffer(LIZ_sf, dist = buffers)
  
  area.extent <- extent(mask.sp)
  r <- raster(area.extent)
  projection(r) <- CRS(crs)
  res(r) <- c(3,3)
  
  LIZ.layer <- fasterize(buffers.sp, r, background=0)
  
  return(LIZ.layer)
}

LIZ.objects <- function(LIZ_sf) {

  for (m in 1:nrow(LIZ_sf)) {
    LIZ_sf$age[m] <- floor(runif(1,1,3.9999))
    rand <- runif(1,0,1)
    if (rand < 0.826590) { #Probability of zebra
      LIZ_sf$species[m] <- 1
      LIZ_sf$size[m] <- rnorm(1, 350, (350/8))
      size.factor <- log(LIZ_sf$size[m])/log(350)
      LIZ_sf$buffer[m] <- size.factor * (rnorm(1, (3 * (4 - LIZ_sf$age[m])), ((3 * (4 - LIZ_sf$age[m]))/8)))
    } else if (rand > 0.826589 && rand < 0.974950) {
      LIZ_sf$species[m] <- 2
      LIZ_sf$size[m] <- rnorm(1, 35, (35/8))
      size.factor <- log(LIZ_sf$size[m])/log(350)
      LIZ_sf$buffer[m] <- size.factor * (rnorm(1, (3 * (4 - LIZ_sf$age[m])), ((3 * (4 - LIZ_sf$age[m]))/8)))
    } else {
      LIZ_sf$species[m] <- 3
      LIZ_sf$size[m] <- rnorm(1, 3500, (3500/8))
      size.factor <- log(LIZ_sf$size[m])/log(350)
      LIZ_sf$buffer[m] <- size.factor * (rnorm(1, (3 * (4 - LIZ_sf$age[m])), ((3 * (4 - LIZ_sf$age[m]))/8)))
    }
  }
  
  return(LIZ_sf)
}

####################################################################

for (i in 2:16) {
  orig <- readRDS(paste0('Final_Runs_Markov/ABM_Final_Markov_k',i,'.rds'))
  new <- readRDS(paste0('Final_Runs_Markov/ABM_Final_Markov_k',i,'_part2.rds'))
  full <- list()
  for (j in 1:500) {
    full[[j]] <- orig[[j]]
  }
  for (k in 1:500) {
    full[[(k+500)]] <- new[[k]]
  }
  saveRDS(full, paste0('Final_Runs_Markov/ABM_Final_Markov_k',i,'_Full.rds'))
  print(i)
}

####################################################################

##### Pattern Matching Step Length Distirbution #####

true_foraging_mu <- c(193.6, 78.3, 68, 188.1, 173.4, 205.6,
                      187, 110.9, 154.1, 97.1, 108.4)
true_directed_mu <- c(853.8, 405.8, 386.4, 758, 907.7, 818.9,
                      852, 600.1, 785.5, 560.1, 581.8)

all.step.means <- data.frame(matrix(0,0,3))
colnames(all.step.means) <- c('foraging.mu', 'directed.mu', 'k')
t.test.results <- data.frame(matrix(0,9,2))

for (k in 2:16) {

  all_inds <- readRDS(paste0('Final_Runs_Markov/ABM_Final_Markov_k',k,'_Full.rds'))
  rand.ids <- sample(x=seq(1,1000,1), size=11, replace=FALSE)
  
  foraging.mu <- c()
  directed.mu <- c()
  for (i in 1:length(rand.ids)) {
    temp_path <- all_inds[[rand.ids[i]]]
    temp_HMM <- HMM(temp_path)
    foraging.mu[i] <- temp_HMM$mle$stepPar[1,2]
    directed.mu[i] <- temp_HMM$mle$stepPar[1,3]
    print(i)
  }
  
  t.test.results[k,1] <- t.test(true_foraging_mu, foraging.mu)$p.value
  t.test.results[k,2] <- t.test(true_directed_mu, directed.mu)$p.value
  
  temp <- data.frame(cbind(foraging.mu, directed.mu, rep(k,11)))
  colnames(temp) <- c('foraging.mu', 'directed.mu', 'k')
  all.step.means <- rbind(all.step.means, temp)
}

all.step.means$k <- as.factor(all.step.means$k)

ggplot() + geom_histogram(aes(x=rnorm(1000,mean(true_foraging_mu),sd(true_foraging_mu))), bins=20) +
  geom_histogram(data=all.step.means, aes(x=foraging.mu, fill=k), bins=20) +
  xlab('Mean Step Length (Foraging)') + ylab('Frequency')

ggplot() + geom_histogram(aes(x=rnorm(1000,mean(true_directed_mu),sd(true_directed_mu))), bins=20) +
  geom_histogram(data=all.step.means, aes(x=directed.mu, fill=k), bins=20) +
  xlab('Mean Step Length (Directed)') + ylab('Frequency')

####################################################################

setwd("~/Box Sync/Dissertation/Behavioral_SSF")

zebra09 <- read_csv("Zebra_Anthrax_2009_Cleaned.csv") %>%
  dplyr::select(x,y,date,ID) %>%
  dplyr::filter(!is.na(x)) %>%
  mutate(row =  seq(1,nrow(.),1), dists = 0, weights = 1) %>%
  st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

zebra10 <- read_csv("Zebra_Anthrax_2010_Cleaned.csv") %>%
  dplyr::select(x,y,date,ID) %>%
  dplyr::filter(!is.na(x)) %>%
  mutate(row =  seq(1,nrow(.),1), dists = 0, weights = 1) %>%
  st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

#all_zebra <- rbind(zebra09, zebra10)

zebra.ext <- extent(c(extent(zebra10)@xmin - 5000,
                      extent(zebra10)@xmax + 5000,
                      extent(zebra09)@ymin - 5000,
                      extent(zebra10)@ymax + 5000))

zebra_ext <- as(zebra.ext, 'SpatialPolygons')
crs(zebra_ext) <- '+init=epsg:32733'
zebra_ext <- st_as_sf(zebra_ext)

ENP <- read_sf('ENP_Predictors/enp fence poly.shp')
ENP_crs <- st_crs(ENP, '+proj=longlat')
st_crs(ENP) <- ENP_crs
ENP <- st_transform(ENP, '+init=epsg:32733') %>%
  st_union(.) %>% st_intersection(zebra_ext)

ENP.sp <- ENP %>%
  as(., 'Spatial')

#####################################################

setwd('~/Desktop/GitHub/New_Anthrax_Sim')

Risk <- raster('Layers/Mean_Risk.tif') %>% 
  mask(ENP.sp)

LIZs <- weighted_sampling(input.raster = Risk, N = 282)
LIZs <- st_as_sf(LIZs, crs=37233)
#LIZ_sf <- LIZ.objects(LIZs)
#st_write(LIZ_sf, "LIZ_buffers.csv", 
#         layer_options = "GEOMETRY=AS_XY")
LIZ_sf <- read_csv('/Users/ericdougherty/Desktop/GitHub/New_Anthrax_Sim/LIZ_buffers.csv') %>%
  st_as_sf(coords = 1:2, crs=32733)
LIZ_buffers <- st_buffer(LIZ_sf, dist = LIZ_sf$buffer)

all_inds <- readRDS('/Users/ericdougherty/Desktop/GitHub/New_Anthrax_Sim/Final_Runs_Markov/ABM_Final_Markov_k1_Full.rds')

all_contacts <- list()
for (i in 1:1000) {
  temp_df <- data.frame(matrix(0,1,2))
  temp_path <- all_inds[[i]]
  temp_path$ID <- 1
  temp_buff <- temp_path %>%
    st_as_sf(coords=1:2, crs=32733) %>%
    group_by(ID) %>%
    summarize(m = mean(state), do_union=FALSE) %>%
    st_cast('LINESTRING')
    #st_buffer(15)
  
  all_lines <- temp_path[1:2,] %>%
    st_as_sf(coords=1:2, crs=32733) %>%
    group_by(ID) %>%
    summarize(m = mean(state), do_union=FALSE) %>%
    st_cast('LINESTRING')
  
  for (j in 2:nrow(temp_path)-1) {
    temp_line <- temp_path[j:(j+1),] %>%
      st_as_sf(coords=1:2, crs=32733) %>%
      group_by(ID) %>%
      summarize(m = mean(state), do_union=FALSE) %>%
      st_cast('LINESTRING')
    all_lines[j,] <- temp_line
    print(j)
  }

  all_lines$state <- temp_path$state[2:nrow(temp_path)]
  all_lines$datetime <- temp_path$datetime[2:nrow(temp_path)]
  forage_lines <- all_lines %>%
    dplyr::select(-ID, -m) %>%
    filter(state == 2)
  forage_buff <- forage_lines %>%
    st_buffer(15)
  
  overlaps <- st_intersection(forage_buff, LIZ_buffers)
  temp_df[1,1] <- nrow(forage_buff)
  temp_df[1,2] <- length(overlaps)
  colnames(temp_df) <- c('foraging.steps', 'contact.count')
  out_list <- list(temp_df, overlaps)
  all_contacts[[i]] <- out_list
}

#####################################################

stack <- function(date_list, mask.sp) {
  
  #Green <- raster('Layers/Mean_Greenness.tif') %>% 
  #  mask(mask.sp)
  #Wet <- raster('Layers/Mean_Wetness.tif') %>% 
  #  mask(mask.sp)
  
  layer.stats <- data.frame(matrix(0,11,3))
  for (i in 1:length(date_list)) {
    Green <- raster(paste0('Layers/Greenness_', date_list[i], '.tif')) %>% 
      mask(mask.sp)
    Wet <- raster(paste0('Layers/Wetness_', date_list[i], '.tif')) %>% 
      mask(mask.sp)
      
    
    green_mean <- cellStats(Green, stat='mean', na.rm=TRUE)
    wet_mean <- cellStats(Wet, stat='mean', na.rm=TRUE)
    layer.stats[i,1] <- date_list[i]
    layer.stats[i,2] <- green_mean
    layer.stats[i,3] <- wet_mean
    print(i)
  }
  colnames(layer.stats) <- c('Date', 'Green', 'Wet')
  return(layer.stats)
}

########################################################

layers <- read.csv('LayerStats.csv')

for (j in 1:9) {

  all_contacts <- readRDS(paste0('ABM_Contact_List_k',j,'.rds'))

  temp_list <- list()
  for (i in 1:1000) {
    temp_list[[i]] <- all_contacts[[i]][[1]]
  }
  
  contact_rates <- data.frame(matrix(unlist(temp_list), nrow=1000, byrow=TRUE))
  
  print(table(contact_rates[,2]))
  
  # for (i in 1:1000) {
  #   contact_rates[i,2] <- nrow(all_contacts[[i]][[2]])
  #   all_contacts[[i]][[1]][1,2] <- contact_rates[i,2]
  # }
  # 
  # layers$mean[j] <- mean(contact_rates[,2])
  # layers$sd[j] <- sd(contact_rates[,2])
  # 
  # saveRDS(all_contacts, paste0('ABM_Contact_List_k',j,'.rds'))
  
}

###########################################################

library(lattice)

my.palette <- c("#FFF5F0", "#FEE0D2", "#FEE0D2", 
                "#FCBBA1", "#FCBBA1", "#FC9272", 
                "#FC9272", "#FB6A4A", "#FB6A4A", 
                "#EF3B2C", "#EF3B2C", "#CB181D", 
                "#CB181D", "#A50F15", "#A50F15", "#67000D")

level <- levelplot(layers$mean ~ layers$Wet * layers$Green, 
                   xlab="Wet", ylab="Green", zlab="mean", 
                   screen = list(z = -30, x=-60), 
                   regions=TRUE, cuts=15, 
                   col.regions=my.palette)
wire <- wireframe(layers$mean ~ layers$Wet * layers$Green, 
                  xlab="Wet", ylab="Green", zlab="mean",
                  screen = list(z = -30, x=-60), 
                  drape=TRUE, cuts=15, 
                  col.regions=my.palette)


