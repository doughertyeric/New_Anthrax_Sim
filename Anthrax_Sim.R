library(raster)
library(survival)
library(tidyverse)
library(moveHMM)
library(rgeos)

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

initialize.states <- function(N) {
  state.vector <- data.frame(matrix(0,N,1))
  for (i in 1:N) {
    rand <- runif(1,0,1)
    if (rand < 0.08007182) {
      state.vector[i,1] <- 1
    } else if (rand > 0.08007181 && rand < 0.6216637) {
      state.vector[i,1] <- 2
    } else {
      state.vector[i,1] <- 3
    }
  }
  colnames(state.vector) <- c('behav.state')
  return(state.vector)
}

initialize.agents <- function(N, sp.mean, SSFs, input.stack, pred_df) {
  
  traits <- data.frame(matrix(0,N,6)) #traitframe
  state.vector <- initialize.states(N)
  forage.maps <- list()
  directed.maps <- list()
  
  # Loop over individuals of a given species
  for (j in 1:N) {
    
    # Derive individual preference maps based on selection function
    foraging.coeff <- SSFs[[1]]$coefficients
    directed.coeff <- SSFs[[2]]$coefficients
    
    new.foraging.coeff <- foraging.coeff + rnorm(4,0,0.05)
    new.directed.coeff <- directed.coeff + rnorm(4,0,0.05)
    
    SSFs[[1]]$coefficients <- new.foraging.coeff
    SSFs[[2]]$coefficients <- new.directed.coeff
  
    pred <- raster(input.stack@layers[[1]])
    pred[] <- 0
    
    forage_pred <- predict(object=SSFs[[1]], newdata=pred_df, type='risk')
    directed_pred <- predict(object=SSFs[[2]], newdata=pred_df, type='risk')
    
    pred@data@values <- forage_pred
    pred <- pred/(1 + pred)
    v_forage <- velox(pred)
    forage.maps[[j]] <- v_forage
    
    pred@data@values <- directed_pred
    pred <- pred/(1 + pred)
    v_directed <- velox(pred)
    directed.maps[[j]] <- v_directed
    
    current.state <- state.vector[j,1]
    
    # Use individual selection raster to probabilistically set starting location location
    if (current.state < 3) {
      selection.raster <- v_forage$as.RasterLayer(band=1)
    } else if (current.state == 3) {
      selection.raster <- v_directed$as.RasterLayer(band=1)
    }
    x <- getValues(selection.raster)
    x[is.na(x)] <- 0
    cellID <- sample(nrow(selection.raster)*ncol(selection.raster), size=1, prob=x)
    rast.temp <- raster(selection.raster)
    rast.temp[cellID] <- 1
    
    points <- data.frame(rasterToPoints(rast.temp, fun=function(x){x == 1}))[,1:2]
    points[1,1] <- points[1,1] + runif(1,-14.99,14.99)
    points[1,2] <- points[1,2] + runif(1,-14.99,14.99)
    #points <- SpatialPoints(points)
  
    size <- rnorm(1,sp.mean,sp.mean/8) 
    infected <- 0

    traits[j,1] <- j
    traits[j,2] <- size
    traits[j,3] <- infected
    traits[j,4] <- current.state
    traits[j,5] <- points[1,1]
    traits[j,6] <- points[1,2]
  }
  colnames(traits) <- c('ID', 'body.mass', 'infected', 'behav.state', 'x', 'y')
  ind.out <- list(traits, forage.maps, directed.maps)
  
  return(ind.out)
}

state.shift <- function(N, prev.vector, trans.mat) {
  state.vector <- data.frame(matrix(0,N,1))
  
  for (j in 1:N) {
    if (prev.vector[j] == 1) {
      probs = trans.mat[1,]
    } else if (prev.vector[j] == 2) {
      probs = trans.mat[2,]
    } else if (prev.vector[j] == 3) {
      probs = trans.mat[3,]
    }
    
    state.vector[j,1] <- sample(c(1,2,3), size=1, prob=probs)
  }
  colnames(state.vector) <- 'behav.state'
  return(state.vector)
}

extract.ranges <- function(HMM.out, states = 3) {
  ranges <- data.frame(matrix(0,states,3))
  for (i in 1:states) {
    step.mean <- HMM.out$mle[[1]][1,i]
    step.sd <- HMM.out$mle[[1]][2,i]
    step.shape <- (step.mean^2)/(step.sd^2)
    step.rate <- step.mean/(step.sd^2)
    radius <- qgamma(0.975, shape=step.shape, rate=step.rate)
    
    ranges[i,1] <- step.shape
    ranges[i,2] <- step.rate
    ranges[i,3] <- radius
  }
  colnames(ranges) <- c('gamma.shape', 'gamma.rate', 'radius')
  return(ranges)
}

anglefun <- function(xx, yy, bearing=TRUE) {
  ## calculates the compass bearing of the line between two points
  ## xx and yy are the differences in x and y coordinates between two points
  ## Options:
  ## bearing = FALSE returns +/- pi instead of 0:2*pi
  
  b<-sign(yy)
  b[b==0]<-1  #corrects for the fact that sign(0) == 0
  tempangle = b*(xx<0)*pi+atan(yy/xx)
  if(bearing){
    #return a compass bearing 0 to 2pi
    #if bearing==FALSE then a heading (+/- pi) is returned
    tempangle[tempangle<0]<-tempangle[tempangle<0]+2*pi
  }
  return(tempangle)
}

movement <- function(xy, step, heading) {
  
  pi = 3.141593
  x_init <- xy[1,1]
  y_init <- xy[1,2]
  
  if (heading < 0) {
    heading <- abs(heading) + pi
  }
  
  #rad_y <- angle*0.0174532925
  y_change <- sin(heading)*step
  y_new <- y_init + y_change
  
  # Use cosine to determine the movement in x (longitude)
  #rad_x <- angle*0.0174532925
  x_change <- cos(heading)*step
  x_new <- x_init + x_change
  
  x_init <- x_new
  y_init <- y_new
  
  move.temp <- as.data.frame(matrix(0,1,2))
  move.temp[1,1] <- x_new
  move.temp[1,2] <- y_new
  
  return(move.temp)
}

move.func <- function(N, inds, sp.mean, ranges, crs, v_bg) {
  
  moves <- data.frame(matrix(0,N,2))
  state.vector <- inds[[1]][,'behav.state']
  body.size <- inds[[1]][,'body.mass']
  
  for (j in 1:N) {
    curr.pos <- data.frame(inds[[1]][j,c('x','y')])
    curr.pos <- st_as_sf(curr.pos, coords = 1:2, crs = crs)
    curr.state <- state.vector[j]
    radius <- ranges[curr.state,3]
    percep.range <- (body.size[j]/sp.mean) * radius
    perception <- st_buffer(curr.pos, dist=percep.range)
    
    if (curr.state == 2) {
      v_selection <- inds[[2]][[j]]$copy()
    } else if (curr.state == 3) {
      v_selection <- inds[[3]][[j]]$copy()
    } else {
      v_selection <- v_bg$copy()
    }
    
    # Agents move probablistically according to their selection map within their perceptual range
    v_selection$crop(perception)
    coords <- data.frame(v_selection$getCoordinates()) %>%
      st_as_sf(coords=1:2, crs=crs) %>%
      st_intersection(., perception)
    
    nearby <- v_selection$extract(perception, df=TRUE)
    temp.pt <- sample(1:nrow(nearby), size=1, prob = nearby[,2])
    cell.center <- coords[temp.pt,] + runif(2, -14.99, 14.99)
    
    # Agents could move to the next cell with some error
    #diff_x <- st_coordinates(cell.center)[1] - st_coordinates(curr.pos)[1]
    #diff_y <- st_coordinates(cell.center)[2] - st_coordinates(curr.pos)[2]
    #heading <- anglefun(diff_x, diff_y)
    #dist <- sqrt(diff_x^2 + diff_y^2) + rnorm(1,0,30)
    #move <- movement(st_coordinates(curr.pos), dist, heading)
    
    moves[j,] <- st_coordinates(cell.center$geometry)
  }
  colnames(moves) <- c('x', 'y')
  return(moves)
}

####################################################################

name_list = c('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009',
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010',
              'AG253_2010', 'AG255_2010', 'AG256_2010')

date_list = c('20160213', '20160325', '20160426', '20160528', '20160629')

# Run general HMM across all 11 zebra tracks during the anthrax season
zebra09 <- read_csv("Zebra_Data/Zebra_Anthrax_2009_Cleaned.csv") %>% 
  dplyr::select(x,y,date,ID) %>% 
  dplyr::filter(!is.na(x)) %>%
  st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

zebra10 <- read_csv("Zebra_Data/Zebra_Anthrax_2010_Cleaned.csv") %>% 
  dplyr::select(x,y,date,ID) %>% 
  dplyr::filter(!is.na(x)) %>%
  st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

all_data <- rbind(zebra09, zebra10)

zebra.ext <- extent(c(extent(zebra10)@xmin - 5000,
                      extent(zebra10)@xmax + 5000,
                      extent(zebra09)@ymin - 5000,
                      extent(zebra10)@ymax + 5000))

zebra_ext <- as(zebra.ext, 'SpatialPolygons')
crs(zebra_ext) <- '+init=epsg:32733'
zebra_ext <- st_as_sf(zebra_ext)

ENP <- read_sf('Layers/enp fence poly.shp')
ENP_crs <- st_crs(ENP, '+proj=longlat')
st_crs(ENP) <- ENP_crs
ENP <- st_transform(ENP, '+init=epsg:32733') %>%
  st_union(.) %>% st_intersection(zebra_ext)
ENP.sp <- as(ENP, 'Spatial')

# Import predictor layers, create raster stack, normalize variables and create dataframe
Risk <- raster('Layers/Mean_Risk.tif') %>% 
  mask(ENP.sp)
road_dens <- raster('Layers/Road_Density.tif') %>% 
  mask(ENP.sp)
Green <- raster(paste0('Layers/Greenness_', date_list[1], ".tif")) %>% 
  mask(ENP.sp)
Wet <- raster(paste0('Layers/Wetness_', date_list[1], ".tif")) %>% 
  mask(ENP.sp)

pred_stack <- stack(Green, Wet, road_dens, Risk)

stack_means <- cellStats(pred_stack, stat='mean', na.rm=TRUE)
stack_sd <- cellStats(pred_stack, stat='sd', na.rm=TRUE)

norm_stack <- stack((Green - stack_means) / stack_sd[1],
                    (Wet - stack_means[2]) / stack_sd[2],
                    (road_dens - stack_means[3])/ stack_sd[3],
                    (Risk - stack_means[4])/ stack_sd[4])

norm_stack <- stack(norm_stack@layers[[1]],
                    norm_stack@layers[[2]],
                    norm_stack@layers[[3]],
                    norm_stack@layers[[4]])
names(norm_stack) <- c("Green_Norm", "Wet_Norm", "Road_Dens_Norm", "Risk_Norm")

pred_df <- as.data.frame(norm_stack)
pred_df$fix <- 1
pred_df$ID <- "AG068_2009"

pred_bg <- raster(norm_stack@layers[[1]])
pred_bg[] <- runif(1,0,1)
v_bg <- velox(pred_bg)

#### Data Recorders ####

behav.states <- list()
all.steps <- list()

#### Initialization ####

HMMs <- HMM(data = all_data, states = 3)
trans.mat <- HMMs$mle$gamma
ranges <- extract.ranges(HMMs, states = 3)

LIZs <- weighted_sampling(input.raster = Risk, N = 200)
SSFs <- selection_functions(name.list = name_list)
strt <- Sys.time()
agents <- initialize.agents(N=20, sp.mean=350, SSFs, norm_stack, pred_df)
print(Sys.time() - strt)

#### Model Implementation ###

N = 20
days = 90
t = days*24*3
crs <- "+proj=utm +south +zone=33 +ellps=WGS84"
init.state <- agents[[1]]$behav.state
behav.states[[1]] <- init.state
init.step <- agents[[1]][,c('x','y')]
all.steps[[1]] <- init.step

for (z in 2:t) {
  # Assign new states and add results to behav.states tracker
  prev.state <- behav.states[[z-1]]
  new.states <- state.shift(N, prev.state, trans.mat)
  behav.states[[z]] <- new.states[,1]
  # Find new positions and add results to all.steps tracker
  strt <- Sys.time()
  step <- move.func(N, agents, sp.mean=350, ranges, crs, v_bg)
  print(Sys.time() - strt)
  all.steps[[z]] <- step
  # Update agents object with new states and positions
  agents[[1]]$behav.state <- new.states
  agents[[1]]$x <- step[,1]
  agents[[1]]$y <- step[,2]
}

#################################################################




