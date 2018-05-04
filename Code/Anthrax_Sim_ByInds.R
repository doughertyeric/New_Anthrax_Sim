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

build.stack <- function(date_list, mask.sp) {
  
  stack_list <- list()
  
  # Import predictor layers, create raster stack, normalize variables and create dataframe
  Risk <- raster('Layers/Mean_Risk_v2.tif') %>% 
    mask(mask.sp)
  road_dens <- raster('Layers/Road_Density.tif') %>% 
    mask(mask.sp)
  Green <- raster('Layers/Mean_Greenness.tif') %>% 
    mask(mask.sp)
  Wet <- raster('Layers/Mean_Wetness.tif') %>% 
    mask(mask.sp)
  
  original_stack <- raster::stack(Green, Wet, road_dens, Risk)
  
  original_means <- cellStats(original_stack, stat='mean', na.rm=TRUE)
  original_sd <- cellStats(original_stack, stat='sd', na.rm=TRUE)
  
  k = 1
  for (i in 1:length(date_list)) {
    for (j in 1:length(date_list)) {
      Green <- raster(paste0('Layers/Greenness_', date_list[i], '.tif')) %>% 
        mask(mask.sp)
      Wet <- raster(paste0('Layers/Wetness_', date_list[j], '.tif')) %>% 
        mask(mask.sp)
      
      new_stack <- raster::stack(Green, Wet, road_dens, Risk)
      new_means <- cellStats(new_stack, stat='mean', na.rm=TRUE)
      new_sd <- cellStats(new_stack, stat='sd', na.rm=TRUE)
      
      norm_stack <- raster::stack((Green - original_means[1]) / original_sd[1],
                          (Wet - original_means[2]) / original_sd[2],
                          (road_dens - original_means[3])/ original_sd[3],
                          (Risk - original_means[4])/ original_sd[4])
      
      norm_stack <- raster::stack(norm_stack@layers[[1]],
                          norm_stack@layers[[2]],
                          norm_stack@layers[[3]],
                          norm_stack@layers[[4]])
      names(norm_stack) <- c("Green_Norm", "Wet_Norm", "Road_Dens_Norm", "Risk_Norm")
      assign(paste0("stack_Green0",i,"_Wet0",j), norm_stack)
      
      stack_list[[k]] <- get(paste0("stack_Green0",i,"_Wet0",j))
      print(k)
      k = k + 1
    }
  }
  return(stack_list)
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

initialize.maps <- function(SSFs, input.stack) {
  
  forage.map <- list()
  directed.map <- list()
  
  forage_clogit <- SSFs[[1]]
  directed_clogit <- SSFs[[2]]
  
  pred <- raster(input.stack@layers[[1]])
  pred[] <- 0
  
  pred_df <- as.data.frame(input.stack)
  pred_df$fix <- 1
  pred_df$ID <- "AG068_2009"
  
  forage_pred <- predict(object=forage_clogit, newdata=pred_df, type='risk')
  directed_pred <- predict(object=directed_clogit, newdata=pred_df, type='risk')
  
  pred@data@values <- forage_pred
  pred <- pred/(1 + pred)
  v_forage <- velox(pred)
  forage.map[[1]] <- v_forage
  
  pred@data@values <- directed_pred
  pred <- pred/(1 + pred)
  v_directed <- velox(pred)
  directed.map[[1]] <- v_directed
  
  out_list <- list(forage.map, directed.map)
  return(out_list)
}

initial.position <- function(state.vector, maps, add.mask) {

  current.state <- state.vector[1,1]
  
  # Use individual selection raster to probabilistically set starting location
  if (current.state < 3) {
    selection.raster <- maps[[1]][[1]]$as.RasterLayer(band=1)
  } else if (current.state == 3) {
    selection.raster <- maps[[2]][[1]]$as.RasterLayer(band=1)
  }
  #selection.raster <- mask(selection.raster, add.mask)
  x <- getValues(selection.raster)
  x[is.na(x)] <- 0
  cellID <- sample(nrow(selection.raster)*ncol(selection.raster), size=1, prob=x)
  rast.temp <- raster(selection.raster)
  rast.temp[cellID] <- 1
  
  points <- data.frame(rasterToPoints(rast.temp, fun=function(x){x == 1}))[,1:2]
  points[1,1] <- points[1,1] + runif(1,-14.99,14.99)
  points[1,2] <- points[1,2] + runif(1,-14.99,14.99)
  points <- st_as_sf(points, coords=1:2)
  st_crs(points) <- 32733
  
  while (nrow(st_intersection(points, add.mask)) != 0) {
    cellID <- sample(nrow(selection.raster)*ncol(selection.raster), size=1, prob=x)
    rast.temp <- raster(selection.raster)
    rast.temp[cellID] <- 1
    
    points <- data.frame(rasterToPoints(rast.temp, fun=function(x){x == 1}))[,1:2]
    points[1,1] <- points[1,1] + runif(1,-14.99,14.99)
    points[1,2] <- points[1,2] + runif(1,-14.99,14.99)
    points <- st_as_sf(points, coords=1:2)
    st_crs(points) <- 32733
  }
  
  points <- data.frame(st_coordinates(points$geometry))
  colnames(points) <- c('x', 'y')
  return(points)
}

designate.size <- function(sp.mean) {
  return(rnorm(1, sp.mean, sp.mean/8))
}

state.shift <- function(prev.vector, trans.mat) {
  state.vector <- data.frame(matrix(0,1,1))
  if (prev.vector[1] == 1) {
    probs <- trans.mat[1,]
  } else if (prev.vector[1] == 2) {
    probs <- trans.mat[2,]
  } else if (prev.vector[1] == 3) {
    probs <- trans.mat[3,]
  }
  state.vector[1,1] <- sample(c(1,2,3), size=1, prob=probs)
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

move.func <- function(t, start.state, trans.mat, size, sp.mean, start.pos, small, med, large, crs, maps, v_bg) {
    
    moves <- data.frame(matrix(0,t,3))
    colnames(moves) <- c('x', 'y', 'behav.state')
    moves[1,c('x','y')] <- start.pos
    moves[1,'behav.state'] <- start.state
    body.size <- size
    
    for (j in 2:t) {
      prev.pos <- moves[(j-1), c('x','y')]
      prev.pos <- st_as_sf(prev.pos, coords = 1:2, crs = crs)
      prev.state <- moves[(j-1), 3]
      new.state <- state.shift(prev.state, trans.mat)
      rand <- runif(1,0,1)
      if (rand <= 0.68) {
        radius <- small[new.state[1,1]]
      } else if (rand > 0.68 && rand <= 0.95) {
        radius <- med[new.state[1,1]]
      } else {
        radius <- large[new.state[1,1]]
      }
      
      percep.range <- (body.size/sp.mean) * radius
      perception <- st_buffer(prev.pos, dist=percep.range)
      st_crs(perception) <- 32733

      if (new.state[1,1] == 2) {
        v_selection <- maps[[1]][[1]]$copy()
      } else if (new.state[1,1] == 3) {
        v_selection <- maps[[1]][[1]]$copy()
      } else {
        v_selection <- v_bg$copy()
      }
      
      # Agents move probablistically according to their selection map within their perceptual range
      v_selection$crop(perception)
      coords <- data.frame(v_selection$getCoordinates()) %>%
        st_as_sf(coords=1:2)
      st_crs(coords) <- 32733
      coords <- st_intersection(coords, perception)

      if (nrow(coords) != 0) {
        nearby <- v_selection$extract_points(coords)
        nearby[is.na(nearby)] <- 0
        temp.pt <- sample(1:nrow(nearby), size=1, prob = nearby[,1])
        cell.center <- coords[temp.pt,] + runif(2, -14.99, 14.99)
      } else {
        temp.coords <- st_coordinates(prev.pos$geometry)
        cell.center <- temp.coords + runif(2, -14.99, 14.99)
        cell.center <- data.frame(cell.center)
        cell.center <- st_as_sf(cell.center, coords=1:2)
        st_crs(cell.center) <- 32733
      }

      moves[j,c('x','y')] <- st_coordinates(cell.center$geometry)
      moves[j,'behav.state'] <- new.state
    }
    return(moves)
  }

save.paths <- function(all.steps, N, t) {
  all.paths <- list()
  start.day <- as.Date("2010/02/01", format="%Y/%m/%d")
  for (l in 1:N) {
    ind.steps <- data.frame(matrix(0,t,4))
    ind.steps[,1:3] <- all.steps[[l]]
    for (m in 1:t) {
      day <- m %/% 72
      date <- start.day + day
      hour <- ((m - ((day-1)*72)) %/% 3) - 24
      hour.char <- if (hour < 10) {paste0("0",hour)} else {as.character(hour)}
      min <- (m - ((day-1)*72)) %% 3
      ind.steps[m,4] <- paste0(date, " ", hour.char, ":", if (min==0) {"00"} else {(min*20)})
    }
    colnames(ind.steps) <- c('x', 'y', 'state', 'datetime')
    all.paths[[l]] <- data.frame(ind.steps)
  }
  return(all.paths)
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
  path.dens <- data.frame(matrix(0,0,2))
  for (i in 1:(nrow(xy) - 1)) {
    dist.x <- abs(xy[i,1] - xy[(i+1),1])
    dist.y <- abs(xy[i,2] - xy[(i+1),2])
    current.dist <- sqrt(dist.x^2 + dist.y^2)
    
    v_max <- rtruncnorm(n=1,mean=1.1,sd=0.05,a=1.01)*current.dist
    path.dens.temp <- sim.paths(xy[i,], xy[(i+1),], v_max, 21, crs)
    path.dens <- data.frame(rbind(path.dens, path.dens.temp))
    print(i)
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

####################################################################

name_list = c('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009',
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010',
              'AG253_2010', 'AG255_2010', 'AG256_2010')

#date_list = c('20160309', '20160410', '20160512', '20160613')
#date_list = c('20120203', '20120306', '20120407', '20120509')
#date_list = c('20120306', '20120407', '20120509')
date_list = c('20100205', '20100528', '20090525', '20090322')

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

pans <- read_sf('Layers/enp pans.shp')
pans_crs <- st_crs(pans, '+proj=longlat')
st_crs(pans) <- pans_crs
pans <- st_transform(pans, '+init=epsg:32733') %>%
  st_union(.) %>% st_intersection(zebra_ext)
pans.sp <- as(pans, 'Spatial')

ENP_nopans <- st_difference(ENP, pans) %>% 
  st_union %>% 
  st_sf() %>%
  as('Spatial')

#### Overarching Objects/Parameters ####

#HMMs <- HMM(data = all_data, states = 3)
HMMs <- readRDS('Population_Hidden_Markov_Model.rds')
trans.mat <- HMMs$mle$gamma
ranges <- extract.ranges(HMMs, states = 3)
radii <- as.vector(ranges$radius)

#env_covariates <- build.stack(date_list, ENP.sp)
env_covariates <- readRDS('Covariate_Stacks_Reduced.rds')

#### Initialization ####

k = 1
N = 20
days = 1
t = days*24*3
crs <- "+proj=utm +south +zone=33 +ellps=WGS84"
norm_stack <- env_covariates[[k]]
Risk <- raster('Layers/Mean_Risk.tif') %>% 
  mask(ENP.sp)

#### Model Implementation ###

LIZs <- weighted_sampling(input.raster = Risk, N = 200)
SSFs <- selection_functions(name.list = name_list)
maps <- initialize.maps(SSFs, norm_stack)
pred_bg <- raster(norm_stack@layers[[1]])
pred_bg[] <- runif(1,0,1)
v_bg <- velox(pred_bg)

strt <- Sys.time()
all_inds <- list()
for (z in 1:N) {
  start.state <- initialize.states(1)
  start.pos <- initial.position(start.state, maps, pans)
  size <- designate.size(sp.mean=350)
  moves <- move.func(t, start.state, trans.mat, size, sp.mean=350,
                     start.pos, radii, crs, maps, v_bg)
  all_inds[[z]] <- moves
}
print(Sys.time() - strt)

all_paths <- save.paths(all_inds, N, t)

#################################################################
