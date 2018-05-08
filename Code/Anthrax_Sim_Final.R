library(raster)
library(survival)
library(moveHMM)
library(dplyr)
library(magrittr)
library(readr)
library(rgeos)
library(foreach)
library(doParallel)
library(sf)
library(velox)

name_list = c('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009',
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010',
              'AG253_2010', 'AG255_2010', 'AG256_2010')

date_list = c('20120306', '20120407', '20120509')

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

####################################################################
############                Functions                 ##############
####################################################################

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

save.paths <- function(all.steps, N, t) {
  all.paths <- list()
  start.day <- as.Date("2010/02/01", format="%Y/%m/%d")
  for (l in 1:length(all.steps)) {
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

####################################################################

env_covariates <- readRDS('Covariates_Reduced_New.rds')
HMMs <- readRDS('Population_Hidden_Markov_Model.rds')

extract.ranges <- function(HMM.out, states = 3) {
  ranges <- data.frame(matrix(0,states,4))
  for (i in 1:states) {
    step.mean <- HMM.out$mle[[1]][1,i]
    step.sd <- HMM.out$mle[[1]][2,i]
    step.shape <- (step.mean^2)/(step.sd^2)
    step.rate <- step.mean/(step.sd^2)
    #radius.small <- qgamma(0.841, shape=step.shape, rate=step.rate)
    radius.small <- qgamma(0.682, shape=step.shape, rate=step.rate)
    radius.med <- qgamma(0.95, shape=step.shape, rate=step.rate)
    radius.large <- qgamma(0.997, shape=step.shape, rate=step.rate)
    
    ranges[i,1] <- step.shape
    ranges[i,2] <- step.rate
    ranges[i,3] <- radius.small
    ranges[i,4] <- radius.med
    ranges[i,5] <- radius.large
  }
  colnames(ranges) <- c('gamma.shape', 'gamma.rate', 'radius.small', 'radius.med', 'radius.large')
  return(ranges)
}

trans.mat <- HMMs$mle$gamma
ranges <- extract.ranges(HMMs, states = 3)
radii.small <- as.vector(ranges$radius.small)
radii.med <- as.vector(ranges$radius.med)
radii.large <- as.vector(ranges$radius.large)

#### Initialization ####

k = 1
N = 505
days = 150
t = days*24*3
crs <- "+init=epsg:32733"
norm_stack <- env_covariates[[k]]
#Risk <- raster('Layers/Mean_Risk.tif') %>% 
#  mask(ENP.sp)

#### Model Implementation ###

SSFs <- selection_functions(name.list = name_list)
maps <- initialize.maps(SSFs, norm_stack)
pred_bg <- raster(norm_stack@layers[[1]])
pred_bg[] <- runif(1,0,1)
v_bg <- velox(pred_bg)

num_cores <- detectCores()
cl <- makeCluster(24, outfile = '')
registerDoParallel(cl)

loop_out <- foreach(z = 1:N, .packages=c('raster', 'dplyr', 'magrittr', 'sf', 'rgeos', 'velox'), .errorhandling = 'remove') %dopar% {

  strt <- Sys.time()
  
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
  
  initial.position <- function(state.vector, maps, add.mask) {
    
    current.state <- state.vector[1,1]
    
    # Use individual selection raster to probabilistically set starting location
    if (current.state < 3) {
      selection.raster <- maps[[1]][[1]]$as.RasterLayer(band=1)
    } else if (current.state == 3) {
      selection.raster <- maps[[2]][[1]]$as.RasterLayer(band=1)
    }

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
    st_crs(add.mask) <- 32733

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
  
  move.func <- function(t, start.state, trans.mat, start.pos, small, med, large, crs, maps, v_bg) {
    
    moves <- data.frame(matrix(0,t,3))
    colnames(moves) <- c('x', 'y', 'behav.state')
    moves[1,c('x','y')] <- start.pos
    moves[1,'behav.state'] <- start.state

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
      
      perception <- st_buffer(prev.pos, dist=radius)
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
  
  start.state <- initialize.states(1)
  start.pos <- initial.position(start.state, maps, pans)
  print(paste('initialized.position:',z))

  print('starting.moves')
  moves <- move.func(t, start.state, trans.mat, start.pos, radii.small, 
                     radii.med, radii.large, crs, maps, v_bg)
  
  print(paste(Sys.time() - strt, ":", z))
  return(moves)
}

stopCluster(cl)

all_paths <- save.paths(loop_out, N, t)
saveRDS(all_paths, paste0('ABM_Final_Markov_k',k,'.rds'))

#################################################################
