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

env_covariates <- readRDS('Covariate_Stacks_Reduced.rds')
HMMs <- readRDS('Population_Hidden_Markov_Model.rds')

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

trans.mat <- HMMs$mle$gamma
ranges <- extract.ranges(HMMs, states = 3)
radii <- as.vector(ranges$radius)


####################################################################
###########              Parallel Loop                 #############
####################################################################

num_cores <- detectCores()
cl <- makeCluster(9, outfile = '')
registerDoParallel(cl)

# Begin parallel loop to create movement paths for N agents over each of the k stacks
loop_out <- foreach(k = 1:9, .packages=c('raster', 'survival', 'dplyr', 'magrittr', 'readr', 'sf', 'rgeos', 'velox')) %dopar% {
  
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
  
  # Need to consider how to accurately model individual variation in selection
  initialize.agents <- function(N, sp.mean, SSFs, input.stack, add.mask) {
    
    traits <- data.frame(matrix(0,N,6)) #traitframe
    state.vector <- initialize.states(N)
    forage.maps <- list()
    directed.maps <- list()
    
    # Loop over individuals of a given species
    for (j in 1:N) {
      
      forage_clogit <- SSFs[[1]]
      directed_clogit <- SSFs[[2]]
      
      # Derive individual preference maps based on selection function
      foraging.coeff <- SSFs[[1]]$coefficients
      directed.coeff <- SSFs[[2]]$coefficients
      
      new.foraging.coeff <- foraging.coeff + rnorm(4,0,0.05) ####
      new.directed.coeff <- directed.coeff + rnorm(4,0,0.05) ####
      
      forage_clogit$coefficients <- new.foraging.coeff
      directed_clogit$coefficients <- new.directed.coeff
      
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
      forage.maps[[j]] <- v_forage
      
      pred@data@values <- directed_pred
      pred <- pred/(1 + pred)
      v_directed <- velox(pred)
      directed.maps[[j]] <- v_directed
      
      current.state <- state.vector[j,1]
      
      # Use individual selection raster to probabilistically set starting location
      if (current.state < 3) {
        selection.raster <- v_forage$as.RasterLayer(band=1)
      } else if (current.state == 3) {
        selection.raster <- v_directed$as.RasterLayer(band=1)
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
  
  move.func <- function(N, inds, sp.mean, radii, crs, v_bg) {
    
    moves <- data.frame(matrix(0,N,2))
    state.vector <- inds[[1]][,'behav.state']
    body.size <- inds[[1]][,'body.mass']
    
    for (j in 1:N) {
      curr.pos <- data.frame(inds[[1]][j,c('x','y')])
      curr.pos <- st_as_sf(curr.pos, coords = 1:2, crs = crs)
      curr.state <- state.vector[j]
      radius <- radii[curr.state]
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
      
      nearby <- v_selection$extract_points(coords)
      nearby[is.na(nearby)] <- 0
      temp.pt <- sample(1:nrow(nearby), size=1, prob = nearby[,1])
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
  
  save.paths <- function(all.steps, N, t) {
    all.paths <- list()
    for (l in 1:N) {
      ind.steps <- data.frame(matrix(0,length(all.steps),4))
      start.day <- as.Date("2010/02/01", format="%Y/%m/%d")
      for (m in 1:length(all.steps)) {
        pos <- all.steps[[m]]
        state <- behav.states[[m]]
        ind.steps[m,1] <- pos[l,1]
        ind.steps[m,2] <- pos[l,2]
        ind.steps[m,3] <- state[l]
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
  
  #### Initialization ####
  
  norm_stack <- env_covariates[[k]]
  Risk <- raster('Layers/Mean_Risk.tif') %>% 
    mask(ENP.sp)
  SSFs <- selection_functions(name.list = name_list)
  strt <- Sys.time()
  agents <- initialize.agents(N=10, sp.mean=350, SSFs, norm_stack, pans)
  print(Sys.time() - strt)
  
  #### Data Recorders ####
  
  behav.states <- list()
  all.steps <- list()
  
  #### Model Implementation ###
  
  N = 10
  days = 1
  t = days*24*3
  crs <- "+proj=utm +south +zone=33 +ellps=WGS84"
  init.state <- agents[[1]]$behav.state
  behav.states[[1]] <- init.state
  init.step <- agents[[1]][,c('x','y')]
  all.steps[[1]] <- init.step
  norm_stack <- env_covariates[[k]]
  
  pred_bg <- raster(norm_stack@layers[[1]])
  pred_bg[] <- runif(1,0,1)
  v_bg <- velox(pred_bg)
  
  for (z in 2:t) {
    # Assign new states and add results to behav.states tracker
    prev.state <- behav.states[[z-1]]
    new.states <- state.shift(N, prev.state, trans.mat)
    behav.states[[z]] <- new.states[,1]
    # Find new positions and add results to all.steps tracker
    strt <- Sys.time()
    step <- move.func(N, agents, sp.mean=350, radii, crs, v_bg)
    #step <- move.func.par(N, agents, sp.mean=350, radii, crs, v_bg)
    print(paste(Sys.time() - strt,':',z))
    all.steps[[z]] <- step
    # Update agents object with new states and positions
    agents[[1]]$behav.state <- new.states[,1]
    agents[[1]]$x <- step[,1]
    agents[[1]]$y <- step[,2]
  }
  
  all_paths <- save.paths(all.steps, N, t)
  return(all_paths)
}

saveRDS(loop_out, 'ABM_Test.rds')




