library(ggplot2)
library(ggpubr)

########################################################
########     Perceptual Range Schematic         ########
########################################################

env_covariates <- readRDS('Covariates_Reduced_New.rds')
norm_stack <- env_covariates[[1]]

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

forage_clogit <- SSFs[[1]]

pred <- raster(input.stack@layers[[1]])
pred[] <- 0
pred_df <- as.data.frame(norm_stack)
pred_df$fix <- 1
pred_df$ID <- "AG068_2009"

forage_pred <- predict(object=forage_clogit, newdata=pred_df, type='risk')
pred@data@values <- forage_pred
pred <- pred/(1 + pred)
v_temp <- velox(pred)
v_selection <- v_temp$copy()

prev.pos <- data.frame(matrix(c(638111.219451494, 7863348.8914231),1,2))
prev.pos <- st_as_sf(prev.pos, coords = 1:2, crs = 32733)

percep.range <- 154.7842
perception <- st_buffer(prev.pos, dist=percep.range)
st_crs(perception) <- 32733

v_selection$crop(perception)
coords <- data.frame(v_selection$getCoordinates()) %>%
  st_as_sf(coords=1:2)
st_crs(coords) <- 32733
coords <- st_intersection(coords, perception)

pred.temp <- crop(pred, extent(as(perception,'Spatial')))

pdf('Figures/Perceptual_Range_Schematic.pdf', width=8, height=6)
plot(pred.temp)
plot(perception$geometry, add=TRUE)
plot(coords$geometry, pch=19, cex=0.3, add=TRUE)
plot(prev.pos$geometry, pch=19, col='red', add=TRUE)
dev.off()

###############################################################
#######           Step Distribution Fig               #########
###############################################################

HMMs <- readRDS('Population_Hidden_Markov_Model.rds')

extract.ranges <- function(HMM.out, states = 3) {
  ranges <- data.frame(matrix(0,states,4))
  for (i in 1:states) {
    step.mean <- HMM.out$mle[[1]][1,i]
    step.sd <- HMM.out$mle[[1]][2,i]
    step.shape <- (step.mean^2)/(step.sd^2)
    step.rate <- step.mean/(step.sd^2)
    
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

for (i in 1:3) {
  step.mean <- HMMs$mle[[1]][1,i]
  step.sd <- HMMs$mle[[1]][2,i]
  step.shape <- (step.mean^2)/(step.sd^2)
  step.rate <- step.mean/(step.sd^2)
  
  sim.steps <- data.frame(rgamma(10000, shape=step.shape, rate=step.rate))
  for (j in 1:nrow(sim.steps)) {
    if (sim.steps[j,1] < ranges$radius.small[i]) {
      sim.steps$size[j] <- 'small'
    } else if (sim.steps[j,1] > ranges$radius.small[i] && sim.steps[j,1] < ranges$radius.med[i]) {
      sim.steps$size[j] <- 'medium'
    } else {
      sim.steps$size[j] <- 'large'
    }
  }
  colnames(sim.steps) <- c('steps', 'size')
  
  p <- ggplot() + geom_histogram(data=sim.steps, aes(x=steps, fill=size), bins=25) +
    ylab(NULL) + xlab('Step Size (m)') + guides(fill=guide_legend(title=NULL))
  assign(paste0('p',i), p)
}

pdf('Figures/Step_Size_Distributions.pdf', width=8, height=8)
ggarrange(p1, p2, p3,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)
dev.off()

#####################################################
########         POM Mean Steps              ########
#####################################################

true_foraging_mu <- c(193.6, 78.3, 68, 188.1, 173.4, 205.6,
                      187, 110.9, 154.1, 97.1, 108.4)
true_directed_mu <- c(853.8, 405.8, 386.4, 758, 907.7, 818.9,
                      852, 600.1, 785.5, 560.1, 581.8)

all.step.means <- read.csv('Markov_Step_Means_POM.csv')

all.step.means$k <- as.factor(all.step.means$k)

pdf('Step_Means_POM_Foraging.pdf', height=8, width=8)
ggplot() + geom_histogram(aes(x=rnorm(1000,mean(true_foraging_mu),sd(true_foraging_mu))), bins=20) +
   geom_histogram(data=all.step.means, aes(x=foraging.mu, fill=k), bins=20) +
   xlab('Mean Step Length (Foraging)') + ylab('Frequency')
dev.off()

pdf('Step_Means_POM_Directed.pdf', height=8, width=8)
ggplot() + geom_histogram(aes(x=rnorm(1000,mean(true_directed_mu),sd(true_directed_mu))), bins=20) +
   geom_histogram(data=all.step.means, aes(x=directed.mu, fill=k), bins=20) +
   xlab('Mean Step Length (Directed)') + ylab('Frequency')
dev.off()
