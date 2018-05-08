library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)

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

all.step.means <- read.csv('Final_Runs_Markov/Markov_NoPan_Step_Means_POM.csv')

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

#######################################################

layers <- read.csv('Layer_Stats.csv')

pdf('Figures/Environmental_Covariate_Grid.pdf', width=8, height=6)
ggplot(data=layers, aes(x=Green, y=Wet)) + 
  geom_point() + 
  geom_text(aes(label=ID), nudge_x = 1) + 
  geom_point(aes(x=-21.00, y=-36.89, size=1), col='red') + 
  theme_bw() + 
  theme(legend.position = 'none') + 
  xlab('Mean Greenness') + 
  ylab('Mean Wetness')
dev.off()

########################################################

ENP <- read_sf('Layers/enp fence poly.shp')
ENP_crs <- st_crs(ENP, '+proj=longlat')
st_crs(ENP) <- ENP_crs
ENP <- st_transform(ENP, '+init=epsg:32733') %>%
     st_union(.) %>% st_intersection(zebra_ext)

pans <- read_sf('Layers/enp pans.shp')
pans_crs <- st_crs(pans, '+proj=longlat')
st_crs(pans) <- pans_crs
pans <- st_transform(pans, '+init=epsg:32733') %>%
  st_union(.) %>% st_intersection(zebra_ext)

ENP_nopans <- st_difference(ENP, pans) %>% 
  st_union %>% 
  st_sf()

LIZ_sf <- read_csv('LIZ_buffers.csv') %>%
  st_as_sf(coords = 1:2, crs=32733)

pdf('Figures/LIZ_Sites.pdf', width=8, height=6)
ggplot() + 
  geom_sf(data=ENP_nopans) + 
  geom_sf(data=LIZ_sf, aes(color=factor(species,
                                        levels=c(1,2,3),
                                        labels=c('Zebra', 'Springbok', 'Elephant')),
                           stroke=0)) + 
  theme_bw() + 
  theme(legend.position = 'none') +
  #theme(legend.position = 'bottom') + 
  #guides(color = guide_legend(title="Species")) +
  coord_sf(crs = st_crs(ENP_nopans), datum = NA)
dev.off()

#########################################################
###########             Results                 #########
#########################################################

quad.lm <- lm(mean ~ Wet + Wet2, data=layers)
wetvals <- seq(-80,-20,1)
pred.list <- list(Wet=wetvals, Wet2=wetvals^2)
predicted <- predict(quad.lm, pred.list)

pred <- data.frame(cbind(Wet=wetvals, Wet2=wetvals^2,
                         ContactRate = predicted))

p1 <- ggplot() + 
  geom_point(data=layers, aes(x=-Wet, y=mean)) +
  geom_line(data=pred, aes(x=-Wet, y=ContactRate)) + 
  scale_x_continuous(labels = c('Wet', 'Moderately Wet',
                                'Moderately Dry', 'Dry')) +
  ylab('Mean Contact Rate') + xlab('Mean Wetness') +
  theme_bw()

pdf('Figures/Contact_Rate_Trend.pdf', width=8, height=4)
p1
dev.off()

carcass <- read_csv('Carcass_Data.csv') %>%
  filter(BAPos == 'POS') %>%
  dplyr::select(edod, Species, BAPos) %>%
  mutate(Month = lubridate::month(edod),
         Year = lubridate::year(edod)) %>%
  filter(Month > 1) %>%
  filter(Month < 7) %>%
  filter(Species == "EB")

carcass$Month <- factor(carcass$Month,
                        levels=c(2,3,4,5,6),
                        labels = c('Feb', 'Mar', 'Apr',
                                   'May', 'Jun'))

p2 <- ggplot() +
  geom_histogram(data=carcass, aes(x=Month), stat='count') +
  ylab('Confirmed Anthrax Mortalities') + xlab('Month')

carcass <- read_csv('ENP_carcasses.csv') %>%
  dplyr::select(DATE, SPECIES) %>%
  mutate(Date = lubridate::mdy(DATE),
         Month = lubridate::month(Date),
         Year = lubridate::year(Date)) %>%
  filter(Month > 1) %>%
  filter(Month < 7) %>%
  filter(SPECIES == "EB") %>%
  filter(Year > 2003) %>%
  filter(Year < 2010)

carcass$Month <- factor(carcass$Month,
                        levels=c(2,3,4,5,6),
                        labels = c('Feb', 'Mar', 'Apr',
                                   'May', 'Jun'))

p3 <- ggplot() +
  geom_histogram(data=carcass, aes(x=Month), stat='count') +
  ylab('Confirmed Anthrax Mortalities (1996-2009)') + xlab('Month')

pdf('Figures/Contact_Rate_Trend_v2.pdf', width=8, height=8)
grid.arrange(p1, p3, ncol=1)
dev.off()

rainfall <- read_csv('ENP_Monthly_Rain.csv') %>%
  dplyr::select(Year, Month, Okaukuejo) %>%
  filter(Month %in% c('feb', 'mar', 'apr',
                      'may', 'jun')) %>%
  filter(Year > 1995) %>%
  filter(Year < 2010) %>%
  group_by(Month) %>%
  summarise(avg.rain = mean(Okaukuejo),
            sd.rain = sd(Okaukuejo),
            se.rain = sd.rain/n())

rainfall$Month <- factor(rainfall$Month,
                         levels=c('feb', 'mar', 'apr',
                                  'may', 'jun'),
                         labels = c('Feb', 'Mar', 'Apr',
                                    'May', 'Jun'))

p4 <- ggplot(data=rainfall, aes(x=Month, y=avg.rain)) + 
  geom_errorbar(aes(ymin=avg.rain-se.rain, ymax=avg.rain+se.rain), width=.2) +
  geom_point() +
  geom_line(aes(group=1)) +
  ylab('Mean Rainfall (cm)')

pdf('Figures/Contact_Rate_Trend_v3.pdf', width=8, height=12)
grid.arrange(p1, p3, p4, ncol=1)
dev.off()

vp <- viewport(width = 0.35, height = 0.3, 
               x = 0.8, y = 0.75)
full <- function() {
  print(p3)
  theme_set(theme_bw(base_size = 8))
  print(p4, vp = vp)
  theme_set(theme_bw())
}

pdf('Figures/Mortality_Rainfall.pdf', width=8, height=4)
full()
dev.off()
