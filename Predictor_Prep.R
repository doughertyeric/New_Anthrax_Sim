risk09 <- raster('Layers/Anthrax_Risk_2009.tif')
risk10 <- raster('Layers/Anthrax_Risk_2010.tif')
risk_mean <- mean(risk09, risk10)
writeRaster(risk_mean, 'Layers/Mean_Risk.tif', format='GTiff')

wet09 <- raster('Layers/Mean_Wetness_2009.tif')
wet10 <- raster('Layers/Mean_Wetness_2010.tif')
wet_mean <- mean(wet09, wet10)
writeRaster(wet_mean, 'Layers/Mean_Wetness.tif', format='GTiff')

green09 <- raster('Layers/Mean_Greenness_2009.tif')
green10 <- raster('Layers/Mean_Greenness_2010.tif')
green_mean <- mean(green09, green10)
writeRaster(green_mean, 'Layers/Mean_Greenness.tif', format='GTiff')

#################################################################

#date_list = c('20160309', '20160410', '20160512', '20160613')
date_list = c('20120203', '20120306', '20120407', '20120509')

zebra09 <- read_csv("Zebra_Data/Zebra_Anthrax_2009_Cleaned.csv") %>% 
  dplyr::select(x,y,date,ID) %>% 
  dplyr::filter(!is.na(x)) %>%
  st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

zebra10 <- read_csv("Zebra_Data/Zebra_Anthrax_2010_Cleaned.csv") %>% 
  dplyr::select(x,y,date,ID) %>% 
  dplyr::filter(!is.na(x)) %>%
  st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

zebra.ext <- extent(c(extent(zebra10)@xmin - 5000,
                      extent(zebra10)@xmax + 5000,
                      extent(zebra09)@ymin - 5000,
                      extent(zebra10)@ymax + 5000))

for (i in 1:length(date_list)) {
  RED <- raster(paste0("Raw_Layers/", date_list[i], "/B4.TIF")) %>%
    projectRaster(road_dens) %>% crop(zebra.ext) %>%
    resample(road_dens, method='bilinear')
  NIR <- raster(paste0("Raw_Layers/", date_list[i], "/B5.TIF")) %>%
    projectRaster(road_dens) %>% crop(zebra.ext) %>%
    resample(road_dens)
  NDVI <- (NIR - RED) / (NIR + RED)
  writeRaster(NDVI, paste0("Layers/NDVI_", date_list[i]), format="GTiff", overwrite=TRUE)
  print(i)
}

#Brightness.vals <- c(0.3037, 0.2793, 0.4343, 0.5585, 0.5082, 0.1863)
Greenness.vals <-	c(-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800)
Wetness.vals <-	c(0.1509, 0.1793, 0.3299, 0.3406, -0.7112, -0.4572)

for (i in 1:length(date_list)) {
  CH1 <- raster(paste0("Raw_Layers/", date_list[i], "/B1.TIF")) %>%
    projectRaster(road_dens) %>% crop(zebra.ext) %>%
    resample(road_dens, method='bilinear')
  CH2 <- raster(paste0("Raw_Layers/", date_list[i], "/B2.TIF")) %>%
    projectRaster(road_dens) %>% crop(zebra.ext) %>%
    resample(road_dens, method='bilinear')
  CH3 <- raster(paste0("Raw_Layers/", date_list[i], "/B3.TIF")) %>%
    projectRaster(road_dens) %>% crop(zebra.ext) %>%
    resample(road_dens, method='bilinear')
  CH4 <- raster(paste0("Raw_Layers/", date_list[i], "/B4.TIF")) %>%
    projectRaster(road_dens) %>% crop(zebra.ext) %>%
    resample(road_dens, method='bilinear')
  CH5 <- raster(paste0("Raw_Layers/", date_list[i], "/B5.TIF")) %>%
    projectRaster(road_dens) %>% crop(zebra.ext) %>%
    resample(road_dens, method='bilinear')
  CH7 <- raster(paste0("Raw_Layers/", date_list[i], "/B7.TIF")) %>%
    projectRaster(road_dens) %>% crop(zebra.ext) %>%
    resample(road_dens, method='bilinear')
  Greenness <- (Greenness.vals[1] * CH1) + (Greenness.vals[2] * CH2) + (Greenness.vals[3] * CH3) + 
    (Greenness.vals[4] * CH4) + (Greenness.vals[5] * CH5) + (Greenness.vals[6] * CH7)
  writeRaster(Greenness, paste0("Layers/Greenness_",date_list[i]), format="GTiff", overwrite=TRUE)
  Wetness <- (Wetness.vals[1] * CH1) + (Wetness.vals[2] * CH2) + (Wetness.vals[3] * CH3) + 
    (Wetness.vals[4] * CH4) + (Wetness.vals[5] * CH5) + (Wetness.vals[6] * CH7)
  writeRaster(Wetness, paste0("Layers/Wetness_",date_list[i]), format="GTiff", overwrite=TRUE)
  print(i)
}
