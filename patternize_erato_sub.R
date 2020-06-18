library(patternize)
library(viridis)

# H. erato
IDList_era   <- c('BC0147','BC0148','BC0149','BC0154','BC0163','BC0198','BC0200','BC0327','BC0340','BC0351')
IDList_hydFG <- c('BC0004','BC0049','BC0050','BC0061','BC0071','BC0076','BC0077','BC0079','BC0082','BC0125')
IDList_not   <- c('CAM016057','CAM016058','CAM016060','CAM016067','CAM016797','CAM016894','CAM016900','CAM016915','CAM017174','CAM017177')
IDList_lat   <- c('CAM016583','CAM016586','CAM016973','CAM016999','CAM017013','CAM017027','CAM017124','CAM017170','CAM017395','CAM017405')
IDList_emm   <- c('BC2563','BC2578','BC2579','BC2580','BC2611','BC2612','BC2620','BC2624','10429077','10429078')
IDList_ety   <- c('BC3000','BC3001','BC3002','BC3003','BC3006','BC3007','BC3009','BC3010','10429062','10429064')
IDList_fav   <- c('BC2634','BC2635','BC2636','BC2637','BC2638','BC2640','BC2643','BC2646','10428999')
IDList_hydP  <- c('IMG_1855','IMG_1861','IMG_1859','IMG_1863','IMG_1869','IMG_1871','IMG_1877','IMG_1885','IMG_1899','IMG_1993')
IDList_dem   <- c('IMG_1960','IMG_1972','IMG_1974','IMG_1978','IMG_1980','IMG_1982','IMG_1987','IMG_2049','IMG_2125','IMG_2135','IMG_2141')
IDList_phy   <- c('6_MARA','10_MARA','35_SAG','49_SAG','52_SAG','85_IG','112_BOQ2','114_BOQ2','10428488','10428490')
IDList_cyr   <- c('CAM040321','CAM040354','CAM040355','CAM040356','CAM040357','CAM040361','CAM040362','CAM040363','CAM040364','CAM040365')
IDList_ven   <- c('CS000046','CS000047','CS000265','CS000275','CS000278','CS000284','CS000286','CS000288','CS003656','CS003659')
IDList_amal  <- c('BC2108','BC2146','BC2168','BC2181','BC2184','BC2187','BC2193','BC2223','BC2227','BC2351')
IDList_mic   <- c('10428928','10428929','10428930','10428931','10428932','10428934','10428935','10428936','10428937','10428941')

# Landmarks
landList_hydP <- makeList(IDList_hydP, 'landmark', 'landmarks/H.e.hydaraP', '_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_hydFG <- makeList(IDList_hydFG, 'landmark', 'landmarks/H.e.hydaraFG', '_landmarks_LFW.txt', skipLandmark = c(2:5,7:9))
landList_amal <- makeList(IDList_amal, 'landmark', 'landmarks/H.e.amalfreda', '-D.txt', skipLandmark = c(2:5,7:9))
landList_cyr <- makeList(IDList_cyr, 'landmark', 'landmarks/H.e.cyrbia', '_d_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_not <- makeList(IDList_not, 'landmark', 'landmarks/H.e.notabilis', '_d_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_ety <- makeList(IDList_ety, 'landmark', 'landmarks/H.e.etylus', '_d_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_fav <- makeList(IDList_fav, 'landmark', 'landmarks/H.e.favorinus', '_d_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_ven <- makeList(IDList_ven, 'landmark', 'landmarks/H.e.venus', '_d_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_era <- makeList(IDList_era, 'landmark', 'landmarks/H.e.erato', '_landmarks_LFW.txt', skipLandmark = c(2:5,7:9))
landList_emm <- makeList(IDList_emm, 'landmark', 'landmarks/H.e.emma', '_d_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_phy <- makeList(IDList_phy, 'landmark', 'landmarks/H.e.phyllis', '_calibrated_D_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_dem <- makeList(IDList_dem, 'landmark', 'landmarks/H.e.demophoon', '_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_lat <- makeList(IDList_lat, 'landmark', 'landmarks/H.e.lativitta', '_d_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_mic <- makeList(IDList_mic, 'landmark', 'landmarks/H.e.microclea', '_D_butterfly_landmarks.txt', skipLandmark = c(2:5,7:9))

# images
imageList_hydP <- makeList(IDList_hydP, 'image', 'images/H.e.hydaraP', '-D.JPG')
imageList_hydFG <- makeList(IDList_hydFG, 'image', 'images/H.e.hydaraFG', '-D.JPG')
imageList_amal <- makeList(IDList_amal, 'image', 'images/H.e.amalfreda', '-D.JPG')
imageList_cyr <- makeList(IDList_cyr, 'image', 'images/H.e.cyrbia', '_d.JPG')
imageList_not <- makeList(IDList_not, 'image', 'images/H.e.notabilis', '_d.JPG')
imageList_ety <- makeList(IDList_ety, 'image', 'images/H.e.etylus', '_d.JPG')
imageList_fav <- makeList(IDList_fav, 'image', 'images/H.e.favorinus', '_d.JPG')
imageList_ven <- makeList(IDList_ven, 'image', 'images/H.e.venus', '_d.JPG')
imageList_era <- makeList(IDList_era, 'image', 'images/H.e.erato', '-D.JPG')
imageList_emm <- makeList(IDList_emm, 'image', 'images/H.e.emma', '_d.JPG')
imageList_phy <- makeList(IDList_phy, 'image', 'images/H.e.phyllis', '_calibrated_D.JPG')
imageList_dem <- makeList(IDList_dem, 'image', 'images/H.e.demophoon', '.JPG')
imageList_lat <- makeList(IDList_lat, 'image', 'images/H.e.lativitta', '_d.JPG')
imageList_mic <- makeList(IDList_mic, 'image', 'images/H.e.microclea', '_D_butterfly.JPEG')

########################################################################

# choose target image
target <- as.matrix(read.table('cartoon/BC0004_landmarks_LFW.txt',h = F))
target <- target[c(-2:-5,-7:-9),]

# cartoon
outline_BC0004 <- read.table('cartoon/BC0004_outline.txt', h = F)
lines_BC0004 <- list.files(path ='cartoon', pattern ='BC0004_vein', full.names = T)

########################################################################
# run alignment of color patterns

# erato
RGB <- c(231, 240, 173) 
rasterList_era_sub <- patLanRGB(imageList_era, landList_era, RGB, transformRef = target, resampleFactor = 3, colOffset= 0.2,
                       adjustCoords = TRUE, res = 200, plot = 'stack')
save(rasterList_era_sub, file = 'aligned_rasterLists/rasterList_era_sub.rda')

# hydara FG
RGB <- c(223, 106, 13) 
rasterList_hydFG_sub <- patLanRGB(imageList_hydFG, landList_hydFG, RGB, transformRef = target, resampleFactor = 3, colOffset= 0.15,
                           adjustCoords = TRUE, res = 200, plot = 'stack')
save(rasterList_hydFG_sub, file = 'aligned_rasterLists/rasterList_hydFG_sub.rda')

# hydara P
RGB <- c(207, 43, 45) 
rasterList_hydP_sub <- patLanRGB(imageList_hydP, landList_hydP, RGB, transformRef = target, resampleFactor = 2, colOffset= 0.15,
                          adjustCoords = TRUE, res = 200, plot = 'stack', iterations = 3)
save(rasterList_hydP_sub, file = 'aligned_rasterLists/rasterList_hydP_sub.rda')

# demophoon
RGB <- c(207, 43, 45) 
rasterList_dem_sub <- patLanRGB(imageList_dem, landList_dem, RGB, transformRef = target, resampleFactor = 2, colOffset= 0.15,
                            adjustCoords = TRUE, res = 200, plot = 'stack', iterations = 3)
save(rasterList_dem_sub, file = 'aligned_rasterLists/rasterList_dem_sub.rda')

# notabilis
RGB <- c(240, 241, 200) 
rasterList_not_white_sub <- patLanRGB(imageList_not, landList_not, RGB, transformRef = target, resampleFactor = 3, colOffset= 0.3,
                            adjustCoords = TRUE, res = 200, plot = 'stack')

RGB <- c(253, 98, 16) 
rasterList_not_red_sub <- patLanRGB(imageList_not, landList_not, RGB, transformRef = target, resampleFactor = 3, colOffset= 0.25,
                                  adjustCoords = TRUE, res = 200, plot = 'stack', iterations = 3)

save(rasterList_not_white_sub, file = 'aligned_rasterLists/rasterList_not_white_sub.rda')
save(rasterList_not_red_sub, file = 'aligned_rasterLists/rasterList_not_red_sub.rda')

# combine red and white
rasterList_not_sub <- list()
for(e in IDList_not){
  print(e)
  subRasterList <- list(rasterList_not_white_sub[[e]], rasterList_not_red_sub[[e]])

  names(subRasterList) <- NULL
  subRasterList$fun <- sum
  subRasterList$na.rm <- TRUE
  summedRaster <- do.call(raster::mosaic,subRasterList)
  summedRaster[summedRaster == 2] <- 1
  rasterList_not_sub[[e]] <- summedRaster

}
save(rasterList_not_sub, file = 'aligned_rasterLists/rasterList_not_sub.rda')

# lativitta
RGB <- c(248, 253, 217) 
rasterList_lat_sub <- patLanRGB(imageList_lat, landList_lat, RGB, transformRef = target, resampleFactor = 3, colOffset= 0.2,
                            adjustCoords = TRUE, res = 200, plot = 'stack')
save(rasterList_lat_sub, file = 'aligned_rasterLists/rasterList_lat_sub.rda')

# amalfreda
RGB <- c(233, 236, 155)
rasterList_amal_sub <- patLanRGB(imageList_amal, landList_amal, RGB, transformRef = target, resampleFactor = 2, colOffset= 0.2,
                                 adjustCoords = TRUE, res = 200, plot = 'stack', iterations = 3)
save(rasterList_amal_sub, file = 'aligned_rasterLists/rasterList_amal_sub.rda')

# emm
RGB <- c(146, 199, 122) 
rasterList_emm1_sub <- patLanRGB(imageList_emm[c(1:8)], landList_emm[c(1:8)], RGB, transformRef = target, resampleFactor = 3, colOffset= 0.2,
                             adjustCoords = TRUE, res = 200, plot = 'stack')
RGB <- c(226, 235, 154)
rasterList_emm2_sub <- patLanRGB(imageList_emm[c(9:10)], landList_emm[c(9:10)], RGB, transformRef = target, resampleFactor = 3, colOffset= 0.2,
                            adjustCoords = TRUE, res = 200, plot = 'stack')
rasterList_emm_sub <- c(rasterList_emm1_sub, rasterList_emm2_sub)
save(rasterList_emm_sub, file = 'aligned_rasterLists/rasterList_emm_sub.rda')

# etylus
RGB <- c(184, 226, 144) 
rasterList_ety1_sub <- patLanRGB(imageList_ety[c(1:8)], landList_ety[c(1:8)], RGB, transformRef = target, resampleFactor = 3, colOffset= 0.2,
                                  adjustCoords = TRUE, res = 200, plot = 'stack')
RGB <- c(220, 225, 143) 
rasterList_ety2_sub <- patLanRGB(imageList_ety[c(9:10)], landList_ety[c(9:10)], RGB, transformRef = target, resampleFactor = 3, colOffset= 0.2,
                             adjustCoords = TRUE, res = 200, plot = 'stack')
rasterList_ety_sub <- c(rasterList_ety1_sub, rasterList_ety2_sub)
save(rasterList_ety_sub, file = 'aligned_rasterLists/rasterList_ety_sub.rda')

# favorinus
RGB <- c(149, 67, 0) 
rasterList_fav1_sub <- patLanRGB(imageList_fav[c(1:8)], landList_fav[c(1:8)], RGB, transformRef = target, resampleFactor = 3, colOffset= 0.15,
                                   adjustCoords = TRUE, res = 200, plot = 'stack')
RGB <- c(218, 108, 10) 
rasterList_fav2_sub <- patLanRGB(imageList_fav[c(9)], landList_fav[c(9)], RGB, transformRef = target, resampleFactor = 3, colOffset= 0.15,
                            adjustCoords = TRUE, res = 200, plot = 'stack')
rasterList_fav_sub <- c(rasterList_fav1_sub, rasterList_fav2_sub)
save(rasterList_fav_sub, file = 'aligned_rasterLists/rasterList_fav_sub.rda')

# phyllis
RGB <- c(186, 57, 25) 
rasterList_phy1_sub <- patLanRGB(imageList_phy[c(1:8)], landList_phy[c(1:8)], RGB, transformRef = target, resampleFactor = 1, colOffset= 0.15,
                                   adjustCoords = TRUE, res = 200, plot = 'stack')
RGB <- c(216, 121, 15) 
rasterList_phy2_sub <- patLanRGB(imageList_phy[c(9:10)], landList_phy[c(9:10)], RGB, transformRef = target, resampleFactor = 1, colOffset= 0.15,
                             adjustCoords = TRUE, res = 200, plot = 'stack')
rasterList_phy_sub <- c(rasterList_phy1_sub, rasterList_phy2_sub)
save(rasterList_phy_sub, file = 'aligned_rasterLists/rasterList_phy_sub.rda')

# cyrbia
RGB <- c(229, 59, 71) 
rasterList_cyr_sub <- patLanRGB(imageList_cyr, landList_cyr, RGB, transformRef = target, resampleFactor = 3, colOffset= 0.15,
                                adjustCoords = TRUE, res = 200, plot = 'stack', iterations = 3)
save(rasterList_cyr_sub, file = 'aligned_rasterLists/rasterList_cyr_sub.rda')

# venus
RGB <- c(229, 59, 71) 
rasterList_ven_sub <- patLanRGB(imageList_ven, landList_ven, RGB, transformRef = target, resampleFactor = 3, colOffset= 0.15,
                                adjustCoords = TRUE, res = 200, plot = 'stack', iterations = 3)
save(rasterList_ven_sub, file = 'aligned_rasterLists/rasterList_ven_sub.rda')

# microclea
RGB <- c(216, 122, 38)
rasterList_mic_sub <- patLanRGB(imageList_mic, landList_mic, RGB, transformRef = target, resampleFactor = 2, colOffset= 0.15,
                             adjustCoords = TRUE, res = 200, plot = 'stack', iterations = 3)
save(rasterList_mic_sub, file = 'aligned_rasterLists/rasterList_mic_sub.rda')


#####
load('aligned_rasterLists/rasterList_era_sub.rda')
load('aligned_rasterLists/rasterList_hydFG_sub.rda')
load('aligned_rasterLists/rasterList_not_sub.rda')
load('aligned_rasterLists/rasterList_lat_sub.rda')
load('aligned_rasterLists/rasterList_emm_sub.rda')
load('aligned_rasterLists/rasterList_ety_sub.rda')
load('aligned_rasterLists/rasterList_fav_sub.rda')
load('aligned_rasterLists/rasterList_hydP_sub.rda')
load('aligned_rasterLists/rasterList_dem_sub.rda')
load('aligned_rasterLists/rasterList_phy_sub.rda')
load('aligned_rasterLists/rasterList_cyr_sub.rda')
load('aligned_rasterLists/rasterList_ven_sub.rda')
load('aligned_rasterLists/rasterList_amal_sub.rda')
load('aligned_rasterLists/rasterList_mic_sub.rda')
#####

# remove noise by masking
# mask
mask_BC0004 <- read.table('mask/BC0004_mask1.txt', h = F)

rasterList_era_M_sub <-list()
for(e in IDList_era){
  rasterList_era_M_sub[[e]] <- maskOutline(rasterList_era_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_era)
}

rasterList_hydFG_M_sub <-list()
for(e in IDList_hydFG){
  rasterList_hydFG_M_sub[[e]] <- maskOutline(rasterList_hydFG_sub[[e]], mask_BC0004, refShape = 'target', 
                                         imageList = imageList_hydFG)
}

rasterList_not_M_sub <-list()
for(e in IDList_not){
  rasterList_not_M_sub[[e]] <- maskOutline(rasterList_not_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_not)
}

rasterList_lat_M_sub <-list()
for(e in IDList_lat){
  rasterList_lat_M_sub[[e]] <- maskOutline(rasterList_lat_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_lat)
}

rasterList_emm_M_sub <-list()
for(e in IDList_emm){
  rasterList_emm_M_sub[[e]] <- maskOutline(rasterList_emm_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_emm)
}

rasterList_ety_M_sub <-list()
for(e in IDList_ety){
  rasterList_ety_M_sub[[e]] <- maskOutline(rasterList_ety_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_ety)
}

rasterList_fav_M_sub <-list()
for(e in IDList_fav){
  rasterList_fav_M_sub[[e]] <- maskOutline(rasterList_fav_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_fav)
}

rasterList_hydP_M_sub <-list()
for(e in IDList_hydP){
  rasterList_hydP_M_sub[[e]] <- maskOutline(rasterList_hydP_sub[[e]], mask_BC0004, refShape = 'target', 
                                        imageList = imageList_hydP)
}

rasterList_dem_M_sub <-list()
for(e in IDList_dem){
  rasterList_dem_M_sub[[e]] <- maskOutline(rasterList_dem_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_dem)
}

rasterList_phy_M_sub <-list()
for(e in IDList_phy){
  rasterList_phy_M_sub[[e]] <- maskOutline(rasterList_phy_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_phy)
}

rasterList_cyr_M_sub <-list()
for(e in IDList_cyr){
  rasterList_cyr_M_sub[[e]] <- maskOutline(rasterList_cyr_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_cyr)
}

rasterList_ven_M_sub <-list()
for(e in IDList_ven){
  rasterList_ven_M_sub[[e]] <- maskOutline(rasterList_ven_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_ven)
}

rasterList_mic_M_sub <-list()
for(e in IDList_mic){
  rasterList_mic_M_sub[[e]] <- maskOutline(rasterList_mic_sub[[e]], mask_BC0004, refShape = 'target', 
                                       imageList = imageList_mic)
}

# setMask(summedRaster_amal, IDListAmal, cartoonID = 'BC0004', 'mask/mask_amal.txt')

mask_amal <- read.table("mask/mask_amal.txt", h=T)
rasterList_amal_M_sub <-list()
for(e in IDList_amal){
  rasterList_amal_M_sub[[e]] <- maskOutline(rasterList_amal_sub[[e]], mask_amal, refShape = 'target', 
                                       imageList = imageList_amal)
}

# setMask(summedRaster_lat_M_sub, IDList_lat, cartoonID = 'BC0004', 'mask/mask_lat.txt', imageList = imageList_hydFG)

mask_lat <- read.table("mask/mask_lat.txt", h=T)
rasterList_lat_M_sub <-list()
for(e in IDList_lat){
  rasterList_lat_M_sub[[e]] <- maskOutline(rasterList_lat_sub[[e]], mask_lat, refShape = 'target', 
                                            imageList = imageList_lat)
}

# setMask(summedRaster_emm_M_sub, IDList_emm, cartoonID = 'BC0004', 'mask/mask_emm.txt', imageList = imageList_hydFG)

mask_emm <- read.table("mask/mask_emm.txt", h=T)
rasterList_emm_M_sub <-list()
for(e in IDList_emm){
  rasterList_emm_M_sub[[e]] <- maskOutline(rasterList_emm_sub[[e]], mask_emm, refShape = 'target', 
                                           imageList = imageList_emm)
}

# setMask(summedRaster_ety_M_sub, IDList_ety, cartoonID = 'BC0004', 'mask/mask_ety.txt', imageList = imageList_hydFG)

mask_ety <- read.table("mask/mask_ety.txt", h=T)
rasterList_ety_M_sub <-list()
for(e in IDList_ety){
  rasterList_ety_M_sub[[e]] <- maskOutline(rasterList_ety_sub[[e]], mask_ety, refShape = 'target', 
                                           imageList = imageList_ety)
}

# setMask(summedRaster_not_M_sub, IDList_not, cartoonID = 'BC0004', 'mask/mask_not.txt', imageList = imageList_hydFG)

mask_not <- read.table("mask/mask_not.txt", h=T)
rasterList_not_M_sub <-list()
for(e in IDList_not){
  rasterList_not_M_sub[[e]] <- maskOutline(rasterList_not_sub[[e]], mask_not, refShape = 'target', 
                                           imageList = imageList_not)
}

###
# Plot masked (final)
###

# save(rasterList_era_M_sub, file = 'aligned_rasterLists/rasterList_era_M_sub.rda')
# save(rasterList_not_M_sub, file = 'aligned_rasterLists/rasterList_not_M_sub.rda')
# save(rasterList_lat_M_sub, file = 'aligned_rasterLists/rasterList_lat_M_sub.rda')
# save(rasterList_emm_M_sub, file = 'aligned_rasterLists/rasterList_emm_M_sub.rda')
# save(rasterList_ety_M_sub, file = 'aligned_rasterLists/rasterList_ety_M_sub.rda')
# save(rasterList_phy_M_sub, file = 'aligned_rasterLists/rasterList_phy_M_sub.rda')
# save(rasterList_hydFG_M_sub, file = 'aligned_rasterLists/rasterList_hydFG_M_sub.rda')
# save(rasterList_hydP_M_sub, file = 'aligned_rasterLists/rasterList_hydP_M_sub.rda')
# save(rasterList_fav_M_sub, file = 'aligned_rasterLists/rasterList_fav_M_sub.rda')
# save(rasterList_dem_M_sub, file = 'aligned_rasterLists/rasterList_dem_M_sub.rda')
# save(rasterList_ven_M_sub, file = 'aligned_rasterLists/rasterList_ven_M_sub.rda')
# save(rasterList_cyr_M_sub, file = 'aligned_rasterLists/rasterList_cyr_M_sub.rda')
# save(rasterList_amal_M_sub, file = 'aligned_rasterLists/rasterList_amal_M_sub.rda')
# save(rasterList_mic_M_sub, file = 'aligned_rasterLists/rasterList_mic_M_sub.rda')

load('aligned_rasterLists/rasterList_era_M_sub.rda')
load('aligned_rasterLists/rasterList_hydFG_M_sub.rda')
load('aligned_rasterLists/rasterList_not_M_sub.rda')
load('aligned_rasterLists/rasterList_lat_M_sub.rda')
load('aligned_rasterLists/rasterList_emm_M_sub.rda')
load('aligned_rasterLists/rasterList_ety_M_sub.rda')
load('aligned_rasterLists/rasterList_fav_M_sub.rda')
load('aligned_rasterLists/rasterList_hydP_M_sub.rda')
load('aligned_rasterLists/rasterList_dem_M_sub.rda')
load('aligned_rasterLists/rasterList_phy_M_sub.rda')
load('aligned_rasterLists/rasterList_cyr_M_sub.rda')
load('aligned_rasterLists/rasterList_ven_M_sub.rda')
load('aligned_rasterLists/rasterList_amal_M_sub.rda')
load('aligned_rasterLists/rasterList_mic_M_sub.rda')

summedRaster_era_M_sub <- sumRaster(rasterList_era_M_sub, IDList_era, type = 'RGB')
summedRaster_not_M_sub <- sumRaster(rasterList_not_M_sub, IDList_not, type = 'RGB')
summedRaster_lat_M_sub <- sumRaster(rasterList_lat_M_sub, IDList_lat, type = 'RGB')
summedRaster_emm_M_sub <- sumRaster(rasterList_emm_M_sub, IDList_emm, type = 'RGB')
summedRaster_ety_M_sub <- sumRaster(rasterList_ety_M_sub, IDList_ety, type = 'RGB')
summedRaster_phy_M_sub <- sumRaster(rasterList_phy_M_sub, IDList_phy, type = 'RGB')
summedRaster_hydFG_M_sub <- sumRaster(rasterList_hydFG_M_sub, IDList_hydFG, type = 'RGB')
summedRaster_hydP_M_sub  <- sumRaster(rasterList_hydP_M_sub, IDList_hydP, type = 'RGB')
summedRaster_fav_M_sub   <- sumRaster(rasterList_fav_M_sub, IDList_fav, type = 'RGB')
summedRaster_dem_M_sub   <- sumRaster(rasterList_dem_M_sub, IDList_dem, type = 'RGB')
summedRaster_ven_M_sub   <- sumRaster(rasterList_ven_M_sub, IDList_ven, type = 'RGB')
summedRaster_cyr_M_sub   <- sumRaster(rasterList_cyr_M_sub, IDList_cyr, type = 'RGB')
summedRaster_amal_M_sub  <- sumRaster(rasterList_amal_M_sub, IDList_amal, type = 'RGB')
summedRaster_mic_M_sub   <- sumRaster(rasterList_mic_M_sub, IDList_mic, type = 'RGB')


colfunc <- inferno(100)

plotHeat(summedRaster_era_M_sub, IDList_era, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_era, adjustCoords = TRUE, 
         imageList = imageList_era, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_hydFG_M_sub, IDList_hydFG, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_hydFG, adjustCoords = TRUE, 
         imageList = imageList_hydFG, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_hydP_M_sub, IDList_hydP, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_hydP, adjustCoords = TRUE, 
         imageList = imageList_hydP, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_not_M_sub, IDList_not, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004,
         lines = lines_BC0004, landList = landmarkList_not, adjustCoords = TRUE,
         imageList = imageList_not, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_lat_M_sub, IDList_lat, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_lat, adjustCoords = TRUE, 
         imageList = imageList_lat, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_emm_M_sub, IDList_emm, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_emm, adjustCoords = TRUE, 
         imageList = imageList_emm, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_fav_M_sub, IDList_fav, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_fav, adjustCoords = TRUE, 
         imageList = imageList_fav, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_dem_M_sub, IDList_dem, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_dem, adjustCoords = TRUE, 
         imageList = imageList_dem, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_ven_M_sub, IDList_ven, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_ven, adjustCoords = TRUE, 
         imageList = imageList_ven, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_cyr_M_sub, IDList_cyr, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_cyr, adjustCoords = TRUE, 
         imageList = imageList_cyr, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_ety_M_sub, IDList_ety, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_ety, adjustCoords = TRUE, 
         imageList = imageList_ety, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_phy_M_sub, IDList_phy, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_phy, adjustCoords = TRUE, 
         imageList = imageList_phy, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_amal_M_sub, IDList_amal, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_amal, adjustCoords = TRUE, 
         imageList = imageList_amal, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

plotHeat(summedRaster_mic_M_sub, IDList_mic, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_mic, adjustCoords = TRUE, 
         imageList = imageList_mic, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',colpalette = colfunc)

