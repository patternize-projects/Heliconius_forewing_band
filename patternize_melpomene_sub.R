library(patternize)
library(viridis)

# H. melpomene
IDList_the  <- c('10428329','10428331','10428332','10428333','10428334','10428335','10428336','10428337','10428338','10428339','10428340','10428341','10428342','10428343')
IDList_mal  <- c('CAM016144','CAM016267','CAM016609','CAM016610','CAM016611','CAM017049','CAM017064','CAM017070','CAM016224','CAM016606','CAM016772')#,'10428198','10428199') #'CAM017030','CAM017120'
IDList_ecu  <- c('CAM009120','CAM009141','CAM009112','CAM009113','CAM009114','CAM009115','CAM009119','10428219','10428220') # these are faded, but we can extract them with different color values
IDList_amar  <- c('MJ12.3137','MJ12.3371','MJ12.3392','MJ12.3393','MJ12.3396','MJ12.3414','MJ12.3416','MJ12.3417','MJ12.3442','10428116','BC2639') # 'MJ12.3306'
IDList_cyt  <- c('15N005','15N006','15N009','15N020','15N022','15N023','CAM008510_d','CAM040383_d','CAM040459_d','CAM040474_d')
IDList_melFG<- c('CAM001422','CAM008171','CAM008215','CAM008218','CAM009316','CAM009317','CAM000413','10428369','10428371') #'CAM009315'
IDList_melP <- c('CAM008863','CAM008887','CAM008888','CAM008954','CAM008955','CAM008956','CAM008957','CAM008979')
IDList_mer  <- c('13715_H_m_meriana','melp_14-103','melp_14-108','melp_14-110','melp_14-111','melp_14-122','melp_14-133','melp_14-138')# ,'melp_14-138' bad landmarks
IDList_agl  <- c('CAM008689', 'CAM008702','10428229','10428230','10428231','10428232','10428233','10428234','10428235','10428238','10428240','10428266','10428267')
IDList_nan  <- c('LMCI_105-62','LMCI_105-63','LMCI_105-64','LMCI_183-10','LMCI_183-11','LMCI_183-14','LMCI_183-15','LMCI_183-16','LMCI_183-18','LMCI_183-19')
IDList_ple  <- c('CAM016347','CAM016349','CAM016354','CAM016355','CAM016378','CAM016810','CAM017185','CAM017186','CAM017187','CAM017614')
IDList_ros  <- c('CAM000903','CAM000947','CAM001015','CAM001027','CAM001067', 'CAM001137','CAM001391','CAM002901','CAM008052','CAM009554')
IDList_vul  <- c('CAM000058','CAM000059','CAM000060','CAM000061','CAM000062','CAM000063','CAM000064','CAM000129','CAM000132','CAM000134')
IDList_xen  <- c('MJ12.3606','MJ12.3608','MJ12.3636','MJ12.3638','MJ12.3647','MJ12.3648','MJ12.3651','MJ12.3653')

# All 18 landmarks
landList_melP <- makeList(IDList_melP, 'landmark', 'landmarks/H.m.melpomeneP', '_d.txt', skipLandmark = c(2:5,7:9))
landList_melFG <- makeList(IDList_melFG, 'landmark', 'landmarks/H.m.melpomeneFG', '_d.txt', skipLandmark = c(2:5,7:9))
landList_mer <- makeList(IDList_mer, 'landmark', 'landmarks/H.m.meriana', '.txt', skipLandmark = c(2:5,7:9))
landList_cyt <- makeList(IDList_cyt, 'landmark', 'landmarks/H.m.cythera', '.txt', skipLandmark = c(2:5,7:9))
landList_ple <- makeList(IDList_ple, 'landmark', 'landmarks/H.m.plesseni', '_d.txt', skipLandmark = c(2:5,7:9))
landList_ecu <- makeList(IDList_ecu, 'landmark', 'landmarks/H.m.ecuadorensis', '_d.txt', skipLandmark = c(2:5,7:9))
landList_amar <- makeList(IDList_amar, 'landmark', 'landmarks/H.m.amaryllis', '-d.txt', skipLandmark = c(2:5,7:9))
landList_vul <- makeList(IDList_vul, 'landmark', 'landmarks/H.m.vulcanus', '_d.txt', skipLandmark = c(2:5,7:9))
landList_the <- makeList(IDList_the, 'landmark', 'landmarks/H.m.thelxiopeia', '_D_butterfly.txt', skipLandmark = c(2:5,7:9))
landList_mal <- makeList(IDList_mal, 'landmark', 'landmarks/H.m.malleti', '_d.txt', skipLandmark = c(2:5,7:9))
landList_nan <- makeList(IDList_nan, 'landmark', 'landmarks/H.m.nanna', '_Hmnanna_upper.txt', skipLandmark = c(2:5,7:9))
landList_ros <- makeList(IDList_ros, 'landmark', 'landmarks/H.m.rosina', '_d.txt', skipLandmark = c(2:5,7:9))
landList_agl <- makeList(IDList_agl, 'landmark', 'landmarks/H.m.aglaope', '_D_butterfly_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_xen <- makeList(IDList_xen, 'landmark', 'landmarks/H.m.xenoclea', '-d.txt', skipLandmark = c(2:5,7:9))

# images
imageList_melP <- makeList(IDList_melP, 'image', 'images/H.m.melpomeneP', '_d.JPG')
imageList_melFG <- makeList(IDList_melFG, 'image', 'images/H.m.melpomeneFG', '_d.JPG')
imageList_mer <- makeList(IDList_mer, 'image', 'images/H.m.meriana', '.JPG')
imageList_cyt <- makeList(IDList_cyt, 'image', 'images/H.m.cythera', '.JPG')
imageList_ple <- makeList(IDList_ple, 'image', 'images/H.m.plesseni', '_d.JPG')
imageList_ecu <- makeList(IDList_ecu, 'image', 'images/H.m.ecuadorensis', '_d.JPG')
imageList_amar <- makeList(IDList_amar, 'image', 'images/H.m.amaryllis', '-d.JPEG')
imageList_vul <- makeList(IDList_vul, 'image', 'images/H.m.vulcanus', '_d.JPG')
imageList_the <- makeList(IDList_the, 'image', 'images/H.m.thelxiopeia', '_D_butterfly.JPEG')
imageList_mal <- makeList(IDList_mal, 'image', 'images/H.m.malleti', '_d.JPG')
imageList_nan <- makeList(IDList_nan, 'image', 'images/H.m.nanna', '_Hmnanna_upper.JPG')
imageList_ros <- makeList(IDList_ros, 'image', 'images/H.m.rosina', '_d.JPG')
imageList_agl <- makeList(IDList_agl, 'image', 'images/H.m.aglaope', '_d.JPG')
imageList_xen <- makeList(IDList_xen, 'image', 'images/H.m.xenoclea', '-d.JPEG')

########################################################################

# choose target image
target <- as.matrix(read.table('cartoon/BC0004_landmarks_LFW.txt',h = F))
target <- target[c(-2:-5,-7:-9),]

# cartoon
outline_BC0004 <- read.table('cartoon/BC0004_outline.txt', h = F)
lines_BC0004 <- list.files(path ='cartoon', pattern ='BC0004_vein', full.names = T)

########################################################################
# run alignment of color patterns

# malleti
RGB <- c(225,229,184) 
rasterList_mal1_sub <- patLanRGB(imageList_mal[c(1:8)], landList_mal[c(1:8)], RGB, transformRef = target, resampleFactor = 3,
                                     colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
RGB <- c(196,182,107) 
rasterList_mal2_sub <- patLanRGB(imageList_mal[c(9:11)], landList_mal[c(9:11)], RGB, transformRef = target, resampleFactor = 3,
                                 colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
rasterList_mal_sub <- c(rasterList_mal1_sub, rasterList_mal2_sub)
save(rasterList_mal_sub, file = 'aligned_rasterLists/rasterList_mal_sub.rda')

# ecuadorensis
RGB <- c(243,249,215) 
rasterList_ecu1_sub <- patLanRGB(imageList_ecu[c(1:7)], landList_ecu[c(1:7)], RGB, transformRef = target, resampleFactor = 3,
                                    colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
RGB <- c(236,230,182) 
rasterList_ecu2_sub <- patLanRGB(imageList_ecu[c(8:9)], landList_ecu[c(8:9)], RGB, transformRef = target, resampleFactor = 3,
                             colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')

rasterList_ecu_sub <- c(rasterList_ecu1_sub, rasterList_ecu2_sub)
save(rasterList_ecu_sub, file = 'aligned_rasterLists/rasterList_ecu_sub.rda')

# amaryllis
RGB <- c(231,71,33) 
rasterList_amar1_sub <- patLanRGB(imageList_amar[c(1:9)], landList_amar[c(1:9)], RGB, transformRef = target, resampleFactor = 3,
                                         colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
RGB <- c(246,146,6) 
rasterList_amar2_sub <- patLanRGB(imageList_amar[c(10)], landList_amar[c(10)], RGB, transformRef = target, resampleFactor = 3,
                              colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
RGB <- c(149,49,6) 
rasterList_amar3_sub<- patLanRGB(imageList_amar[c(11)], landList_amar[c(11)], RGB, transformRef = target, resampleFactor = 3,
                              colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
rasterList_amar_sub <- c(rasterList_amar1_sub, rasterList_amar2_sub,rasterList_amar3_sub)
save(rasterList_amar_sub, file = 'aligned_rasterLists/rasterList_amar_sub.rda')

# cythera
RGB <- c(252,126,120) 
rasterList_cyt_redA_sub <- patLanRGB(imageList_cyt[c(1:6)], landList_cyt[c(1:6)], RGB, transformRef = target, resampleFactor = 2,
                                       colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack', iterations = 3)
RGB <- c(242,112,36) 
rasterList_cyt_redB_sub <- patLanRGB(imageList_cyt[c(7:10)], landList_cyt[c(7:10)], RGB, transformRef = target, resampleFactor = 2,
                                       colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack', iterations = 3)
rasterList_cyt_red_sub <- c(rasterList_cyt_redA_sub, rasterList_cyt_redB_sub)

RGB <- c(241,250,254) 
rasterList_cyt_white_sub <- patLanRGB(imageList_cyt, landList_cyt, RGB, transformRef = target, resampleFactor = 2,
                                       colOffset = 0.3, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')

# combine red and white
rasterList_cyt_sub <- list()
for(e in IDList_cyt){
  print(e)
  subRasterList <- list(rasterList_cyt_red_sub[[e]], rasterList_cyt_white_sub[[e]])

  names(subRasterList) <- NULL
  subRasterList$fun <- sum
  subRasterList$na.rm <- TRUE
  summedRaster <- do.call(raster::mosaic,subRasterList)
  summedRaster[summedRaster == 2] <- 1
  rasterList_cyt_sub[[e]] <- summedRaster

}
save(rasterList_cyt_sub, file = 'aligned_rasterLists/rasterList_cyt_sub.rda')

# melpomene FG
RGB <- c(222,93,19) 
rasterList_melFG_sub <- patLanRGB(imageList_melFG, landList_melFG, RGB, transformRef = target, resampleFactor = 3,
                                     colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack', iterations = 3)
save(rasterList_melFG_sub, file = 'aligned_rasterLists/rasterList_melFG_sub.rda')

# melpomene P
RGB <- c(231,71,33) 
rasterList_melP_sub <- patLanRGB(imageList_melP, landList_melP, RGB, transformRef = target, resampleFactor = 3,
                                    colOffset = 0.18, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack', iterations = 3)
save(rasterList_melP_sub, file = 'aligned_rasterLists/rasterList_melP_sub.rda')

# meriana
RGB <- c(167,138,94)
rasterList_merA_sub <- patLanRGB(imageList_mer[c(1)], landList_mer[c(1)], RGB, transformRef = target, resampleFactor = 2,
                                        colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
RGB <- c(247,222,180)
rasterList_merB_sub <- patLanRGB(imageList_mer[c(2:8)], landList_mer[c(2:8)], RGB, transformRef = target, resampleFactor = 1,
                                        colOffset = 0.3, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
rasterList_mer_sub <- c(rasterList_merA_sub, rasterList_merB_sub)

save(rasterList_mer_sub, file = 'aligned_rasterLists/rasterList_mer_sub.rda')

# aglaope
RGB <- c(235,250,213) 
rasterList_agl_sub <- patLanRGB(imageList_agl, landList_agl, RGB, transformRef = target, resampleFactor = 3,
                                       colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
save(rasterList_agl_sub, file = 'aligned_rasterLists/rasterList_agl_sub.rda')

# nanna
RGB <- c(178,25,14) 
rasterList_nan_sub <- patLanRGB(imageList_nan, landList_nan, RGB, transformRef = target, resampleFactor = 3,
                                     colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack', iterations = 3)
save(rasterList_nan_sub, file = 'aligned_rasterLists/rasterList_nan_sub.rda')

# plesseni
RGB <- c(231,234,237) # white
rasterList_ple_white_sub <- patLanRGB(imageList_ple, landList_ple, RGB, transformRef = target, resampleFactor = 3,
                                        colOffset = 0.3, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
RGB <- c(248,106,38) # red
rasterList_ple_red_sub <- patLanRGB(imageList_ple, landList_ple, RGB, transformRef = target, resampleFactor = 3,
                                              colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack', iteration = 3)

# combine red and white
rasterList_ple_sub <- list()
for(e in IDList_ple){
  print(e)
  subRasterList <- list(rasterList_ple_white_sub[[e]], rasterList_ple_red_sub[[e]])

  names(subRasterList) <- NULL
  subRasterList$fun <- sum
  subRasterList$na.rm <- TRUE
  summedRaster <- do.call(raster::mosaic,subRasterList)
  summedRaster[summedRaster == 2] <- 1
  rasterList_ple_sub[[e]] <- summedRaster

}
save(rasterList_ple_sub, file = 'aligned_rasterLists/rasterList_ple_sub.rda')

# rosina
RGB <- c(224,102,37) 
rasterList_ros_sub <- patLanRGB(imageList_ros, landList_ros, RGB, transformRef = target, resampleFactor = 3,
                                      colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack',
                                      iteration = 3)
save(rasterList_ros_sub, file = 'aligned_rasterLists/rasterList_ros_sub.rda')

# vulcanus
RGB <- c(240,110,38) 
rasterList_vul_sub <- patLanRGB(imageList_vul, landList_vul, RGB, transformRef = target, resampleFactor = 3,
                                        colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack',
                                        iteration = 3)
save(rasterList_vul_sub, file = 'output/rasterList_vul_sub.rda')

# thelxiopea
RGB <- c(248,249,191) 
rasterList_the_sub <- patLanRGB(imageList_the, landList_the, RGB, transformRef = target, resampleFactor = 3,
                                         colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')
save(rasterList_the_sub, file = 'aligned_rasterLists/rasterList_the_sub.rda')

# xenoclea
RGB <- c(191,67,31) 
rasterList_xen_sub <- patLanRGB(imageList_xen, landList_xen, RGB, transformRef = target, resampleFactor = 3,
                                         colOffset = 0.2, crop = TRUE, res = 200, adjustCoords = TRUE, plot = 'stack', iterations = TRUE)
save(rasterList_xen_sub, file = 'aligned_rasterLists/rasterList_xen_sub.rda')

# load 
load('aligned_rasterLists/rasterList_mal_sub.rda')
load('aligned_rasterLists/rasterList_ecu_sub.rda')
load('aligned_rasterLists/rasterList_amar_sub.rda')
load('aligned_rasterLists/rasterList_cyt_sub.rda')
load('aligned_rasterLists/rasterList_melFG_sub.rda')
load('aligned_rasterLists/rasterList_melP_sub.rda')
load('aligned_rasterLists/rasterList_mer_sub.rda')
load('aligned_rasterLists/rasterList_agl_sub.rda')
load('aligned_rasterLists/rasterList_nan_sub.rda')
load('aligned_rasterLists/rasterList_ple_sub.rda')
load('aligned_rasterLists/rasterList_ros_sub.rda')
load('aligned_rasterLists/rasterList_vul_sub.rda')
load('aligned_rasterLists/rasterList_the_sub.rda')
load('aligned_rasterLists/rasterList_xen_sub.rda')


# setMask(summedRaster_mal_M, IDList_mal ,refImage = target, 'mask/mask_mal.txt')
# setMask(summedRaster_ecu, IDlist_ecu, 'mask/mask_ecu.txt')
# setMask(summedRaster_cythera, IDlist_cyt, 'mask/mask_cyt.txt')
# setMask(summedRaster_melFG, IDlist_melFG, 'mask/mask_melFG.txt')
# setMask(summedRaster_meriana, IDlist_meri, 'mask/mask_mer.txt')
# setMask(summedRaster_aglaope, IDlist_agla, 'mask/mask_agl.txt')
# setMask(summedRaster_nanna, IDlist_na, 'mask/mask_nan.txt')
# setMask(summedRaster_plesseni, IDlist_ple, 'mask/mask_ple.txt')
# setMask(summedRaster_thelxixopeia, IDlist_thel,refImage = target, 'mask/mask_thel.txt')
# setMask(summedRaster_xenoclea, IDlist_xeno, refImage = target, 'mask/mask_xeno.txt')


# mask
mask_mal <- read.table("mask/mask_mal.txt", h=T)
rasterList_mal_M_sub <-list()
for(e in 1:length(rasterList_mal_sub)){
  ID <- names(rasterList_mal_sub)[[e]]
  rasterList_mal_M_sub[[ID]] <- maskOutline(rasterList_mal_sub[[ID]], IDlist = IDList_mal, mask_mal, 
                                               refShape = 'target', imageList = imageList_mal)
}


mask_ecu <- read.table("mask/mask_ecu.txt", h=T)
rasterList_ecu_M_sub <-list()
for(e in 1:length(rasterList_ecu_sub)){
  ID <- names(rasterList_ecu_sub)[[e]]
  rasterList_ecu_M_sub[[ID]] <- maskOutline(rasterList_ecu_sub[[ID]], IDlist = IDList_ecu, mask_ecu, 
                                               refShape = 'target', imageList = imageList_ecu)
}

mask_cyt <- read.table("mask/mask_cyt.txt", h=T)
rasterList_cyt_M_sub <-list()
for(e in 1:length(rasterList_cyt_sub)){
  ID <- names(rasterList_cyt_sub)[[e]]
  rasterList_cyt_M_sub[[ID]] <- maskOutline(rasterList_cyt_sub[[ID]], IDlist = IDList_cyt, mask_cyt, 
                                                   refShape = 'target', imageList = imageList_cyt)
}

mask_melFG <- read.table("mask/mask_melFG.txt", h=T)
rasterList_melFG_M_sub <-list()
for(e in 1:length(rasterList_melFG_sub)){
  ID <- names(rasterList_melFG_sub)[[e]]
  rasterList_melFG_M_sub[[ID]] <- maskOutline(rasterList_melFG_sub[[ID]], IDlist = IDList_melFG, mask_melFG, 
                                                 refShape = 'target', imageList = imageList_melFG)
}

mask_mer <- read.table("mask/mask_mer.txt", h=T)
rasterList_mer_M_sub<-list()
for(e in 1:length(rasterList_mer_sub)){
  ID <- names(rasterList_mer_sub)[[e]]
  rasterList_mer_M_sub[[ID]] <- maskOutline(rasterList_mer_sub[[ID]], IDlist = IDList_mer, mask_mer, 
                                                   refShape = 'target', imageList = imageList_mer)
}

mask_agl <- read.table("mask/mask_agl.txt", h=T)
rasterList_agl_M_sub<-list()
for(e in 1:length(rasterList_agl_sub)){
  ID <- names(rasterList_agl_sub)[[e]]
  rasterList_agl_M_sub[[ID]] <- maskOutline(rasterList_agl_sub[[ID]], IDlist = IDList_agl, mask_agl, 
                                                   refShape = 'target', imageList = imageList_agl)
}

mask_nan <- read.table("mask/mask_nan.txt", h=T)
rasterList_nan_M_sub<-list()
for(e in 1:length(rasterList_nan_sub)){
  ID <- names(rasterList_nan_sub)[[e]]
  rasterList_nan_M_sub[[ID]] <- maskOutline(rasterList_nan_sub[[ID]], IDlist = IDList_nan, mask_nan, 
                                                 refShape = 'target', imageList = imageList_nan)
}

mask_ple <- read.table("mask/mask_ple.txt", h=T)
rasterList_ple_M_sub<-list()
for(e in 1:length(rasterList_ple_sub)){
  ID <- names(rasterList_ple_sub)[[e]]
  rasterList_ple_M_sub[[ID]] <- maskOutline(rasterList_ple_sub[[ID]], IDlist = IDList_ple, mask_ple, 
                                                    refShape = 'target', imageList = imageList_ple)
}

mask_the <- read.table("mask/mask_thel.txt", h=T)
rasterList_the_M_sub<-list()
for(e in 1:length(rasterList_the_sub)){
  ID <- names(rasterList_the_sub)[[e]]
  rasterList_the_M_sub[[ID]] <- maskOutline(rasterList_the_sub[[ID]], IDlist = IDList_the, mask_the, 
                                                       refShape = 'target', imageList = imageList_the)
}

mask_xen <- read.table("mask/mask_xeno.txt", h=T)
rasterList_xen_M_sub<-list()
for(e in 1:length(rasterList_xen_sub)){
  ID <- names(rasterList_xen_sub)[[e]]
  rasterList_xen_M_sub[[ID]] <- maskOutline(rasterList_xen_sub[[ID]], IDlist = IDList_xen, mask_xen, 
                                                    refShape = 'target', imageList = imageList_xen)
}

# rasterList_amar_M_sub <- rasterList_amar_sub
# rasterList_melP_M_sub <- rasterList_melP_sub
# rasterList_ros_M_sub <- rasterList_ros_sub
# rasterList_vul_M_sub <- rasterList_vul_sub
# 
# save(rasterList_mal_M_sub, file = 'aligned_rasterLists/rasterList_mal_M_sub.rda')
# save(rasterList_ecu_M_sub, file = 'aligned_rasterLists/rasterList_ecu_M_sub.rda')
# save(rasterList_amar_M_sub, file = 'aligned_rasterLists/rasterList_amar_M_sub.rda')
# save(rasterList_cyt_M_sub, file = 'aligned_rasterLists/rasterList_cyt_M_sub.rda')
# save(rasterList_melFG_M_sub, file = 'aligned_rasterLists/rasterList_melFG_M_sub.rda')
# save(rasterList_melP_M_sub, file = 'aligned_rasterLists/rasterList_melP_M_sub.rda')
# save(rasterList_mer_M_sub, file = 'aligned_rasterLists/rasterList_mer_M_sub.rda')
# save(rasterList_agl_M_sub, file = 'aligned_rasterLists/rasterList_agl_M_sub.rda')
# save(rasterList_nan_M_sub, file = 'aligned_rasterLists/rasterList_nan_M_sub.rda')
# save(rasterList_ple_M_sub, file = 'aligned_rasterLists/rasterList_ple_M_sub.rda')
# save(rasterList_ros_M_sub, file = 'aligned_rasterLists/rasterList_ros_M_sub.rda')
# save(rasterList_vul_M_sub, file = 'aligned_rasterLists/rasterList_vul_M_sub.rda')
# save(rasterList_the_M_sub, file = 'aligned_rasterLists/rasterList_the_M_sub.rda')
# save(rasterList_xen_M_sub, file = 'aligned_rasterLists/rasterList_xen_M_sub.rda')

load('aligned_rasterLists/rasterList_mal_M_sub.rda')
load('aligned_rasterLists/rasterList_ecu_M_sub.rda')
load('aligned_rasterLists/rasterList_amar_M_sub.rda')
load('aligned_rasterLists/rasterList_cyt_M_sub.rda')
load('aligned_rasterLists/rasterList_melFG_M_sub.rda')
load('aligned_rasterLists/rasterList_melP_M_sub.rda')
load('aligned_rasterLists/rasterList_mer_M_sub.rda')
load('aligned_rasterLists/rasterList_agl_M_sub.rda')
load('aligned_rasterLists/rasterList_nan_M_sub.rda')
load('aligned_rasterLists/rasterList_ple_M_sub.rda')
load('aligned_rasterLists/rasterList_ros_M_sub.rda')
load('aligned_rasterLists/rasterList_vul_M_sub.rda')
load('aligned_rasterLists/rasterList_the_M_sub.rda')
load('aligned_rasterLists/rasterList_xen_M_sub.rda')

# sum the colorpatterns
summedRaster_mal_M_sub <- sumRaster(rasterList_mal_M_sub, IDList_mal, type = 'RGB')
summedRaster_ecu_M_sub <- sumRaster(rasterList_ecu_M_sub, IDList_ecu, type = 'RGB')
summedRaster_amar_M_sub <- sumRaster(rasterList_amar_M_sub, IDList_amar, type = 'RGB')
summedRaster_cyt_M_sub <- sumRaster(rasterList_cyt_M_sub, IDList_cyt, type = 'RGB')
summedRaster_melFG_M_sub <- sumRaster(rasterList_melFG_M_sub, IDList_melFG, type = 'RGB')
summedRaster_melP_M_sub <- sumRaster(rasterList_melP_M_sub, IDList_melP, type = 'RGB')
summedRaster_mer_M_sub <- sumRaster(rasterList_mer_M_sub, IDList_mer, type = 'RGB')
summedRaster_agl_M_sub <- sumRaster(rasterList_agl_M_sub, IDList_agl, type = 'RGB')
summedRaster_nan_M_sub <- sumRaster(rasterList_nan_M_sub, IDList_nan, type = 'RGB')
summedRaster_ple_M_sub <- sumRaster(rasterList_ple_M_sub, IDList_ple, type = 'RGB')
summedRaster_ros_M_sub <- sumRaster(rasterList_ros_M_sub, IDList_ros, type = 'RGB')
summedRaster_vul_M_sub <- sumRaster(rasterList_vul_M_sub, IDList_vul, type = 'RGB')
summedRaster_the_M_sub <- sumRaster(rasterList_the_M_sub, IDList_the, type = 'RGB')
summedRaster_xen_M_sub <- sumRaster(rasterList_xen_M_sub, IDList_xen, type = 'RGB')


colfunc <- inferno(100)
plotHeat(summedRaster_mal_M_sub, IDList_mal, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mal, adjustCoords = TRUE, imageList = imageList_mal, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_ecu_M_sub, IDList_ecu, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_ecu, adjustCoords = TRUE, imageList = imageList_ecu, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_amar_M_sub, IDList_amar, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_amar, adjustCoords = TRUE, imageList = imageList_amar, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_cyt_M_sub, IDList_cyt, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_cyt, adjustCoords = TRUE, imageList = imageList_cyt, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_melFG_M_sub, IDList_melFG, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_melFG, adjustCoords = TRUE, imageList = imageList_melFG, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_melP_M_sub, IDList_melP, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_melP, adjustCoords = TRUE, imageList = imageList_melP, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_mer_M_sub, IDList_mer, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mer, adjustCoords = TRUE, imageList = imageList_mer, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_agl_M_sub, IDList_agl, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_agl, adjustCoords = TRUE, imageList = imageList_agl, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_nan_M_sub, IDList_nan, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_nan, adjustCoords = TRUE,  imageList = imageList_nan, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_ple_M_sub, IDList_ple, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_ple, adjustCoords = TRUE, imageList = imageList, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_ros_M_sub, IDList_ros, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_ros, adjustCoords = TRUE, imageList = imageList_ros, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_vul_M, IDList_vul, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_vul, adjustCoords = TRUE, imageList = imageList_vul, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_the_M_sub, IDList_the, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_the, adjustCoords = TRUE, imageList = imageList_the, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

plotHeat(summedRaster_xen_M_sub, IDList_xen, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_xen, adjustCoords = TRUE, imageList = imageList_xen, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

