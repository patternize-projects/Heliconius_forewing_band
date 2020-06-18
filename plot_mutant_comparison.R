library(patternize)
library(viridis)
library(RColorBrewer)

###
# Load some necessary data
###

IDList_hydFG <- c('BC0004','BC0049','BC0050','BC0061','BC0071','BC0076','BC0077','BC0079','BC0082','BC0125')
IDList_dem   <- c('IMG_1960','IMG_1972','IMG_1974','IMG_1978','IMG_1980','IMG_1982','IMG_1987','IMG_2049','IMG_2125','IMG_2135','IMG_2141')
IDList_ros  <- c('CAM000903','CAM000947','CAM001015','CAM001027','CAM001067', 'CAM001137','CAM001391','CAM002901','CAM008052','CAM009554')

landList_dem <- makeList(IDList_dem, 'landmark', 'landmarks/H.e.demophoon', '_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_hydFG <- makeList(IDList_hydFG, 'landmark', 'landmarks/H.e.hydaraFG', '_landmarks_LFW.txt', skipLandmark = c(2:5,7:9))
landList_ros <- makeList(IDList_ros, 'landmark', 'landmarks/H.m.rosina', '_d.txt', skipLandmark = c(2:5,7:9))

imageList_dem <- makeList(IDList_dem, 'image', 'images/H.e.demophoon', '.JPG')
imageList_hydFG <- makeList(IDList_hydFG, 'image', 'images/H.e.hydaraFG', '-D.JPG')
imageList_ros <- makeList(IDList_ros, 'image', 'images/H.m.rosina', '_d.JPG')

load('aligned_rasterLists/rasterList_dem_M_sub.rda')
load('aligned_rasterLists/rasterList_ros_M_sub.rda')

summedRaster_dem_sub_m <- sumRaster(rasterList_dem_M_sub, IDList_dem, type = 'RGB')
summedRaster_ros_sub_m <- sumRaster(rasterList_ros_M_sub, IDList_ros, type = 'RGB')



###
# Run patternize for mutants
###

IDList_mutDem  <- c('EratoDem_mut9_LDFW','EratoDem_mut9_RDFW','EratoDem_mut11_LDFW','EratoDem_mut11_RDFW','MU4_DLFW','MU4_DRFW','MU5_DLFW','MU5_DRFW','MU6_DLFW','MU6_DRFW')
IDList_mutRos <- c('Hmr_mut1_DRFW','Hmr_mut2_DLFW','Hmr_mut2_DRFW','Hmr_mut3_DLFW','Hmr_mut3_DRFW','Hmr_mut4_DLFW','Hmr_mut4_DRFW','rosina_mut_dfw','Rosina-mutant-FK1-1','Rosina-mutant-FK1-2')


IDlist <- c('WOM7998_RDFW')
library("jpeg")
library("tiff")
for(e in 1:length(IDlist)){
  img <- readTIFF(paste(IDlist[e], ".tif", sep=""), native=TRUE)
  writeJPEG(img, target = paste(IDlist[e], ".jpg", sep=""), quality = 1)
}


landList_mutDem <- makeList(IDList_mutDem, 'landmark', 'landmarks/mutants_demophoon/', '.txt', skipLandmark = c(2:5,7:9))
landList_mutRos <- makeList(IDList_mutRos, 'landmark', 'landmarks/mutants_rosina/', '.txt', skipLandmark = c(2:5,7:9))

imageList_mutDem <- makeList(IDList_mutDem, 'image', 'images/mutants_demophoon', '.jpg')
imageList_mutRos <- makeList(IDList_mutRos, 'image', 'images/mutants_rosina', '.jpg')

# target
target <- as.matrix(read.table('cartoon/BC0004_landmarks_LFW.txt',h = F))
target <- target[c(-2:-5,-7:-9),]

# cartoon
outline_BC0004 <- read.table('cartoon/BC0004_outline.txt', h = F)
lines_BC0004 <- list.files(path ='cartoon', pattern ='BC0004_vein', full.names = T)

# sampleRGB(imageList_mutDem[[1]])
RGB <- c(203,43,49)
rasterList_mutDem_sub <- patLanRGB(imageList_mutDem, landList_mutDem, RGB, transformRef = target, resampleFactor = 1,
                                 colOffset = 0.2, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack', iterations = 3)
save(rasterList_mutDem_sub, file = 'aligned_rasterLists/rasterList_mutDem_sub.rda')

# sampleRGB(imageList_mutRos[[1]])
RGB <- c(214,105,84)
rasterList_mutRos_red_sub <- patLanRGB(imageList_mutRos, landList_mutRos, RGB, transformRef = target, resampleFactor = 1,
                                 colOffset = 0.15, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack', iterations = 3)

# sampleRGB(imageList_mutRos[[8]], type='area')
RGB <- c(221,226,176)
rasterList_mutRos_yellow_sub <- patLanRGB(imageList_mutRos, landList_mutRos, RGB, transformRef = target, resampleFactor = 1,
                                       colOffset = 0.15, crop = FALSE, res = 200, adjustCoords = TRUE, plot = 'stack')


# combine red and yellow
rasterList_mutRos_sub <- list()
for(e in IDList_mutRos){
  print(e)
  subRasterList <- list(rasterList_mutRos_red_sub[[e]], rasterList_mutRos_yellow_sub[[e]])

  names(subRasterList) <- NULL
  subRasterList$fun <- sum
  subRasterList$na.rm <- TRUE
  summedRaster <- do.call(raster::mosaic,subRasterList)
  summedRaster[summedRaster == 2] <- 1
  rasterList_mutRos_sub[[e]] <- summedRaster
}

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
summedRaster_mutRos_sub <- sumRaster(rasterList_mutRos_sub, IDList_mutRos, type = 'RGB')
setMask(summedRaster_mutRos_sub, IDList_mutRos, refImage = target, 'mask/mask_mut_ros.txt')

mask_mut_ros <- read.table("mask/mask_mut_ros.txt", h=T)
rasterList_mutRos_sub_M<-list()
for(e in 1:length(rasterList_mutRos_sub)){
  ID <- names(rasterList_mutRos_sub)[[e]]
  rasterList_mutRos_sub_M[[ID]] <- maskOutline(rasterList_mutRos_sub[[ID]], IDlist = IDList_mutRos, mask_mut_ros, 
                                        refShape = 'target', imageList = imageList_mutRos)
}

save(rasterList_mutRos_sub_M, file = 'aligned_rasterLists/rasterList_mutRos_sub_M.rda')


###
# Plot heatmaps mutants
###

load('aligned_rasterLists/rasterList_mutDem_sub.rda')
load('aligned_rasterLists/rasterList_mutRos_sub_M.rda')

summedRaster_mutDem_sub <- sumRaster(rasterList_mutDem_sub, IDList_mutDem, type = 'RGB')
summedRaster_mutRos_sub_M <- sumRaster(rasterList_mutRos_sub_M, IDList_mutRos, type = 'RGB')

summedRaster_mutDem_sub_int <- summedRaster_mutDem_sub
summedRaster_mutDem_sub_int[summedRaster_mutDem_sub_int < 5] <- 0
summedRaster_mutDem_sub_int[summedRaster_mutDem_sub_int >= 5] <- 1

summedRaster_mutRos_sub_int <- summedRaster_mutRos_sub_M
summedRaster_mutRos_sub_int[summedRaster_mutRos_sub_int < 5] <- 0
summedRaster_mutRos_sub_int[summedRaster_mutRos_sub_int >= 5] <- 1

# colfunc <- inferno(100)
colfunc <- colorRampPalette(brewer.pal(9,'Blues'))(100)[100:1]
plotHeat(summedRaster_mutDem_sub, IDList_mutDem, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutDem, adjustCoords = TRUE, imageList = imageList_hydFG, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, legend = F)

raster::contour(summedRaster_mutRos_sub_int, add=T, col = 'yellow', lwd=3, maxpixels=5000, nlevels=1)


plotHeat(summedRaster_dem_sub_m, IDList_dem, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004,
         lines = lines_BC0004, landList = landList_dem, adjustCoords = TRUE, imageList = imageList_hydFG,
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, legend = F)

colfunc <- c("blue","lightblue","black","burlywood1","orange")
raster_diff <- summedRaster_mutDem_sub/length(IDList_mutDem) - summedRaster_dem_sub_m/length(IDList_dem)

plotHeat(raster_diff, IDList_mutDem, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutDem, adjustCoords = TRUE, imageList = imageList_hydFG, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, normalized = T, 
         zlim = c(-1,1), legendTitle = "", legend = F)

raster_diff_min0 <- raster::calc(raster_diff, function(x){x[x<0]<-0; return(x)})

plotHeat(raster_diff_min0, IDList_mutDem, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutDem, adjustCoords = TRUE, imageList = imageList_hydFG, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, normalized = T, 
         zlim = c(-1,1), legendTitle = "", legend = F)


colfunc <- colorRampPalette(brewer.pal(9,'Blues'))(100)[100:1]
plotHeat(summedRaster_mutRos_sub_M, IDList_mutRos, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutRos, adjustCoords = TRUE, imageList = imageList_hydFG, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, legend = F)

raster::contour(summedRaster_mutDem_sub_int, add=T, col = 'green', lwd=3, maxpixels=5000, nlevels=1)


plotHeat(summedRaster_ros_sub_m, IDList_ros, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004,
         lines = lines_BC0004, landList = landList_ros, adjustCoords = TRUE, imageList = imageList_hydFG,
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, legend = F)

colfunc <- c("blue","lightblue","black","burlywood1","orange")
raster_diff <- summedRaster_mutRos_sub_M/length(IDList_mutRos) - summedRaster_ros_sub_m/length(IDList_ros)

plotHeat(raster_diff, IDList_mutRos, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutRos, adjustCoords = TRUE, imageList = imageList_hydFG, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, normalized = T, 
         zlim = c(-1,1), legendTitle = "", legend = F)

raster_diff_min0 <- raster::calc(raster_diff, function(x){x[x<0]<-0; return(x)})

### Calculation of area of wing impacted
# setMask(raster_diff, IDList_mutRos, colpalette = colfunc, zlim = c(-1,1), normalized = T, 'mask/mask_mut1.txt')
# setMask(raster_diff, IDList_mutRos, colpalette = colfunc, zlim = c(-1,1), normalized = T, 'mask/mask_mut2.txt')

mask_mut2 <- read.table("mask/mask_mut2.txt", h=T)

raster_diff_M <- maskOutline(raster_diff, mask_mut2, refShape = 'target', 
                                            imageList = imageList_hydFG)

plotHeat(raster_diff_M, IDList_mutRos, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutRos, adjustCoords = TRUE, imageList = imageList_hydFG, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, normalized = T, 
         zlim = c(-1,1), legendTitle = "", legend = F)

raster_diff_min0 <- raster::calc(raster_diff_M, function(x){x[x<0]<-0; return(x)})



##
# calculate area outline
rasterEx <- raster::extent(raster_diff)
rRe <- raster::raster(nrow=200,ncol=200)
raster::extent(rRe) <- rasterEx
newRaster <- raster::resample(raster_diff_min0, rRe)

poly <- sp::Polygons(list(sp::Polygon(outline_BC0004)),paste("r"))

polyList  <- c(poly)
polyNames <- c(paste("r"))
sr=sp::SpatialPolygons(polyList)
srdf=sp::SpatialPolygonsDataFrame(sr, data.frame(1:length(polyNames), row.names=polyNames))

r <- raster::raster(rasterEx, nrow=dim(newRaster)[1], ncol=dim(newRaster)[2])
rr <-raster::rasterize(srdf, r)

nrCellsOutline  <- raster::freq(rr, value=1)

# calculate area difference (normalized for average difference)
df <- abs(raster::as.data.frame(raster_diff_min0))

df2 <- subset(df, df$layer != 0)

AvDiff <- sum(df$layer)/nrCellsOutline
sdDiff <- (sd(df$layer)*nrow(df2))/nrCellsOutline
##


####
summedRaster_mutDem_sub_min <- summedRaster_mutDem_sub
summedRaster_mutDem_sub_min[summedRaster_mutDem_sub_min >= 1] <- 1

plotHeat(summedRaster_mutDem_sub_min, IDList_mutDem[1], plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutDem, adjustCoords = TRUE, imageList = imageList_mutDem, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

summedRaster_mutDem_sub_int <- summedRaster_mutDem_sub
summedRaster_mutDem_sub_int[summedRaster_mutDem_sub_int < 2] <- 0
summedRaster_mutDem_sub_int[summedRaster_mutDem_sub_int >= 2] <- 1

colfunc <- c('green')
plotHeat(summedRaster_mutDem_sub_int, IDList_mutDem[1], plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutDem, adjustCoords = TRUE, imageList = imageList_mutDem, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, legend = F)

summedRaster_mutDem_sub_max <- summedRaster_mutDem_sub
summedRaster_mutDem_sub_max[summedRaster_mutDem_sub_max < 10] <- 0
summedRaster_mutDem_sub_max[summedRaster_mutDem_sub_max == 10] <- 1


plotHeat(summedRaster_mutDem_sub_max, IDList_mutDem[1], plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutDem, adjustCoords = TRUE, imageList = imageList_mutDem, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)
####
summedRaster_mutRos_sub_min <- summedRaster_mutRos_sub_M
summedRaster_mutRos_sub_min[summedRaster_mutRos_sub_min >= 1] <- 1

plotHeat(summedRaster_mutRos_sub_min, IDList_mutRos[1], plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutRos, adjustCoords = TRUE, imageList = imageList_mutRos, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)

summedRaster_mutRos_sub_int <- summedRaster_mutRos_sub_M
summedRaster_mutRos_sub_int[summedRaster_mutRos_sub_int < 2] <- 0
summedRaster_mutRos_sub_int[summedRaster_mutRos_sub_int >= 2] <- 1

colfunc <- c('yellow')
plotHeat(summedRaster_mutRos_sub_int, IDList_mutRos[1], plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutRos, adjustCoords = TRUE, imageList = imageList_mutRos, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, legend =F)

summedRaster_mutRos_sub_max <- summedRaster_mutRos_sub_M
summedRaster_mutRos_sub_max[summedRaster_mutRos_sub_max < 10] <- 0
summedRaster_mutRos_sub_max[summedRaster_mutRos_sub_max == 10] <- 1

plotHeat(summedRaster_mutRos_sub_max, IDList_mutRos[1], plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutRos, adjustCoords = TRUE, imageList = imageList_mutRos, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc)


###
# Compare mutants to demophoon/rosina
###

load('aligned_rasterLists/rasterList_dem_M_sub.rda')
load('aligned_rasterLists/rasterList_ros_M_sub.rda')

IDList_dem   <- c('IMG_1960','IMG_1972','IMG_1974','IMG_1978','IMG_1980','IMG_1982','IMG_1987','IMG_2049','IMG_2125','IMG_2135','IMG_2141')
IDList_ros  <- c('CAM000903','CAM000947','CAM001015','CAM001027','CAM001067', 'CAM001137','CAM001391','CAM002901','CAM008052','CAM009554')

summedRaster_dem_M_sub <- sumRaster(rasterList_dem_M_sub, IDList_dem, type = 'RGB')
summedRaster_ros_M_sub <- sumRaster(rasterList_ros_M_sub, IDList_ros, type = 'RGB')

# save(summedRaster_dem_M_sub, file = 'aligned_rasterLists/summedRaster_dem_M_sub.rda')
# save(summedRaster_ros_M_sub, file = 'aligned_rasterLists/summedRaster_ros_M_sub.rda')


colfunc <- c("blue","lightblue","black","burlywood1","orange")
raster_dem_ros <- summedRaster_dem_M_sub/length(IDList_dem) - summedRaster_ros_M_sub/length(IDList_ros)
plotHeat(raster_dem_ros, IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList, adjustCoords = TRUE, 
         imageList = imageList, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, normalized = T, zlim = c(-1,1), legendTitle = "", legend=F)

# raster_ven_vul <- summedRaster_ven_M_sub/length(IDList_ven) - summedRaster_vul_M_sub/length(IDList_vul)
# plotHeat(raster_ven_vul, IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
#          lines = lines_BC0004, landList = landList, adjustCoords = TRUE, 
#          imageList = imageList, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
#          colpalette = colfunc, normalized = T, zlim = c(-1,1), legendTitle = "", legend=F)
# 
# colfunc <- c("blue","lightblue","black","burlywood1","orange")
# raster_hydP_melP <- summedRaster_hydP_M_sub/length(IDList_hydP) - summedRaster_melP_M_sub/length(IDList_melP)
# plotHeat(raster_hydP_melP, IDList_hydP, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
#          lines = lines_BC0004, landList = landList, adjustCoords = TRUE, 
#          imageList = imageList, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
#          colpalette = colfunc, normalized = T, zlim = c(-1,1), legendTitle = "", legend=F)

raster::contour(summedRaster_mutDem_sub_int, add=T, col = 'green', lwd=3, maxpixels=5000, nlevels=1)
# raster::contour(summedRaster_mutDem_sub_int, add=T, col = 'green', lwd=1, maxpixels=5000, nlevels=1)
# raster::contour(summedRaster_mutDem_sub_max, add=T, col = 'green', lwd=0.5, maxpixels=5000, nlevels=1)
raster::contour(summedRaster_mutRos_sub_int, add=T, col = 'yellow', lwd=3, maxpixels=5000, nlevels=1)


### Calculation of area of wing impacted
setMask(raster_dem_ros, IDlist, colpalette = colfunc, zlim = c(-1,1), normalized = T, 'mask/mask_mut3.txt')

mask_mut3 <- read.table("mask/mask_mut3.txt", h=T)

raster_dem_ros_M <- maskOutline(raster_dem_ros, mask_mut3, refShape = 'target', 
                             imageList = imageList_hydFG)

plotHeat(raster_dem_ros_M, IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutRos, adjustCoords = TRUE, imageList = imageList_hydFG, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, normalized = T, 
         zlim = c(-1,1), legendTitle = "", legend = F)

raster_diff_min0 <- raster::calc(raster_dem_ros_M, function(x){x[x<0]<-0; return(x)})



##
# calculate area outline
rasterEx <- raster::extent(raster_diff)
rRe <- raster::raster(nrow=200,ncol=200)
raster::extent(rRe) <- rasterEx
newRaster <- raster::resample(raster_diff_min0, rRe)

poly <- sp::Polygons(list(sp::Polygon(outline_BC0004)),paste("r"))

polyList  <- c(poly)
polyNames <- c(paste("r"))
sr=sp::SpatialPolygons(polyList)
srdf=sp::SpatialPolygonsDataFrame(sr, data.frame(1:length(polyNames), row.names=polyNames))

r <- raster::raster(rasterEx, nrow=dim(newRaster)[1], ncol=dim(newRaster)[2])
rr <-raster::rasterize(srdf, r)

nrCellsOutline  <- raster::freq(rr, value=1)

# calculate area difference (normalized for average difference)
df <- abs(raster::as.data.frame(raster_diff_min0))

df2 <- subset(df, df$layer != 0)

AvDiff <- sum(df$layer)/nrCellsOutline
sdDiff <- (sd(df$layer)*nrow(df2))/nrCellsOutline
##


###
# Compare postman pattern to mutants
###

load('aligned_rasterLists/pcaOut_postman_sub.rda')

pcdata <- pcaOut[[3]]$x
rotation <- pcaOut[[3]]$rotation

PCx <- 1
PCy <- 2

PCxmin <- min(pcdata[,PCx])
PCxmax <- max(pcdata[,PCx])

PCymin <- min(pcdata[,PCy])
PCymax <- max(pcdata[,PCy])


pc.vecMix <- rep(0, dim(pcdata)[1])
pc.vecMix[PCx] <- PCxmin

pc.vecMax <- rep(0, dim(pcdata)[1])
pc.vecMax[PCx] <- PCxmax

pc.vecMiy <- rep(0, dim(pcdata)[1])
pc.vecMiy[PCy] <- PCymin

pc.vecMay <- rep(0, dim(pcdata)[1])
pc.vecMay[PCy] <- PCymax


xMi <- pc.vecMix %*%  t(rotation)
xMa <- pc.vecMax %*%  t(rotation)

x2Mi <- t(matrix(xMi,ncol = dim(raster_dem_ros)[1], nrow = dim(raster_dem_ros)[2]))
x2Ma <- t(matrix(xMa,ncol = dim(raster_dem_ros)[1], nrow = dim(raster_dem_ros)[2]))

yMi <- pc.vecMiy %*%  t(rotation)
yMa <- pc.vecMay %*%  t(rotation)

y2Mi <- t(matrix(yMi,ncol = dim(raster_dem_ros)[1], nrow = dim(raster_dem_ros)[2]))
y2Ma <- t(matrix(yMa,ncol = dim(raster_dem_ros)[1], nrow = dim(raster_dem_ros)[2]))


mapMix <-raster::raster(x2Mi)
mapMax <-raster::raster(x2Ma)

mapMiy <-raster::raster(y2Mi)
mapMay <-raster::raster(y2Ma)

raster::extent(mapMix) <- raster::extent(raster_dem_ros)
raster::extent(mapMax) <- raster::extent(raster_dem_ros)

raster::extent(mapMiy) <- raster::extent(raster_dem_ros)
raster::extent(mapMay) <- raster::extent(raster_dem_ros)

colfunc <- c("blue","lightblue","black","indianred1","firebrick1")
plotHeat(mapMiy/max(abs(xMi)), IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList, adjustCoords = TRUE, 
         imageList = imageList, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, normalized = T, zlim = c(-1,1), legendTitle = "", legend=F)

raster::contour(summedRaster_mutDem_sub_int, add=T, col = 'green', lwd=3, maxpixels=5000, nlevels=1)
# raster::contour(summedRaster_mutDem_sub_int, add=T, col = 'green', lwd=1, maxpixels=5000, nlevels=1)
# raster::contour(summedRaster_mutDem_sub_max, add=T, col = 'green', lwd=0.5, maxpixels=5000, nlevels=1)
raster::contour(summedRaster_mutRos_sub_int, add=T, col = 'yellow', lwd=3, maxpixels=5000, nlevels=1)



# empty <- raster()
# extent(empty) <- extent(summedRaster_mutDem_sub_int)
# plotHeat(empty, IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
#          lines = lines_BC0004, landList = landList, adjustCoords = TRUE, 
#          imageList = imageList, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
#          colpalette = colfunc, normalized = T, zlim = c(-1,1), legendTitle = "", legend=F)


# setMask(mapMiy/max(abs(xMi)), IDlist, colpalette = colfunc, zlim = c(-1,1), normalized = T, 'mask/mask_mut4.txt')
# setMask(mapMiy/max(abs(xMi)), IDlist, colpalette = colfunc, zlim = c(-1,1), normalized = T, 'mask/mask_mut5.txt')

mask_mut5 <- read.table("mask/mask_mut5.txt", h=T)

mapMiy_M <- maskOutline(mapMiy/max(abs(xMi)), mask_mut5, refShape = 'target', 
                                imageList = imageList_hydFG)

plotHeat(mapMiy_M, IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landList_mutRos, adjustCoords = TRUE, imageList = imageList_hydFG, 
         cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, normalized = T, 
         zlim = c(-1,1), legendTitle = "", legend = F)

raster_diff_min0 <- raster::calc(mapMiy_M, function(x){x[x>0]<-0; return(x)})



##
# calculate area outline
rasterEx <- raster::extent(raster_diff)
rRe <- raster::raster(nrow=200,ncol=200)
raster::extent(rRe) <- rasterEx
newRaster <- raster::resample(raster_diff_min0, rRe)

poly <- sp::Polygons(list(sp::Polygon(outline_BC0004)),paste("r"))

polyList  <- c(poly)
polyNames <- c(paste("r"))
sr=sp::SpatialPolygons(polyList)
srdf=sp::SpatialPolygonsDataFrame(sr, data.frame(1:length(polyNames), row.names=polyNames))

r <- raster::raster(rasterEx, nrow=dim(newRaster)[1], ncol=dim(newRaster)[2])
rr <-raster::rasterize(srdf, r)

nrCellsOutline  <- raster::freq(rr, value=1)

# calculate area difference (normalized for average difference)
df <- abs(raster::as.data.frame(raster_diff_min0))

df2 <- subset(df, df$layer != 0)

AvDiff <- sum(df$layer)/nrCellsOutline
sdDiff <- (sd(df$layer)*nrow(df2))/nrCellsOutline
##

