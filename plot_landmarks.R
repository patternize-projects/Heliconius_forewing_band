library(patternize)
library(Morpho)
library(Momocs)

############################
# Modified Momocs functions
############################
tps_iso2 <- function(fr, to, amp = 1,
                     grid = FALSE, outline = NULL, over = 1.2, palette = col_spring,
                     iso.nb = 1000, iso.levels = 12,
                     cont = TRUE, cont.col = "black",
                     poly = TRUE, shp = TRUE, shp.border = col_qual(2),
                     shp.lwd = c(2, 2), shp.lty = c(1, 1),
                     legend = TRUE, legend.text, maxArr = NULL, ...) {
  fr.n <- substitute(fr)
  to.n <- substitute(to)  # otherwise problems with substitute in legend below
  if (!missing(amp))
    to <- to + (to - fr) * amp
  if (grid) {
    grid0 <- .grid.sample(fr, to, nside = round(sqrt(iso.nb)), over = over)
  } else {
    grid0 <- sp::spsample(sp::Polygon(coo_close(fr)), iso.nb, type='regular')@coords
  }
  if (!is.null(outline)){
    grid0 <- sp::spsample(sp::Polygon(coo_close(outline)), iso.nb, type='regular')@coords
  }
  grid1 <- tps2d(grid0, fr, to)
  def <- edm(grid0, grid1)
  print(max(edm(grid0, grid1)))
  x1 <- length(unique(grid0[, 1]))
  y1 <- length(unique(grid0[, 2]))
  im <- matrix(NA, x1, y1)
  xind <- (1:x1)[as.factor(rank(grid0[, 1]))]
  yind <- (1:y1)[as.factor(rank(grid0[, 2]))]
  n <- length(xind)
  for (i in 1:n) im[xind[i], yind[i]] <- def[i]
  if(!is.null(maxArr)){
    iso.cols <- palette(iso.levels)[1:as.integer(iso.levels*(max(edm(grid0, grid1))/maxArr))]
  }
  else{
    iso.cols <- palette(iso.levels)
  }
  x <- sort(unique(grid0[, 1]))
  y <- sort(unique(grid0[, 2]))
  op <- par(mar = rep(1, 4))
  on.exit(par(op))
  image(x, y, im, col = iso.cols, asp = 1,
        xlim = range(grid0[, 1])*over,
        ylim = range(grid0[, 2])*over,
        axes = FALSE, frame = FALSE,
        ann = FALSE)
  if (cont) {
    contour(x, y, im, nlevels = iso.levels,
            add = TRUE, drawlabels = FALSE,
            col = cont.col, lty = 2)
  }
  if (shp) {
    points <- ifelse(poly, FALSE, TRUE)
    coo_draw(fr, border = shp.border[1], col = NA, lwd = shp.lwd[1],
             lty = shp.lty[1], points = points, first.point = FALSE,
             centroid = FALSE, ...)
    coo_draw(to, border = shp.border[2], col = NA, lwd = shp.lwd[2],
             lty = shp.lty[2], points = points, first.point = FALSE,
             centroid = FALSE, ...)
    if (legend | !missing(legend.text)) {
      if (missing(legend.text)) legend.text <- c(fr.n, to.n)
      legend("topright", legend = legend.text, col = shp.border,
             lwd = shp.lwd, bty = "n")
    }
  }
}

tps_arr2 <- function(fr, to, amp = 1,
                     grid = FALSE, outline = NULL, over = 1.2, palette = col_summer,
                     arr.nb = 200, arr.levels = 100, arr.len = 0.1, arr.ang = 20,
                     arr.lwd = 0.75, arr.col = "grey50", poly = TRUE,
                     shp = TRUE, shp.col = rep(NA, 2),
                     shp.border = col_qual(2),
                     shp.lwd = c(2, 2), shp.lty = c(1, 1),
                     legend = TRUE, legend.text, maxArr = NULL, ...) {
  fr.n <- substitute(fr)
  to.n <- substitute(to)  # otherwise problems with substitute in legend below
  if (!missing(amp))
    to <- to + (to - fr) * amp
  if (grid){
    grid0 <- .grid.sample(fr, to, nside = round(sqrt(arr.nb)), over = over)
  }
  else {
    grid0 <- sp::spsample(sp::Polygon(coo_close(fr)), arr.nb, type='regular')@coords
  }
  if (!is.null(outline)){
    grid0 <- sp::spsample(sp::Polygon(coo_close(outline)), arr.nb, type='regular')@coords
  }
  grid1 <- tps2d(grid0, fr, to)
  # grille simple, on affiche d'abord les deux courbes
  op <- par(mar = rep(0, 4))
  on.exit(par(op))
  plot(NA, xlim = range(grid0[, 1])*over, ylim = range(grid0[, 2])*over,
       asp = 1, axes = FALSE, ann = FALSE, mar = rep(0, 4))
  if (missing(arr.levels)) {
    arr.levels = arr.nb
  }
  if (!missing(palette)) {
    print(max(edm(grid0, grid1)))
    if(!is.null(maxArr)){
      q.lev <- (edm(grid0, grid1)/maxArr)*arr.nb
    }
    else{
      q.lev <- cut(edm(grid0, grid1), breaks = arr.levels,
                   labels = FALSE)
    }
    
    
    arr.cols <- palette(arr.levels)[q.lev]
  } else {
    arr.cols <- rep(arr.col, nrow(grid0))
  }
  arrows(grid0[, 1], grid0[, 2], grid1[, 1], grid1[, 2], length = arr.len,
         angle = arr.ang, lwd = arr.lwd, col = arr.cols)
  if (shp) {
    points <- ifelse(poly, FALSE, TRUE)
    coo_draw(fr, border = shp.border[1], col = NA, lwd = shp.lwd[1],
             lty = shp.lty[1], points = points, first.point = FALSE,
             centroid = FALSE, ...)
    coo_draw(to, border = shp.border[2], col = NA, lwd = shp.lwd[2],
             lty = shp.lty[2], points = points, first.point = FALSE,
             centroid = FALSE, ...)
    if (legend | !missing(legend.text)) {
      if (missing(legend.text)) legend.text <- c(fr.n, to.n)
      legend("topright", legend = legend.text, col = shp.border,
             lwd = shp.lwd, bty = "n")
    }
  }
}

############################
# Sample and landmark lists
############################

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
landList_hydP <- makeList(IDList_hydP, 'landmark', 'landmarks/H.e.hydaraP', '_LFW_landmarks.txt')
landList_hydFG <- makeList(IDList_hydFG, 'landmark', 'landmarks/H.e.hydaraFG', '_landmarks_LFW.txt')
landList_amal <- makeList(IDList_amal, 'landmark', 'landmarks/H.e.amalfreda', '-D.txt')
landList_cyr <- makeList(IDList_cyr, 'landmark', 'landmarks/H.e.cyrbia', '_d_LFW_landmarks.txt')
landList_not <- makeList(IDList_not, 'landmark', 'landmarks/H.e.notabilis', '_d_landmarks.txt')
landList_ety <- makeList(IDList_ety, 'landmark', 'landmarks/H.e.etylus', '_d_LFW_landmarks.txt')
landList_fav <- makeList(IDList_fav, 'landmark', 'landmarks/H.e.favorinus', '_d_LFW_landmarks.txt')
landList_ven <- makeList(IDList_ven, 'landmark', 'landmarks/H.e.venus', '_d_LFW_landmarks.txt')
landList_era <- makeList(IDList_era, 'landmark', 'landmarks/H.e.erato', '_landmarks_LFW.txt')
landList_emm <- makeList(IDList_emm, 'landmark', 'landmarks/H.e.emma', '_d_LFW_landmarks.txt')
landList_phy <- makeList(IDList_phy, 'landmark', 'landmarks/H.e.phyllis', '_calibrated_D_LFW_landmarks.txt')
landList_dem <- makeList(IDList_dem, 'landmark', 'landmarks/H.e.demophoon', '_LFW_landmarks.txt')
landList_lat <- makeList(IDList_lat, 'landmark', 'landmarks/H.e.lativitta', '_d_landmarks.txt')
landList_mic <- makeList(IDList_mic, 'landmark', 'landmarks/H.e.microclea', '_D_butterfly_landmarks.txt')

landList_melP <- makeList(IDList_melP, 'landmark', 'landmarks/H.m.melpomeneP', '_d.txt')
landList_melFG <- makeList(IDList_melFG, 'landmark', 'landmarks/H.m.melpomeneFG', '_d.txt')
landList_mer <- makeList(IDList_mer, 'landmark', 'landmarks/H.m.meriana', '.txt')
landList_cyt <- makeList(IDList_cyt, 'landmark', 'landmarks/H.m.cythera', '.txt')
landList_ple <- makeList(IDList_ple, 'landmark', 'landmarks/H.m.plesseni', '_d.txt')
landList_ecu <- makeList(IDList_ecu, 'landmark', 'landmarks/H.m.ecuadorensis', '_d.txt')
landList_amar <- makeList(IDList_amar, 'landmark', 'landmarks/H.m.amaryllis', '-d.txt')
landList_vul <- makeList(IDList_vul, 'landmark', 'landmarks/H.m.vulcanus', '_d.txt')
landList_the <- makeList(IDList_the, 'landmark', 'landmarks/H.m.thelxiopeia', '_D_butterfly.txt')
landList_mal <- makeList(IDList_mal, 'landmark', 'landmarks/H.m.malleti', '_d.txt')
landList_nan <- makeList(IDList_nan, 'landmark', 'landmarks/H.m.nanna', '_Hmnanna_upper.txt')
landList_ros <- makeList(IDList_ros, 'landmark', 'landmarks/H.m.rosina', '_d.txt')
landList_agl <- makeList(IDList_agl, 'landmark', 'landmarks/H.m.aglaope', '_D_butterfly_landmarks.txt')
landList_xen <- makeList(IDList_xen, 'landmark', 'landmarks/H.m.xenoclea', '-d.txt')

# Landmark subset
landList_hydP_sub  <- makeList(IDList_hydP, 'landmark', 'landmarks/H.e.hydaraP', '_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_hydFG_sub <- makeList(IDList_hydFG, 'landmark', 'landmarks/H.e.hydaraFG', '_landmarks_LFW.txt', skipLandmark = c(2:5,7:9))
landList_amal_sub  <- makeList(IDList_amal, 'landmark', 'landmarks/H.e.amalfreda', '-D.txt', skipLandmark = c(2:5,7:9))
landList_cyr_sub   <- makeList(IDList_cyr, 'landmark', 'landmarks/H.e.cyrbia', '_d_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_not_sub   <- makeList(IDList_not, 'landmark', 'landmarks/H.e.notabilis', '_d_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_ety_sub   <- makeList(IDList_ety, 'landmark', 'landmarks/H.e.etylus', '_d_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_fav_sub   <- makeList(IDList_fav, 'landmark', 'landmarks/H.e.favorinus', '_d_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_ven_sub   <- makeList(IDList_ven, 'landmark', 'landmarks/H.e.venus', '_d_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_era_sub   <- makeList(IDList_era, 'landmark', 'landmarks/H.e.erato', '_landmarks_LFW.txt', skipLandmark = c(2:5,7:9))
landList_emm_sub   <- makeList(IDList_emm, 'landmark', 'landmarks/H.e.emma', '_d_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_phy_sub   <- makeList(IDList_phy, 'landmark', 'landmarks/H.e.phyllis', '_calibrated_D_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_dem_sub   <- makeList(IDList_dem, 'landmark', 'landmarks/H.e.demophoon', '_LFW_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_lat_sub   <- makeList(IDList_lat, 'landmark', 'landmarks/H.e.lativitta', '_d_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_mic_sub   <- makeList(IDList_mic, 'landmark', 'landmarks/H.e.microclea', '_D_butterfly_landmarks.txt', skipLandmark = c(2:5,7:9))

landList_melP_sub  <- makeList(IDList_melP, 'landmark', 'landmarks/H.m.melpomeneP', '_d.txt', skipLandmark = c(2:5,7:9))
landList_melFG_sub <- makeList(IDList_melFG, 'landmark', 'landmarks/H.m.melpomeneFG', '_d.txt', skipLandmark = c(2:5,7:9))
landList_mer_sub   <- makeList(IDList_mer, 'landmark', 'landmarks/H.m.meriana', '.txt', skipLandmark = c(2:5,7:9))
landList_cyt_sub   <- makeList(IDList_cyt, 'landmark', 'landmarks/H.m.cythera', '.txt', skipLandmark = c(2:5,7:9))
landList_ple_sub   <- makeList(IDList_ple, 'landmark', 'landmarks/H.m.plesseni', '_d.txt', skipLandmark = c(2:5,7:9))
landList_ecu_sub   <- makeList(IDList_ecu, 'landmark', 'landmarks/H.m.ecuadorensis', '_d.txt', skipLandmark = c(2:5,7:9))
landList_amar_sub  <- makeList(IDList_amar, 'landmark', 'landmarks/H.m.amaryllis', '-d.txt', skipLandmark = c(2:5,7:9))
landList_vul_sub   <- makeList(IDList_vul, 'landmark', 'landmarks/H.m.vulcanus', '_d.txt', skipLandmark = c(2:5,7:9))
landList_the_sub   <- makeList(IDList_the, 'landmark', 'landmarks/H.m.thelxiopeia', '_D_butterfly.txt', skipLandmark = c(2:5,7:9))
landList_mal_sub   <- makeList(IDList_mal, 'landmark', 'landmarks/H.m.malleti', '_d.txt', skipLandmark = c(2:5,7:9))
landList_nan_sub   <- makeList(IDList_nan, 'landmark', 'landmarks/H.m.nanna', '_Hmnanna_upper.txt', skipLandmark = c(2:5,7:9))
landList_ros_sub   <- makeList(IDList_ros, 'landmark', 'landmarks/H.m.rosina', '_d.txt', skipLandmark = c(2:5,7:9))
landList_agl_sub   <- makeList(IDList_agl, 'landmark', 'landmarks/H.m.aglaope', '_D_butterfly_landmarks.txt', skipLandmark = c(2:5,7:9))
landList_xen_sub   <- makeList(IDList_xen, 'landmark', 'landmarks/H.m.xenoclea', '-d.txt', skipLandmark = c(2:5,7:9))

# ID list and color vector for plotting
IDlist_ALL <- c(IDList_hydP, IDList_hydFG, IDList_amal, IDList_cyr, IDList_not, IDList_ety, IDList_fav, IDList_ven, 
                IDList_era, IDList_emm, IDList_phy, IDList_dem, IDList_lat, IDList_mic,
                
                IDList_melP, IDList_melFG, IDList_mer, IDList_cyt, IDList_ple, IDList_ecu, IDList_amar, IDList_vul,
                IDList_the, IDList_mal, IDList_nan, IDList_ros, IDList_agl, IDList_xen)

colVec <- c(rep('orange',length(landList_hydP)),rep('orange',length(landList_hydFG)),rep('orange',length(landList_amal)),
            rep('orange',length(landList_cyr)),rep('orange',length(landList_not)),rep('orange',length(landList_ety)),
            rep('orange',length(landList_fav)),rep('orange',length(landList_ven)),rep('orange',length(landList_era)),
            rep('orange',length(landList_emm)),rep('orange',length(landList_phy)),rep('orange',length(landList_dem)),
            rep('orange',length(landList_lat)),rep('orange',length(landList_mic)),
            
            rep('royalblue1',length(landList_melP)),rep('royalblue1',length(landList_melFG)),rep('royalblue1',length(landList_mer)),
            rep('royalblue1',length(landList_cyt)),rep('royalblue1',length(landList_ple)),rep('royalblue1',length(landList_ecu)),
            rep('royalblue1',length(landList_amar)),rep('royalblue1',length(landList_vul)),rep('royalblue1',length(landList_the)),
            rep('royalblue1',length(landList_mal)),rep('royalblue1',length(landList_nan)),rep('royalblue1',length(landList_ros)),
            rep('royalblue1',length(landList_agl)),rep('royalblue1',length(landList_xen)))

# Landmark array all 18 landmarks
landArray <- lanArray(c(landList_hydP, landList_hydFG, landList_amal, landList_cyr, landList_not, landList_ety, landList_fav, landList_ven, 
                        landList_era, landList_emm, landList_phy, landList_dem, landList_lat, landList_mic,
                        
                        landList_melP, landList_melFG, landList_mer, landList_cyt, landList_ple, landList_ecu, landList_amar, landList_vul,
                        landList_the, landList_mal, landList_nan, landList_ros, landList_agl, landList_xen))

landArray_E <- lanArray(c(landList_hydP, landList_hydFG, landList_amal, landList_cyr, landList_not, landList_ety, landList_fav, landList_ven, 
                        landList_era, landList_emm, landList_phy, landList_dem, landList_lat, landList_mic))

landArray_M <- lanArray(c(landList_melP, landList_melFG, landList_mer, landList_cyt, landList_ple, landList_ecu, landList_amar, landList_vul,
                        landList_the, landList_mal, landList_nan, landList_ros, landList_agl, landList_xen))
                        
# Landmark array subset landmarks
landArray_sub <- lanArray(c(landList_hydP_sub, landList_hydFG_sub, landList_amal_sub, landList_cyr_sub, landList_not_sub, landList_ety_sub, landList_fav_sub, landList_ven_sub, 
                        landList_era_sub, landList_emm_sub, landList_phy_sub, landList_dem_sub, landList_lat_sub, landList_mic_sub,
                        
                        landList_melP_sub, landList_melFG_sub, landList_mer_sub, landList_cyt_sub, landList_ple_sub, landList_ecu_sub, landList_amar_sub, landList_vul_sub,
                        landList_the_sub, landList_mal_sub, landList_nan_sub, landList_ros_sub, landList_agl_sub, landList_xen_sub))

landArray_E_sub <- lanArray(c(landList_hydP_sub, landList_hydFG_sub, landList_amal_sub, landList_cyr_sub, landList_not_sub, landList_ety_sub, landList_fav_sub, landList_ven_sub, 
                          landList_era_sub, landList_emm_sub, landList_phy_sub, landList_dem_sub, landList_lat_sub, landList_mic_sub))

landArray_M_sub <- lanArray(c(landList_melP_sub, landList_melFG_sub, landList_mer_sub, landList_cyt_sub, landList_ple_sub, landList_ecu_sub, landList_amar_sub, landList_vul_sub,
                          landList_the_sub, landList_mal_sub, landList_nan_sub, landList_ros_sub, landList_agl_sub, landList_xen_sub))

# Calculate transformation
transformed      <- Morpho::procSym(landArray)
transformed_E    <- Morpho::procSym(landArray_E)
transformed_M    <- Morpho::procSym(landArray_M)
transformed_sub  <- Morpho::procSym(landArray_sub)
transformed_sub_E<- Morpho::procSym(landArray_E_sub)
transformed_sub_M<- Morpho::procSym(landArray_M_sub)

############################
# Transform cartoon data to mean shape
############################
outlineLan <- as.matrix(read.table('landmarks/cartoon/BC0004_landmarks_LFW.txt'))
outline <- as.matrix(read.table('landmarks/cartoon/BC0004_outline.txt'))
lines <- list.files(path ='landmarks/cartoon', pattern ='BC0004_vein', full.names = T)

cartoonLandTrans <- Morpho::computeTransform(transformed$mshape, outlineLan, type="tps")
cartoonLandTrans_sub <- Morpho::computeTransform(transformed_sub$mshape, outlineLan[c(-2:-5,-7:-9),], type="tps")

outlineTrans <- Morpho::applyTransform(as.matrix(outline), cartoonLandTrans)
outlineTrans_sub <- Morpho::applyTransform(as.matrix(outline), cartoonLandTrans_sub)

XLIM <- c(min(outlineTrans[,1]),max(outlineTrans[,1]))
YLIM <- c(min(outlineTrans[,2]),max(outlineTrans[,2]))

XLIM_sub <- c(min(outlineTrans_sub[,1]),max(outlineTrans_sub[,1]))
YLIM_sub <- c(min(outlineTrans_sub[,2]),max(outlineTrans_sub[,2]))

##
lineList <- list()
for(e in 1:length(lines)){
  lineList[[e]] <- read.table(lines[e], header = FALSE)
}

cartoonLinesTrans <- list()
cartoonLinesTrans_sub <- list()
for(e in 1:length(lines)){
  cartoonLinesTrans[[e]]     <- Morpho::applyTransform(as.matrix(lineList[[e]]), cartoonLandTrans)
  cartoonLinesTrans_sub[[e]] <- Morpho::applyTransform(as.matrix(lineList[[e]]), cartoonLandTrans_sub)
}

############################
# Plot mean landmarks
############################
pdf('Landmark_locations.pdf',width=15,height=10)
par(mar=c(0,2,0,0), oma=c(0,1,0,0), pty='m')
plot(NULL, xlim=XLIM, ylim=YLIM, asp = 1, xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
polygon(outlineTrans, col='white', border='gray60', xlim = XLIM, ylim= YLIM, asp=1)
for(e in 1:length(lineList)){
  lines(cartoonLinesTrans[[e]], col='gray60', xlim = XLIM, ylim= YLIM, asp=1)
}
par(new=T)
plot(transformed$mshape[,1], transformed$mshape[,2], pch=c(19,1,1,1,1,1,1,1,1,1,19,19,19,19,19,19,19,19), cex=2, xlim=XLIM, ylim=YLIM, asp=1, col = adjustcolor('red',alpha=0.3), xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
text(transformed$mshape[,1], transformed$mshape[,2], labels = c(1:18), cex = 2, pos = 2)
dev.off()

############################
# Plot points all 18 landmarks
############################
# pdf('Landmark_all18.pdf',width=15,height=10)
png('Landmark_all18.png',width=1500,height=1000)
par(mar=c(0,2,0,0), oma=c(0,1,0,0), pty='m')
plot(NULL, xlim=XLIM, ylim=YLIM, asp = 1, xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
polygon(outlineTrans, col='gray30', border='gray', xlim = XLIM, ylim= YLIM, asp=1)
for(e in 1:length(lineList)){
  lines(cartoonLinesTrans[[e]], col='gray60', xlim = XLIM, ylim= YLIM, asp=1)
}

for(e in 1:281){
  par(new=T)
  plot(transformed$rotated[,,e][,1], transformed$rotated[,,e][,2], pch=19, cex=1, xlim=XLIM, ylim=YLIM, asp=1, 
       col = adjustcolor(colVec[e], alpha=1), xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
}
dev.off()

############################
# Plot points subset landmarks
############################
# pdf('Landmark_allsubset.pdf',width=15,height=10)
png('Landmark_allsubset.png',width=1500,height=1000)
par(mar=c(0,2,0,0), oma=c(0,1,0,0), pty='m')
plot(NULL, xlim=XLIM_sub, ylim=YLIM_sub, asp = 1, xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
polygon(outlineTrans_sub, col='gray30', border='gray', xlim = XLIM_sub, ylim= YLIM_sub, asp=1)
for(e in 1:length(lineList)){
  lines(cartoonLinesTrans_sub[[e]], col='gray60', xlim = XLIM_sub, ylim= YLIM_sub, asp=1)
}

for(e in 1:281){
  par(new=T)
  plot(transformed_sub$rotated[,,e][,1], transformed_sub$rotated[,,e][,2], pch=19, cex=1, xlim=XLIM_sub, ylim=YLIM_sub, asp=1, 
       col = adjustcolor(colVec[e], alpha=1), xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
}
dev.off()

############################
# plot PCA
############################

# All 18 landmarks
par(mar=c(4,2,1,1), oma=c(1,1,1,1), pty='s')
pdf('Landmark_PCA_all18.pdf',width=10,height=10)
plot(transformed$PCscores[,1], transformed$PCscores[,2], pch = 19, col = colVec, 
     xlab = paste('PC',' (', round((transformed$eigenvalues[1]/sum(transformed$eigenvalues))*100, 1), ' %)'), 
     ylab = paste('PC',' (', round((transformed$eigenvalues[2]/sum(transformed$eigenvalues))*100, 1), ' %)'))
dev.off()

# Subset landmarks
pdf('Landmark_PCA_allsubset.pdf',width=10,height=10)
plot(transformed_sub$PCscores[,1], transformed_sub$PCscores[,2], pch = 19, col = colVec, 
     xlab = paste('PC',' (', round((transformed_sub$eigenvalues[1]/sum(transformed_sub$eigenvalues))*100, 1), ' %)'), 
     ylab = paste('PC',' (', round((transformed_sub$eigenvalues[2]/sum(transformed_sub$eigenvalues))*100, 1), ' %)'))
dev.off()

############################
# Calculate mean transformation for erato and melpomene
############################

# All 18 landmarks
transformed_mean_E <- transformed$rotated[,,1]
transformed_mean_M <- transformed$rotated[,,141]
for(i in 1:18){
  for(e in 2:140){
    transformed_mean_E[i,1] <- transformed_mean_E[i,1] + transformed$rotated[,,e][i,1]
    transformed_mean_E[i,2] <- transformed_mean_E[i,2] + transformed$rotated[,,e][i,2]
  }
  for(m in 142:281){
    transformed_mean_M[i,1] <- transformed_mean_M[i,1] + transformed$rotated[,,m][i,1]
    transformed_mean_M[i,2] <- transformed_mean_M[i,2] + transformed$rotated[,,m][i,2]
  }
}

transformed_mean_E <- transformed_mean_E/140
transformed_mean_M <- transformed_mean_M/141

# Subset landmarks
transformed_mean_E_sub <- transformed_sub$rotated[,,1]
transformed_mean_M_sub <- transformed_sub$rotated[,,141]
for(i in 1:9){
  for(e in 2:140){
    transformed_mean_E_sub[i,1] <- transformed_mean_E_sub[i,1] + transformed_sub$rotated[,,e][i,1]
    transformed_mean_E_sub[i,2] <- transformed_mean_E_sub[i,2] + transformed_sub$rotated[,,e][i,2]
  }
  for(m in 142:281){
    transformed_mean_M_sub[i,1] <- transformed_mean_M_sub[i,1] + transformed_sub$rotated[,,m][i,1]
    transformed_mean_M_sub[i,2] <- transformed_mean_M_sub[i,2] + transformed_sub$rotated[,,m][i,2]
  }
}

transformed_mean_E_sub <- transformed_mean_E_sub/140
transformed_mean_M_sub <- transformed_mean_M_sub/141


############################
# plot transformation melpomene to erato
############################

# All 18 landmarks
pdf('Landmark_transformation_melpomene_erato_all18.pdf',width=15,height=10)
par(mar=c(0,2,0,0), oma=c(0,1,0,0), pty='m')
plot(NULL, xlim=XLIM, ylim=YLIM, asp = 1, xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
polygon(outlineTrans, col='white', border='gray60', xlim = XLIM, ylim= YLIM, asp=1)
for(e in 1:length(lineList)){
  lines(cartoonLinesTrans[[e]], col='gray60', xlim = XLIM, ylim= YLIM, asp=1)
}
par(new=T)
plot(transformed_mean_M[,1], transformed_mean_M[,2], pch=19, cex=2, xlim=XLIM, ylim=YLIM, asp=1, col = adjustcolor('blue',alpha=1), xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
par(new=T)
plot(transformed_mean_E[,1], transformed_mean_E[,2], pch=19, cex=2, xlim=XLIM, ylim=YLIM, asp=1, col = adjustcolor('orange',alpha=1), xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
dev.off()

# Subset landmarks
pdf('Landmark_transformation_melpomene_erato_subset.pdf',width=15,height=10)
par(mar=c(0,2,0,0), oma=c(0,1,0,0), pty='m')
plot(NULL, xlim=XLIM_sub, ylim=YLIM_sub, asp = 1, xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
polygon(outlineTrans_sub, col='white', border='gray60', xlim = XLIM, ylim= YLIM, asp=1)
for(e in 1:length(lineList)){
  lines(cartoonLinesTrans_sub[[e]], col='gray60', xlim = XLIM_sub, ylim= YLIM_sub, asp=1)
}
par(new=T)
plot(transformed_mean_M_sub[,1], transformed_mean_M_sub[,2], pch=19, cex=2, xlim=XLIM_sub, ylim=YLIM_sub, asp=1, col = adjustcolor('blue',alpha=1), xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
par(new=T)
plot(transformed_mean_E_sub[,1], transformed_mean_E_sub[,2], pch=19, cex=2, xlim=XLIM_sub, ylim=YLIM_sub, asp=1, col = adjustcolor('orange',alpha=1), xlab='', ylab='', xaxt='n', yaxt='n', axes=F)
dev.off()

############################
# Plot tension grid
############################

# #All 18 landmarks
# plot(NULL, xlim=XLIM, ylim=YLIM, xlab='x', ylab='y')
# polygon(outlineTrans, col='black', border='gray', xlim = XLIM, ylim= YLIM, asp=1)
# for(e in 1:length(lineList)){
# 
#   lines(cartoonLinesTrans[[e]], col='gray60', xlim = XLIM, ylim= YLIM, asp=1)
# 
# }
# deformGrid2d(transformed_mean_M, transformed_mean_E, ngrid=30, pch=19, add=T, col1 = 'royalblue1', col2 = 'orange', xlim = XLIM, ylim= YLIM, asp=1)
# 
# # Subset landmarks
# plot(NULL, xlim=XLIM_sub, ylim=YLIM_sub, xlab='x', ylab='y')
# polygon(outlineTrans_sub, col='black', border='gray', xlim = XLIM, ylim= YLIM, asp=1)
# for(e in 1:length(lineList)){
#   
#   lines(cartoonLinesTrans_sub[[e]], col='gray60', xlim = XLIM, ylim= YLIM, asp=1)
#   
# }
# deformGrid2d(transformed_mean_M_sub, transformed_mean_E_sub, ngrid=30, pch=19, add=T, col1 = 'royalblue1', col2 = 'orange', xlim = XLIM_sub, ylim= YLIM_sub, asp=1)
# 


# # plot tension map erato to mean
# plot(NULL, xlim=XLIM, ylim=YLIM, xlab='x', ylab='y')
# polygon(outlineTrans, col='black', border='gray', xlim = XLIM, ylim= YLIM, asp=1)
# for(e in 1:length(lineList)){
# 
#   lines(cartoonLinesTrans[[e]], col='gray60', xlim = XLIM, ylim= YLIM, asp=1)
# 
# }
# deformGrid2d(transformed_mean_E, transformed$mshape, ngrid=30, pch=19, add=T, col1= 'orange')
# 
# # plot tension map melpomene to mean
# plot(NULL, xlim=XLIM, ylim=YLIM, xlab='x', ylab='y')
# polygon(outlineTrans, col='black', border='gray', xlim = XLIM, ylim= YLIM, asp=1)
# for(e in 1:length(lineList)){
# 
#   lines(cartoonLinesTrans[[e]], col='gray60', xlim = XLIM, ylim= YLIM, asp=1)
# 
# }
# transformed_M$mshape[,1] <- -transformed_M$mshape[,1]
# transformed_M$mshape[,2] <- -transformed_M$mshape[,2]
# deformGrid2d(transformed_mean_M, transformed$mshape, ngrid=30, pch=19, add=T, col1 = 'royalblue1')


############################
# Plot tension map arrows
############################

# all 18 landmarks
pdf('Landmark_tension_melpomene_erato_all18.pdf',width=15,height=10)
# png('Landmark_tension_melpomene_erato_all18.png',width=1500,height=1000)
tps_arr2(transformed_mean_M, transformed_mean_E, palette = col_summer, arr.nb= 500, shp = F, arr.lwd = 4,
         outline = outlineTrans, arr.len = 0.15, maxArr = 0.019)
polygon(outlineTrans, col=NA, border='gray', xlim = XLIM, ylim= YLIM, asp=1)
for(e in 1:length(lineList)){
  lines(cartoonLinesTrans[[e]], col='gray60', xlim = XLIM, ylim= YLIM, asp=1)
}
dev.off()

# all 18 landmarks
pdf('Landmark_tension_melpomene_erato_subset.pdf',width=15,height=10)
# png('Landmark_tension_melpomene_erato_subset.png',width=1500,height=1000)
tps_arr2(transformed_mean_M_sub, transformed_mean_E_sub, palette = col_summer, arr.nb= 500, shp = F, arr.lwd = 4,
         outline = outlineTrans_sub, arr.len = 0.15, maxArr = 0.019)
polygon(outlineTrans_sub, col=NA, border='gray', xlim = XLIM_sub, ylim= YLIM_sub, asp=1)
for(e in 1:length(lineList)){
  lines(cartoonLinesTrans_sub[[e]], col='gray60', xlim = XLIM_sub, ylim= YLIM_sub, asp=1)
}
dev.off()

############################
# Plot tension map colors
############################

# all 18 landmarks
# pdf('Landmark_tension2_melpomene_erato_all18.pdf',width=15,height=10)
png('Landmark_tension2_melpomene_erato_all18.png',width=1500,height=1000)
tps_iso2(transformed_mean_M, transformed_mean_E, grid =F, outline = outlineTrans, iso.nb =50000, 
        legend = F, shp = F, iso.levels = 1000, cont = F, poly = T, xlim = XLIM, ylim= YLIM, maxArr = 0.019)
polygon(outlineTrans, col=NA, border='gray', xlim = XLIM, ylim= YLIM, asp=1)
for(e in 1:length(lineList)){
  lines(cartoonLinesTrans[[e]], col='gray60', xlim = XLIM, ylim= YLIM, asp=1)
}
dev.off()

# subset landmarks
# pdf('Landmark_tension2_melpomene_erato_subset.pdf',width=15,height=10)
png('Landmark_tension2_melpomene_erato_subset.png',width=1500,height=1000)
tps_iso2(transformed_mean_M_sub, transformed_mean_E_sub, grid =F, outline = outlineTrans_sub, iso.nb =50000, 
         legend = F, shp = F, iso.levels = 1000, cont = F, poly = T, xlim = XLIM_sub, ylim= YLIM_sub, maxArr = 0.019)
polygon(outlineTrans_sub, col=NA, border='gray', xlim = XLIM_sub, ylim= YLIM_sub, asp=1)
for(e in 1:length(lineList)){
  lines(cartoonLinesTrans_sub[[e]], col='gray60', xlim = XLIM_sub, ylim= YLIM_sub, asp=1)
}
dev.off()


# plot legend

# Define color vector and transform to RGB colors (should have the same length as number of populations)
colVec <- c("#c851b5", "#5cb248", "#7f66d2", "#b6b444", "#6984c8", 
            "#80da8b", "#4cbcd1", "#cd4e34", "#54a978", "#c54b6c", 
            "#737a32", "#ba71a9", "#c27e5a", "#c851b5", "#000000")

colFunc <- colorRampPalette(col_spring(1000))

colVec <- colFunc(1000)

# # Show color scheme
# plot(0, xlim=c(0,1000), ylim=c(0,1), pch = '', xaxt='n', yaxt = 'n', axes=F)
# for(e in 1:length(colVec)){
#   rect(e-1, 0, e, 1, col = colVec[e], border=colVec[e])
# }
         
############################
# Males versus females
############################       

sexTable <- read.table('landmarks/sex_table.txt', h=T)
head(sexTable) 
head(as.data.frame(IDlist_ALL))
colnames(sexTable) <- c('IDlist_ALL','Sex')

table_sex <- merge(as.data.frame(IDlist_ALL), sexTable, by = 'IDlist_ALL', sort = F)
table_sex$Sex <- gsub('m', 21, table_sex$Sex)
table_sex$Sex <- gsub('f', 19, table_sex$Sex)
table_sex$Sex <- as.numeric(table_sex$Sex)

colVec <- c(rep('orange',length(landList_hydP)),rep('orange',length(landList_hydFG)),rep('orange',length(landList_amal)),
            rep('orange',length(landList_cyr)),rep('orange',length(landList_not)),rep('orange',length(landList_ety)),
            rep('orange',length(landList_fav)),rep('orange',length(landList_ven)),rep('orange',length(landList_era)),
            rep('orange',length(landList_emm)),rep('orange',length(landList_phy)),rep('orange',length(landList_dem)),
            rep('orange',length(landList_lat)),rep('orange',length(landList_mic)),
            
            rep('royalblue1',length(landList_melP)),rep('royalblue1',length(landList_melFG)),rep('royalblue1',length(landList_mer)),
            rep('royalblue1',length(landList_cyt)),rep('royalblue1',length(landList_ple)),rep('royalblue1',length(landList_ecu)),
            rep('royalblue1',length(landList_amar)),rep('royalblue1',length(landList_vul)),rep('royalblue1',length(landList_the)),
            rep('royalblue1',length(landList_mal)),rep('royalblue1',length(landList_nan)),rep('royalblue1',length(landList_ros)),
            rep('royalblue1',length(landList_agl)),rep('royalblue1',length(landList_xen)))

############################
# plot PCA
############################

# All 18 landmarks
par(mar=c(4,2,1,1), oma=c(1,1,1,1), pty='s')
pdf('Landmark_PCA_all18_sex.pdf',width=10,height=10)
plot(transformed$PCscores[,1], transformed$PCscores[,2], col = colVec, pch = table_sex$Sex,
     xlab = paste('PC',' (', round((transformed$eigenvalues[1]/sum(transformed$eigenvalues))*100, 1), ' %)'), 
     ylab = paste('PC',' (', round((transformed$eigenvalues[2]/sum(transformed$eigenvalues))*100, 1), ' %)'))
  dev.off()

pdf('Landmark_PCA_all18_sex_erato.pdf',width=10,height=10)
plot(transformed_E$PCscores[,1], transformed_E$PCscores[,2], col = colVec[1:142], pch = table_sex$Sex[1:142],
     xlab = paste('PC',' (', round((transformed_E$eigenvalues[1]/sum(transformed_E$eigenvalues))*100, 1), ' %)'), 
     ylab = paste('PC',' (', round((transformed_E$eigenvalues[2]/sum(transformed_E$eigenvalues))*100, 1), ' %)'))
dev.off()

pdf('Landmark_PCA_all18_sex_melpomene.pdf',width=10,height=10)
plot(transformed_M$PCscores[,1], transformed_M$PCscores[,2], col = colVec[143:280], pch = table_sex$Sex[143:280],
     xlab = paste('PC',' (', round((transformed_M$eigenvalues[1]/sum(transformed_M$eigenvalues))*100, 1), ' %)'), 
     ylab = paste('PC',' (', round((transformed_M$eigenvalues[2]/sum(transformed_M$eigenvalues))*100, 1), ' %)'))
dev.off()

# Subset landmarks
pdf('Landmark_PCA_allsubset_sex.pdf',width=10,height=10)
plot(transformed_sub$PCscores[,1], transformed_sub$PCscores[,2], col = colVec, pch = table_sex$Sex,
     xlab = paste('PC',' (', round((transformed_sub$eigenvalues[1]/sum(transformed_sub$eigenvalues))*100, 1), ' %)'), 
     ylab = paste('PC',' (', round((transformed_sub$eigenvalues[2]/sum(transformed_sub$eigenvalues))*100, 1), ' %)'))
dev.off()
