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


#load
load('aligned_rasterLists/rasterList_era_M.rda')
load('aligned_rasterLists/rasterList_hydFG_M.rda')
load('aligned_rasterLists/rasterList_not_M.rda')
load('aligned_rasterLists/rasterList_lat_M.rda')
load('aligned_rasterLists/rasterList_emm_M.rda')
load('aligned_rasterLists/rasterList_ety_M.rda')
load('aligned_rasterLists/rasterList_fav_M.rda')
load('aligned_rasterLists/rasterList_hydP_M.rda')
load('aligned_rasterLists/rasterList_dem_M.rda')
load('aligned_rasterLists/rasterList_phy_M.rda')
load('aligned_rasterLists/rasterList_cyr_M.rda')
load('aligned_rasterLists/rasterList_ven_M.rda')
load('aligned_rasterLists/rasterList_amal_M.rda')
load('aligned_rasterLists/rasterList_mic_M.rda')

load('aligned_rasterLists/rasterList_mal_M.rda')
load('aligned_rasterLists/rasterList_ecu_M.rda')
load('aligned_rasterLists/rasterList_amar_M.rda')
load('aligned_rasterLists/rasterList_cyt_M.rda')
load('aligned_rasterLists/rasterList_melFG_M.rda')
load('aligned_rasterLists/rasterList_melP_M.rda')
load('aligned_rasterLists/rasterList_mer_M.rda')
load('aligned_rasterLists/rasterList_agl_M.rda')
load('aligned_rasterLists/rasterList_nan_M.rda')
load('aligned_rasterLists/rasterList_ple_M.rda')
load('aligned_rasterLists/rasterList_ros_M.rda')
load('aligned_rasterLists/rasterList_vul_M.rda')
load('aligned_rasterLists/rasterList_the_M.rda')
load('aligned_rasterLists/rasterList_xen_M.rda')

IDList_list_E <- list(IDList_hydFG, IDList_hydP, IDList_dem, IDList_ven, IDList_fav,  IDList_phy, IDList_cyr, IDList_lat, IDList_emm, IDList_ety, IDList_not, IDList_mic, IDList_amal, IDList_era)
IDList_list_M <- list(IDList_melFG, IDList_melP, IDList_ros, IDList_vul, IDList_amar, IDList_nan, IDList_cyt, IDList_mal, IDList_agl, IDList_ecu, IDList_ple, IDList_xen, IDList_mer, IDList_the)


landList_list_E <- list(landList_hydFG, landList_hydP, landList_dem, landList_ven, landList_fav,  landList_phy, landList_cyr, landList_lat, landList_emm, landList_ety, landList_not, landList_mic, landList_amal, landList_era)
landList_list_M <- list(landList_melFG, landList_melP, landList_ros, landList_vul, landList_amar, landList_nan, landList_cyt, landList_mal, landList_agl, landList_ecu, landList_ple, landList_xen, landList_mer, landList_the)


imageList_list_E <- list(imageList_hydFG, imageList_hydP, imageList_dem, imageList_ven, imageList_fav,  imageList_phy, imageList_cyr, imageList_lat, imageList_emm, imageList_ety, imageList_not, imageList_mic, imageList_amal, imageList_era)
imageList_list_M <- list(imageList_melFG, imageList_melP, imageList_ros, imageList_vul, imageList_amar, imageList_nan, imageList_cyt, imageList_mal, imageList_agl, imageList_ecu, imageList_ple, imageList_xen, imageList_mer, imageList_the)

rasterList_list_E <- list(rasterList_hydFG_M, rasterList_hydP_M, rasterList_dem_M, rasterList_ven_M, rasterList_fav_M,  rasterList_phy_M, rasterList_cyr_M, rasterList_lat_M, rasterList_emm_M, rasterList_ety_M, rasterList_not_M, rasterList_mic_M, rasterList_amal_M, rasterList_era_M)
rasterList_list_M <- list(rasterList_melFG_M, rasterList_melP_M, rasterList_ros_M, rasterList_vul_M, rasterList_amar_M, rasterList_nan_M, rasterList_cyt_M, rasterList_mal_M, rasterList_agl_M, rasterList_ecu_M, rasterList_ple_M, rasterList_xen_M, rasterList_mer_M, rasterList_the_M)


#figure creation 
png('comparison_noNames.png', width=4000, height=11200)
layout(matrix(c(1:42), nrow=14, byrow=TRUE))
# layout.show(n=21)
par(mar=c(0,0,0,0), oma=c(0,3,0,3))


for(e in 1:length(IDList_list_E)){
  
  summedRaster_E_M_m <- sumRaster(rasterList_list_E[[e]], IDList_list_E[[e]], type = 'RGB')
  
  summedRaster_M_M_m <- sumRaster(rasterList_list_M[[e]], IDList_list_M[[e]], type = 'RGB')
  
  colfunc <- inferno(100)
  plotHeat(summedRaster_E_M_m, IDList_list_E[[e]], plotCartoon = TRUE, refShape = 'target', 
           outline = outline_BC0004, lines = lines_BC0004, landList = landList_list_E[[e]], adjustCoords = TRUE, 
           imageList = imageList_list_E[[e]], cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',
           colpalette = colfunc, legend = F)
  # text(1787.333,1998, substitute(paste(italic(nn), ' French Guyana'), list(nn='H. e. hydara')), cex=1, adj=0)
  
  plotHeat(summedRaster_M_M_m, IDList_list_M[[e]], plotCartoon = TRUE, refShape = 'target', 
           outline = outline_BC0004, lines = lines_BC0004, landList = landList_list_M[[e]], adjustCoords = TRUE, 
           imageList = imageList_list_M[[e]], cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under',
           colpalette = colfunc, legend = F)
  # text(1787.333,1998, substitute(paste(italic(nn), ' French Guyana'), list(nn='H. m. melpomene')), cex=1, adj=0)
  
  colfunc <- c("blue","lightblue","black","burlywood1","orange")
  raster_diff <- summedRaster_E_M_m/length(IDList_list_E[[e]]) - summedRaster_M_M_m/length(IDList_list_M[[e]])
  plotHeat(raster_diff, IDList_list_E[[e]], plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
           lines = lines_BC0004, landList = landList_list_E[[e]], adjustCoords = TRUE, imageList = imageList_list_E[[e]], 
           cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', colpalette = colfunc, normalized = T, 
           zlim = c(-1,1), legendTitle = "", legend = F)
}
dev.off()