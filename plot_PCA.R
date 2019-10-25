library(RColorBrewer)

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

########################################################
### load of rasters
#erato
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

#melpomene
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

########################################################
#color palette
colbli_palette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", 
                    "#FF7F00", "#1B9E77", "#6A3D9A", "#67001F","#B15928", "#F0027F","#000000") 
#colfunc <- inferno(100)
colf <- colorRampPalette(c("blue","lightblue","black","indianred1","firebrick1"))
colfunc <- colf(100)

## target image
target <- as.matrix(read.table('cartoon/BC0004_landmarks_LFW.txt',h = F))

# cartoon
outline_BC0004 <- read.table('cartoon/BC0004_outline.txt', h = F)
lines_BC0004 <- list.files(path ='cartoon', pattern ='BC0004_vein', full.names = T)


######################################################
##PCA for H. erato 
######################################################
# Make population and color list
popList_era <- list(IDList_hydP,IDList_hydFG, IDList_dem, IDList_ven, IDList_cyr, IDList_lat, IDList_emm,
                    IDList_not,IDList_ety, IDList_fav, IDList_phy, IDList_amal, IDList_era, IDList_mic)

symbolList_era <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1)

TotalList_era <- c(rasterList_hydP_M, rasterList_hydFG_M, rasterList_dem_M, rasterList_ven_M, rasterList_cyr_M, 
                   rasterList_lat_M, rasterList_emm_M, rasterList_not_M, rasterList_ety_M, rasterList_fav_M,
                   rasterList_phy_M, rasterList_amal_M, rasterList_era_M, rasterList_mic_M)

# plot PCA

pcaOut <- patPCA(TotalList_era, popList_era, colbli_palette, symbolList = symbolList_era, plot = TRUE, plotType = 'points', 
                 plotChanges = TRUE, PCx = 1, PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black', 
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)


########################################################
#PCA for H. melpomene 
########################################################
# Make population and color list
popList_mel <- list(IDList_melP,IDList_melFG, IDList_ros, IDList_vul, IDList_cyt, IDList_mal, IDList_agl, 
                    IDList_ple, IDList_ecu, IDList_nan, IDList_amar, IDList_mer, IDList_the, IDList_xen)

symbolList_mel <- c(16,16,16,16,16,16,16,16,16,16,16,16,16,16)


TotalList_mel <- c(rasterList_melP_M, rasterList_melFG_M, rasterList_ros_M, rasterList_vul_M, rasterList_cyt_M,
                  rasterList_mal_M, rasterList_agl_M, rasterList_ple_M, rasterList_ecu_M, rasterList_nan_M, 
                  rasterList_amar_M, rasterList_mer_M, rasterList_the_M, rasterList_xen_M)
                  

# plot PCA  #par(mar=c(6,6,2,2))

#colfunc <- inferno(100)
pcaOut <- patPCA(TotalList_mel, popList_mel, colbli_palette, symbolList = symbolList_mel, plot = TRUE, plotType = 'points', 
                 plotChanges = TRUE, PCx =1 , PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)
 

########################################################
###PCA for H. erato and H. melpomene 
########################################################
# population and color list
popList_mel_era <- c(popList_era, popList_mel)

colList_mel_era <- c(colbli_palette, colbli_palette)

symbolList_mel_era <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,16,16,16,16,16,16,16,16,16,16,16,16,16,16)

TotalList_mel_era <- c(TotalList_era, TotalList_mel)

png('PCA_all18.png',width=1000,height=1000)
pcaOut <- patPCA(TotalList_mel_era, popList_mel_era, colList_mel_era, symbolList = symbolList_mel_era, plot = TRUE, plotType = 'points',
                 plotChanges = TRUE, PCx = 1, PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004,
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                  cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)
dev.off()
pdf('PCA_all18.pdf',width=10,height=10)
pcaOut <- patPCA(TotalList_mel_era, popList_mel_era, colList_mel_era, symbolList = symbolList_mel_era, plot = TRUE, plotType = 'points',
                 plotChanges = FALSE, PCx = 1, PCy = 2, plotCartoon = FALSE, refShape = 'target', outline = outline_BC0004,
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)
dev.off()

# ###############################################
# # Postman PCA for H.erato
# ###############################################

postman_palette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#6A3D9A", "#67001F") 



######################################################
## Postmans PCA for H. erato 
######################################################
# Make population and color list
popPostmanList_era <- list(IDList_hydP,IDList_hydFG, IDList_dem, IDList_ven, 
                           IDList_fav, IDList_phy)

symbolPostmanList_era <- c(1,1,1,1,1,1)

TotalPostmanList_era <- c(rasterList_hydP_M, rasterList_hydFG_M, rasterList_dem_M, rasterList_ven_M,  
                          rasterList_fav_M,rasterList_phy_M)

# plot PCA
pcaOut <- patPCA(TotalPostmanList_era, popPostmanList_era, postman_palette, symbolList = symbolPostmanList_era, plot = TRUE, plotType = 'points', 
                 plotChanges = TRUE, PCx = 1, PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black', 
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)



########################################################
# Postmans PCA for H. melpomene 
########################################################
# Make population and color list
popPostmanList_mel <- list(IDList_melP,IDList_melFG, IDList_ros,  IDList_vul, 
                           IDList_amar, IDList_nan)

symbolPostmanList_mel <- c(16,16,16,16,16,16)


TotalPostmanList_mel <- c(rasterList_melP_M, rasterList_melFG_M, rasterList_ros_M, rasterList_vul_M,
                          rasterList_amar_M, rasterList_nan_M)


# plot PCA  #par(mar=c(6,6,2,2))

#colfunc <- inferno(100)
pcaOut <- patPCA(TotalPostmanList_mel, popPostmanList_mel, postman_palette, symbolList = symbolPostmanList_mel, plot = TRUE, plotType = 'points', 
                 plotChanges = TRUE, PCx =1 , PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)


########################################################
###Postman PCA for H. erato and H. melpomene 
########################################################
# population and color list
popPostmanList_mel_era <- c(popPostmanList_era, popPostmanList_mel)

PostmanPalette_mel_era <- c(postman_palette, postman_palette)

symbolList_mel_era <- c(symbolPostmanList_era, symbolPostmanList_mel)

TotalPostmanList_mel_era <- c(TotalPostmanList_era, TotalPostmanList_mel)

png('PCA_all18_postman.png',width=1000,height=1000)
pcaOut <- patPCA(TotalPostmanList_mel_era, popPostmanList_mel_era, PostmanPalette_mel_era, symbolList = symbolList_mel_era, plot = TRUE, plotType = 'points',
                 plotChanges = TRUE, PCx = 1, PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004,
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc, cex=3)
dev.off()
pdf('PCA_all18_postman.pdf',width=10,height=10)
pcaOut <- patPCA(TotalPostmanList_mel_era, popPostmanList_mel_era, PostmanPalette_mel_era, symbolList = symbolList_mel_era, plot = TRUE, plotType = 'points',
                 plotChanges = FALSE, PCx = 1, PCy = 2, plotCartoon = FALSE, refShape = 'target', outline = outline_BC0004,
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc, cex=3)
dev.off()

save(pcaOut, file = 'aligned_rasterLists/pca_postman.rda')

# ###############################################
# ###############################################