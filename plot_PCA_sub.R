library(RColorBrewer)
library(patternize)

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

imageList <- c(imageList_hydP,
               imageList_hydFG,
               imageList_amal,
               imageList_cyr,
               imageList_not,
               imageList_ety,
               imageList_fav,
               imageList_ven,
               imageList_era,
               imageList_emm,
               imageList_phy,
               imageList_dem,
               imageList_lat,
               imageList_mic)

########################################################
### load of rasters
#erato
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

#melpomene
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

########################################################
sexTable <- read.table('landmarks/sex_table.txt', h=T)
sexTable$Image.ID <- gsub('_calibrated', '', sexTable$Image.ID)
head(sexTable) 
########################################################

par(mar=c(4,4,1,1), oma=c(0,0,0,0))
######################################################
##PCA for H. erato 
######################################################
# Make population and color list
popList_era <- list(IDList_hydP,IDList_hydFG, IDList_dem, IDList_ven, IDList_cyr, IDList_lat, IDList_emm,
                    IDList_not, IDList_ety, IDList_fav, IDList_phy, IDList_amal, IDList_era, IDList_mic)

symbolList_era <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1)

TotalList_era <- c(rasterList_hydP_M_sub, rasterList_hydFG_M_sub, rasterList_dem_M_sub, rasterList_ven_M_sub, rasterList_cyr_M_sub, 
                   rasterList_lat_M_sub, rasterList_emm_M_sub, rasterList_not_M_sub, rasterList_ety_M_sub, rasterList_fav_M_sub,
                   rasterList_phy_M_sub, rasterList_amal_M_sub, rasterList_era_M_sub, rasterList_mic_M_sub)

# plot PCA

pcaOut <- patPCA(TotalList_era, popList_era, colbli_palette, symbolList = symbolList_era, plot = TRUE, plotType = 'points', 
                 plotChanges = TRUE, PCx = 1, PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black', 
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)

### Stat analysis
# Run permutation to evaluate significance of PC axes
library(jackstraw)
matr <- t(pcaOut[[1]])
out <- pcaOut[[3]]
plot(out$x[,c(1,2)])
# permutationPA(matr, B=100)

pca_sign <- out$x[,c(1:10)]

summ <- summary(out)
summ$importance[2,c(1:10)]

# Create groups list
popList <- popList_era
pop <- c('hydP', 'hydFG', 'dem', 'ven', 'cyr', 'lat', 'emm', 'not', 'ety', 'fav', 'phy', 'amal', 'era', 'mic')

group <- c()
for(p in 1:length(popList)){
  for(ind in 1:length(popList[[p]])){
    group <- rbind(group, c(pop[p], as.character(sexTable[match(popList[[p]][ind], sexTable$Image.ID),][1,2])))
  }
}

group <- as.data.frame(group)
colnames(group) <- c('pop','sex')

# Run MANOVA
res.man <- manova(pca_sign ~ group$pop*group$sex)
summ1 <- summary(res.man)
summ2 <- summary.aov(res.man)

aov_table <- c()
for(e in 1:length(summ2)){
  aov_table <- rbind(aov_table, c(summ2[[e]]$`F value`[2], summ2[[e]]$`Pr(>F)`[2]))
}

aov_table

# Run LDA and posterior classification
library(MASS)
ldaOut <- lda(x = pca_sign, grouping = as.character(group$pop), CV = TRUE)

class_mat <- ldaOut$posterior >=.5

class_table <- c()
for(e in 1:length(pop)){
  class_vec <- class_mat[, pop[e]]
  class_true <- class_vec[which(as.character(group$pop) %in% pop[e])]
  class_rate <- sum(class_true)/length(class_true)*100
  class_table <- rbind(class_table, c(pop[e], round(class_rate,2)))
}
as.data.frame(class_table)


ldaOut <- lda(x = pca_sign, grouping = as.character(group$sex), CV = TRUE, prior=c(0.5,0.5))

class_mat <- ldaOut$posterior >=.5

sexes <- c('m','f')
class_table <- c()
for(e in 1:length(sexes)){
  class_vec <- class_mat[, sexes[e]]
  class_true <- class_vec[which(as.character(group$sex) %in% sexes[e])]
  class_rate <- sum(class_true)/length(class_true)*100
  class_table <- rbind(class_table, c(sexes[e], round(class_rate,2)))
}
as.data.frame(class_table)

# permutation of radom sample set
mresult <-c()
fresult <- c()
for(x in c(1:100)){
  ldaOut <- lda(x = pca_sign, grouping = sample(sample(c(rep('m',70),rep('f',70)))), CV = TRUE, prior=c(0.5,0.5))
  
  class_mat <- ldaOut$posterior >=.5
  
  sexes <- c('m','f')
  class_table <- c()
  for(e in 1:length(sexes)){
    class_vec <- class_mat[, sexes[e]]
    class_true <- class_vec[which(as.character(group$sex) %in% sexes[e])]
    class_rate <- sum(class_true)/length(class_true)*100
    class_table <- rbind(class_table, c(sexes[e], round(class_rate,2)))
  }
  mresult <- c(mresult, as.numeric(as.character(as.data.frame(class_table)$V2[1])))
  fresult <- c(fresult, as.numeric(as.character(as.data.frame(class_table)$V2[2])))
}
mean(mresult)
mean(fresult)

sd(mresult)
sd(fresult)

########################################################
# PCA for H. melpomene 
########################################################
# Make population and color list
popList_mel <- list(IDList_melP,IDList_melFG, IDList_ros,  IDList_vul, IDList_cyt, IDList_mal, IDList_agl, 
                    IDList_ple, IDList_ecu,  IDList_amar, IDList_nan, IDList_mer, IDList_the, IDList_xen)

symbolList_mel <- c(16,16,16,16,16,16,16,16,16,16,16,16,16,16)


TotalList_mel <- c(rasterList_melP_M_sub, rasterList_melFG_M_sub, rasterList_ros_M_sub, rasterList_vul_M_sub, rasterList_cyt_M_sub,
                  rasterList_mal_M_sub, rasterList_agl_M_sub, rasterList_ple_M_sub, rasterList_ecu_M_sub, 
                  rasterList_amar_M_sub, rasterList_nan_M_sub, rasterList_mer_M_sub, rasterList_the_M_sub, rasterList_xen_M_sub)
                  

# plot PCA  #par(mar=c(6,6,2,2))

#colfunc <- inferno(100)
pcaOut <- patPCA(TotalList_mel, popList_mel, colbli_palette, symbolList = symbolList_mel, plot = TRUE, plotType = 'points', 
                 plotChanges = TRUE, PCx =1 , PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)
 
### Stat analysis
# Run permutation to evaluate significance of PC axes
library(jackstraw)
matr <- t(pcaOut[[1]])
out <- pcaOut[[3]]
plot(out$x[,c(1,2)])
# permutationPA(matr, B=100)

pca_sign <- out$x[,c(1:11)]

summ <- summary(out)
summ$importance[2,c(1:11)]

# Create groups list
popList <- popList_mel
pop <- c('melP', 'melFG', 'ros', 'vul', 'cyt', 'mal', 'agl', 'ple', 'ecu',  'amar', 'nan', 'mer', 'the', 'xen')

group <- c()
for(p in 1:length(popList)){
  for(ind in 1:length(popList[[p]])){
    group <- rbind(group, c(pop[p], as.character(sexTable[match(popList[[p]][ind], sexTable$Image.ID),][1,2])))
  }
}

group <- as.data.frame(group)
colnames(group) <- c('pop','sex')

# Run MANOVA
res.man <- manova(pca_sign ~ group$pop*group$sex)
summ1 <- summary(res.man)
summ2 <- summary.aov(res.man)

aov_table <- c()
for(e in 1:length(summ2)){
  aov_table <- rbind(aov_table, c(summ2[[e]]$`F value`[2], summ2[[e]]$`Pr(>F)`[2]))
}

aov_table

# Run LDA and posterior classification
library(MASS)
ldaOut <- lda(x = pca_sign, grouping = as.character(group$pop), CV = TRUE)

class_mat <- ldaOut$posterior >=.5

class_table <- c()
for(e in 1:length(pop)){
  class_vec <- class_mat[, pop[e]]
  class_true <- class_vec[which(as.character(group$pop) %in% pop[e])]
  class_rate <- sum(class_true)/length(class_true)*100
  class_table <- rbind(class_table, c(pop[e], round(class_rate,2)))
}
as.data.frame(class_table)


ldaOut <- lda(x = pca_sign, grouping = as.character(group$sex), CV = TRUE, prior=c(0.5,0.5))

class_mat <- ldaOut$posterior >=.5

sexes <- c('m','f')
class_table <- c()
for(e in 1:length(sexes)){
  class_vec <- class_mat[, sexes[e]]
  class_true <- class_vec[which(as.character(group$sex) %in% sexes[e])]
  class_rate <- sum(class_true)/length(class_true)*100
  class_table <- rbind(class_table, c(sexes[e], round(class_rate,2)))
}
as.data.frame(class_table)

# permutation of radom sample set
mresult <-c()
fresult <- c()
for(x in c(1:100)){
  ldaOut <- lda(x = pca_sign, grouping = sample(c(rep('m',70),rep('f',71))), CV = TRUE, prior=c(0.5,0.5))
  
  class_mat <- ldaOut$posterior >=.5
  
  sexes <- c('m','f')
  class_table <- c()
  for(e in 1:length(sexes)){
    class_vec <- class_mat[, sexes[e]]
    class_true <- class_vec[which(as.character(group$sex) %in% sexes[e])]
    class_rate <- sum(class_true)/length(class_true)*100
    class_table <- rbind(class_table, c(sexes[e], round(class_rate,2)))
  }
  mresult <- c(mresult, as.numeric(as.character(as.data.frame(class_table)$V2[1])))
  fresult <- c(fresult, as.numeric(as.character(as.data.frame(class_table)$V2[2])))
}
mean(mresult)
mean(fresult)

sd(mresult)
sd(fresult)

########################################################
###PCA for H. erato and H. melpomene 
########################################################
# population and color list
popList_mel_era <- c(popList_era, popList_mel)

colList_mel_era <- c(colbli_palette, colbli_palette)

symbolList_mel_era <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,16,16,16,16,16,16,16,16,16,16,16,16,16,16)

TotalList_mel_era <- c(TotalList_era, TotalList_mel)

png('PCA_subset.png',width=1000,height=1000)
pcaOut <- patPCA(TotalList_mel_era, popList_mel_era, colList_mel_era, symbolList = symbolList_mel_era, plot = TRUE, plotType = 'points',
                 plotChanges = TRUE, PCx = 1, PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004,
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                  cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)
dev.off()
pdf('PCA_subset.pdf',width=10,height=10)
pcaOut <- patPCA(TotalList_mel_era, popList_mel_era, colList_mel_era, symbolList = symbolList_mel_era, plot = TRUE, plotType = 'points',
                 plotChanges = FALSE, PCx = 1, PCy = 2, plotCartoon = FALSE, refShape = 'target', outline = outline_BC0004,
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)
dev.off()

### Stat analysis
# Run permutation to evaluate significance of PC axes
library(jackstraw)
matr <- t(pcaOut[[1]])
out <- pcaOut[[3]]
plot(out$x[,c(1,2)])
# permutationPA(matr, B=100)

pca_sign <- out$x[,c(1:18)]

summ <- summary(out)
summ$importance[2,c(1:18)]

# Create groups list
popList1 <- popList_era
popList2 <- popList_mel

pop1 <- c('hydP', 'hydFG', 'dem', 'ven', 'cyr', 'lat', 'emm', 'not', 'ety', 'fav', 'phy', 'amal', 'era', 'mic')
pop2 <- c('melP', 'melFG', 'ros', 'vul', 'cyt', 'mal', 'agl', 'ple', 'ecu',  'amar', 'nan', 'mer', 'the', 'xen')
pop <- c(pop1, pop2)

group <- c()
for(p in 1:length(popList1)){
  for(ind in 1:length(popList1[[p]])){
    group <- rbind(group, c(pop1[p], "erato"))
  }
}
for(p in 1:length(popList2)){
  for(ind in 1:length(popList2[[p]])){
    group <- rbind(group, c(pop2[p], "melp"))
  }
}

group <- as.data.frame(group)
colnames(group) <- c('pop', 'spec')

# Run MANOVA
dat <- data.frame(pca_sign, group)
res.man <- manova(as.matrix(dat[,c(1:18)]) ~ dat$spec*dat$pop)
summ1 <- summary(res.man)
summ2 <- summary.aov(res.man)

aov_table <- c()
for(e in 1:length(summ2)){
  aov_table <- rbind(aov_table, c(summ2[[e]]$`F value`[1], summ2[[e]]$`Pr(>F)`[1]))
}

aov_table

# Run LDA and posterior classification
library(MASS)
ldaOut <- lda(x = pca_sign, grouping = as.character(group$spec), CV = TRUE)

class_mat <- ldaOut$posterior >=.5

species <- c('erato', 'melp')
class_table <- c()
for(e in 1:length(species)){
  class_vec <- class_mat[, species[e]]
  class_true <- class_vec[which(as.character(group$spec) %in% species[e])]
  class_rate <- sum(class_true)/length(class_true)*100
  class_table <- rbind(class_table, c(species[e], round(class_rate,2)))
}
as.data.frame(class_table)

# ###############################################
# # Postman PCA for H.erato
# ###############################################

######################################################
##PCA for H. erato 
######################################################
# Make population and color list

colbli_palette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#6A3D9A", "#67001F") 


popList_era <- list(IDList_hydP, IDList_hydFG, IDList_dem, IDList_ven, IDList_fav, IDList_phy)

symbolList_era <- c(1,1,1,1,1,1)

TotalList_era <- c(rasterList_hydP_M_sub, rasterList_hydFG_M_sub, rasterList_dem_M_sub, rasterList_ven_M_sub, rasterList_fav_M_sub,
                   rasterList_phy_M_sub)

# plot PCA

pcaOut <- patPCA(TotalList_era, popList_era, colbli_palette, symbolList = symbolList_era, plot = TRUE, plotType = 'points', 
                 plotChanges = TRUE, PCx = 1, PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black', 
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)

########################################################
#PCA for H. melpomene 
########################################################
# Make population and color list
popList_mel <- list(IDList_melP, IDList_melFG, IDList_ros, IDList_vul, IDList_amar, IDList_nan)

symbolList_mel <- c(16,16,16,16,16,16)


TotalList_mel <- c(rasterList_melP_M_sub, rasterList_melFG_M_sub, rasterList_ros_M_sub, rasterList_vul_M_sub, rasterList_amar_M_sub, rasterList_nan_M_sub)


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

symbolList_mel_era <- c(1,1,1,1,1,1,16,16,16,16,16,16)

TotalList_mel_era <- c(TotalList_era, TotalList_mel)

png('PCA_subset_postman.png',width=1000,height=1000)
pcaOut <- patPCA(TotalList_mel_era, popList_mel_era, colList_mel_era, symbolList = symbolList_mel_era, plot = TRUE, plotType = 'points',
                 plotChanges = TRUE, PCx = 1, PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004,
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc, cex=2)
dev.off()
pdf('PCA_subset_postman.pdf',width=10,height=10)
pcaOut <- patPCA(TotalList_mel_era, popList_mel_era, colList_mel_era, symbolList = symbolList_mel_era, plot = TRUE, plotType = 'points',
                 plotChanges = FALSE, PCx = 1, PCy = 2, plotCartoon = FALSE, refShape = 'target', outline = outline_BC0004,
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc, cex=2)
dev.off()

save(pcaOut, file = "aligned_rasterLists/pcaOut_postman_sub.rda")

### Stat analysis
# Run permutation to evaluate significance of PC axes
library(jackstraw)
matr <- t(pcaOut[[1]])
out <- pcaOut[[3]]
plot(out$x[,c(1,2)])
# permutationPA(matr, B=100)

pca_sign <- out$x[,c(1:15)]

summ <- summary(out)
summ$importance[2,c(1:15)]

# Create groups list
popList1 <- popList_era
popList2 <- popList_mel

pop1 <- c('hydP', 'hydFG', 'dem', 'ven', 'fav', 'phy')
pop2 <- c('melP', 'melFG', 'ros', 'vul', 'amar', 'nan')
pop <- c(pop1, pop2)

group <- c()
for(p in 1:length(popList1)){
  for(ind in 1:length(popList1[[p]])){
    group <- rbind(group, c(pop1[p], "erato"))
  }
}
for(p in 1:length(popList2)){
  for(ind in 1:length(popList2[[p]])){
    group <- rbind(group, c(pop2[p], "melp"))
  }
}

group <- as.data.frame(group)
colnames(group) <- c('pop', 'spec')

# Run MANOVA
dat <- data.frame(pca_sign, group)
res.man <- manova(as.matrix(dat[,c(1:15)]) ~ dat$spec*dat$pop)
summ1 <- summary(res.man)
summ2 <- summary.aov(res.man)

aov_table <- c()
for(e in 1:length(summ2)){
  aov_table <- rbind(aov_table, c(summ2[[e]]$`F value`[1], summ2[[e]]$`Pr(>F)`[1]))
}

aov_table

# Run LDA and posterior classification
library(MASS)
ldaOut <- lda(x = pca_sign, grouping = as.character(group$spec), CV = TRUE)

class_mat <- ldaOut$posterior >=.5

species <- c('erato', 'melp')
class_table <- c()
for(e in 1:length(species)){
  class_vec <- class_mat[, species[e]]
  class_true <- class_vec[which(as.character(group$spec) %in% species[e])]
  class_rate <- sum(class_true)/length(class_true)*100
  class_table <- rbind(class_table, c(species[e], round(class_rate,2)))
}
as.data.frame(class_table)




########################################################
###PCA for H. erato and H. melpomene - non postman
########################################################
# population and color list

popList_era <- list(IDList_cyr, IDList_lat, IDList_emm,
                    IDList_not, IDList_ety, IDList_amal, IDList_era, IDList_mic)

popList_mel <- list(IDList_cyt, IDList_mal, IDList_agl, 
                    IDList_ple, IDList_ecu, IDList_mer, IDList_the, IDList_xen)

popList_mel_era <- c(popList_era, popList_mel)

colList_mel_era <- c(colbli_palette, colbli_palette)

symbolList_mel_era <- c(1,1,1,1,1,1,1,1,16,16,16,16,16,16,16,16)

TotalList_era <- c(rasterList_cyr_M_sub, 
                   rasterList_lat_M_sub, rasterList_emm_M_sub, rasterList_not_M_sub, rasterList_ety_M_sub, 
                   rasterList_amal_M_sub, rasterList_era_M_sub, rasterList_mic_M_sub)

TotalList_mel <- c(rasterList_cyt_M_sub,
                   rasterList_mal_M_sub, rasterList_agl_M_sub, rasterList_ple_M_sub, rasterList_ecu_M_sub, 
                   rasterList_mer_M_sub, rasterList_the_M_sub, rasterList_xen_M_sub)



TotalList_mel_era <- c(TotalList_era, TotalList_mel)

png('PCA_subset.png',width=1000,height=1000)
pcaOut <- patPCA(TotalList_mel_era, popList_mel_era, colList_mel_era, symbolList = symbolList_mel_era, plot = TRUE, plotType = 'points',
                 plotChanges = TRUE, PCx = 1, PCy = 2, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004,
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)
dev.off()
pdf('PCA_subset.pdf',width=10,height=10)
pcaOut <- patPCA(TotalList_mel_era, popList_mel_era, colList_mel_era, symbolList = symbolList_mel_era, plot = TRUE, plotType = 'points',
                 plotChanges = FALSE, PCx = 1, PCy = 2, plotCartoon = FALSE, refShape = 'target', outline = outline_BC0004,
                 imageList = imageList, cartoonID = 'BC0004', normalized = TRUE, cartoonFill = 'black',
                 cartoonOrder = 'under', legendTitle = 'Predicted', colpalette = colfunc)
dev.off()

### Stat analysis
# Run permutation to evaluate significance of PC axes
library(jackstraw)
matr <- t(pcaOut[[1]])
out <- pcaOut[[3]]
plot(out$x[,c(1,2)])
permutationPA(matr, B=100)

pca_sign <- out$x[,c(1:13)]

summ <- summary(out)
summ$importance[2,c(1:13)]

# Create groups list
popList1 <- popList_era
popList2 <- popList_mel

pop1 <- c('cyr', 'lat', 'emm', 'not', 'ety', 'amal', 'era', 'mic')
pop2 <- c('cyt', 'mal', 'agl', 'ple', 'ecu', 'mer', 'the', 'xen')
pop <- c(pop1, pop2)

group <- c()
for(p in 1:length(popList1)){
  for(ind in 1:length(popList1[[p]])){
    group <- rbind(group, c(pop1[p], "erato"))
  }
}
for(p in 1:length(popList2)){
  for(ind in 1:length(popList2[[p]])){
    group <- rbind(group, c(pop2[p], "melp"))
  }
}

group <- as.data.frame(group)
colnames(group) <- c('pop', 'spec')

# Run MANOVA
dat <- data.frame(pca_sign, group)
res.man <- manova(as.matrix(dat[,c(1:18)]) ~ dat$spec*dat$pop)
summ1 <- summary(res.man)
summ2 <- summary.aov(res.man)

aov_table <- c()
for(e in 1:length(summ2)){
  aov_table <- rbind(aov_table, c(summ2[[e]]$`F value`[1], summ2[[e]]$`Pr(>F)`[1]))
}

aov_table

# Run LDA and posterior classification
library(MASS)
ldaOut <- lda(x = pca_sign, grouping = as.character(group$spec), CV = TRUE)

class_mat <- ldaOut$posterior >=.5

species <- c('erato', 'melp')
class_table <- c()
for(e in 1:length(species)){
  class_vec <- class_mat[, species[e]]
  class_true <- class_vec[which(as.character(group$spec) %in% species[e])]
  class_rate <- sum(class_true)/length(class_true)*100
  class_table <- rbind(class_table, c(species[e], round(class_rate,2)))
}
as.data.frame(class_table)

