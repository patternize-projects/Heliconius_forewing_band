##############################################################
# load dependencies
##############################################################
library(ggplot2)
library(maps)
library(mapdata)
library(plyr)
library(alphahull)
library(G1DBN)
library(geosphere)
library(RColorBrewer)

##############################################################
# distribution function
##############################################################
distribution <- function(data,
                         info,
                         title = NULL,
                         xlim = NULL,
                         ylim = NULL,
                         axes = FALSE,
                         border = "gray",
                         bg = "white",
                         polygon = TRUE,
                         points = TRUE, 
                         alphaFill = 0.5,
                         alphaPoints = 0.5,
                         alphaShape = 1,
                         runoption = 'A'){
  
  map("worldHires", xlim = xlim, ylim = ylim, col="gray70", border = border, fill = TRUE, bg = bg, lforce="e")
  
  title(main = title, line = 1)
  mtext('longitude', side = 1, line = 2.5, cex = 0.8)
  mtext('latitude', side = 2, line = 2.5, cex = 0.8)
  
  if(axes){map.axes()}
  
  e <- 1
  alphaS <- alphaShape
  areas <- c()
  
  while(e <= nrow(info)){
    
    Loc2 <- merge(data, info, by = 'taxon_nombre')
    
    sub <- subset(Loc2, as.numeric(as.character(Loc2$order)) == e)
    
    LocUnique <- unique(sub[,c("lon_dec", "lat_dec")])
    
    colS <- as.vector(col2rgb(unique(sub$color)))/255
    
    
    if(polygon){
      
      n10 <- ashape(LocUnique, alpha = alphaS)
      n10g = graph.edgelist(cbind(as.character(n10$edges[, "ind1"]), as.character(n10$edges[, "ind2"])), 
                            directed = FALSE)
      
      go <- 'NoGo'
      if(runoption == 'A'){
        if(!is.connected(n10g) || (clusters(n10g)$no > 1)) {  # any(degree(n10g) > 2) ||
          
          alphaS <- alphaS + 0.1
          print(paste('increasing alpha with 0.1 for', info$taxon_nombre[e], sep = ' '))
        }
        else{
          go <- 'ok'
        }
      }
      if(runoption == 'B'){
        if(!is.connected(n10g) || any(degree(n10g) > 2) || (clusters(n10g)$no > 1)) {  #
          
          alphaS <- alphaS + 0.1
          print(paste('increasing alpha with 0.1 for', info$taxon_nombre[e], sep = ' '))
        }
        else{
          go <- 'ok'
        }
      }
      
      if(go == 'ok'){
        
        cutg = n10g - E(n10g)[1]
        # find chain end points
        ends = names(which(degree(cutg) == 1))

        if(identical(ends, character(0))){
          alphaS <- alphaS + 0.1
          print(paste('increasing alpha with 0.1 for', info$taxon_nombre[e], sep = ' '))
          next
        }
        
        path <- tryCatch(get.shortest.paths(cutg, ends[1], ends[2])[[1]],
                         warning = function(err) {
                           print(paste('increasing alpha', info$taxon_nombre[e], sep = ' '))
                           return(NULL)
                         })
        
        if(is.null(path)){
          alphaS <- alphaS + 0.1
          next
        }
        
        # this is an index into the points
        pathX = as.numeric(V(n10g)[path[[1]]]$name)
        # join the ends
        pathX = c(pathX, pathX[1])
        
        ashapem <- as.matrix(n10$x[pathX, ])
        
        ashapem[,2][ashapem[,2] < ylim[1]] <- ylim[1]
        ashapem[,2][ashapem[,2] > ylim[2]] <- ylim[2]
        ashapem[,1][ashapem[,1] < xlim[1]] <- xlim[1]
        ashapem[,1][ashapem[,1] > xlim[2]] <- xlim[2]
        
        polygon(ashapem, col = rgb(colS[1], colS[2], colS[3], alpha = alphaFill), border = NA)

        areas <- rbind(areas, c(as.character(info$taxon_nombre[e]), areaPolygon(ashapem)/1000000))
        
        if(points){
          
          points(x = LocUnique$lon_dec, y = LocUnique$lat_dec, col = rgb(colS[1], colS[2], colS[3], alpha = alphaPoints), 
                 pch = 19, cex = 0.5)
        }
        
        e <- e + 1
        alphaS <- alphaShape
      }
    }
    
    if(points == TRUE && polygon == FALSE){
      
      points(x = LocUnique$lon_dec, y = LocUnique$lat_dec, col = rgb(colS[1], colS[2], colS[3], alpha = alphaPoints), 
             pch = 17, cex = 1)
      
      e <- e + 1
    }
    
  }
  colnames(areas) <- c('taxon_nombre','area_km2')
  return(as.data.frame(areas))
}



##############################################################
# Run example data
##############################################################
layout(matrix(c(1:2), nrow=1, byrow=TRUE))
layout.show(n=2)

par(mar=c(1,1,1,1), oma=c(1,2,1,1))

###
# Heliconius erato
###

# Load locality files
# files <- list.files(path='Localities/H_erato', pattern='.csv', full.names = T)

files <- c('H_e_hydara','H_e_erato','H_e_demophoon', 'H_e_venus','H_e_cyrbia', 'H_e_lativitta','H_e_emma',
           'H_e_notabilis','H_e_etylus','H_e_favorinus','H_e_phyllis','H_e_amalfreda', 'H_e_microclea')

# Combine files
Loc <- c()
for(f in 1:length(files)){
  table <- read.csv(paste('Localities/H_erato/',files[f],'.csv',sep=''), h=T)
  Loc <- rbind(Loc,table)
}

# Remove NA records
Loc <- Loc[!is.na(Loc$lon_dec),]
Loc <- Loc[!is.na(Loc$lat_dec),]

# Extract species/race names
LocSpec <- unique(Loc[,c("taxon_nombre")])

# Define color vector and transform to RGB colors (should have the same length as number of populations)
#colVec <- c("#c851b5", "#5cb248", "#7f66d2", "#b6b444", "#6984c8", 
 #           "#80da8b", "#4cbcd1", "#cd4e34", "#54a978", "#c54b6c", 
  #          "#737a32")#, "#ba71a9", "#c27e5a", "#c851b5", "#000000")
colVec <- c("#1F78B4","#F0027F", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
           "#FDBF6F", "#FF7F00",  "#1B9E77", "#6A3D9A", "#67001F",
            "#B15928",'#000000')#,"#003C30")#first is blue for hydP"#A6CEE3", )"#FFFF99"
# # Show color scheme
# plot(0, xlim=c(0,15), ylim=c(0,1), pch = '')
# for(e in 1:length(colVec)){
#   rect(e-1, 0, e, 1, col = colVec[e])
# }

# Define order for plotting (should have the same length as number of populations)
#order <- c(5,3,7,13,2,1,4,6,12,11,10,9,8,14,15)
order <- c(1,12,2,3,4,5,6,7,8,9,10,11,13)#,12,13,14)

# Combine colors and species/races
Spec <- as.data.frame(cbind(as.character(LocSpec), colVec, order))
colnames(Spec) <- c('taxon_nombre', 'color', 'order')

# Plot distribution
distribution(data = Loc, info = Spec, axes = TRUE,
             title = substitute(paste(italic("Heliconius erato "), "", sep=" ")),
             xlim = c(-95,-30), ylim = c(-40, 20), 
             polygon = TRUE, points = TRUE, alphaFill = 0.6, alphaShape = 3, runoption = 'A')


###
# Heliconius melpomene
###

# Load locality files
files <- list.files(path='Localities/H_melpomene', pattern='.csv', full.names = T)

files <- c('H_m_melpomene','H_m_thelxiopeia','H_m_rosina','H_m_vulcanus', 'H_m_cythera','H_m_malleti','H_m_aglaope',
           'H_m_plesseni','H_m_ecuadoriensis','H_m_amaryllis', 'H_m_nanna','H_m_meriana','H_m_xenoclea')

# Combine files
Loc <- c()
for(f in 1:length(files)){
  table <- read.csv(paste('Localities/H_melpomene/',files[f],'.csv',sep=''), h=T)
  Loc <-rbind(Loc,table)
}

# Remove NA records
Loc <- Loc[!is.na(Loc$lon_dec),]
Loc <- Loc[!is.na(Loc$lat_dec),]

# Extract species/race names
LocSpec <- unique(Loc[,c("taxon_nombre")])

# Define color vector and transform to RGB colors (should have the same length as number of populations)
#colVec <- c("#c851b5", "#5cb248", "#7f66d2", "#b6b444", "#6984c8", 
 #           "#80da8b", "#4cbcd1", "#cd4e34", "#54a978", "#c54b6c", 
  #          "#737a32")#, "#ba71a9", "#c27e5a", "#c851b5", "#000000")
# '''
# colVec <- c("#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
#             "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99","#B15928")#first is blue for melP"#A6CEE3", )
# '''
# # Show color scheme
# # Show color scheme
# plot(0, xlim=c(0,15), ylim=c(0,1), pch = '')
# for(e in 1:length(colVec)){
#   rect(e-1, 0, e, 1, col = colVec[e])
# }

# Define order for plotting (should have the same length as number of populations)
#order <- c(5,3,7,13,2,1,4,6,12,11,10,9,8,14,15)
order <- c(1,12,2,3,4,5,6,7,8,9,10,11,13)#,12,13,14)

# Combine colors and species/races
Spec <- as.data.frame(cbind(as.character(LocSpec), colVec, order))
colnames(Spec) <- c('taxon_nombre', 'color', 'order')

# Plot distribution
distribution(data = Loc, info = Spec, axes = TRUE,
             title = substitute(paste(italic("Heliconius melpomene "), "", sep=" ")),
             xlim = c(-95,-30), ylim = c(-40, 20), 
             polygon = TRUE, points = TRUE, alphaFill = 0.6, alphaShape = 3, runoption = 'B')


