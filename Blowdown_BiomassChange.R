##### Read in lidar data #####

  # 5 m resolution canopy height rasters from 2009 and 2019
    lidarCHM09 <- raster::raster("CHM09_AOI_res500.tif")
    lidarCHM19 <- raster::raster("CHM19_AOI_res500.tif")

  # Plot canopy height from each year  
    pal <- colorRampPalette(c("blue","lightyellow","red"))
    htBrks <- seq(-10,70,5)
    par(mar=c(1,1,1,1),oma=c(1,1,1,1))
    
    raster::plot(lidarCHM09, breaks=htBrks, col=pal(length(htBrks)))
    raster::plot(lidarCHM19, breaks=htBrks, col=pal(length(htBrks)))

##### Mask based on land use types #####
    
  # Read in land use class shapefile  
    LaSelvaClass <- rgdal::readOGR("LandUseShapefile/LU2000.shp")
  
  # Make a polygon of all forested area
    # Define classes to keep
      ForestClasses <- c("Old-growth Forests","Ecological Reserve","Forested Swamps","Secondary Forests","Abandoned Agroforestry","Abandoned Plantation", "Selectively-logged Forests")
    # Subset forest classes  
      ForestPoly <- LaSelvaClass[LaSelvaClass$DESCRIPTIO %in% ForestClasses,]
    # Remove forest area added in or after 1990
      ForestPoly <- ForestPoly[!(ForestPoly$YEAR >= 1990),]
    # Remove tiny isolated fragment of abandoned plantation
      ForestPoly <- ForestPoly[!(ForestPoly$LU00_UTM0_==68),]

  # Make polygon of old growth forest area
    OldGrowthPoly <- ForestPoly[ForestPoly$DESCRIPTIO %in% c("Old-growth Forests","Ecological Reserve","Forested Swamps"),]
    
  # Make polygon of secondary (and previously disturbed) forest area
    SecondaryPoly <- ForestPoly[ForestPoly$DESCRIPTIO %in% c("Secondary Forests","Abandoned Agroforestry","Abandoned Plantation", "Selectively-logged Forests"),]
    
  # Make canopy height raster of all forest area
    CHM09Forest <- raster::mask(lidarCHM09, ForestPoly)
    CHM19Forest <- raster::mask(lidarCHM19, ForestPoly)  
  
  # Make canopy height raster of old growth forest area
    CHM09OldGrowth <- raster::mask(lidarCHM09, OldGrowthPoly)
    CHM19OldGrowth <- raster::mask(lidarCHM19, OldGrowthPoly)
  
  # Make canopy height raster of secondary forest area    
    CHM09Secondary <- raster::mask(lidarCHM09, SecondaryPoly)
    CHM19Secondary <- raster::mask(lidarCHM19, SecondaryPoly)
  
  # Plot forest area
  raster::plot(CHM09Forest, breaks=htBrks, col=pal(length(htBrks)))
  raster::plot(CHM19Forest, breaks=htBrks, col=pal(length(htBrks)))

##### Make aboveground biomass (AGB) rasters with MC uncertainty #####

  # Make an empty raster with correct resolution for AGB model (half hectare)
    AGB_template <- raster::raster(ext=raster::extent(ForestPoly), resolution=(5000^0.5), crs=sp::proj4string(ForestPoly))
  
  # Resample CHM rasters to lower AGB model resolution
    CHM09_low <- raster::resample(x=CHM09Forest, y=AGB_template, method="bilinear")
    CHM19_low <- raster::resample(x=CHM19Forest, y=AGB_template, method="bilinear")
  
  # Predict AGB values using Monte Carlo simulation to get uncertainty
    
    # Make data frame to store AGB results
    AGBresults <- data.frame(n = 1:1000,
                             allForest09 = NA,
                             allForestSD09 = NA,
                             oldGrowth09 = NA,
                             oldGrowthSD09 = NA,
                             secondaryForest09 = NA,
                             secondaryForestSD09 = NA,
                             
                             allForest19 = NA,
                             allForestSD19 = NA,
                             oldGrowth19 = NA,
                             oldGrowthSD19 = NA,
                             secondaryForest19 = NA,
                             secondaryForestSD19 = NA,
                             
                             forestAGBchange = NA,
                             oldGrowthAGBchange = NA,
                             secondaryForestAGBchange = NA)
    
    set.seed(1) # make random sampling repeatable
    for(i in 1:1000){
      AGB09 <- CHM09_low
      AGB09@data@values[!is.na(AGB09@data@values)] <- 0.47*((8.881)*CHM09_low@data@values[!is.na(CHM09_low@data@values)]^(1.021)
                                                       + rnorm(n=length(AGB09@data@values[!is.na(AGB09@data@values)]), mean = 0, sd = 19.174093))
      AGB19 <- CHM19_low
      AGB19@data@values[!is.na(AGB19@data@values)] <-  0.47*((8.881)*CHM19_low@data@values[!is.na(CHM19_low@data@values)]^(1.021)
                                                        + rnorm(n=length(AGB19@data@values[!is.na(AGB19@data@values)]), mean = 0, sd = 19.174093))

      # Verify location of NA pixels is the same
      # which(is.na(AGB09@data@values)!=is.na(AGB19@data@values))
      
      # Find mean and SD of AGB values in all forest
      AGBresults$allForest09[i] <- mean(raster::mask(AGB09,ForestPoly)@data@values,na.rm=T)
      AGBresults$allForestSD09[i] <- sd(raster::mask(AGB09,ForestPoly)@data@values,na.rm=T)
      AGBresults$allForest19[i] <- mean(raster::mask(AGB19,ForestPoly)@data@values,na.rm=T)
      AGBresults$allForestSD19[i] <- sd(raster::mask(AGB19,ForestPoly)@data@values,na.rm=T)
      
      # Find mean and SD of AGB in old-growth forest
      AGBresults$oldGrowth09[i] <- mean(raster::mask(AGB09,OldGrowthPoly)@data@values,na.rm=T)
      AGBresults$oldGrowthSD09[i] <- sd(raster::mask(AGB09,OldGrowthPoly)@data@values,na.rm=T)
      AGBresults$oldGrowth19[i] <- mean(raster::mask(AGB19,OldGrowthPoly)@data@values,na.rm=T)
      AGBresults$oldGrowthSD19[i] <- sd(raster::mask(AGB19,OldGrowthPoly)@data@values,na.rm=T)
      
      # Find mean and SD of AGB in secondary forest
      AGBresults$secondaryForest09[i] <- mean(raster::mask(AGB09,SecondaryPoly)@data@values,na.rm=T)
      AGBresults$secondaryForestSD09[i] <- sd(raster::mask(AGB09,SecondaryPoly)@data@values,na.rm=T)
      AGBresults$secondaryForest19[i] <- mean(raster::mask(AGB19,SecondaryPoly)@data@values,na.rm=T)
      AGBresults$secondaryForestSD19[i] <- sd(raster::mask(AGB19,SecondaryPoly)@data@values,na.rm=T)
      
      # Find the mean absolute biomass change 
      AGBresults$forestAGBchange[i] <- mean(raster::mask(AGB09,ForestPoly)@data@values - raster::mask(AGB19,ForestPoly)@data@values,na.rm=T)/AGBresults$allForest09[i]*100
      AGBresults$oldGrowthAGBchange[i] <- mean(raster::mask(AGB09,OldGrowthPoly)@data@values - raster::mask(AGB19,OldGrowthPoly)@data@values,na.rm=T)/AGBresults$oldGrowth09[i]*100
      AGBresults$secondaryForestAGBchange[i] <- mean(raster::mask(AGB09,SecondaryPoly)@data@values - raster::mask(AGB19,SecondaryPoly)@data@values,na.rm=T)/AGBresults$secondaryForest09[i]*100
    }
    
    round(mean(AGBresults$allForest09),1)
    paste0("[",round(quantile(AGBresults$allForest09,0.025),1),", ",round(quantile(AGBresults$allForest09,0.975),1),"]")
    round(mean(AGBresults$allForest19),1)
    paste0("[",round(quantile(AGBresults$allForest19,0.025),1),", ",round(quantile(AGBresults$allForest19,0.975),1),"]")

    round(mean(AGBresults$oldGrowth09),1)
    paste0("[",round(quantile(AGBresults$oldGrowth09,0.025),1),", ",round(quantile(AGBresults$oldGrowth09,0.975),1),"]")
    round(mean(AGBresults$oldGrowth19),1)
    paste0("[",round(quantile(AGBresults$oldGrowth19,0.025),1),", ",round(quantile(AGBresults$oldGrowth19,0.975),1),"]")
    
    round(mean(AGBresults$secondaryForest09),1)
    paste0("[",round(quantile(AGBresults$secondaryForest09,0.025),1),", ",round(quantile(AGBresults$secondaryForest09,0.975),1),"]")
    round(mean(AGBresults$secondaryForest19),1)
    paste0("[",round(quantile(AGBresults$secondaryForest19,0.025),1),", ",round(quantile(AGBresults$secondaryForest19,0.975),1),"]")
    
    round(mean(AGBresults$allForestSD09),1)
    paste0("[",round(quantile(AGBresults$allForestSD09,0.025),1),", ",round(quantile(AGBresults$allForestSD09,0.975),1),"]")
    round(mean(AGBresults$allForestSD19),1)
    paste0("[",round(quantile(AGBresults$allForestSD19,0.025),1),", ",round(quantile(AGBresults$allForestSD19,0.975),1),"]")

    round(mean(AGBresults$oldGrowthSD09),1)
    paste0("[",round(quantile(AGBresults$oldGrowthSD09,0.025),1),", ",round(quantile(AGBresults$oldGrowthSD09,0.975),1),"]")
    round(mean(AGBresults$oldGrowthSD19),1)
    paste0("[",round(quantile(AGBresults$oldGrowthSD19,0.025),1),", ",round(quantile(AGBresults$oldGrowthSD19,0.975),1),"]")
    
    round(mean(AGBresults$secondaryForestSD09),1)
    paste0("[",round(quantile(AGBresults$secondaryForestSD09,0.025),1),", ",round(quantile(AGBresults$secondaryForestSD09,0.975),1),"]")
    round(mean(AGBresults$secondaryForestSD19),1)
    paste0("[",round(quantile(AGBresults$secondaryForestSD19,0.025),1),", ",round(quantile(AGBresults$secondaryForestSD19,0.975),1),"]")
    
    round(mean(AGBresults$forestAGBchange),2)
    paste0("[",round(quantile(AGBresults$forestAGBchange,0.025),1),", ",round(quantile(AGBresults$forestAGBchange,0.975),1),"]")
    
    round(mean(AGBresults$oldGrowthAGBchange),2)
    paste0("[",round(quantile(AGBresults$oldGrowthAGBchange,0.025),1),", ",round(quantile(AGBresults$oldGrowthAGBchange,0.975),1),"]")
    
    round(mean(AGBresults$secondaryForestAGBchange),2)
    paste0("[",round(quantile(AGBresults$secondaryForestAGBchange,0.025),1),", ",round(quantile(AGBresults$secondaryForestAGBchange,0.975),1),"]")
    
    
  # Calculate area in each forest type
    length(raster::mask(AGB09,ForestPoly)@data@values[!is.na(raster::mask(AGB09,ForestPoly)@data@values)])/2 # Divide by 2 to convert from 0.5 ha pixels to area in ha
    length(raster::mask(AGB09,OldGrowthPoly)@data@values[!is.na(raster::mask(AGB09,OldGrowthPoly)@data@values)])/2
    length(raster::mask(AGB09,SecondaryPoly)@data@values[!is.na(raster::mask(AGB09,SecondaryPoly)@data@values)])/2
    
    length(raster::mask(AGB19,ForestPoly)@data@values[!is.na(raster::mask(AGB19,ForestPoly)@data@values)])/2 # Divide by 2 to convert from 0.5 ha pixels to area in ha
    length(raster::mask(AGB19,OldGrowthPoly)@data@values[!is.na(raster::mask(AGB19,OldGrowthPoly)@data@values)])/2
    length(raster::mask(AGB19,SecondaryPoly)@data@values[!is.na(raster::mask(AGB19,SecondaryPoly)@data@values)])/2
    
#### Make raster plot of mean biomass loss data ####
    
    # Make a raster AGB plots with height-AGB model mean values
      AGB09 <- CHM09_low
      AGB09@data@values[!is.na(AGB09@data@values)] <- 0.47*(8.881)*CHM09_low@data@values[!is.na(CHM09_low@data@values)]^(1.021)
                                                      
      AGB19 <- CHM19_low
      AGB19@data@values[!is.na(AGB19@data@values)] <- 0.47*(8.881)*CHM19_low@data@values[!is.na(CHM19_low@data@values)]^(1.021)
                                                        
      AGBchange <- AGB09
      AGBchange@data@values <- (AGB19@data@values-AGB09@data@values)
      
        pal <- colorRampPalette(c("red4","red3","red","tomato","white","lightblue","blue"))
        agbBrks <- seq(-60,20,5)
        raster::plot(AGBchange, breaks=agbBrks, col=pal(length(agbBrks)),colNA="yellow")
        
        raster::writeRaster(AGBchange,'Output/AGBchange.tif',options=c('TFW=YES'))

    
##### Make plot of delta NPV versus delta AGB for 2009-2019 #####
  
  # Make raster of change in AGB  
    deltaAGB <- AGB19-AGB09
    
  # Scale to proportion of biomass lost
    deltaAGBProp <- deltaAGB
    deltaAGBProp@data@values <- (deltaAGB@data@values)/AGB09@data@values*100
    
  # Load delta NPV data
    deltaNPV <- raster::raster(rgdal::readGDAL("SMA_Delta_NPV.tif"))
    deltaNPV_resampled <- raster::resample(x=deltaNPV, y=deltaAGB , method="bilinear")
    deltaNPV_forest <- raster::mask(deltaNPV_resampled, ForestPoly)
    
  mean(deltaNPV_forest@data@values[!is.na(deltaAGB@data@values)],na.rm=T)
  sd(deltaNPV_forest@data@values[!is.na(deltaAGB@data@values)],na.rm=T)
  min(deltaNPV_forest@data@values[!is.na(deltaAGB@data@values)],na.rm=T)
  max(deltaNPV_forest@data@values[!is.na(deltaAGB@data@values)],na.rm=T)
  
  summary(lm(deltaAGB@data@values[!is.na(deltaAGB@data@values)]~deltaNPV_forest@data@values[!is.na(deltaAGB@data@values)]))  
  
  par(mfrow=c(2,1),mar=c(3,4,0,1),oma=c(2,2,2,0))  
  plot(x=deltaNPV_forest@data@values[!is.na(deltaAGBProp@data@values)],
       y=deltaAGBProp@data@values[!is.na(deltaAGBProp@data@values)],
       xlab=NA,ylab=NA,
       pch=20, cex.axis=1)
  abline(h=0, lwd=2, col="grey")
  abline(v=0, lwd=2, col="grey")
  points(x=deltaNPV_forest@data@values[!is.na(deltaAGB@data@values)],
          y=deltaAGBProp@data@values[!is.na(deltaAGBProp@data@values)],
          pch=20)
  text("A", x=-0.055, y=60)
  mtext(expression(ACD~change~"(%)"),side=2,line=3, cex=1)


  plot(x=deltaNPV_forest@data@values[!is.na(deltaAGB@data@values)],
       y=deltaAGB@data@values[!is.na(deltaAGB@data@values)],
       xlab=NA,ylab=NA,
       pch=20, cex.axis=1)
  abline(h=0, lwd=2, col="grey")
  abline(v=0, lwd=2, col="grey")
  points(x=deltaNPV_forest@data@values[!is.na(deltaAGB@data@values)],
          y=deltaAGB@data@values[!is.na(deltaAGB@data@values)],
          pch=20)
  text("B", x=-0.055, y=40)
  mtext(expression(Delta~NPV~fraction),side=1,line=3, cex=1)
  mtext(expression(Delta~ACD~"(Mg C"~ha^{-1}~")"),side=2,line=3, cex=1)
  
  
  summary(lm(deltaAGBProp@data@values[!is.na(deltaAGBProp@data@values)]~deltaNPV_forest@data@values[!is.na(deltaAGBProp@data@values)]))  
  

##### Plot changes in AGB 2009-2019 ##### 
  usedOG <- sp::spTransform(usedOG, sp::proj4string(lidarCHM19))
  usedSF <- sp::spTransform(usedSF, sp::proj4string(lidarCHM19))
  
  usedAll <- rgeos::gUnion(usedOG,usedSF)
  
  pal <- colorRampPalette(c("red","orange","yellow","white","blue"))
  htBrks <- seq(-200,80,10)
  
  raster::plot(deltaAGB,ext=raster::extent(usedSF), lab.breaks=seq(-200,100,50), col=pal(length(htBrks)),
               axes=F,box=F)
  raster::plot(usedOG,add=T, lwd=3)
  raster::plot(usedSF, add=T, lwd=3,  lty=2)
  
##### Calculate annual AGBD change in CARBONO plots 1997 - 2017 #####
  # Get Wood production data per plot
    Wood <- read.csv("Wood.csv")
    
  # Make a plot variable
    Wood$plot <- substr(Wood$plot_treeid,start=1,stop=2)
    
  # Calculate AGB from allometric equation for each stem
    
  # Define function (agb.allometry) to estimate aboveground biomass from wood specific
  # gravity (wsg), tree diameter in centimeters (dbh).
    agb.allometry <- function(E,wsg,dbh){exp(-1.803-0.976*E+0.976*log(wsg) + 2.673*log(dbh) - 0.0299*(log(dbh)^2))}
    
  # Retrieve values of E for each plot using database from Chave et al. 2014
  # Call in source info:
  # source("http://chave.ups-tlse.fr/pantropical_allometry/readlayers.r")
  # Define coordinates of La Selva:
  # LS.coord <- cbind(-84.00,10.43); E.LaSelva <- retrieve_raster("E",LS.coord)
    E.LaSelva <- -0.06340053
    
  # Merge with other datasets to assign wood density to each stem
    TreeIDs <- read.csv("TreeIDs.csv") #Tree species from ID for each stem
    Species <- read.csv("Species.csv") #Wood density from species
    
    Wood <- merge(x=Wood, y=TreeIDs[,c("plot_treeid","genspcode")])
    Wood <- merge(x=Wood, y=Species[,c("genspcode","family","wooddens")], all.x=T)
    
  # First, fill in missing wood density values with family-level average
    FamilyWSG <- aggregate(Species$wooddens, by=list(Species$family), FUN="mean", na.rm=T)
    names(FamilyWSG) <- c("family","wooddens.family")
    
    Wood <- merge(x=Wood, y=FamilyWSG, all.x=T)
    
    Wood[is.na(Wood$wooddens),"wooddens"] <- Wood[is.na(Wood$wooddens),"wooddens.family"]
    
  # Second, fill in remaining missing wood density values with site-level average
    SiteWSG <- mean(Wood[!duplicated(Wood$plot_treeid),'wooddens'], na.rm=T)
    Wood[is.na(Wood$wooddens),"wooddens"] <- SiteWSG
    
  # Calculate AGB (convert diameter to cm)
    Wood$AGB <- agb.allometry(E.LaSelva, Wood$wooddens, Wood$dia_calc/10)
    
  # REMOVE one tree with error in DBH record?
    Wood <- Wood[!(Wood$plot_treeid=="L3268"),]
    
    plotIDs <- unique(Wood$plot)
    years <- range(Wood$dia_year)[1]:range(Wood$dia_year)[2]
    
    carbonoAGBD <- data.frame(Plot = rep(plotIDs,length(years)),
                              Year = rep(years,each=length(plotIDs)),
                              AGBD = NA)
    
    for(i in 1:length(plotIDs)){
      for(j in 1:length(years)){
        carbonoAGBD[carbonoAGBD$Plot==plotIDs[i] & carbonoAGBD$Year == years[j],"AGBD"] <- 
          sum(Wood[Wood$plot==plotIDs[i] & Wood$dia_year==years[j] & !(Wood$EBMkg_SB==-999), 'AGB'],na.rm=T)
        }
    }
          
    # Convert from kg to Mg C/ha
    carbonoAGBD$AGBD <- 0.47*2*(carbonoAGBD$AGBD)/1000
    
    # Average annual biomass gain between 1997 and 2016 
    avgAnnualChange <- (mean(carbonoAGBD[carbonoAGBD$Year==2016,"AGBD"]) - mean(carbonoAGBD[carbonoAGBD$Year==1997,"AGBD"]))/(2016-1997)
    
    avgAnnualRecent <- (mean(carbonoAGBD[carbonoAGBD$Year==2016,"AGBD"]) - mean(carbonoAGBD[carbonoAGBD$Year==2009,"AGBD"]))/(2016-2009)
    
    pctAvgAnnualRecent <- 100*((mean(carbonoAGBD[carbonoAGBD$Year==2016,"AGBD"])/mean(carbonoAGBD[carbonoAGBD$Year==2009,"AGBD"]))^(1/(2016-2009))-1)

    # Find time of regrowth
    
    # Original estimate
    regrowthTime <- (mean(AGBresults$allForest09) - mean(AGBresults$allForest19))/avgAnnualChange
    
    # AGBD loss including recent growth
    highLoss <- (mean(AGBresults$allForest09) + avgAnnualRecent*(2018-2009) - mean(AGBresults$allForest19))/mean(AGBresults$allForest09)*100
    
    # Regrowth estimate including recent growth
    regrowthTimeb <- (mean(AGBresults$allForest09) + avgAnnualRecent*(2018-2009) - mean(AGBresults$allForest19))/avgAnnualChange
    
    # Find limits by using maximum and minimum plot-level change
    changeAGBD <- data.frame(Plot = rep(plotIDs,each = length(years)-1),
                             Year1 = rep(years[1:length(years)-1],length(plotIDs)),
                             Year2 = rep(years[2:length(years)],length(plotIDs)),
                             dAGBD = NA)
    
    for(i in 1:length(plotIDs)){
      for(j in 1:length(years)-1){
        changeAGBD[changeAGBD$Plot==plotIDs[i] & changeAGBD$Year1 == years[j],"dAGBD"] <- 
          carbonoAGBD[carbonoAGBD$Plot==plotIDs[i] & carbonoAGBD$Year == years[j]+1, "AGBD"]-carbonoAGBD[carbonoAGBD$Plot==plotIDs[i] & carbonoAGBD$Year == years[j], "AGBD"]
        }
    }
    
    # Extreme change occurences
    
    changeAGBD[which(changeAGBD$dAGBD < -25*0.47),]
    
    recoveryEst <- data.frame(yr = 0:9,
                              dAGBD = NA)
    for(i in 1:10){
      recoveryEst$dAGBD[i] <- mean(changeAGBD[which(changeAGBD$dAGBD < -25*0.47)+recoveryEst$yr[i],"dAGBD"])
    }

    par(mar=c(4,4,1,1),oma=c(0,0,0,0))
    plot(dAGBD~yr, data = recoveryEst,
         pch=20,
         xlab = NA,
         ylab = NA)
    mtext(expression(Time~after~disturbance~"(yr)"),side=1,line=2.5, cex=1)
    mtext(expression(Delta~ACD~"(Mg C"~ha^{-1}~yr^{-1}~")"),side=2,line=2.5, cex=1)
    abline(h=0,col="black", lwd=2)
    abline(h=avgAnnualChange, lty=2, lwd=2, col="red")
    
    regrowthTime2 <- (mean(AGBresults$allForest09) - mean(AGBresults$allForest19) - sum(recoveryEst$dAGBD[2:6]))/avgAnnualChange + 6
      
      # Including recent growth
    regrowthTime2b <- (mean(AGBresults$allForest09) + avgAnnualRecent*(2018-2009) - mean(AGBresults$allForest19) - sum(recoveryEst$dAGBD[2:6]))/avgAnnualChange + 6

##### Revision: look at AGBD loss on trails, and in higher vs. lower forests #####
      
      # Must run sections "Make raster plot of mean biomass loss data" and 
      # "Make plot of delta NPV versus delta AGB for 2009-2019" first
      
      # Read in trail file
      trails <- rgdal::readOGR("Figure S1 info/trailsV3_utm16N.kml")
      trails <- sp::spTransform(trails, sp::proj4string(deltaAGB))
      
      # Divide pixels based on trail intersection or not
      trailPix <- raster::mask(deltaAGB,trails)
      noTrailPix <- raster::mask(deltaAGB,trails, inverse = T)
      
      trailPixProp <- raster::mask(deltaAGBProp,trails)
      noTrailPixProp <- raster::mask(deltaAGBProp,trails, inverse = T)
      
      # Get trail vs. no trail pixel values (Proportion ACD lost)
      trailVals <- trailPix@data@values[!is.na(trailPix@data@values)]
      noTrailVals <- noTrailPix@data@values[!is.na(noTrailPix@data@values)]
      
      trailValsProp <- trailPixProp@data@values[!is.na(trailPixProp@data@values)]
      noTrailValsProp <- noTrailPixProp@data@values[!is.na(noTrailPixProp@data@values)]
      
      # Not significanly different by t-test or K-S test
      t.test(trailVals, noTrailVals)
      ks.test(trailVals, noTrailVals)

      trailPixDens <- density(trailVals)
      noTrailPixDens <- density(noTrailVals)
      
      pdf("FigureS6_TrailEffects.PDF", width=5, height=4)
        par(mfrow=c(1,1), mar=c(3,3,1,1), oma=c(1,2,0,0), las=1)
        
        plot(noTrailPixDens, main=NA, lwd=2,
             ylim=range(c(trailPixDens$y,noTrailPixDens$y))+c(0,0.008))
        lines(trailPixDens, lty=2, lwd=2)
        abline(v=mean(trailVals),col="red",lwd=2,lty=2)
        abline(v=mean(noTrailVals),col="red",lwd=2,lty=1)
  
        legend(x="topleft",y=NULL,
               c("Off trail (n = 151 pixels)",
                 "On trail (n = 92 pixels)"),
               lwd=2,lty=c(1,2),
               bty="n", cex = 0.8)
        mtext(expression(Delta~ACD~"(Mg C"~ha^{-1}~")"), side=1, outer=F, line=3)
        mtext("Frequency", side=2, outer=F, las=0, line=3.5)
      dev.off()
      
      
      htData <- data.frame(AGB09 = AGB09@data@values,
                           deltaAGB = deltaAGB@data@values,
                           deltaAGBProp = deltaAGBProp@data@values)
      hta <- lm(deltaAGB~AGB09, data = htData)
      htb <- lm(deltaAGBProp~AGB09, data = htData)
      
      newx <- seq(min(AGB09@data@values,na.rm=T),
                  max(AGB09@data@values,na.rm=T),length.out = 100)
      
      conf_a <- predict(hta, newdata = data.frame(AGB09=newx), interval="confidence",
                        level = 0.95)
      conf_b <- predict(htb, newdata = data.frame(AGB09=newx),interval="confidence",
                        level = 0.95)

      
      pdf("FigureS7_OringinalForestHeightEffects.PDF", width=5, height=7)
        par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(1,2,0,0), las=1)
        
        plot(deltaAGB~AGB09, data = htData,
             pch=20,
             xlab= NA,ylab=NA)
        abline(hta)
        text("A", x=50,y=16)
        mtext(expression(Delta~ACD~"(Mg C"~ha^{-1}~")"),side=2,line=2.5, cex=1, las=0)
        
        plot(deltaAGBProp~AGB09, data = htData,
             pch=20,
             xlab= NA,ylab=NA)
        abline(htb)
        text("B", x=50,y=20)
  
        mtext(expression(Delta~ACD~"(%)"),side=2,line=2.5, cex=1, las=0)
        mtext(expression("2009"~ACD~"(Mg C"~ha^{-1}~")"),side=1,line=2.5, cex=1, las=0)
      dev.off()

 



      