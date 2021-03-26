#### Define new funtions used below ####

# Make binary gap raster
  makeGapRaster <- function(CHMraster, gapHeight=2){
    gapRaster <- CHMraster
    gapRaster@data@values[!is.na(CHMraster@data@values) & (CHMraster@data@values <= gapHeight)] <- 1 
    gapRaster@data@values[!is.na(CHMraster@data@values) & (CHMraster@data@values > gapHeight)] <- 0 
    return(gapRaster)
  }

# Get vector of gap sizes     
  getGapSizes <- function(gapRaster){
    clumpRaster <- raster::clump(gapRaster, directions=8, gaps=F)
    gapFreq <- raster::freq(clumpRaster)
    gapSizes <- gapFreq[!is.na(gapFreq[,1]),]
    gapSizeVector <- gapSizes[,2]
    return(gapSizeVector)
  }

# Define a Metropolis-Hastings MCMC algorithm to esimates the posterior distribution of lambda    
# Based on recommendations from: https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
  
  # Data: vector of individual gap sizes (in # of pixels)
  
  # Likelihood function: 
  
    dpl <- function(param, data){
      sum(log(data^-param/VGAM::zeta(param)))
    }
  
  # Define a prior for lambda
    prior <- function(param){
      lambdaPrior <- dunif(param, min=1.001, max=5, log=T)
      #lambdaPrior <- dgamma(param-1, 0.001, 0.001, log=T)
      return(lambdaPrior)
    }
    
  # Define the posterior, the sum of the likelihood function and the prior (because log-transformed)
    
    posterior <- function(data, param){
      return(dpl(param,data) + prior(param))
    }
    
  # Define proposal function (adjust to get acceptance rate ~ 20%)
    proposalFunction <- function(param){
      return(rnorm(1,  mean=param, sd=0.11))
      #return(rgamma(1, shape=2500, rate=2500/(param-1))+1)
    }
    

  # Define MCMC routine
    
    runMCMC <- function(startvalue, iterations, data){
      
      chain = array(dim = c(iterations+1,1))
      chain[1] = startvalue
      accept = array(dim= c(iterations,1))
      
      for (i in 1:iterations){
        proposal = proposalFunction(chain[i])
        
        # Proposal must be greater than 1
        while(proposal < 1){
          proposal = proposalFunction(chain[i])
        }
        
        q1 <- dgamma(chain[i]-1, shape=2500, rate=2500/(proposal-1), log=T)
        q2 <- dgamma((proposal-1), shape=2500, rate=2500/(chain[i]-1), log=T)
        
        probabNumerator = exp(posterior(data, param = proposal) - posterior(data, param = chain[i]))
        probabDenominator = exp(q2-q1)
        #probab <- probabNumerator/probabDenominator
        probab <- probabNumerator
        
        if (runif(1) < probab){
          chain[i+1,] = proposal
          accept[i] = 1
        }else{
          chain[i+1,] = chain[i,]
          accept[i] = 0
        }
      }
      
      return(list(chain,accept))
      
    }    
    
    # # Test MCMC routine on simulated data
    #   lambda <- 2.5
    #   gapValues <- c(1:10000)
    #   gapProbs <- (gapValues^-lambda)/VGAM::zeta(lambda)
    #   gapSims <- sample(size=1000, x=gapValues, prob=gapProbs, replace=T)
    #   
    #   fitdp <- optimize(dpl, data=gapSims, lower = 1.0001, upper = 20, maximum = T)
    # 
    #   fitSim <- runMCMC(startvalue = fitdp$maximum, iterations = 5000, data = gapSims)
    #   hist(fitSim[[1]][1001:5000]);abline(v=lambda,col="red")
    #   sum(fitSim[[2]]);mean(fitSim[[1]][1001:5000])-lambda
    
#### Read in data ####

  # 1.25 m resolution canopy height rasters from 2009 and 2019
    lidarCHM09 <- raster::raster("CHM09_AOI_res125.tif")
    lidarCHM19 <- raster::raster("CHM19_AOI_res125.tif")

    # Assign negative values to 0  
      lidarCHM09@data@values[lidarCHM09@data@values < 0 & !is.na(lidarCHM09@data@values)] <- 0
      lidarCHM19@data@values[lidarCHM19@data@values < 0 & !is.na(lidarCHM19@data@values)] <- 0

    # Make forest type polygons
      # Read shapefile
      LaSelvaClass <- rgdal::readOGR("LandUseShapefile/LU2000.shp")
      
      # All forested area
      ForestClasses <- c("Old-growth Forests","Ecological Reserve","Forested Swamps","Secondary Forests","Abandoned Agroforestry","Abandoned Plantation", "Selectively-logged Forests")
      ForestPoly <- LaSelvaClass[LaSelvaClass$DESCRIPTIO %in% ForestClasses,]
      # Remove forest area added in or after 1996
       ForestPoly <- ForestPoly[!(ForestPoly$YEAR >= 1990),]
      # Remove tiny fragment of abandoned plantation
      ForestPoly <- ForestPoly[!(ForestPoly$LU00_UTM0_==68),]
    
    # Old growth forest area
    OldGrowthPoly <- ForestPoly[ForestPoly$DESCRIPTIO %in% c("Old-growth Forests","Ecological Reserve","Forested Swamps"),]
    
    # Secondary forest area
    SecondaryPoly <- ForestPoly[ForestPoly$DESCRIPTIO %in% c("Secondary Forests","Abandoned Agroforestry","Abandoned Plantation", "Selectively-logged Forests"),]
    
  # Make rasters for all forst, old growth forest, and secondary forest  
    Forest09 <- raster::mask(lidarCHM09,ForestPoly)
    Forest19 <- raster::mask(lidarCHM19,ForestPoly)
  
    OldGrowth09 <- raster::mask(Forest09, OldGrowthPoly)
    OldGrowth19 <- raster::mask(Forest19, OldGrowthPoly)
    
    Secondary09 <- raster::mask(Forest09, SecondaryPoly)
    Secondary19 <- raster::mask(Forest19, SecondaryPoly)

##### Cycle through multiple gap height thresholds, calculating size-frequency distribution values #####  

  # Use gap height thresholds from 2 m to 20 m, in increments of 2  
    gapHts <- seq(2,20,2)  
  
  # Make lists to hold vectors of gap sizes
  # Calculate separately for all, old growth, and secondary forests  
    listGaps09 <- list()
    listGaps19 <- list()
    listGaps09_OldGrowth <- list()
    listGaps19_OldGrowth <- list()
    listGaps09_Secondary <- list()
    listGaps19_Secondary <- list()

  # Make lists to hold gap size-freq distribution parameter from each MCMC step
    listFit09 <- list()
    listFit19 <- list()
    listFit09_OldGrowth <- list()
    listFit19_OldGrowth <- list()
    listFit09_Secondary <- list()
    listFit19_Secondary <- list()

# Loop through height thresholds and estimate gap-size frequency distributions
# NOTE: this takes ~2 hrs to run    
set.seed(1)
for(j in 1:length(gapHts)){
  
  # Make a raster of gap pixels
    gapRaster09 <- makeGapRaster(CHMraster = Forest09, gapHeight = gapHts[j])
    gapRaster19 <- makeGapRaster(CHMraster = Forest19, gapHeight = gapHts[j])
    
    gapRaster09_OldGrowth <- makeGapRaster(CHMraster = OldGrowth09, gapHeight = gapHts[j])
    gapRaster19_OldGrowth <- makeGapRaster(CHMraster = OldGrowth19, gapHeight = gapHts[j])
    
    gapRaster09_Secondary <- makeGapRaster(CHMraster = Secondary09, gapHeight = gapHts[j])
    gapRaster19_Secondary <- makeGapRaster(CHMraster = Secondary19, gapHeight = gapHts[j])
    
  # Get a vector of all gap sizes (one entry per gap)
    gaps09 <- getGapSizes(gapRaster09)
    gaps19 <- getGapSizes(gapRaster19)
    
    gaps09_OldGrowth <- getGapSizes(gapRaster09_OldGrowth)
    gaps19_OldGrowth <- getGapSizes(gapRaster19_OldGrowth)
    
    gaps09_Secondary <- getGapSizes(gapRaster09_Secondary)
    gaps19_Secondary <- getGapSizes(gapRaster19_Secondary)

  # Run MCMC to estimate gap size frequency distribution parameter  
    sims <- 100000
    fit09 <- runMCMC(startvalue = optimize(dpl, data=gaps09, lower = 1.0001, upper = 20, maximum = T)$maximum, iterations = sims, data = gaps09)
    fit19 <- runMCMC(startvalue = optimize(dpl, data=gaps19, lower = 1.0001, upper = 20, maximum = T)$maximum, iterations = sims, data = gaps19)
    
    fit09_OldGrowth <- runMCMC(startvalue = optimize(dpl, data=gaps09_OldGrowth, lower = 1.0001, upper = 20, maximum = T)$maximum, iterations = sims, data = gaps09_OldGrowth)
    fit19_OldGrowth <- runMCMC(startvalue = optimize(dpl, data=gaps19_OldGrowth, lower = 1.0001, upper = 20, maximum = T)$maximum, iterations = sims, data = gaps19_OldGrowth)
    
    fit09_Secondary <- runMCMC(startvalue = optimize(dpl, data=gaps09_Secondary, lower = 1.0001, upper = 20, maximum = T)$maximum, iterations = sims, data = gaps09_Secondary)
    fit19_Secondary <- runMCMC(startvalue = optimize(dpl, data=gaps19_Secondary, lower = 1.0001, upper = 20, maximum = T)$maximum, iterations = sims, data = gaps19_Secondary) 
  
  # Save results in lists
    listGaps09[[j]] <- gaps09
    listGaps19[[j]] <- gaps19
    listGaps09_OldGrowth[[j]] <- gaps09_OldGrowth
    listGaps19_OldGrowth[[j]] <- gaps19_OldGrowth
    listGaps09_Secondary[[j]] <- gaps09_Secondary
    listGaps19_Secondary[[j]] <- gaps19_Secondary
    
    listFit09[[j]] <- fit09
    listFit19[[j]] <- fit19
    listFit09_OldGrowth[[j]] <- fit09_OldGrowth
    listFit19_OldGrowth[[j]] <- fit19_OldGrowth
    listFit09_Secondary[[j]] <- fit09_Secondary
    listFit19_Secondary[[j]] <- fit19_Secondary
}

#### Consolidate results for different height threshholds ####

# Make a data frame to store results (this is table )
  gapAreaResults <- data.frame(forestType = rep(rep(c("All","OldGrowth","Secondary"),each=length(gapHts)),2),
                               gapHt = rep(rep(gapHts,3),2),
                               Year = rep(c(2009,2019),each=length(rep(gapHts,3))),
                               Area = NA,
                               LambdaMin = NA,
                               LambdaMed = NA,
                               LambdaMax = NA)
# Define burn-in period and sampling periodicity for estimating gap size-frequency 
# parameter from MCMC results
    burnIn <- 5001
    skipN <- 25

# Cycle through gap thresholds to calculate results    
for(i in 1:length(gapHts)){
  
  # Add total area in gaps
    area09 <- round(sum(listGaps09[[i]])*(1.25^2)/10000,2)
    area19 <- round(sum(listGaps19[[i]])*(1.25^2)/10000,2)
    
    old09 <- round(sum(listGaps09_OldGrowth[[i]])*(1.25^2)/10000,2)
    old19 <- round(sum(listGaps19_OldGrowth[[i]])*(1.25^2)/10000,2)
    
    sec09 <- round(sum(listGaps09_Secondary[[i]])*(1.25^2)/10000,2)
    sec19 <- round(sum(listGaps19_Secondary[[i]])*(1.25^2)/10000,2)
    
  # Estimate lambda (gap size frequency parameter)
    lambdaMin09 <- quantile(listFit09[[i]][[1]][seq(burnIn,sims,skipN)], 0.025)
    lambdaMed09 <- quantile(listFit09[[i]][[1]][seq(burnIn,sims,skipN)], 0.50)
    lambdaMax09 <- quantile(listFit09[[i]][[1]][seq(burnIn,sims,skipN)], 0.975)
    
    lambdaMin19 <- quantile(listFit19[[i]][[1]][seq(burnIn,sims,skipN)], 0.025)
    lambdaMed19 <- quantile(listFit19[[i]][[1]][seq(burnIn,sims,skipN)], 0.50)
    lambdaMax19 <- quantile(listFit19[[i]][[1]][seq(burnIn,sims,skipN)], 0.975)
    
    lambdaMin09old <- quantile(listFit09_OldGrowth[[i]][[1]][seq(burnIn,sims,skipN)], 0.025)
    lambdaMed09old <- quantile(listFit09_OldGrowth[[i]][[1]][seq(burnIn,sims,skipN)], 0.50)
    lambdaMax09old <- quantile(listFit09_OldGrowth[[i]][[1]][seq(burnIn,sims,skipN)], 0.975)
    
    lambdaMin19old <- quantile(listFit19_OldGrowth[[i]][[1]][seq(burnIn,sims,skipN)], 0.025)
    lambdaMed19old <- quantile(listFit19_OldGrowth[[i]][[1]][seq(burnIn,sims,skipN)], 0.50)
    lambdaMax19old <- quantile(listFit19_OldGrowth[[i]][[1]][seq(burnIn,sims,skipN)], 0.975)
    
    lambdaMin09sec <- quantile(listFit09_Secondary[[i]][[1]][seq(burnIn,sims,skipN)], 0.025)
    lambdaMed09sec <- quantile(listFit09_Secondary[[i]][[1]][seq(burnIn,sims,skipN)], 0.50)
    lambdaMax09sec <- quantile(listFit09_Secondary[[i]][[1]][seq(burnIn,sims,skipN)], 0.975)
    
    lambdaMin19sec <- quantile(listFit19_Secondary[[i]][[1]][seq(burnIn,sims,skipN)], 0.025)
    lambdaMed19sec <- quantile(listFit19_Secondary[[i]][[1]][seq(burnIn,sims,skipN)], 0.50)
    lambdaMax19sec <- quantile(listFit19_Secondary[[i]][[1]][seq(burnIn,sims,skipN)], 0.975)
  
  # Store results in data frame  
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "All","Area"] <- area09
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "All","Area"] <- area19
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "OldGrowth","Area"] <- old09
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "OldGrowth","Area"] <- old19
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "Secondary","Area"] <- sec09
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "Secondary","Area"] <- sec19
    
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "All","LambdaMin"] <- lambdaMin09
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "All","LambdaMin"] <- lambdaMin19
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "OldGrowth","LambdaMin"] <- lambdaMin09old
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "OldGrowth","LambdaMin"] <- lambdaMin19old
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "Secondary","LambdaMin"] <- lambdaMin09sec
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "Secondary","LambdaMin"] <- lambdaMin19sec
    
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "All","LambdaMed"] <- lambdaMed09
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "All","LambdaMed"] <- lambdaMed19
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "OldGrowth","LambdaMed"] <- lambdaMed09old
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "OldGrowth","LambdaMed"] <- lambdaMed19old
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "Secondary","LambdaMed"] <- lambdaMed09sec
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "Secondary","LambdaMed"] <- lambdaMed19sec
    
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "All","LambdaMax"] <- lambdaMax09
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "All","LambdaMax"] <- lambdaMax19
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "OldGrowth","LambdaMax"] <- lambdaMax09old
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "OldGrowth","LambdaMax"] <- lambdaMax19old
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2009" & gapAreaResults$forestType == "Secondary","LambdaMax"] <- lambdaMax09sec
    gapAreaResults[gapAreaResults$gapHt == gapHts[i] & gapAreaResults$Year == "2019" & gapAreaResults$forestType == "Secondary","LambdaMax"] <- lambdaMax19sec
}

# Quantify change in gap area for results
    
  # Calculate total area in all, old growth, and secondary forests  
    haTot <- length(Forest09@data@values[!is.na(Forest09@data@values)])*(1.25^2)/10000
    haOld <- length(OldGrowth09@data@values[!is.na(OldGrowth09@data@values)])*(1.25^2)/10000
    Sec09Values <- raster::getValues(Secondary09)
    haSec <- length(Sec09Values[!is.na(Sec09Values)])*(1.25^2)/10000
  
  # Calculate gap area as a percent of total area for each forest class  
    gapAreaResults$pctArea <- NA
    gapAreaResults[gapAreaResults$forestType=="All","pctArea"] <- gapAreaResults[gapAreaResults$forestType=="All","Area"]/haTot*100
    gapAreaResults[gapAreaResults$forestType=="OldGrowth","pctArea"] <- gapAreaResults[gapAreaResults$forestType=="OldGrowth","Area"]/haOld*100
    gapAreaResults[gapAreaResults$forestType=="Secondary","pctArea"] <- gapAreaResults[gapAreaResults$forestType=="Secondary","Area"]/haSec*100
        
  # Calculate percent increase in gap area
    pctIncrease <- round(gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2019,"pctArea"]-gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2009,"pctArea"],1)
    propIncresase <- pctIncrease/gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2009,"pctArea"]*100      
  
    pctIncreaseO <- round(gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2019,"pctArea"]-gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2009,"pctArea"],1)
    propIncreaseO <- pctIncrease/gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2009,"pctArea"]*100      
  
    pctIncreaseS <- round(gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2019,"pctArea"]-gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2009,"pctArea"],1)
    propIncreaseS <- pctIncrease/gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2009,"pctArea"]*100      

    mean(propIncreaseO)
    mean(propIncreaseS)
  
#### Plot Fig 2: gap area and scaling parameter vs. height ####

     par(mfrow=c(1,2), mar=c(2,2,1,1), oma=c(2,2,0,0)) 
     
      plot(gapHt ~ pctArea, data = gapAreaResults[gapAreaResults$forestType=="All",],
      xlim=range(c(gapAreaResults$pctArea,gapAreaResults$pctArea)),
      type = "n", 
      xlab = NA, ylab = NA,
      las = 1)
      
      lines(gapHt ~ pctArea, data = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2009,],
             col = "black",
             lwd=2)
      lines(gapHt ~ pctArea, data = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2019,],
             col = "red",
             lwd=2)
      
       text(x=1.5, y=20, "a", cex=1.2)
       
             legend(x=22,y=8,bty="n",
             c("Before","After"),
             seg.len = 0.5,x.intersp = 0.5,
             col=c("black","red"),
             lwd=3, cex=1.2)
       
      mtext(side = 1, line = 2.5,expression(Gap~area~("%")), cex=1.2)
      mtext(side = 2, line = 2, expression(Gap~height~("m")), cex=1.2)
      
     plot(gapHt ~ LambdaMin, data = gapAreaResults[gapAreaResults$forestType=="All",],
          xlim=range(c(gapAreaResults$LambdaMin,gapAreaResults$LambdaMax)),
          type = "n", yaxt = "n",
          xlab = NA, ylab = NA)
     axis(2, labels=F)
     
     arrows(x0 = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2009,"LambdaMin"],
            x1 = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2009,"LambdaMax"],
            y0 = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2009,"gapHt"],
            y1 = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2009,"gapHt"],
            col = adjustcolor("black",alpha.f = 0.8),
            angle = 90, lwd=1.5, length = 0.05, code=3)
     
      arrows(x0 = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2019,"LambdaMin"],
            x1 = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2019,"LambdaMax"],
            y0 = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2019,"gapHt"],
            y1 = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2019,"gapHt"],
            col = adjustcolor("red",alpha.f = 0.8),
            angle = 90, lwd=1.5, length = 0.05, code=3)
      
      text(x=1.355, y=20, "b", cex=1.2)
      
      
      points(gapHt ~ LambdaMed, data = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2009,],
             col = "black",
             pch = 20, cex = 1.5)
      points(gapHt ~ LambdaMed, data = gapAreaResults[gapAreaResults$forestType=="All" & gapAreaResults$Year==2019,],
             col = "red",
             pch = 20, cex = 1.5)
      
      mtext(side = 1, line = 2.3,expression(lambda), cex=1.2)
      
      
#### Plot Fig S5: plots of size frequency ####
    
      par(mfrow=c(3,2), mar=c(2,2,1,1), oma=c(2,3,0,0))
      
        # all forest
        gapFreq2 <- aggregate(listGaps09[[1]], by=list(listGaps09[[1]]), FUN="length")
        gapFreq6 <- aggregate(listGaps09[[3]], by=list(listGaps09[[3]]), FUN="length")
        gapFreq10 <- aggregate(listGaps09[[5]], by=list(listGaps09[[5]]), FUN="length")
        gapFreq20 <- aggregate(listGaps09[[10]], by=list(listGaps09[[10]]), FUN="length")
        
        xmax <- max(c(gapFreq2[,1],gapFreq6[,1],gapFreq10[,1],gapFreq20[,1]))*(1.25^2)
        ymax <- max(c(gapFreq2[,2],gapFreq6[,2],gapFreq10[,2],gapFreq20[,2]))
        
        plot(y=gapFreq2[,2], x=gapFreq2[,1]*(1.25^2), log="xy", pch=20,
             cex=2, col=adjustcolor("black",alpha.f = 0.5),
             xlim=c(1,xmax),ylim=c(1,ymax), cex.axis=1.3,
             main = "2009")
        points(y=gapFreq6[,2], x=gapFreq6[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("blue",alpha.f = 0.5))
        points(y=gapFreq10[,2], x=gapFreq10[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("brown",alpha.f = 0.5))
        points(y=gapFreq20[,2], x=gapFreq20[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("pink",alpha.f = 0.5))
        text("A", x=0.9,y=ymax, cex=1)
        
        legend(x=2000, y=ymax+1000,
               c("2 m","6 m","10 m","20 m"),
               col=adjustcolor(c("black","blue","brown","pink"),alpha.f = 0.5),
               y.intersp = 0.8,
               pch=20,cex=1.5,
               bty="n")

        gapFreq2 <- aggregate(listGaps19[[1]], by=list(listGaps19[[1]]), FUN="length")
        gapFreq6 <- aggregate(listGaps19[[2]], by=list(listGaps19[[2]]), FUN="length")
        gapFreq10 <- aggregate(listGaps19[[3]], by=list(listGaps19[[3]]), FUN="length")
        gapFreq20 <- aggregate(listGaps19[[4]], by=list(listGaps19[[4]]), FUN="length")
        
        xmax <- max(c(gapFreq2[,1],gapFreq6[,1],gapFreq10[,1],gapFreq20[,1]))*(1.25^2)
        ymax <- max(c(gapFreq2[,2],gapFreq6[,2],gapFreq10[,2],gapFreq20[,2]))
        
        plot(y=gapFreq2[,2], x=gapFreq2[,1]*(1.25^2), log="xy", pch=20,
             cex=2, col=adjustcolor("black",alpha.f = 0.5),
             xlim=c(1,xmax),ylim=c(1,ymax), cex.axis=1.3,
             main = "2019")
        points(y=gapFreq6[,2], x=gapFreq6[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("blue",alpha.f = 0.5))
        points(y=gapFreq10[,2], x=gapFreq10[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("brown",alpha.f = 0.5))
        points(y=gapFreq20[,2], x=gapFreq20[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("pink",alpha.f = 0.5))
        text("B", x=0.9,y=ymax, cex=1)
      
        # Old growth forest
        gapFreq2 <- aggregate(listGaps09_OldGrowth[[1]], by=list(listGaps09_OldGrowth[[1]]), FUN="length")
        gapFreq6 <- aggregate(listGaps09_OldGrowth[[3]], by=list(listGaps09_OldGrowth[[3]]), FUN="length")
        gapFreq10 <- aggregate(listGaps09_OldGrowth[[5]], by=list(listGaps09_OldGrowth[[5]]), FUN="length")
        gapFreq20 <- aggregate(listGaps09_OldGrowth[[10]], by=list(listGaps09_OldGrowth[[10]]), FUN="length")
        
        xmax <- max(c(gapFreq2[,1],gapFreq6[,1],gapFreq10[,1],gapFreq20[,1]))*(1.25^2)
        ymax <- max(c(gapFreq2[,2],gapFreq6[,2],gapFreq10[,2],gapFreq20[,2]))
        
        plot(y=gapFreq2[,2], x=gapFreq2[,1]*(1.25^2), log="xy", pch=20,
             cex=2, col=adjustcolor("black",alpha.f = 0.5),
             xlim=c(1,xmax),ylim=c(1,ymax), cex.axis=1.3)
        points(y=gapFreq6[,2], x=gapFreq6[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("blue",alpha.f = 0.5))
        points(y=gapFreq10[,2], x=gapFreq10[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("brown",alpha.f = 0.5))
        points(y=gapFreq20[,2], x=gapFreq20[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("pink",alpha.f = 0.5))
        text("C", x=0.9,y=ymax, cex=1)
        
        gapFreq2 <- aggregate(listGaps19_OldGrowth[[1]], by=list(listGaps19_OldGrowth[[1]]), FUN="length")
        gapFreq6 <- aggregate(listGaps19_OldGrowth[[2]], by=list(listGaps19_OldGrowth[[2]]), FUN="length")
        gapFreq10 <- aggregate(listGaps19_OldGrowth[[3]], by=list(listGaps19_OldGrowth[[3]]), FUN="length")
        gapFreq20 <- aggregate(listGaps19_OldGrowth[[4]], by=list(listGaps19_OldGrowth[[4]]), FUN="length")
        
        xmax <- max(c(gapFreq2[,1],gapFreq6[,1],gapFreq10[,1],gapFreq20[,1]))*(1.25^2)
        ymax <- max(c(gapFreq2[,2],gapFreq6[,2],gapFreq10[,2],gapFreq20[,2]))
        
        plot(y=gapFreq2[,2], x=gapFreq2[,1]*(1.25^2), log="xy", pch=20,
             cex=2, col=adjustcolor("black",alpha.f = 0.5),
             xlim=c(1,xmax),ylim=c(1,ymax), cex.axis=1.3)
        points(y=gapFreq6[,2], x=gapFreq6[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("blue",alpha.f = 0.5))
        points(y=gapFreq10[,2], x=gapFreq10[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("brown",alpha.f = 0.5))
        points(y=gapFreq20[,2], x=gapFreq20[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("pink",alpha.f = 0.5))
        text("D", x=0.9,y=ymax, cex=1)
      
        # Secondary forest
        gapFreq2 <- aggregate(listGaps09_Secondary[[1]], by=list(listGaps09_Secondary[[1]]), FUN="length")
        gapFreq6 <- aggregate(listGaps09_Secondary[[3]], by=list(listGaps09_Secondary[[3]]), FUN="length")
        gapFreq10 <- aggregate(listGaps09_Secondary[[5]], by=list(listGaps09_Secondary[[5]]), FUN="length")
        gapFreq20 <- aggregate(listGaps09_Secondary[[10]], by=list(listGaps09_Secondary[[10]]), FUN="length")
        
        xmax <- max(c(gapFreq2[,1],gapFreq6[,1],gapFreq10[,1],gapFreq20[,1]))*(1.25^2)
        ymax <- max(c(gapFreq2[,2],gapFreq6[,2],gapFreq10[,2],gapFreq20[,2]))
        
        plot(y=gapFreq2[,2], x=gapFreq2[,1]*(1.25^2), log="xy", pch=20,
             cex=2, col=adjustcolor("black",alpha.f = 0.5),
             xlim=c(1,xmax),ylim=c(1,ymax), cex.axis=1.3)
        points(y=gapFreq6[,2], x=gapFreq6[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("blue",alpha.f = 0.5))
        points(y=gapFreq10[,2], x=gapFreq10[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("brown",alpha.f = 0.5))
        points(y=gapFreq20[,2], x=gapFreq20[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("pink",alpha.f = 0.5))
        text("E", x=0.9,y=ymax, cex=1)
  
        gapFreq2 <- aggregate(listGaps19_Secondary[[1]], by=list(listGaps19_Secondary[[1]]), FUN="length")
        gapFreq6 <- aggregate(listGaps19_Secondary[[2]], by=list(listGaps19_Secondary[[2]]), FUN="length")
        gapFreq10 <- aggregate(listGaps19_Secondary[[3]], by=list(listGaps19_Secondary[[3]]), FUN="length")
        gapFreq20 <- aggregate(listGaps19_Secondary[[4]], by=list(listGaps19_Secondary[[4]]), FUN="length")
        
        xmax <- max(c(gapFreq2[,1],gapFreq6[,1],gapFreq10[,1],gapFreq20[,1]))*(1.25^2)
        ymax <- max(c(gapFreq2[,2],gapFreq6[,2],gapFreq10[,2],gapFreq20[,2]))
        
        plot(y=gapFreq2[,2], x=gapFreq2[,1]*(1.25^2), log="xy", pch=20,
             cex=2, col=adjustcolor("black",alpha.f = 0.5),
             xlim=c(1,xmax),ylim=c(1,ymax), cex.axis=1.3)
        points(y=gapFreq6[,2], x=gapFreq6[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("blue",alpha.f = 0.5))
        points(y=gapFreq10[,2], x=gapFreq10[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("brown",alpha.f = 0.5))
        points(y=gapFreq20[,2], x=gapFreq20[,1]*(1.25^2), pch=20,
               cex=2, col=adjustcolor("pink",alpha.f = 0.5))
        text("F", x=0.9,y=ymax, cex=1)
      
      mtext(expression(Gap~size~(m^{2})),side=1,outer=T, cex=1.3, line=1)
      mtext(expression(Frequency),side=2,outer=T, cex=1.3, line=0.5) 
      
#### Plot Fig. S6: gap area and scaling parameter vs. height by forest type ####
      par(mfrow=c(2,2), mar=c(2,2,1,1), oma=c(2,2,0,0)) 
     
      plot(gapHt ~ pctArea, data = gapAreaResults[gapAreaResults$forestType=="OldGrowth",],
      xlim=range(c(gapAreaResults$pctArea,gapAreaResults$pctArea)),
      type = "n", 
      xlab = NA, ylab = NA)
      
      lines(gapHt ~ pctArea, data = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2009,],
             col = "black",
             lwd=2)
      lines(gapHt ~ pctArea, data = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2019,],
             col = "red",
             lwd=2)
      
       text(x=1.2, y=20, "A", cex=1.1)
       
             legend(x=15,y=10,bty="n",
             c(2009,2019),
             seg.len = 0.5,x.intersp = 0.5,
             col=c("black","red"),
             lwd=3, cex=1.3)
       

      mtext(side = 2, line = 0, expression(Gap~height~("m")), outer = T)
      
     plot(gapHt ~ LambdaMin, data = gapAreaResults[gapAreaResults$forestType=="OldGrowth",],
          xlim=range(c(gapAreaResults$LambdaMin,gapAreaResults$LambdaMax))+c(-0.02,0),
          type = "n", 
          xlab = NA, ylab = NA)
     
     arrows(x0 = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2009,"LambdaMin"],
            x1 = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2009,"LambdaMax"],
            y0 = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2009,"gapHt"],
            y1 = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2009,"gapHt"],
            col = adjustcolor("black",alpha.f = 0.8),
            angle = 90, lwd=1.5, length = 0.05, code=3)
     
      arrows(x0 = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2019,"LambdaMin"],
            x1 = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2019,"LambdaMax"],
            y0 = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2019,"gapHt"],
            y1 = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2019,"gapHt"],
            col = adjustcolor("red",alpha.f = 0.8),
            angle = 90, lwd=1.5, length = 0.05, code=3)
      
      text(x=1.33, y=20, "B", cex=1.1)
      
      
      points(gapHt ~ LambdaMed, data = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2009,],
             col = "black",
             pch = 20, cex = 1.5)
      points(gapHt ~ LambdaMed, data = gapAreaResults[gapAreaResults$forestType=="OldGrowth" & gapAreaResults$Year==2019,],
             col = "red",
             pch = 20, cex = 1.5)
      
      
      # SECONDARY
      plot(gapHt ~ pctArea, data = gapAreaResults[gapAreaResults$forestType=="Secondary",],
      xlim=range(c(gapAreaResults$pctArea,gapAreaResults$pctArea)),
      type = "n", 
      xlab = NA, ylab = NA)
      
      lines(gapHt ~ pctArea, data = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2009,],
             col = "black",
             lwd=2)
      lines(gapHt ~ pctArea, data = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2019,],
             col = "red",
             lwd=2)
      
       text(x=1.2, y=20, "C", cex=1.1)
       
      mtext(side = 1, line = 2.3,expression(Gap~area~("%")))
      
     plot(gapHt ~ LambdaMin, data = gapAreaResults[gapAreaResults$forestType=="Secondary",],
          xlim=range(c(gapAreaResults$LambdaMin,gapAreaResults$LambdaMax))+c(-0.02,0),
          type = "n", 
          xlab = NA, ylab = NA)
     
     arrows(x0 = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2009,"LambdaMin"],
            x1 = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2009,"LambdaMax"],
            y0 = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2009,"gapHt"],
            y1 = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2009,"gapHt"],
            col = adjustcolor("black",alpha.f = 0.8),
            angle = 90, lwd=1.5, length = 0.05, code=3)
     
      arrows(x0 = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2019,"LambdaMin"],
            x1 = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2019,"LambdaMax"],
            y0 = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2019,"gapHt"],
            y1 = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2019,"gapHt"],
            col = adjustcolor("red",alpha.f = 0.8),
            angle = 90, lwd=1.5, length = 0.05, code=3)
      
      text(x=1.33, y=20, "D", cex=1.1)
      
      
      points(gapHt ~ LambdaMed, data = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2009,],
             col = "black",
             pch = 20, cex = 1.5)
      points(gapHt ~ LambdaMed, data = gapAreaResults[gapAreaResults$forestType=="Secondary" & gapAreaResults$Year==2019,],
             col = "red",
             pch = 20, cex = 1.5)
      
      mtext(side = 1, line = 2,expression(lambda))
      
