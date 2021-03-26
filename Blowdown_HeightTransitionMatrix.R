#### Define relevant functions ####

# An function to generate a transition matrix
# data is the file that contains observations from which the transition matrix will be assembled
# In the current version, the function assumes that the first column of the data is the "from"
# variable, and the second column is the "to" variable.  Here, column 1 is initial canopy height
# and column 2 is final canopy height
# max.height is the maximum height for which transitions will be calculated.

transition.matrix <- function(data, max.height){

transition <<- matrix(NA, nrow=max.height, ncol=max.height)

for (i in 1:max.height){

	temp <- subset(data, data[,1] < i & data[,1] >= i - 1)

	for (j in 1:max.height){

		temp2 <- subset(temp, temp[,2] < j & temp[,2] >= j - 1)

		transition[j,i] <<- dim(temp2)[1]
	}
}

return(transition)
}

# An function to create a transition probability matrix from a transition matrix
# transition.matrix is the transition matrix produced by function transition.matrix
transition.probability.matrix <- function(transition.matrix){

transition.p <<- matrix(NA, nrow=dim(transition.matrix)[1], ncol=dim(transition.matrix)[1])

for (i in 1:dim(transition.matrix)[1]){

	for (j in 1:dim(transition.matrix)[1]){

	transition.p[i,j] <<- transition.matrix[i,j]/sum(transition.matrix[,j])

	}
}
transition.p[transition.p[,]=="NaN"] <<- 0

return(transition.p)

}


# An function to simulate the stable stage distribution using non-parametric bootstrap sampling
# n.sim is the number of iterations

simulate <- function(n.sim){

	time <- Sys.time()

	results <<- matrix(NA, nrow=56, ncol=n.sim)

	for (i in 1:n.sim){

		rows <- round(runif(dim(lidar)[1], min=1, max=dim(lidar)[1]), 0)

		temp <- lidar[rows,]

		transition.matrix(temp, 56)

		transition.probability.matrix(transition)

		eigen.results <- popbio::eigen.analysis(transition.p)

		results[,i] <<- eigen.results$stable.stage*sum(transition)

	}

	print(time - Sys.time())

}



#### Read in data ####

  # 5 m resolution canopy height rasters from 2009 and 2019
    lidarCHM09 <- raster::raster("CHM09_AOI_res500.tif")
    lidarCHM19 <- raster::raster("CHM19_AOI_res500.tif")

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
       ForestPoly <- ForestPoly[!(ForestPoly$YEAR >= 1996),]
      # Remove tiny fragment of abandoned plantation
      ForestPoly <- ForestPoly[!(ForestPoly$LU00_UTM0_==68),]
    
    # Old growth forest area
    OldGrowthPoly <- ForestPoly[ForestPoly$DESCRIPTIO %in% c("Old-growth Forests","Ecological Reserve","Forested Swamps"),]
    
    # Secondary forest area
    SecondaryPoly <- ForestPoly[ForestPoly$DESCRIPTIO %in% c("Secondary Forests","Abandoned Agroforestry","Abandoned Plantation", "Selectively-logged Forests"),]
    
    Forest09 <- raster::mask(lidarCHM09,ForestPoly)
    Forest19 <- raster::mask(lidarCHM19,ForestPoly)
        
    OldGrowth09 <- raster::mask(Forest09, OldGrowthPoly)
    OldGrowth19 <- raster::mask(Forest19, OldGrowthPoly)
    
    Secondary09 <- raster::mask(Forest09, SecondaryPoly)
    Secondary19 <- raster::mask(Forest19, SecondaryPoly)
    
    # Make matrices for use in height transition matrix functions (defined above)
    lidar <- matrix(data=c(Forest09@data@values,Forest19@data@values),ncol=2,byrow = F)
    lidar_OldGrowth <- matrix(data=c(OldGrowth09@data@values,OldGrowth19@data@values),ncol=2,byrow = F)
    lidar_Secondary <- matrix(data=c(Secondary09@data@values,Secondary19@data@values),ncol=2,byrow = F)
    
#### Make height transition matrix and transition probability matrix (all forest) ####

  # Find maximum dimension of canopy ehgith transition matrix
    maxHt <- ceiling(max(range(lidar,na.rm=T)))
    
  # Time period 2
    transition <- transition.matrix(lidar, maxHt)
    transition.p <- transition.probability.matrix(transition)

#### Simulate steady state from a Dirichlet distribution (all forest) ####

# Define a function to simulate observations from a Dirichlet distribution

  rdiric <- function(alpha){
    out <- matrix(NA, ncol=length(alpha), nrow=1)
    
    for (i in 1:dim(out)[2]){
      
      out[1,i] <- rgamma(1, alpha[i], 1)
      
    }
    
    out <- out/sum(out)
    return(out)
    
  }
  
  
  sim <- array(data=NA, dim=c(dim(transition)[1],dim(transition)[1],10000))
  
    for(i in 1:dim(sim)[1]){
      
      alpha <- transition[,i]
      
      for(j in 1:dim(sim)[3]){
        
        sim[,i,j] <- rdiric(alpha + 1)
        
      }
      
    }

# Define an array to sample from the posterior distribtuion 10,000 times

  SS <- array(data=NA, dim = c(dim(transition)[1],10000))
  
  for(i in 1:dim(SS)[2]){
    
    transition.p.i <- transition.probability.matrix(sim[,,i])
    
    SS[,i] <- popbio::eigen.analysis(transition.p.i)$stable.stage*sum(transition)
    
  }

# Summarize confidence intervals of posterior distribution 

  SS_summary <- data.frame(index=1:dim(SS)[1],
                                 p2.5=NA,
                                 median=NA,
                                 p97.5=NA)
  
  SS_summary[,"p2.5"] <- apply(X=SS, MARGIN=1, 
                                     FUN=quantile, probs=c(0.025))
  SS_summary[,"median"] <- apply(X=SS, MARGIN=1, 
                                     FUN=quantile, probs=c(0.5))
  SS_summary[,"p97.5"] <- apply(X=SS, MARGIN=1, 
                                     FUN=quantile, probs=c(0.975))

##### Plot Fig. 3 and calculate other results stats #####
par(mfrow=c(1,1), mar=c(3,3,0,0), oma=c(2,2,2,2))
  plot(median~index, data=SS_summary,
       type="l",
       ylim=c(0,2800), col=NA, lty=3, cex.axis=1.2,
       las=1)
  
  polygon(x=c(SS_summary$index,rev(SS_summary$index)),
          y=c(SS_summary$p2.5,rev(SS_summary$p97.5)),
          col=adjustcolor("black",alpha.f = 0.4), border=NA)

  lines(x=1:dim(SS_summary)[1], y=apply(X=transition, MAR=1, FUN=sum), lwd=3, lty=2)
  lines(x=1:dim(SS_summary)[1], y=apply(X=transition, MAR=2, FUN=sum), lwd=3, lty=1)

  legend(x=27,y=2730,bty="n",
         c("Before","After"),
         lty=c(1,2), lwd=3, cex=1)
  legend(x=27,y=2930,bty="n",
         c("Projected steady state"),
         col=adjustcolor("black",alpha.f = 0.4),
         pch=15, cex=1)
  
  mtext("Height (m)",side=1,outer=T,line=-0.5, cex=1.2)
  mtext("Observations",side=2,outer=T,line=0.5, cex=1.2)

  # Find the predicted mean (and confidence intervals) canopy height for steady-state equilibrium
    predMeanHt <- rep(NA,dim(SS)[2])
    for(i in 1:dim(SS)[2]){
      predMeanHt[i] <- weighted.mean(x=0:(dim(SS)[1]-1),w=SS[,i])
    }
  
    mean(predMeanHt)
    quantile(predMeanHt,probs = c(0.025,0.975))

  # Calculate skewness and kurtosis of projected equilibirum 
    momentsSS <- data.frame(index=1:dim(SS)[2], skew = NA, kurtosis = NA)
    
    for(i in 1:dim(SS)[2]){
      htsVector <- c()
      for(j in 1:dim(SS)[1]){
        htsVector <- c(htsVector, rep(j,SS[j,i]))
      }
      momentsSS$skew[i] <- moments::skewness(htsVector)
      momentsSS$kurtosis[i] <- moments::kurtosis(htsVector)
    }
    
    mean(momentsSS$skew)
    quantile(momentsSS$skew,probs = c(0.025,0.975))
    
    mean(momentsSS$kurtosis)
    quantile(momentsSS$kurtosis,probs = c(0.025,0.975))
    
#### Calculate mean canopy height values ####
      mean(Forest09@data@values,na.rm=T)
      mean(Forest19@data@values,na.rm=T)
      
      mean(OldGrowth09@data@values,na.rm=T)
      mean(OldGrowth19@data@values,na.rm=T)
      
      mean(Secondary09@data@values,na.rm=T)
      mean(Secondary19@data@values,na.rm=T)        
#### Make height transition matrix and transition probability matrix (old growth only) ####

  # Find maximum dimension of canopy ehgith transition matrix
    maxHt <- ceiling(max(range(lidar_OldGrowth,na.rm=T)))
    
  # Time period 2
    transition <- transition.matrix(lidar_OldGrowth, maxHt)
    transition.p <- transition.probability.matrix(transition)

#### Simulate steady state from a Dirichlet distribution (old growth only) ####

sim <- array(data=NA, dim=c(dim(transition)[1],dim(transition)[1],10000))

  for(i in 1:dim(sim)[1]){
    
    alpha <- transition[,i]
    
    for(j in 1:dim(sim)[3]){
      
      sim[,i,j] <- rdiric(alpha + 1)
      
    }
    
  }


SS_OG <- array(data=NA, dim = c(dim(transition)[1],10000))

for(i in 1:dim(SS_OG)[2]){
  
  transition.p.i <- transition.probability.matrix(sim[,,i]) # shouldn't do anything
  
  SS_OG[,i] <- popbio::eigen.analysis(transition.p.i)$stable.stage*sum(transition)
  
}


SS_OG_summary <- data.frame(index=1:dim(SS_OG)[1],
                               p2.5=NA,
                               median=NA,
                               p97.5=NA)

SS_OG_summary[,"p2.5"] <- apply(X=SS_OG, MARGIN=1, 
                                   FUN=quantile, probs=c(0.025))
SS_OG_summary[,"median"] <- apply(X=SS_OG, MARGIN=1, 
                                   FUN=quantile, probs=c(0.5))
SS_OG_summary[,"p97.5"] <- apply(X=SS_OG, MARGIN=1, 
                                   FUN=quantile, probs=c(0.975))

#### Plot Fig. S7A and calculate other results stats #####
par(mfrow=c(1,1), mar=c(2,2,0,0), oma=c(2,2,2,0))
  plot(median~index, data=SS_OG_summary,
       type="l",
       ylim=c(0,1000), col=NA, lty=3, cex.axis=1.3)
  
  polygon(x=c(SS_OG_summary$index,rev(SS_OG_summary$index)),
          y=c(SS_OG_summary$p2.5,rev(SS_OG_summary$p97.5)),
          col=adjustcolor("black",alpha.f = 0.4), border=NA)

  lines(x=1:dim(SS_OG_summary)[1], y=apply(X=transition, MAR=1, FUN=sum), lwd=3, lty=2)
  lines(x=1:dim(SS_OG_summary)[1], y=apply(X=transition, MAR=2, FUN=sum), lwd=3, lty=1)

  legend(x=30,y=900,bty="n",
         c(2009,2019),
         lty=c(1,2), lwd=3, cex=1.2)
  legend(x=30,y=1000,bty="n",
         c("Projected steady state"),
         col=adjustcolor("black",alpha.f = 0.4),
         pch=15, cex=1.2)
  text("A. Old growth forest", x=0, y=980, cex=1.5, adj=0)
  
  mtext("Height (m)",side=1,outer=T,line=1, cex=1.3)
  mtext("Observations",side=2,outer=T,line=1, cex=1.3)
  
  # Find the predicted mean (and confidence intervals) canopy height for steady-state equilibrium
    predMeanHt_OG <- rep(NA,dim(SS_OG)[2])
    for(i in 1:dim(SS_OG)[2]){
      predMeanHt_OG[i] <- weighted.mean(x=0:(dim(SS_OG)[1]-1),w=SS_OG[,i])
    }
  
    mean(predMeanHt_OG)
    quantile(predMeanHt_OG,probs = c(0.025,0.975))


#### Make height transition matrix and transition probability matrix (secondary only) ####

  # Find maximum dimension of canopy ehgith transition matrix
    maxHt <- ceiling(max(range(lidar_Secondary,na.rm=T)))
    
  # Time period 2
    transition <- transition.matrix(lidar_Secondary, maxHt)
    transition.p <- transition.probability.matrix(transition)

#### Simulate steady state from a Dirichlet distribution (secondary only) ####

sim <- array(data=NA, dim=c(dim(transition)[1],dim(transition)[1],10000))

  for(i in 1:dim(sim)[1]){
    
    alpha <- transition[,i]
    
    for(j in 1:dim(sim)[3]){
      
      sim[,i,j] <- rdiric(alpha + 1)
      
    }
    
  }


SS_sec <- array(data=NA, dim = c(dim(transition)[1],10000))

for(i in 1:dim(SS_sec)[2]){
  
  transition.p.i <- transition.probability.matrix(sim[,,i]) # shouldn't do anything
  
  SS_sec[,i] <- popbio::eigen.analysis(transition.p.i)$stable.stage*sum(transition)
  
}


SS_sec_summary <- data.frame(index=1:dim(SS_sec)[1],
                               p2.5=NA,
                               median=NA,
                               p97.5=NA)

SS_sec_summary[,"p2.5"] <- apply(X=SS_sec, MARGIN=1, 
                                   FUN=quantile, probs=c(0.025))
SS_sec_summary[,"median"] <- apply(X=SS_sec, MARGIN=1, 
                                   FUN=quantile, probs=c(0.5))
SS_sec_summary[,"p97.5"] <- apply(X=SS_sec, MARGIN=1, 
                                   FUN=quantile, probs=c(0.975))

#### Plot Fig. S7B and calculate other results stats ####
par(mfrow=c(1,1), mar=c(2,2,0,0), oma=c(2,2,2,0))
  plot(median~index, data=SS_sec_summary,
       type="l",
       ylim=c(0,1800), col=NA, lty=3, cex.axis=1.3)
  
  polygon(x=c(SS_sec_summary$index,rev(SS_sec_summary$index)),
          y=c(SS_sec_summary$p2.5,rev(SS_sec_summary$p97.5)),
          col=adjustcolor("black",alpha.f = 0.4), border=NA)

  lines(x=1:dim(SS_sec_summary)[1], y=apply(X=transition, MAR=1, FUN=sum), lwd=3, lty=2)
  lines(x=1:dim(SS_sec_summary)[1], y=apply(X=transition, MAR=2, FUN=sum), lwd=3, lty=1)

  legend(x=30,y=1700,bty="n",
         c(2009,2019),
         lty=c(1,2), lwd=3, cex=1.2)
  legend(x=30,y=1800,bty="n",
         c("Projected steady state"),
         col=adjustcolor("black",alpha.f = 0.4),
         pch=15, cex=1.2)
  text("B. Secondary forest", x=0, y=1750, cex=1.5, adj=0)
  
  mtext("Height (m)",side=1,outer=T,line=1, cex=1.3)
  mtext("Observations",side=2,outer=T,line=1, cex=1.3)

  
  # Find the predicted mean canopy height for steady-state equilibrium
    predMeanHt_Sec <- rep(NA,dim(SS_sec)[2])
    for(i in 1:dim(SS_sec)[2]){
      predMeanHt_Sec[i] <- weighted.mean(x=0:(dim(SS_sec)[1]-1),w=SS_sec[,i])
    }
  
    mean(predMeanHt_Sec)
    quantile(predMeanHt_Sec,probs = c(0.025,0.975))


#### Plot Fig. S11: distribution of canopy height changes in old growth and secondary forests ####
    
    dheight <- density((lidar[,2]-lidar[,1])[!is.na(lidar[,1])])
    dheight_OldGrowth <- density((lidar_OldGrowth[,2]-lidar_OldGrowth[,1])[!is.na(lidar_OldGrowth[,1])])
    dheight_Secondary <- density((lidar_Secondary[,2]-lidar_Secondary[,1])[!is.na(lidar_Secondary[,1])])


    par(mfrow=c(1,1), mar=c(3,3,2,1), oma=c(2,2,0,0))

      plot(x=dheight_OldGrowth$x, y=dheight_OldGrowth$y, type="l",
           lwd=3,
           xlab=NA,ylab=NA, ylim=c(0,0.11),col=adjustcolor("black",alpha.f = 0.6))
      
      lines(x=dheight_Secondary$x, y=dheight_Secondary$y, type="l",
            col=adjustcolor("darkgrey",alpha.f = 0.8), lwd=3)
      abline(v=mean((lidar_OldGrowth[,2]-lidar_OldGrowth[,1]),na.rm=T), lwd=3, lty=3, col=adjustcolor("black",alpha.f = 0.6))
      abline(v=mean((lidar_Secondary[,2]-lidar_Secondary[,1]),na.rm=T), lwd=3, lty=3, col=adjustcolor("darkgrey",alpha.f = 0.8))
      legend(x=-50, y=0.11,
             c("Old-growth forest", "Secondary forest","Mean values"),
             lty=c(1,1,3),lwd=3,
             col=c(adjustcolor("black",alpha.f = 0.6),adjustcolor("darkgrey",alpha.f = 0.8),"black"),
             bty="n")
      mtext("Canopy height change (m)", side=1, line=0, cex=1.2, outer=T)
      mtext("Frequency", side=2, line=0, cex=1.2, outer=T)