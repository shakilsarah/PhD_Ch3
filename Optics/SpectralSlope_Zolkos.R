#===========================================================================================================#
# SpectralSlope.R
# Background: functions for calculating spectral slope from raw absorbance data
# Author: Scott Zolkos (zolkos@ualberta.ca) on March 17, 2018
# *Adapted by: 
  # 1) [include name and email here, if code is edited]
#===========================================================================================================#

  # Set working directory
    #*NOTE: 'wdir' should reference one folder level above raw data; raw data must be in folders names "Aqualog" and "OceanOptics"
    wdir <- "D:/Users/sarah/Dropbox/ThesisDrafts/Chapter3/Optics/DEMO/"
  
  # BACKGROUND FUNCTIONS
    # Horiba Aqualog: function to calculate baseline-corrected Naperian absorption coefficient for absorbance data from 200-800 nm
      Nap_abs_HA <- function(SampleID, rep){
        # Set working directory
          dir <- paste0(wdir,"Aqualog/")
        # Read in data ('absdf' variable = 'spectral slope dataframe')
          absdf <- read.csv(paste0(dir,SampleID,"-",rep,".csv"), header=T) # Read in data
          absdf <- absdf[3:nrow(absdf),] # Subset data, omitting first 2 rows
          absdf <- droplevels(subset(absdf, select=c("Wavelength","Abs"))) # Drop unused factor levels
          colnames(absdf) <- c("l","abs") # Specify column names for dataframe
          rownames(absdf) <- seq(length=nrow(absdf)) # Specify row names for dataframe
          absdf$l <- as.numeric(as.character(absdf$l)) # Store wavelength as numeric variable
          absdf$abs <- as.numeric(as.character(absdf$abs)) # Store absorbance as numeric variable
        
        # Order dataframe by increasing wavelength
          absdf <- absdf[order(absdf$l),]
          
        # Calculate baseline: average absorbance from 700-800 nm, following Helms et al. 2008 L&O p958
          avgbsln <- mean(subset(absdf$abs, absdf$l >= 700 & absdf$l <= 800))
        
        # Calculate Naperian absorption coefficient, using abs corrected for baseline, following Helms et al. 2008 L&O p958
          l <- 0.01 # pathlength (in m); e.g., if 1 cm cuvette, l = 0.01
          absdf$aNap <- ((absdf$abs-avgbsln)*2.303)/l
        
        return(absdf) # Return absorbance values as a variable
      }
        
    # Ocean Optics: function to calculate baseline-corrected Naperian absorption coefficient for absorbance data from 200-800 nm
      Nap_abs_OO <- function(SampleID, rep){
        # Set working directory
          dir <- dir <- paste0(wdir,"OceanOptics/")
        # Read in data with 
          lengths <- read.delim(paste0(dir,"WAVELENGTHS_FLMS012201_09-13-46-184.txt"), skip=14)[1] # Skip first 14 lines of metadata in data files
          names(lengths) <- "l"
        # Read in data ('ssdf' variable = 'spectral slope dataframe')
          ssdf <- read.delim(paste0(dir,SampleID,"-",rep,".txt"), skip=14) # Skip first 14 lines of metadata in data files
          names(ssdf) <- c("l2","abs")
        # Merge absorbance measurements to (non-rounded) Ocean Optics wavelengths
          ssdf <- droplevels(subset(cbind(lengths,ssdf), select=c("l","abs")))
        # Calculate absorbance at each non-fractional wavelength as the average of the fractions within ± 0.5; for instance, abs at 180 = average(abs(179.9), abs(180.2))
          wavelengths <- (seq(179, 881, 1)) # Create list of wavelengths
          absdf <- as.data.frame(matrix(nrow=length(wavelengths), ncol=2))
          colnames(absdf) <- c("l","abs")
          j=1
          for(i in wavelengths){
            absdf[j,1] <- i
            absdf[j,2] <- mean(subset(ssdf$abs, ssdf$l > i-0.5 & ssdf$l < i+0.49))
            j <- j+1
          }
          absdf <- subset(absdf,  absdf$l >= 200 & absdf$l <= 800)
        
        # Calculate baseline: average absorbance from 700-800 nm, following Helms et al. 2008 L&O p958
          avgbsln <- mean(subset(absdf$abs, absdf$l >= 700 & absdf$l <= 800))
        
        # Calculate Naperian absorption coefficient, using abs corrected for baseline, following Helms et al. 2008 L&O p958
          l <- 0.01 # pathlength (in m); e.g., if 1 cm cuvette, l = 0.01
          absdf$aNap <- ((absdf$abs-avgbsln)*2.303)/l
        
        return(absdf) # Return absorbance values as a variable
      }
  
    # Function to calculate baseline-corrected Naperian absorption coefficient for absorbance data from 200-800 nm
      Nap_abs <- function(Instrument, SampleID,rep){
        # Store variables
          SampleID <- SampleID
          rep <- rep
        
        # Calculate absorbance
          if(Instrument=="aqualog"){
            absdf <- Nap_abs_HA(SampleID=SampleID, rep=rep)
            }
          else(
            absdf <- Nap_abs_OO(SampleID=SampleID, rep=rep)
            )
        
        return(absdf)
      }
      
    # Function to calculate mean Naperian abs. coeff. using absorbance data input from two sample replicates
      mean_Nap_abs <- function(Instrument, SampleID){
        # Calculate mean slope ratio
          Instrument <- Instrument
          sr1 <- Nap_abs(Instrument,SampleID,"1")
          sr2 <- Nap_abs(Instrument,SampleID,"2")
          
        # Merge dataframe one and two
          srdf <- merge(x=sr1, y=sr2, by="l", all.x=T)
          names(srdf) <- c("l","abs1","aNap1","abs2","aNap2")
          
        # Calculate average absorbance
          absdf <- as.data.frame(matrix(nrow=nrow(srdf),ncol=2))
          names(absdf) <- c("l","avgabs")
          j=1
          for(i in 1:nrow(absdf)){
            absdf[j,1] <- srdf[i,1]
            absdf[j,2] <- mean(c(srdf[i,2],srdf[i,4]))
            j <- j+1
          }
        
        # Calculate average Naperian absorption coefficient
          aNapdf <- as.data.frame(matrix(nrow=nrow(srdf),ncol=2))
          names(aNapdf) <- c("l","avgaNap")
          j=1
          for(i in 1:nrow(aNapdf)){
            aNapdf[j,1] <- srdf[i,1]
            aNapdf[j,2] <- mean(c(srdf[i,3],srdf[i,5]))
            j <- j+1
          }
        
        # Write dataframe with wavelength, mean absorbance, and mean Naperian absorption coefficient
          MeanSR <- merge(x=absdf, y=aNapdf, by="l", all.x=T)
        
        return(MeanSR)
      }

  # SLOPE RATIO FUNCTIONS (note: 'calc_slope_ratio' can be for-looped for batch SR calculations)
    # Function to calculate slope ratio (SR), given an input SampleID
      calc_slope_ratio <- function(Instrument, SampleID){
        # Read in mean, baseline-corrected absorbance data
          Instrument <- Instrument
          MeanSR <- mean_Nap_abs(Instrument, SampleID)
        
        # Calculate slope ratio
        # Fit a linear model on log-transformed absorption spectra for wavelengths 275-295
          a1 <- subset(MeanSR, MeanSR$l >= 275 & MeanSR$l <= 295)
          s1 <- summary(lm(log(a1$avgaNap) ~ a1$l))$coeff[1] # Extract slope from lm
          a2 <- subset(MeanSR, MeanSR$l >= 350 & MeanSR$l <= 400)
          s2 <- summary(lm(log(a2$avgaNap) ~ a2$l))$coeff[1] # Extract slope from lm
          sr <- s1/s2
          
          return(sr)
        }
  
    # Function to plot slope ratio (SR), given an input SampleID
      plot_slope_ratio <- function(Instrument, SampleID, ylim){
        # Input y-axis limit for graph
          ylim <- ylim
        
        # Read in mean, baseline-corrected absorbance data
          Instrument <- Instrument
          MeanSR <- mean_Nap_abs(Instrument, SampleID)
        
          sr <- round(calc_slope_ratio(Instrument, SampleID),3)
        
        # Plot absorbance data, highlighting absorbance in regions used to calculate slope ratio
          par(mar=c(4.5,4.5,1,1))
          plot(x=MeanSR$l, y=MeanSR$avgabs, ylim=c(0,ylim), lwd=1, ylab="Absorbance", xlab="Wavelength (nm)", main=paste0(SampleID, " ", Instrument))
          lgnd <- bquote(italic(S)[R] == .(round(sr,3)))
          legend("topright", legend=lgnd, bty="n")
          points(x=MeanSR$l[MeanSR$l > 275 & MeanSR$l < 295 ], y=MeanSR$avgabs[MeanSR$l > 275 & MeanSR$l < 295 ], col="orange")
          points(x=MeanSR$l[MeanSR$l > 350 & MeanSR$l < 400 ], y=MeanSR$avgabs[MeanSR$l > 350 & MeanSR$l < 400 ], col="orange")
      }
    
    # Diagnostics
      plot_diagnostics <- function(SampleID){
        # Store data
          absdf_aqualog <- mean_Nap_abs("aqualog", SampleID)
          absdf_OO <- mean_Nap_abs("oceanoptics", SampleID)
        
        # Create plot
          plot(absdf_aqualog$avgabs~absdf_aqualog$l, cex=1.2, ylab=expression(Absorbance~(cm^-1)), xlab="Wavelength (nm)", main="Diagnostics: Aqualog vs. Ocean Optics abs.", type="n")
          points(absdf_aqualog$avgabs~absdf_aqualog$l, cex=1.2)
          points(absdf_OO$avgabs~absdf_OO$l, col="red", cex=0.4)
          legend("topright", legend=c("Aqualog","Ocean Optics"), pch=c(1,1), col=c("black","red"), bty="n")
      }
    
  # DEMO
    # Input parameters
      sample <- "BB01-070116-1101" 
      #sample <- "BB07-062816-1514"
      Rep <- "1" # Replicate # (for calculating 'Nap_abs' for one sample only)
      ylim <- 0.8 # Adjust 'ylim', as needed
    # Calculate baseline-corrected Naperian absorption coefficient for ONE sample
      Nap_abs(Instrument="aqualog", SampleID=sample, rep=Rep)
    # Calculate mean baseline-corrected Naperian absorption coefficient for REPLICATE samples
      mean_Nap_abs(Instrument="aqualog", SampleID=sample)
    # Calculate slope ratio from function "mean_Nap_abs"
      calc_slope_ratio(Instrument="aqualog", SampleID=sample)
      calc_slope_ratio(Instrument="oceanoptics", SampleID=sample)
    # Visualize slope ratio from function "mean_Nap_abs"
      plot_slope_ratio(Instrument="aqualog", SampleID=sample, ylim=ylim)
      plot_slope_ratio(Instrument="oceanoptics", SampleID=sample, ylim=ylim)
    # Compare Aqualog and Ocean Optics absorbance for this sample
      plot_diagnostics(SampleID=sample)
    