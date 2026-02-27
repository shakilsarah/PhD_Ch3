# ===========================================================================================================#
# Streamslope.R
# By: Sarah Shakil
# Contact: shakil@ualberta.ca
# Background: Slope, YSI, and field data
# Last Updated: October 20 2020
#===========================================================================================================#

##### (i) Workspace PREP ==============================================================================

## Clear list 
list=rm(list=ls(all=TRUE))

## Set working directory

df <- "FieldData/data/"

# load libraries
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(vegan) # for check of homogeneity of variances

# load functions
se <- function(x) {
  sqrt(var(x, na.rm=TRUE)/length(x[!is.na(x)])) 
}

## create standard error function for use later
lengthnona <- function(x) {
  length(x[!is.na(x)])
}

## Call book code with necessary functions (e.g. pairs function; Zuur 2009, Mixed Effects Models and Extensions)
#source(paste0(wd,"ThesisDrafts/Statistics/MixedEffectsModels/HighstatLibV10.R"))

##### ========== (1) Slope Data Prep ==========================================================================

### (1.1) Read in data tables ===========

slope <- read_excel(paste0(
  df, "2017_FieldSeasonData-081617.xlsx"), sheet="Slope and Morphometry")

### (1.2) Clean up slope data ===========

slope <- slope %>% 
  select(site=Site, date=Date, time=Time,
         slope_mperm) %>%
  group_by(site, date) %>%
  summarize(slope=mean(slope_mperm, na.rm=TRUE),
            lengthnona(slope_mperm))

## (1.2.1) after making sure no weird positive slopes, make all slopes absolute values

slope$slope <- abs(slope$slope)
slope$`lengthnona(slope_mperm)` <- NULL

## (1.2.2) fix site ids

slope$site <- substr(slope$site, start = 1, stop = 2)
slope$site[slope$site=="0-"] <- "0"
slope$site[slope$site=="1-"] <- "1"
slope$site[slope$site=="3-"] <- "3"
slope$site[slope$site=="5-"] <- "5"
slope$site[slope$site=="6-"] <- "6"
slope$site[slope$site=="7-"] <- "7"
slope$site[slope$site=="8-"] <- "8"
slope$site[slope$site=="9-"] <- "9-alt"
slope$site[slope$site=="12"] <- "12-1"
slope$site[slope$site=="40"] <- "40-alt1"

## (1.2.2) fix dates

slope$date[slope$date==71217] <- "2017-07-12"
slope$date[slope$date==71317] <- "2017-07-13"
slope$date[slope$date==71417] <- "2017-07-14"
slope$date[slope$date==71517] <- "2017-07-15"
slope$date[slope$date==71717] <- "2017-07-17"
slope$date[slope$date==72017] <- "2017-07-20"
slope$date[slope$date==72217] <- "2017-07-22"
slope$date[slope$date==72917] <- "2017-07-29"
slope$date[slope$date==73017] <- "2017-07-30"
slope$date[slope$date==80317] <- "2017-08-03"
slope$date[slope$date==80517] <- "2017-08-05"
slope$date[slope$date==80917] <- "2017-08-09"


##### ========== (2) Slope Data Prep ==========================================================================

### (2.1) Read in data tables ===========

ysi <- read_excel(paste0(
  wd, df, "2017_FieldSeasonData-081617.xlsx"), sheet="YSI")

### (2.2) Make character columns numeric ===========

samplings$Time <- strptime(samplings$Time, "%Y-%m-%d %H:%M:%S") 

ysiav <- ysi %>% 
  select(site, date, time=Time,
         lat=`Latitude (?N)`, long=`Longitude (?W)`, 
         elev_m=`Elevation (m)`, temp_degC=`Temp (?C)`,
         pressure_atm=`Pressure (atm)`, DO_perc=`DO (%)`,
         DO_mgL=`DO (mg/L)`, cond_uscm=`Cond (uS/cm)`,
         pH=pH
        ) %>%
  filter(site!="cal") %>%
  group_by(site, date) %>%
  summarize(time=mean(time, na.rm=TRUE),
            lat=mean(as.numeric(lat), na.rm=TRUE),
            long=mean(as.numeric(long), na.rm=TRUE), 
            melev_m=mean(as.numeric(elev_m), na.rm=TRUE), 
            mtemp_degC=mean(as.numeric(temp_degC), na.rm=TRUE), sdtemp=sd(as.numeric(temp_degC), na.rm=TRUE),
            mDO_perc=mean(as.numeric(DO_perc), na.rm=TRUE), sdDOp=sd(as.numeric(DO_perc), na.rm=TRUE),
            mDO_mgL=mean(as.numeric(DO_mgL), na.rm=TRUE), sdDOmgL=sd(as.numeric(DO_mgL), na.rm=TRUE),
            mcond_uscm=mean(as.numeric(cond_uscm), na.rm=TRUE), sdcond=sd(as.numeric(cond_uscm), na.rm=TRUE),
            mpH=mean(as.numeric(pH), na.rm=TRUE), sdpH=sd(as.numeric(pH), na.rm=TRUE)
           )

### (2.3) Fix 31 and 32 day later sampling for DO ===========

ysiav3132 <- ysiav %>% 
  select(site, date, time=Time,
         lat=`Latitude (?N)`, long=`Longitude (?W)`, 
         elev_m=`Elevation (m)`, temp_degC=`Temp (?C)`,
         pressure_atm=`Pressure (atm)`, DO_perc=`DO (%)`,
         DO_mgL=`DO (mg/L)`, cond_uscm=`Cond (uS/cm)`,
         pH=pH
  ) %>%
  filter(site!="cal" & (site==31|32)) %>%
  group_by(site) %>%
  summarize(lat=mean(as.numeric(lat), na.rm=TRUE), 
            long=mean(as.numeric(long), na.rm=TRUE), 
            elev_m=mean(as.numeric(elev_m), na.rm=TRUE), 
            mtemp_degC=mean(as.numeric(temp_degC), na.rm=TRUE), sdtemp=sd(as.numeric(temp_degC), na.rm=TRUE), cvtemp=sdtemp/mtemp_degC,
            mDOp =mean(as.numeric(DO_perc), na.rm=TRUE), sdDOp=sd(as.numeric(DO_perc), na.rm=TRUE), cvDOp=sdDOp/mDOp,
            mDOmgL=mean(as.numeric(DO_mgL), na.rm=TRUE), sdDOmgL=sd(as.numeric(DO_mgL), na.rm=TRUE), cvDOmgL=sdDOmgL/mDOmgL,
            mcond_uscm=mean(as.numeric(cond_uscm), na.rm=TRUE), sdcond=sd(as.numeric(cond_uscm), na.rm=TRUE), cvcond=sdcond/mcond_uscm,
            mpH=mean(as.numeric(pH), na.rm=TRUE), sdpH=sd(as.numeric(pH), na.rm=TRUE), cvpH=sdpH/mpH
  )
# coefficient of variation for site 31 and 32 less than 1% for temp, less than 3% for conductivity, and less than 9% for pH;
# so we assume this will mean minimal change in DO and use the average water quality values for the two days which means inserting DO for the missing day
# don't merge now, merge before data analysis so that things are organized for data set submission

##### ========== (3) Channel morphometry ==========================================================================

dis <- read_excel(paste0(df, "Discharge Step 6_ Extrapolate Any partial Q.xlsx"))
dis$`Sampling Date` <- as.Date(dis$`Sampling Date`)

dis <- dis %>%
  filter(`Sampling Date` >= "2017-01-01", 
         `Use or Reject` == "U"  ) %>%
  select(site=`Slump Site`,
         date=`Sampling Date`,
         avgvelocity_ms = `AvgOfVelocity (m/s)`,
         Q_m3s = `Full Est Q`,
         wettedwidth_cm = `Wetted Width (cm)`,
         bankfulwidth_cm = `Bankful Width (cm)`,
         maxdepth_m = maxdepth_m,
         avgdepth_m = avgdepth_m
  )

dis$bankfulwidth_m <- dis$bankfulwidth_cm/100
dis$bankfulwidth_cm <- NULL
dis$wettedwidth_m <- dis$wettedwidth_cm/100
dis$wettedwidth_cm <- NULL

##### ========== (4) Merge these supplementary datasets ==========================================================================
slope$date <- as.Date(slope$date)
ysiav$date <- as.Date(ysiav$date)
field2017 <- merge(slope, ysiav, by=c("site", "date"), all=TRUE)
field2017 <- merge(field2017, dis, by=c("site", "date"), all.x=TRUE)

##### ========== (5) Write File ==========================================================================

write.csv(field2017, paste0(df, "field2017datacleaned.csv"))
