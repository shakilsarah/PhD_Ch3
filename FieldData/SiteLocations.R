# ===========================================================================================================#
# SiteLocations.R
# By: Sarah Shakil
# Contact: shakil@ualberta.ca
# Background: Preliminary stats
# Last Updated: October 20 2020
#===========================================================================================================#

##### (i) Workspace PREP ==============================================================================

## Clear list 
list=rm(list=ls(all=TRUE))

## Set working directory
df <- "FieldData/data/"

# load libraries
library(dplyr)
library(readxl)
library(tidyr)
library(data.table)

# load functions
se <- function(x) {
  sqrt(var(x, na.rm=TRUE)/length(x[!is.na(x)])) 
}

## create standard error function for use later
lengthnona <- function(x) {
  length(x[!is.na(x)])
}

##### ========== (1) DATA PREP =========================================================================

# (1.2) trans ====
site <- read_excel(paste0(df, "Daily_Site_Data.xlsx"))
site$year <- format(site$`Sampling Date`, format="%Y")

# average long and lat per site 
#siteav <- site %>% 
 # group_by(`Slump Site`, `Stream Location`, `Transect Location`, year) %>%
#  summarize(lat = mean(Latitude),
 #           latn = lengthnona(Latitude),
  #          latse = se(Latitude),
   #         long = mean(Longitude),
    #        longn = lengthnona(Longitude),
     #       longse = se(Longitude)
      #      ) # think it's better not to

# (1.3) set date to date ====

site$`Sampling Date` <- as.Date(site$`Sampling Date`)

# (1.3) scale to sites of interest ====

ch3sites <- site %>%
            filter((year=="2015" &
                    (`Slump Site`=="SE" | 
                     `Slump Site`=="SD") &
                      `Stream Location`!="HW"&
                      (`Sampling Date`=="2015-07-04" |
                      `Sampling Date`=="2015-07-06" |
                      `Sampling Date`=="2015-07-21" |
                      `Sampling Date`=="2015-07-22" |
                      `Sampling Date`=="2015-08-16" |
                      `Sampling Date`=="2015-08-21" )&
                      `Transect Location`!="TS-0025" &
                      `Transect Location`!="TS-0050") |
                   (year=="2016" &
                    `Stream Location`!="BL" &
                    `Stream Location`!="HW" &
                    `Stream Location`!="SN" &
                    `Stream Location`!="SL" &
                    `Stream Location`!="SMSC" &
                    `Stream Location`!="PP-01" &
                    `Stream Location`!="PP-02" &
                    `Stream Location`!="PP-03" &
                    `Stream Location`!="PP-04" &
                    `Stream Location`!="PP-05" &
                    `Stream Location`!="PP-06" &
                    `Slump Site`!="SA" &
                    `Slump Site`!="SB" &
                    `Slump Site`!="SB" & 
                    `Slump Site`!="SC" &
                    `Slump Site`!="SD")|
                   (year=="2017" &
                    `Stream Location` =="NA"))


##### ========== (2) Create second IDs =======================================================================

setnames(ch3sites, 
         old=c("Slump Site", "Stream Location", "Transect Location",
               "Sampling Date", "Latitude", "Longitude",
               "GPS Data QC Codes (see design view metadata)"),
         new=c("site", "loc", "trans",
               "date", "lat", "long",
               "gpsqc"))

# (2.1) set sampling campaign ====
ch3sites$campaign <- NA
#2015
ch3sites$campaign[ch3sites$year==2015] <- "2015shorttransects"
#2016
ch3sites$campaign[ch3sites$year==2016] <- "2016transect"
#2017
ch3sites$campaign[ch3sites$year==2017] <- "2017synoptic"
ch3sites$campaign[ch3sites$year==2017&
                  (ch3sites$site=="SC Outlet" |
                  ch3sites$site=="SC-1A" |
                  ch3sites$site=="SC-2B" |
                  ch3sites$site=="SC-3B" |
                  ch3sites$site=="SC-4B" |
                  ch3sites$site=="SC-5B" |
                  ch3sites$site=="SC-6B" |
                  ch3sites$site=="SC-7B" |
                  ch3sites$site=="SC-8B" |
                  ch3sites$site=="SC-9A" |
                  ch3sites$site=="DC-3" |
                  ch3sites$site=="DC-4")] <- "2017transect"

# (2.2) set transect number ====

#(2.2.1) 2015
ch3sites$trannum <- 1
ch3sites$trannum[ch3sites$date=="2015-07-04"&ch3sites$site=="SD"] <- 1
ch3sites$trannum[ch3sites$date=="2015-07-22"&ch3sites$site=="SD"] <- 2
ch3sites$trannum[ch3sites$date=="2015-08-16"&ch3sites$site=="SD"] <- 3

ch3sites$trannum[ch3sites$date=="2015-07-06"&ch3sites$site=="SE"] <- 1
ch3sites$trannum[ch3sites$date=="2015-07-21"&ch3sites$site=="SE"] <- 2
ch3sites$trannum[ch3sites$date=="2015-08-21"&ch3sites$site=="SE"] <- 3

#(2.2.2) 2016
unique(ch3sites$date[ch3sites$year==2016])

ch3sites$trannum[ch3sites$date=="2016-07-25"| 
                 ch3sites$date=="2016-07-26"|
                 ch3sites$date=="2016-07-30"] <- 2

#(2.2.3) 2017
unique(ch3sites$date[ch3sites$year==2017])

ch3sites$trannum[ch3sites$date=="2016-08-09"] <- 3

# (2.3) fix formatting to match previous ====

# (2.3.1) 2015 to match 
ch3sites$loc2 <- ch3sites$loc
ch3sites$loc2[ch3sites$year==2015 &
              ch3sites$trans=="TS-0100"] <- "DNTS100m"
ch3sites$loc2[ch3sites$year==2015 &
                ch3sites$trans=="TS-0200"] <- "DNTS200m"
ch3sites$loc2[ch3sites$year==2015 &
                ch3sites$trans=="TS-0300"] <- "DNTS300m"
ch3sites$loc2[ch3sites$year==2015 &
                ch3sites$trans=="TS-0400"] <- "DNTS400m"
ch3sites$loc2[ch3sites$year==2015 &
                ch3sites$trans=="TS-0600"] <- "DNTS600m"
ch3sites$loc2[ch3sites$year==2015 &
                ch3sites$trans=="TS-0800"] <- "DNTS800m"
ch3sites$loc2[ch3sites$year==2015 &
                ch3sites$trans=="TS-1000"] <- "DNTS1000m"

# (2.3.2) 2016 to match Hg

# don't think this is necessary right now so won't do

# (2.4) Remove columns not needed
setnames(ch3sites, old="")
ch3sites <- ch3sites%>%
            select(-gpsqc)

##### ========== (3) Print for ArcGIS =======================================================================

write.csv(ch3sites, paste0(df, "ch3sitecoord.csv"))