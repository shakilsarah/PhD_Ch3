# ===========================================================================================================#
# Endmembertable.R
# By: Sarah Shakil
# Contact: shakil@ualberta.ca
# Background: Preliminary stats
# Last Updated: October 20 2020
#===========================================================================================================#

##### (i) Workspace PREP ==============================================================================

## Clear list 
list=rm(list=ls(all=TRUE))

## Set working directory

df <- "Data/"

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

## Call book code with necessary functions 
# (e.g. pairs function; Zuur 2009, Mixed Effects Models and Extensions)
source("functions/HighstatLibV10.R")

##### ========== (1) DATA PREP ==========================================================================

# (1.1) 2017 streambank and merge ====================

### 13C =====
sed13c <- read_excel(paste0(df, "Shakil_PeelHWSBPO13C.xlsx"))

sed13c <- read_excel(paste0(df, "Sediment_PerPOC.xlsx"))
sed13c$cat <- "SB"
sed13c$cat[sed13c$`Transect Location`=="PLE"] <- "PLE"
sed13c$cat[sed13c$`Transect Location`=="HOL"] <- "HOL"
sed13c$cat[sed13c$`Transect Location`=="OHO"] <- "UAL"
sed13c$cat[sed13c$`Transect Location`=="AHO"|
           sed13c$`Transect Location`=="BHO"] <- "LAL"

sed13c <- sed13c %>%
  filter((`Stream Location`=="NA"|
            `Stream Location`=="erosion"|
            `Stream Location`=="streambank" |
            `Transect Location`=="PLE"|
            `Transect Location`=="HOL"|
            `Transect Location`=="OHO"|
            `Transect Location`=="AHO"|
            `Transect Location`=="BHO")&
           `Sampling Date` > 2017 & 
           `Sample Type`=="Sample") %>%
  select(cat, `Slump Site`, `Stream Location`, 
         `Transect Location`, `Sampling Date`,
         `Sample Type`, d13C)

### %POC ====

sedperpoc <- read_excel(paste0(df, "Sediment_PerPOC.xlsx"))
sedperpoc$cat <- "SB"
sedperpoc$cat[sedperpoc$`Transect Location`=="PLE"] <- "PLE"
sedperpoc$cat[sedperpoc$`Transect Location`=="HOL"] <- "HOL"
sedperpoc$cat[sedperpoc$`Transect Location`=="OHO"] <- "UAL"
sedperpoc$cat[sedperpoc$`Transect Location`=="AHO"|
              sedperpoc$`Transect Location`=="BHO"] <- "LAL"

sedperpoc <- sedperpoc %>%
  filter((`Stream Location`=="NA"|
          `Stream Location`=="erosion"|
          `Stream Location`=="streambank" |
          `Transect Location`=="PLE"|
          `Transect Location`=="HOL"|
          `Transect Location`=="OHO"|
          `Transect Location`=="AHO"|
          `Transect Location`=="BHO")&
          `Sampling Date` > 2017 & 
          `Sample Type`=="Sample") %>%
  select(cat, `Slump Site`, `Stream Location`, 
         `Transect Location`, `Sampling Date`,
         `Sample Type`, PercPOC)

# (1.2) F14C of sed ====================
sed14c <- read_excel(paste0(df, "PO14C.xlsx"))

# Clean, rename, and filter out specific unwanted samples
sed14c <- sed14c %>%
  filter(MatCode=="A" & !is.na(F14C)) %>%
  # NEW: Drop the SL, SM/SC, and SN samples
  filter(!(`Transect Location` %in% c("SL", "SM/SC", "SN"))) %>%
  select(`Slump Site`, `Stream Location`,
         `Transect Location`, `Sampling Date`, 
         F14C_A = F14C, F14Cerror_A = `F14C error`)

# Map the categories exactly like 13C and %POC so they merge perfectly!
sed14c$cat <- "SB"
sed14c$cat[sed14c$`Transect Location`=="PLE"] <- "PLE"
sed14c$cat[sed14c$`Transect Location`=="HOL"] <- "HOL"
sed14c$cat[sed14c$`Transect Location`=="OHO"] <- "UAL"
sed14c$cat[sed14c$`Transect Location`=="AHO"|
             sed14c$`Transect Location`=="BHO"] <- "LAL"

# (1.3) Ch1 headwall and peri samples ====================
c1 <- read.csv(paste0(df, "PeelPlateau_RTSandstream_geochem.csv"))
c1$PO13C[c1$endmembertype=="Periphyton"]
c1$percPOC[c1$endmembertype=="Periphyton"]
c1$cat <- as.character(c1$headwallcat)
c1$cat[c1$endmembertype=="Periphyton"] <- "PERI"
c1 <- c1 %>%
      filter(sampletype=="Endmember") %>%
      select(cat, slumpsite, 
             streamlocation, samplingdate,
             F14C_A,
             F14Cerror_A, F14C_DIC, PO13C, percPOC)

c1$slumpsite <- as.character(c1$slumpsite)
c1$slumpsite[c1$slumpsite=="FM2"] <- "SB"
c1$slumpsite[c1$slumpsite=="FM3"] <- "SC"

c1$streamlocation <- "HW"

# (1.4) Merge files ====================

sed17 <- merge(sed13c, sedperpoc, all=TRUE)
sed <- merge(sed17, sed14c, all=TRUE)

sed$`Sampling Date` <- as.Date(sed$`Sampling Date`)
c1$samplingdate <- as.Date(c1$samplingdate)

em <- merge(c1, sed, 
            by.x=c("cat", "slumpsite", 
                   "streamlocation", "samplingdate",
                   "F14C_A", "F14Cerror_A", "PO13C", "percPOC"), 
            by.y=c("cat", "Slump Site", 
                   "Stream Location", "Sampling Date",
                   "F14C_A", "F14Cerror_A", "d13C", "PercPOC"),
            all=TRUE)

# for some reason PERI is getting exluded after the merge... input the 4 PO13C values manually, they are also reported in Shakil et al. (2020)

write.csv(em, paste0(df, "endmembersummaryall.csv"))

##### ========== (2) Calculate summaries ==========================================================================

# (2.1) Avg, Sd, Range ====================

emsum <- em %>%
  group_by(cat) %>%
  summarize(
    perpocav = mean(PercPOC, na.rm=TRUE),
    perpocsd = sd(PercPOC, na.rm = TRUE),
    perpocse = se(PercPOC),
    # Fix: Split range into explicit min and max
    perpocmin = min(PercPOC, na.rm = TRUE), 
    perpocmax = max(PercPOC, na.rm = TRUE),
    perpoc_n=lengthnona(PercPOC),
    
    d13cav = mean(d13C, na.rm=TRUE),
    d13csd = sd(d13C, na.rm = TRUE),
    d13cse = se(d13C),
    # Fix: Split range into explicit min and max
    d13cmin = min(d13C, na.rm = TRUE),
    d13cmax = max(d13C, na.rm = TRUE),
    d13c_n=lengthnona(d13C),
    
    F14Cav = mean(F14C_A, na.rm=TRUE),
    F14Csd = sd(F14C_A, na.rm = TRUE),
    F14Cse = se(F14C_A),
    # Fix: Split range into explicit min and max
    F14Cmin = min(F14C_A, na.rm = TRUE),
    F14Cmax = max(F14C_A, na.rm = TRUE),
    F14C_n=lengthnona(F14C_A)
  )

write.csv(emsum, paste0("Figures/", "endmembersummaryall.csv"))
