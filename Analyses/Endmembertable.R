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

# MAC: wd <- "~/Dropbox/"

# PC: 
wd <- "D:/Users/sarah/Dropbox/"

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
source(paste0(wd,"ThesisDrafts/Statistics/MixedEffectsModels/HighstatLibV10.R"))

##### ========== (1) DATA PREP ==========================================================================

# (1.1) 2017 streambank and merge ====================

### 13C =====
df <- "ThesisDrafts/Chapter3/Data/"
sed13c <- read_excel(paste0(wd,df, "Shakil_PeelHWSBPO13C.xlsx"))

sed13c <- read_excel(paste0(wd,df, "Sediment_PerPOC.xlsx"))
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

sedperpoc <- read_excel(paste0(wd,df, "Sediment_PerPOC.xlsx"))
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
df2 <- "ThesisDrafts/Chapter3/Data/"
sed14c <- read_excel(paste0(wd, df2, "20162017POCPO14CTrans.xlsx"))

sed14c <- sed14c %>%
          filter(MatCode=="A") %>%
          select(`Slump Site`, `Stream Location`,
                 `Transect Location`, `Sampling Date`, 
                 F14C_A = F14C, F14Cerror_A = `F14C error`)
sed14c$cat <- "SB"
          
# (1.3) Ch1 headwall and peri samples ====================
dfc1 <- "ThesisDrafts/Chapter1/Writing/Database/filestosubmit/"
c1 <- read.csv(paste0(wd, dfc1, "PeelPlateau_RTSandstream_geochem.csv"))
c1$PO13C[c1$endmembertype=="Periphyton"]
c1$cat <- as.character(c1$headwallcat)
c1$cat[c1$endmembertype=="Periphyton"] <- "PERI"
c1 <- c1 %>%
      filter(sampletype=="Endmember" &
             !is.na(F14C_A | F14C_DIC)) %>%
      select(cat, slumpsite, 
             streamlocation, samplingdate,
             F14C_A,
             F14Cerror_A, F14C_DIC)

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
                   "F14C_A", "F14Cerror_A"), 
            by.y=c("cat", "Slump Site", 
                   "Stream Location", "Sampling Date",
                   "F14C_A", "F14Cerror_A"),
            all=TRUE)

write.csv(em,paste0(wd, "ThesisDrafts/Chapter3/Graphs/endmembersummaryall.csv") )
##### ========== (2) Calculate summaries ==========================================================================

# (2.1) Avg, Sd, Range ====================

emsum <- em %>%
      group_by(cat) %>%
      summarize(perpocav=mean(PercPOC, na.rm=TRUE),
                perpocsd=sd(PercPOC, na.rm = TRUE),
                perpocse=se(PercPOC),
                perpocrange=range(PercPOC, na.rm = TRUE),
                
                d13cav=mean(d13C, na.rm=TRUE),
                d13csd=sd(d13C, na.rm = TRUE),
                d13cse=se(d13C),
                d13crange=range(d13C, na.rm = TRUE),
                
                F14Cav=mean(F14C_A, na.rm=TRUE),
                F14Csd=sd(F14C_A, na.rm = TRUE),
                F14Cse=se(F14C_A),
                F14Crange=range(F14C_A, na.rm = TRUE)
       )

write.csv(emsum, paste0(wd, "ThesisDrafts/Chapter3/Graphs/endmembersummary.csv"))
