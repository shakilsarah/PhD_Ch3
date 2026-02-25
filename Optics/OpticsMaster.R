#===========================================================================================================#
# BEPOM_OpticsMaster.R
# Background: Calculate absorption metrics for base-extracted
# Also refer to Optical Parameters Description table for Ch.1
## (Dropbox/Shakil_Thesisdrafts/Chapter1/Writing/Tables.xlsx)
# BEPOM paper: Brym et al. 2014. Optical and chemical characterization of base-extracted particulate organic matter in 
# coastal marine environments. Marine Chemistry, 162, 96 - 113.
# note, the master sheet is in wide form
# Created By: Sarah Shakil, April 18, 2019
# Last Updated: October 29, 2019
#===========================================================================================================#

##### ==============================Section i: Prep Workspace ==============================================

## Clear list 
list=rm(list=ls(all=TRUE))

#load necessary packages
library(readxl) # needed to read in excel sheets
#library(ggplot2) # create the graphing format to feed into ggplotly() function within plotly package
library(tidyr)# needed for separate
library(dplyr) # needed for select() function among others
#library(reshape2) # needed for reshaping data to useable format
#library(data.table) # needed for setting data tables
#library(plyr) #needed for ddply

## Set Working Directory
df <- "Optics/Data/"

## Call book code with necessary functions (e.g. pairs function; Zuur 2009, Mixed Effects Models and Extensions)
source("functions/HighstatLibV10.R")
# change above to functions/
## for calculating summary stats
se <- function(x) sqrt(var(x, na.rm=TRUE)/length(x[!is.na(x)])) 
## create standard error function for use below
lengthnona <- function(x) length(x[!is.na(x)]) 
## create length function that doesn't count na's for use below

##### ==============================Section 1: Optics data from Chapter 1 ===================================================

p <- read.csv(paste0(df, "masteroptics_p.csv"))

p <- p %>%
     select(sampcode,
            slumpsite, streamloc, sampdate, samptype, repnum, dilfac,
            SR_p = SR, atot250450_p = atot250450, a254naperian_p = a254naperian,
            HIX_p = HIX, 
            prcntC1_p = prcntC1, prcntC2_p = prcntC2, prcntC3_p = prcntC3, 
            prcntC4_p = prcntC4, prcntC5_p = prcntC5,
            normC1_p = normC1, normC2_p = normC2, normC3_p = normC3,
            normC4_p = normC4, normC5_p = normC5, 
            pkAB_p = pkAB, pkCB_p = pkCB, pkTC_p = pkTC, pkCM_p = pkCM,
            normA_p = normA,
            normB_p = normB, normC_p = normC, normE_p = normE, 
            normD_p = normD, normM_p = normM, normN_p = normN,
            normT_p = normT) %>%
  filter(!(p$slumpsite==0 & p$repnum==2) &
           p$samptype=="Sample")

p$sampcode <- as.character(p$sampcode)
p$sampcode[p$slumpsite==0] <- "0_NA_SS-0000_2017-07-17_Sample_1"

# for some reason site 0 was measured twice for abs and doesn't have EEMs
# didn't exceed 0.4 so keep undiluted one 
# note that's sample 2 for site 0
# also need need to insert sample code for site0

##### ==============================Section 2: Read in DOM Fluorescence Data ===================================================


# (2.1) Load component Fmax for each sample =====
# fmax of peaks
m3 <- read_excel(paste0(df, "drEEM_PARAFAC_model3.xls"), sheet="Model3Loading")
m3 <- na.omit(subset(m3, select=c("i", "Fmax1", "Fmax2", "Fmax3")))


# (2.2) Load Supplementary Calculations =====
scdom <- read.csv(paste0(df,"SuppCalc.csv"))

# (2.3) Merge by row#s =====
mod3 <- merge(m3, scdom, by=0)
mod3$Row.names <- NULL

# (2.4) Calculate peak ratios =====
attach(mod3)
mod3$`pkAB_d` <- A/B
mod3$`pkCB_d` <- C/B
mod3$`pkTC_d` <- mod3$`T`/C
mod3$`pkCM_d` <- C/M
mod3$`pkCA_d` <- C/A
detach(mod3)

# (2.5) Calculate %F for each component =====
attach(mod3)
sumFmax_d <- Fmax1 + Fmax2 + Fmax3 
mod3$`prcntC1_d` <- Fmax1/sumFmax_d
mod3$`prcntC2_d` <- Fmax2/sumFmax_d
mod3$`prcntC3_d` <- Fmax3/sumFmax_d
detach(mod3)

# (2.6) Calculate Fmax normalized for each component and peak (for correlation matrix table) ======

attach(mod3)
mod3$`normC1_d` <- Fmax1/Fmax
mod3$`normC2_d` <- Fmax2/Fmax
mod3$`normC3_d` <- Fmax3/Fmax
mod3$normA_d <- A/Fmax
mod3$normB_d <- B/Fmax
mod3$normC_d <- C/Fmax
mod3$normE_d <- E/Fmax
mod3$normD_d <- D/Fmax
mod3$normM_d <- M/Fmax
mod3$normN_d <- N/Fmax
mod3$normT_d <- mod3$`T`/Fmax
detach(mod3)

mod3 <- mod3 %>% select(sampcode,
                        FrI_d=FrI, BIX_d=BIX, HIX_d=HIX,
                        pkAB_d, pkCB_d, pkTC_d, pkCM_d,
                        pkCA_d, prcntC1_d, prcntC2_d, prcntC3_d,
                        normC1_d, normC2_d, normC3_d, normA_d, normB_d,
                        normC_d, normE_d, normD_d,
                        normM_d, normN_d, normT_d)

# (2.7) merge with particulate dataset =====
master <- merge(p, mod3, by="sampcode", all=TRUE)


##### ==============================Section 3: Read in DOM Absorbance Data ===================================================

# (3.1) Read in DOM Abs data =====
m4 <- read.csv(paste0(df, "2017DOMAbsIndices.csv"))
m4 <- m4 %>% select(site, date, samptype, SR_d=SR_m, 
                    a254dec_d=a254dec_m, 
                    atot250450_d=atot250450nap_m)

# (3.2) Merge with master =====

master <- merge(m4, master, 
                by.x=c("site", "date", "samptype"),
                by.y=c("slumpsite", "sampdate", "samptype"), all=TRUE)

# (1.8) Write to csv file =====

write.csv(master, paste0(df, "ThesisDrafts/Chapter3/Data/masteroptics.csv"))
