# ===========================================================================================================#
# 2017RDAoptics.R
# By: Sarah Shakil
# Contact: shakil@ualberta.ca
# Background: Preliminary stats
# Last Updated: August 30 2020
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

## Call book code with necessary functions (e.g. pairs function; Zuur 2009, Mixed Effects Models and Extensions)
source(paste0(wd,"ThesisDrafts/Statistics/MixedEffectsModels/HighstatLibV10.R"))

##### ========== (1) DATA PREP ==========================================================================

## (1.1) Read in file ====================
df <- "ThesisDrafts/Chapter3/Data/" 

d <- read.csv(paste0(wd, df, "2017data.csv"))

d <- d %>%
        select(site, date, chla=finalconc_ugL_dlcor, PO13C, POCTSSrat)

plot(d$chla~d$PO13C)

# write the note that, in majority of our sites, Chla concentration tended to be well below DL (1ug/L per 200 mL sample) 
# event at sites where sediment concentrations exceeded xxx mg/L
