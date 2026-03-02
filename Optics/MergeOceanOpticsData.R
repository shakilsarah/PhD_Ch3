#===========================================================================================================#
# MergeOceanOpticsData.R
# Background: functions for calculating spectral slope from ocean optics code
# Author:Scott Zolkos (zolkos@ualberta.ca) on March 17, 2018
# *Adapted by: 
  # 1) Sarah Shakil (shakil@ualberta.ca)
#===========================================================================================================#

library(stringr)
library(plyr)
library(purrr)

df <- "Data/OceanOptics/"


 read_OO <- function(filename){
   read.delim(filename, skip=14)
   }# Skip first 14 lines of metadata in data files


##1. identify files to read in
filenames <- list.files(paste0(wd, df), pattern="*.txt", full.names=TRUE)
filelist <- lapply(filenames, read_OO)
SampleIDs <- stringr::str_remove(str_remove(filenames, paste0(wd,df)), ".txt")
names(filelist) <- SampleIDs
filelist <- mapply(cbind, filelist, SampleIDs, SIMPLIFY=F)
colnames <- c("Wavelength","Absorbance", "SampleIDs")
filelist <- lapply(filelist, setNames, colnames)


abs2017 <- plyr::join_all(filelist, by = c("Wavelength","Absorbance", "SampleIDs"),
               type = "full", match = "all")

write.csv(abs2017, paste0("Optics/data/abs2017.csv"))
write.csv(SampleIDs, paste0("Optics/data/abs2017sampleIDs.csv" ))
