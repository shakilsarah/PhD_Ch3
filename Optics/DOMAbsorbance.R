#===========================================================================================================#
# DOMAbsorbance.R
# Background: Calculate absorption metrics for base-extracted particulate organic matter following 
# Brym et al. 2014. Optical and chemical characterization of base-extracted particulate organic matter in 
# coastal marine environments. Marine Chemistry, 162, 96 - 113.
# Created By: Sarah Shakil, April 18, 2019
# Last Updated: April 18, 2019
#===========================================================================================================#

##### ==============================Section i: Prep Workspace ==============================================

## Clear list 
list=rm(list=ls(all=TRUE))

#load necessary packages
library(readxl) # needed to read in excel sheets
library(reshape2) # needed for reshaping data to useable format
library(data.table) # needed for melt function to get absorbance data into long form rather than wide form
library(tidyverse)
library(zoo) # needed for na.spline to interpolate missing data (BEPOM specific)
library(plotly) # creates html document of graphing file for interactive use of the graph (zoom in, zoom out, click on samples)
library(ggplot2) # create the graphing format to feed into ggplotly() function within plotly package
library(dplyr) # runs calculations on data subsets (e.g. per sample); easier on your computer than R forloop
library(stringr)

## Set Working Directory
#MAC: wd <- "~/Dropbox/"
#PC:
wd <- "D:/Users/sarah/Dropbox/"
datafolder <- "Thesisdrafts/Chapter3/Optics/2017/"

##### ==============================Section 1: Prep Data ===================================================

## === 1.1 Read in Sample Info and merge =====

# ACCESS Database codes (for matching with other datasets)
sampleinfo <- read_excel(paste0(wd, datafolder, "Ch3_OceanOptics.xlsx"), sheet=1)
sampleinfo <- sampleinfo %>% select(-integrationtime_s)
# hold on, maybe merge this in later

# fileids (extracted from filenames in folder)
filenames <- read.csv(paste0(wd, datafolder, "abs2017sampleIDs.csv"))
filenames <- filenames %>% select(filenames=x)

filenames$type <- "Sample"

filenames$type <- ifelse(str_detect(filenames$filenames, "DIblank", negate = FALSE), "Analytical Blank", filenames$type)
filenames$type <- ifelse(str_detect(filenames$filenames, "Blank", negate = FALSE), "Field Blank", filenames$type)

filenames$type[filenames$filenames=="Absorbance_14-39-37-529"|
              filenames$filenames=="Absorbance_14-39-41-452"] <- "Uknown"

## === 1.2 Read in Absorbance Data  =====
# obtain data with raw absorbance from sheet 2 in excel file
a <- read.csv(paste0(wd, datafolder, "abs2017.csv"))
a <- a %>% select(-X)

## === 1.1 Merge with filenames =====

a <- merge(a, filenames, by.x="SampleIDs", by.y="filenames")

## === 1.2 Interpolate absorbance data to obtain fixed wavelength values =====

# summarizing attempt =====

a$Wavelengthnew <- round(a$Wavelength, digits=0)

a <- a%>%
 group_by(SampleIDs, Wavelengthnew, type)

anew <- a %>%
 group_by(SampleIDs, Wavelengthnew, type) %>%
 summarize(abs = mean(Absorbance, na.rm=TRUE)) %>%
 select(SampleIDs, Wavelength=Wavelengthnew, Absorbance=abs, type) %>%
filter(between(Wavelength, 240, 800))

#loess predict method completely changed abs slopes, not good for SR =====
#w_new <- seq(from=240, to=800, by=1)

#test <- a %>%
#  group_by(type, SampleIDs) %>%
#  arrange(SampleIDs, Wavelength) %>%
#  do(mod= loess(Absorbance~Wavelength, .)) %>%
#  summarize(SampleIDs, type, w=w_new, a=predict(mod, w), .groups = "keep")


#anew <- test%>%
#  select(SampleIDs, type, Wavelength=w, Absorbance=a)



#test <- a %>%
#  group_by(SampleIDs) %>%
#  arrange(SampleIDs, Wavelength) %>%
#  summarize(w = w_new, predint = apply_fun(w_new, a$Wavelength, a$Absorbance))


# from Stack Exchange =====

#apply_fun <- function(w_new, w_measured, a_measured) {
#  predfunc <- approxfun(w_measured,a_measured)
#  predfunc(w_new)
#}

#test <- a %>%
#  group_by(SampleIDs) %>%
#  arrange(SampleIDs, Wavelength) %>%
#  summarize(w = w_new, predint = apply_fun(w_new, a$Wavelength, a$Absorbance))



## === 1.2 QC Absorbance Data =====

# plot all absorbance data
ggplotly(ggplot() + 
  geom_line(data=anew[anew$SampleIDs=="45SS03BU-alt2_080317_1a_Absorbance_15-28-53-595" |
                        anew$SampleIDs=="45SS03BU-alt2_080317_1b_Absorbance_15-33-17-583",], aes(x=Wavelength, y=Absorbance, group=SampleIDs, colour=type)) +
  geom_point(data=a[a$SampleIDs=="45SS03BU-alt2_080317_1a_Absorbance_15-28-53-595" |
                      a$SampleIDs=="45SS03BU-alt2_080317_1b_Absorbance_15-33-17-583",], aes(x=Wavelength, y=Absorbance, group=SampleIDs, colour=type)) +
  scale_x_continuous(limits=c(200, 800),
    breaks=c(200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800))
)

a <- anew
## === 1.2.1 Inspecting Blanks and samples visually =====


#plot Field blanks
ggplotly(ggplot() + geom_line(data=a[a$type=="Field Blank",],
                              aes(x=Wavelength, y=Absorbance, group=SampleIDs)) +
           scale_x_continuous(limits=c(200, 800),
                              breaks=c(200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800))
)

# looks like field blanks are lower than lowest sample, good!

## === 1.2.2 Remove optically dense samples ===

# some samples have abs just below 1.5 at 240, think this is okay but may need to 
# double check with Suzanne, code for removal is archived at the end

## === 1.2.3 Dilution Correct ===

# samples weren't diluted
#absallr0.4$absorbanceODdilcor <- (absallr0.4$absorbanceOD)*(1/absallr0.4$dilfac)


#Wavelengthnew <- seq(from=240.000, to=800.000, by=1.000)

#apply_fun <- function(Wavelengthnew, Wmeasured, Absmeasured) {
#  predfunc <- approxfun(Wmeasured,Absmeasured)
#  predfunc(Wavelengthnew)
#}

#test <- a %>%
 # group_by(SampleIDs) %>%
#  arrange(SampleIDs, Wavelength) %>%
#  mutate(predint = apply_fun(Wavelengthnew, a$Wavelength, a$Absorbance))

## === 1.2.5 Baseline adjustment =====

# (a) Baseline adjustment to raw absorbance 
# CDOM absorbance is assumed to be zero above 700 nm
# therefore, avg sample absorbance between 700 and 800 nm was subtracted from 
# the spectrum to correct for offsets due to:
# (a) instrument baseline drift, (b) temperature, (c) scattering, and (d) refractive effects
# See spectral corrections, pg 958 in Helms et al. 2008

# Table with original data (note, groupby can have issues if plyr is masking out dplyr)
absline <- a %>%
     group_by(SampleIDs) %>%
     filter(between(Wavelength, 700, 800)) %>%
     summarize(bsline = mean(Absorbance), sdbsline =sd(Absorbance))

a <- merge(a, absline, by=c("SampleIDs"))
a$absbscor <- (a$Absorbance)-(a$bsline)

# check new data against corrected to make sure no code errors
ggplotly(ggplot() + geom_line(data=a, aes(x=Wavelength, y=Absorbance, group=SampleIDs), colour="blue") +
  geom_line(data=a, aes(x=Wavelength, y=absbscor, group=SampleIDs), colour="red"))

plot(a$absbscor~a$Absorbance) # perfect :) 
# note, slightly shifted upwards

## === 1.2.5 Calculate Naperian Absorption coefficient =====

# (b) Calculate Naperian Absorption coefficients (Helms et al. 2008, pg 958)

l=0.01 # store path length, for these measurements, a path length of 10 mm = 0.01m was used

a$abscoeffnap <- ((a$absbscor)*2.303)/l

# (b) Calculate Decadic asborbance coefficients for SUVA (Weishaar) =====

l=0.01 # store path length, for these measurements, a path length of 10 mm = 0.01m was used

a$abscoeffdec <- ((a$absbscor))/l

## === 1.2 Do corrections for spectral slope as in Helms et al. 2008 =====

# (b) Natural log transform absorbance data (i.e use the linear method for ln transform)

# note, this doesn't need to be done separately because it can be done within dplyr

a$ln.abscoeffnap <- log(a$abscoeffnap) # Table with original data

# (c) remove negative values that are associated with 0 absorbance with blanks

# determine sample IDS with negative values
samples_negabs <- a$SampleIDs[a$Wavelength>=275 & a$Wavelength<=400 & a$abscoeffstd<0]
samples_negabs <- as.data.frame(samples_negabs)
colnames(samples_negabs) <- c("SampleIDs")

# use these codes to remove these groups from table
a <- anti_join(a, samples_negabs, by=c("SampleIDs"))

# look at spectral slopes after log transformed
ggplotly(ggplot() + 
    geom_line(data=a[a$type=="Sample",], aes(x=Wavelength, y=ln.abscoeffnap, group=SampleIDs)) +
    labs(y="Naperian Absorption Coefficient (a, m^-1)", title="zoomin, blue lines are interpolated values")+
    scale_x_continuous(limits=c(200, 600),breaks=c(200, 250, 300, 350, 400, 450, 500, 550, 600)))

##### ==============================Section 2: Calculate Indices ==============================================

# (2.1) calculate slope at 275-295 =====
slope275295 <- a %>%
  group_by(SampleIDs) %>% # keeping all the sample info in :)
  filter(between(Wavelength, 275, 295) & type=="Sample") %>% # here, sample is different from the ACCESS database notations and thus includes samples duplicates (really just differentiates from blanks and less code this way)
  do(model275295=lm(ln.abscoeffnap~Wavelength, data=.)) %>%
  mutate(coef275295=coef(model275295)["Wavelength"],rsq275295=summary(model275295)$r.squared)

# (2.2) calculate slope at 350-400 =====
slope350400 <- a %>%
  group_by(SampleIDs) %>% # keeping all the sample info in :)
  filter(between(Wavelength, 350, 400) & type=="Sample") %>% # here, sample is different from the ACCESS database notations and thus includes samples duplicates (really just differentiates from blanks and less code this way)
  do(model350400=lm(ln.abscoeffnap~Wavelength, data=.)) %>%
  mutate(coef350400=coef(model350400)["Wavelength"],rsq350400=summary(model350400)$r.squared)

# (2.3) calculate Decadic absorption coefficient at 254 nm for SUVA 254 (later calculated when merge with DOC conc data) =====
a254dec <- a %>%
  group_by(SampleIDs) %>% # keeping all the sample info in :)
  filter(Wavelength == 254) %>%
  summarize(a254dec=sum(abscoeffdec))

# (2.4) calculate the total absorbance under the curve (following Helms et al. 2008) =====
atot250450nap <- a %>%
  group_by(SampleIDs) %>%
  filter(between(Wavelength, 250, 450)) %>%
  summarize(atot250450nap=sum(abscoeffnap))

# (2.5) merge =====

MyMerge  <- function(x, y){
  df <- merge(x, y,all.x= TRUE, all.y= TRUE)
  return(df)
}

DOMAbs <- Reduce(MyMerge, list(slope275295, slope350400, a254dec, atot250450nap))

# (2.6) multiply by SR -1 =====
#because Helms uses the mathematical practice for exponential decay curves, but don't use absolute value since directions may differ    

DOMAbs$coef275295 <- DOMAbs$coef275295*-1
DOMAbs$coef350400 <- DOMAbs$coef350400*-1

# (2.7) calculate slope ratios =====
DOMAbs$SR <- DOMAbs$coef275295/DOMAbs$coef350400   

# (2.8) Merge with ACCESS sample info data =====

DOMAbs <- DOMAbs %>% select(-model275295, -model350400)

DOMAbs <- merge(DOMAbs, sampleinfo, by.x = "SampleIDs", by.y="ABS_file")

DOMAbs <- DOMAbs %>%
  select(SampleIDs, SR, a254dec, atot250450nap, site=`Slump Site`,
         loc = `Stream Location`, trans=`Transect Location`,
         date=`Sampling Date`, samptype=`Sample Type`, rep = Rep_num) %>%
  group_by(site, date, samptype) %>%
  summarize(SR_m=mean(SR, na.rm=TRUE), SR_sd=sd(SR, na.rm=TRUE), SR_c = (SR_sd/SR_m)*100,
            a254dec_m=mean(a254dec, na.rm=TRUE), a254dec_sd=sd(a254dec, na.rm=TRUE), a254dec_c = (a254dec_sd/a254dec_m)*100,
            atot250450nap_m=mean(atot250450nap, na.rm=TRUE), 
            atot250450nap_sd=sd(atot250450nap, na.rm=TRUE),
            atot250450nap_c = (atot250450nap_sd/atot250450nap_m)*100)

# note, high error in SR with low absorbance samples ===== 
# may want to note this in thesis, not an issue with 254
plot(DOMAbs$SR_c~DOMAbs$a254dec_m)
plot(DOMAbs$SR_c[DOMAbs$SR_c<100]~DOMAbs$a254dec_m[DOMAbs$SR_c<100]) # only 1 sample is >100, at >800, most likely an error, will remove from analysis

# (2.9) Export data to csv file =====
write.csv(DOMAbs, paste0(wd, datafolder, "2017DOMAbsIndices.csv"))


# Archived code for removing opticall dense samples =====

#samples_to_remove <- absall$samplecode[absall$Wavelength==240 & absall$absorbanceOD>?]
#length(samples_to_remove) # 30 samples need to be removed
#samples_to_remove <- as.data.frame(samples_to_remove)
#colnames(samples_to_remove) <- c("samplecode")

# use remove codes to remove groups from absall
# this can be done using an anti-join

#absallr0.4 <- anti_join(absall, samples_to_remove, by=c("samplecode")) 
# removes any samples in the sample_to_remove table

# check that the difference in length is equal to the length of samples to remove
#((length(absall$absorbanceOD)/(800-239))-(length(absallr0.4$absorbanceOD)/(800-239)))/length(samples_to_remove$samplecode)
# equal to 1, yay!!

# write csv file to save progress
#write.csv(absallr0.4, paste0(wd, datafolder, "absallremove0.4240.csv"))

# view data again
#ggplot() + geom_line(data=absallr0.4[absallr0.4$samptype=="Sample",], 
#                             aes(x=Wavelength, y=absorbanceOD, group=samplecode))


# Archived code for non-linear abs =====

# (b.2) calculate slope at 350-400 using nonlinear methods

#slope350400 <- a %>%
#  group_by(SampleIDs) %>% # keeping all the sample info in :)
#  filter(between(Wavelength, 350, 400) & type=="Sample") %>% # here, sample is different from the ACCESS database notations and thus includes samples duplicates (really just differentiates from blanks and less code this way)
# do(nls(abscoeff ~ a * exp(-b * (Wavelength-440)), start = list(a=0.1, b=0.0012), data=.)) %>% # model350400int=lm(abscoeffint~Wavelength, data=.)) %>%
# mutate(coef350400=coef(model350400)["Wavelength"],rsq350400=summary(model350400)$r.squared,
#         coef350400int=coef(model350400int)["Wavelength"],rsq350400int=summary(model350400int)$r.squared
#  )
