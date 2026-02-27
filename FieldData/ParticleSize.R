# ===========================================================================================================#
# ParticleSize.R
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
source("functions/HighstatLibV10.R")

##### ========== (1) DATA PREP ==========================================================================

ps <- read_excel(paste0(df, "2017_FieldSeasonData-081617.xlsx"), sheet="Bed Material")

### (1.1) remove where pebble counts weren't done ==========

ps <- ps[!is.na(ps$`Hole Size >`) & 
          ps$`Hole Size >`!="na"&
           ps$`Hole Size >`!="NA",]

ps$`Hole Size >`[ps$`Hole Size >`==5] <- "45" # pretty sure 45 entered wrong
ps$`Hole Size >`[ps$`Hole Size >`==6] <- "64" # pretty sure 64 entered wrong
ps$`Hole Size >`[ps$`Hole Size >`==2.6] <- "22.6" # pretty sure 2.6 entered wrong
ps$`Hole Size >`[ps$`Hole Size >`==3.8] <- "2.8" # pretty sure 2.8 entered wrong
ps$`Hole Size >`[ps$`Hole Size >`==54] <- "64" # pretty sure 64 entered wrong
ps$`Hole Size >`[ps$`Hole Size >`==15] <- "16" # pretty sure 16 entered wrong
ps$`Hole Size >`[ps$`Hole Size >`==30] <- "32" # pretty sure 16 entered wrong
ps$`Hole Size >`[ps$`Hole Size >`==226] <- "22.6" # pretty sure 22.6 entered wrong
ps$`Hole Size >`[ps$`Hole Size >`==22] <- "22.6" # pretty sure 22.6 entered wrong

### (1.2) adjust site names ===========
ps$site <- NA
ps$site <- substr(ps$Site, start = 1, stop = 2)
ps$site[ps$site=="9-"] <- "9-alt"
ps$site[ps$site=="6-"] <- "6"
ps$site[ps$site=="0-"] <- "0"
ps$site[ps$site=="3-"] <- "3"
ps$site[ps$site=="7-"] <- "7"
ps$site[ps$site=="12"] <- "12-1"
ps$site[ps$site=="40"] <- "40-alt1"

### (1.3) Check is 100 counts were done per site ===========
pstotcount <- ps%>%
  group_by(site, Date) %>%
  select(site, Date, 
         holesizemm = `Hole Size >`) %>%
  summarize(ncount=lengthnona(holesizemm))

# sites without extra counts purposely have NAs and didn't complete 100 counts,
# not sure why, move one, use counts available (rather than discarding sites)
# but note the sites without 100
# sites: 48, 27, 0, 26

### (1.4) Aggregate counts per bin ===========
psbincount<- ps%>%
  select(site, Date, 
         geomorph = `Geomorph Unit (R=riffle, P=pool, B=bank)`,
         holesizemm = `Hole Size >`) %>%
  group_by(site, Date, holesizemm) %>%
  summarize(ncount=lengthnona(holesizemm))

#psbincount$holesizemm[psbincount$holesizemm=="<2"] <- 0.07/1000
#psbincount$holesizemm <- as.numeric(psbincount$holesizemm)
#psbincount$phi <- -(log2(psbincount$holesizemm))
psbincount$midpoint <- NA
psbincount$midpoint[psbincount$holesizemm=="<2"] <- 1
psbincount$midpoint[psbincount$holesizemm==2] <- (2+2.8)/2
psbincount$midpoint[psbincount$holesizemm==2.8] <- (2.8+4)/2 
psbincount$midpoint[psbincount$holesizemm==4] <- (4+5.6)/2 
psbincount$midpoint[psbincount$holesizemm==5.6] <- (5.6+8)/2 
psbincount$midpoint[psbincount$holesizemm==8] <- (8+11)/2 
psbincount$midpoint[psbincount$holesizemm==11] <- (11+16)/2 
psbincount$midpoint[psbincount$holesizemm==16] <- (16+22.6)/2 
psbincount$midpoint[psbincount$holesizemm==22.6] <- (22.6+32)/2 
psbincount$midpoint[psbincount$holesizemm==32] <- (32+45)/2 
psbincount$midpoint[psbincount$holesizemm==45] <- (45+64)/2 
psbincount$midpoint[psbincount$holesizemm==64] <- (64+90)/2 
psbincount$midpoint[psbincount$holesizemm==90] <- (90+128)/2 
psbincount$midpoint[psbincount$holesizemm==128] <- (128+180)/2 
psbincount$midpoint[psbincount$holesizemm==180] <- (180+256)/2 # though technically this is >256

psbincount <- psbincount[order(psbincount$midpoint),]

psbincount <- psbincount %>%
  group_by(site, Date) %>%
  mutate(cumsum=cumsum(ncount), sum=sum(ncount))

psbincount$cumperc <- (psbincount$cumsum/psbincount$sum)*100

### (1.5) Visualize graphs ===========
library(ggplot2)

ggplot() + geom_line(data=psbincount[psbincount$site==27,], aes(x=midpoint, y=cumperc, colour=site))
# check 48 and 7
### (1.6) Calculate D50 ===========

psd50 <- psbincount %>% 
  group_by(site, Date) %>% 
  mutate(D50 = approx(x = cumperc, y = midpoint, xout = 50)$y) %>% 
  select(site, Date, D50) %>%
  group_by(site, Date) %>%
  summarize(D50=mean(D50))

psand <- psbincount %>% 
  group_by(site, Date) %>% 
  mutate(psand = approx(x = midpoint, y = cumperc, xout = 2)$y) %>% 
  select(site, Date, psand) %>%
  group_by(site, Date) %>%
  summarize(psand=mean(psand))
  
psf <- merge(psd50, psand, all=TRUE)
psf$D50[psf$psand>50] <- 1
psf$psand[is.na(psf$psand)] <- 0

### (1.4) Write to csv ===========

write.csv(psf, paste0(df, "streambedparticlesize.csv"))
