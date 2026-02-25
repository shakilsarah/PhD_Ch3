# ===========================================================================================================#
# datevrain.R
# By: Sarah Shakil
# Contact: shakil@ualberta.ca
# Background: Preliminary stats
# Last Updated: Dec 3 2021
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

##### ========== (1) Read inData ==========================================================================

# (1.1) Rain data ===
df <- "ThesisDrafts/Chapter3/WeatherStationData/"
rain <- read_excel(
  paste0(wd, df, "Peel Plateau DAILY data 2017 Jan 1_Dec 31.xlsx"), 
  sheet="forR")
rain$date <- as.Date(rain$TIMESTAMP)

# (1.2) Peel Q data ===
df <- "ThesisDrafts/Chapter3/peeldataset/loadestfilesformat/"
peelq <- read.delim(paste0(wd, df, 
                         "qloadest.txt"), 
                  header = TRUE, sep = "\t", skip=5,
                  row.names=NULL)
peelq <- peelq[as.Date(peelq$datetime) > "2017-01-01", ]
peelq$date <- as.Date(peelq$datetime)

# (1.3) Sites ===
df <- "ThesisDrafts/Chapter3/Data/"
site <- read.csv(paste0(wd, df, "2017data.csv"))
site <- site[site$campaign=="2017synoptic" |
             site$campaign=="2017transect",
             c("site", "date", "campaign")]
site$date <- as.Date(site$date)
range(site$date[!is.na(site$date)])

##### ========== (2) Graph ==========================================================================

ggplot() +
  geom_bar(data=rain, 
           aes(x=date, y=Rainfall_Tot*100),
           stat="identity", colour="blue") +
  geom_line(data=peelq, 
            aes(x=date, y=(X01_00060_00003/35.3147)),
            size=1) +
  geom_rect( 
    aes(
      xmin=as.Date(c("2017-07-10")),
      xmax=as.Date(c("2017-08-09")),
      ymin=-Inf,
      ymax=Inf),
      fill="grey",alpha=0.5
  ) + 
  scale_y_continuous(
  name = expression("Peel River Discharge (m"^3*"s"^-1*")"),
  sec.axis = sec_axis(trans=~./100, name="Total Rainfall per day")) +
  labs(x="Date") +
  theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line.x = element_line(colour="black"),
                 axis.line.y = element_line(colour="black"),
                 axis.text = element_text(colour="black",size=12),
                 legend.background=element_blank(),
                 text=element_text(size = 12),
                 legend.title=element_blank(),
                 legend.position = c(0,1),
                 legend.direction="horizontal",
                 legend.justification = c(0,1),
                 # aspect.ratio=1
                 legend.key=element_blank(),
                 panel.border = element_rect(colour = "black", fill=NA, size=1),
                 legend.spacing = unit(0, "mm"),
                 legend.key.size = unit(0.001, "cm"),
                 plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
  
rain2 <- rain[rain$date>="2017-07-10",]
rain2 <- rain2[rain2$date<="2017-08-09",]

peelq2 <- peelq[peelq$date>="2017-07-10",]
peelq2 <- peelq2[peelq2$date<="2017-08-09",]

ggplot() +
  geom_bar(data=rain2, 
           aes(x=date, y=Rainfall_Tot*100),
           stat="identity", colour="blue") +
  geom_line(data=peelq2, 
            aes(x=date, y=(X01_00060_00003/35.3147)),
            size=1) +
  geom_rect( 
    aes(
      xmin=as.Date(c("2017-07-10")),
      xmax=as.Date(c("2017-08-09")),
      ymin=-Inf,
      ymax=Inf),
    fill="grey",alpha=0.5
  ) + 
  scale_y_continuous(
    name = expression("Peel River Discharge (m"^3*"s"^-1*")"),
    sec.axis = sec_axis(trans=~./100, name="Total Rainfall per day")) +
  labs(x="Date") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour="black"),
        axis.line.y = element_line(colour="black"),
        axis.text = element_text(colour="black",size=12),
        legend.background=element_blank(),
        text=element_text(size = 12),
        legend.title=element_blank(),
        legend.position = c(0,1),
        legend.direction="horizontal",
        legend.justification = c(0,1),
        # aspect.ratio=1
        legend.key=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.spacing = unit(0, "mm"),
        legend.key.size = unit(0.001, "cm"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

