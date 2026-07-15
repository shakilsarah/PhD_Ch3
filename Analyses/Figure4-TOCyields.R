# ===========================================================================================================#
# Figure4.R
# By: Sarah Shakil
# Contact: shakil@ualberta.ca
# Background: Preliminary stats
# Last Updated: August 30 2020
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
library(patchwork) # Added for publication-ready multi-panel plots

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

## (1.1) Read in file ====================
d <- read.csv(paste0(df, "2017data.csv"))

#dall <- d

#d <- d %>% filter(campaign=="2017synoptic")

#sdist <- read_excel(paste0(wd, df, "slumpstreamdist.xlsx"))
# can't do this because not always a slump present

## (1.2) Select variables ====================
d <- d%>%
  select(site, 
         date,
         campaign,
         tssflux,
         pocflux,
         docflux,
         WatershedArea,
         strahlerstream, 
         percslump17act)

## (1.3) Calculate julian date ====================

d$JDay <- julian(as.Date(d$date), origin = as.Date("2016-12-31"))
# julian() sets origin=0 in jday counts, so start the day before Jan. 1st (i.e Dec. 31st of prev. year)

## (1.4) Insert missing mDO_perc from other day for sites 31 and 32 ====================

d$mDO_perc[d$site=="31"&d$JDay==210] <- d$mDO_perc[d$site=="31"&d$JDay==213] 
d$mDO_perc[d$site=="32"&d$JDay==210] <- d$mDO_perc[d$site=="32"&d$JDay==213] 

## (1.5) Remove extra sites not of interest ====================

d <- d %>%
  filter(d$site!="SC Outlet" & d$JDay!=213)

## (1.6) Fix site codes for RDA plotting ====================

d$site <- as.character(d$site)

d$site[d$site=="40-alt1"] <- "40"
d$site[d$site=="12-1"] <- "12"
d$site[d$site=="9-alt"] <- "9"

## (1.7) Add Y/N for slump impact  ====================

d$slumpYN <- "N"
d$slumpYN[d$percslump17act>0] <- "Y"
d$slumpYN[d$site==29] <- "N"

## (1.8) Calculate POC:TOC and toc flux  =====

d$poctoc <- (d$pocflux/(d$pocflux+d$docflux))*100
d$tocflux <- d$pocflux + d$docflux

## (1.9) Merge chapter 1 with sub-catchments and transects for graph =====
c1 <- read.csv(paste0(df, "PeelPlateau_RTSandstream_geochem.csv"))

c1 <- c1 %>%
  filter(sampletype=="Stream" & streamlocation!="IN") %>%
  select(slumpsite, streamlocation, samplingdate,
         WatershedArea_km2, Discharge_m3s, DOCmgL, POCmgL, TSSavg)

c1$code <- paste0(c1$slumpsite, c1$streamlocation, c1$samplingdate)
c1$tocflux_mgs <- (c1$POCmgL+c1$DOCmgL)*c1$Discharge_m3s*1000
c1$tssflux_mgs <- (c1$TSSavg)*c1$Discharge_m3s*1000
c1$docflux_mgs <- (c1$DOCmgL)*c1$Discharge_m3s*1000
c1$pocflux_mgs <- (c1$POCmgL)*c1$Discharge_m3s*1000
c1$poctoc <- (c1$pocflux_mgs/c1$tocflux_mgs)*100
c1$campaign <- "2015updn"

poctoc <- merge(d, c1, 
                by.x=c("site", "campaign",  "WatershedArea", 
                       "tssflux", "tocflux",
                       "pocflux", "docflux", "poctoc"),
                by.y=c("code", "campaign", "WatershedArea_km2", 
                       "tssflux_mgs", "tocflux_mgs", 
                       "pocflux_mgs", "docflux_mgs", "poctoc"),
                all=TRUE)

poctoc <- poctoc[!is.na(poctoc$campaign),]
poctoc$slumpYN[is.na(poctoc$slumpYN)&poctoc$streamlocation=="DN"] <-"Y"
poctoc$slumpYN[is.na(poctoc$slumpYN)&poctoc$streamlocation=="UP"] <-"N"
poctoc$campaign <- as.character(poctoc$campaign)
poctoc$campaign[poctoc$campaign=="2015updn"&poctoc$streamlocation=="UP"] <-"2015UP"
poctoc$campaign[poctoc$campaign=="2015updn"&poctoc$streamlocation=="DN"] <-"2015DN"
poctoc <- poctoc[!is.na(poctoc$pocflux),]
poctoc$tocyield <- poctoc$tocflux/poctoc$WatershedArea

##### ========== (2) Graph ==========================================================================

## (2.1) Theme ====
options(scipen = 100)

theme <-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.background = element_blank(),axis.line.x = element_line(colour="black"),
              axis.line.y = element_line(colour="black"),
              axis.text = element_text(colour="black",size=14),legend.background=element_blank(),
              text=element_text(size = 16),
            #  legend.title=element_blank(),
              legend.position = c(0,1),
              legend.direction="horizontal",
              legend.justification = c(0,1),
              aspect.ratio=1)

## (2.2) TOC yield ====
g1 <- ggplot() + 
  # Arrange by percentage so the highest impacted (darkest) points are plotted on top
  geom_point(data=poctoc %>% arrange(percslump17act),
             aes(x=WatershedArea, y=tocyield, shape=campaign, fill=percslump17act),
             size=4, alpha=0.9, color="black") + 
  
  scale_shape_manual(limits=c("2017synoptic", "2017transect"),
                     values=c(22, 21),
                     labels=c("subcatch.", "Stony main"),
                     name=NULL) +
  
  # CHANGED 1: Updated the legend label to %RTSactive, using subscript for "active"
  scale_fill_gradient(low="white", high="black", 
                      name=expression("%RTS"["active"]),
                      breaks = c(0, 1, 2, 3)) + 
  
  guides(fill=guide_colorbar(barwidth = unit(0.5, "cm"), 
                             barheight = unit(6, "cm"), 
                             frame.colour = "black",
                             ticks.colour = "black",
                             direction = "vertical"), 
         shape=guide_legend(override.aes = list(size=5))) +
  
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x=expression("Watershed Area (km"^2*")"), y=expression("TOC Yield (mg km"^-2*" s"^-1*")")) +
  
  theme + 
  theme(legend.position="none",
        plot.tag = element_text(size = 18, face = "bold"),
        # CHANGED 2: Added a full black border around the entire plot area
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))


## (2.3) PercPOC ====
g2 <- ggplot() + 
  # Arrange by percentage so the highest impacted (darkest) points are plotted on top
  geom_point(data=poctoc %>% arrange(percslump17act),
             aes(x=WatershedArea, y=poctoc, shape=campaign, fill=percslump17act),
             size=4, alpha=0.9, color="black") + 
  
  scale_shape_manual(limits=c("2017synoptic", "2017transect"),
                     values=c(22, 21),
                     labels=c("subcatch.", "Stony main"),
                     name=NULL) +
  
  # CHANGED 1: Updated the legend label to %RTSactive, using subscript for "active"
  scale_fill_gradient(low="white", high="black", 
                      name=expression("%RTS"["active"]),
                      breaks = c(0, 1, 2, 3)) + 
  
  guides(fill=guide_colorbar(barwidth = unit(0.5, "cm"), 
                             barheight = unit(6, "cm"), 
                             frame.colour = "black",
                             ticks.colour = "black",
                             direction = "vertical"), 
         shape=guide_legend(override.aes = list(size=5))) +
  
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x=expression("Watershed Area (km"^2*")"), y=expression("POC:TOC (%)")) +
  
  theme + 
  theme(legend.position="right", 
        legend.box="vertical", 
        legend.direction="vertical",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key = element_rect(fill = NA, color = NA),
        plot.tag = element_text(size = 18, face = "bold"),
        # CHANGED 2: Added a full black border around the entire plot area
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

##### ============================== Section 5: Export plots ===================================================

# Ensure the Figures directory exists
dfg <- "Figures/"
if(!dir.exists(dfg)) dir.create(dfg)

# We use the patchwork library to combine the toc and perpoc plots row-wise
# guides = "collect" will pull the shared legend to the top.
fig1 <- g1 + g2 + 
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') 

if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Figures/Fig4/")) dir.create("Figures/Fig4/", recursive = TRUE)

# Save Figure 4
ggsave(filename = "Figures/Fig4/Figure4_TOCyields.pdf", 
       plot = fig1, width = 10, height = 5, units = "in", dpi = 300)

ggsave("Figures/Fig4/Figure4_TOCyields.png", plot=fig1, width = 10, height = 5, units="in", dpi = 300)



