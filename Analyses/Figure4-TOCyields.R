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

dall <- d

d$campaign[is.na(d$campaign)] <- "2017synoptic"

d <- d %>% filter(campaign=="2017synoptic")

#sdist <- read_excel(paste0(wd, df, "slumpstreamdist.xlsx"))
# can't do this because not always a slump present

## (1.2) Select variables ====================
d <- d%>%
  select(site, 
         date,
         tssyield,
         tocyield,
         pocyield,
         docyield,
         wateryield,
         Cayield,
         Nayield, 
         Mgyield, 
         Clyield,
         SO4yield,
         Sryield,
         Feyield,
         delta18opermille,
         slope,
         mDO_perc,
         mcond_uscm,
         mpH,
         streampower,
         D50, psand,
         WatershedArea, 
         percshale, # bedrock,
         colluvial_perc, piedmont_perc, # surficial geology
         alluvial_perc, # few sites with alluvial and low percentages when present
         bedrock_perc, # no bedrock perc
         fluvial_perc, # only1 fluvial perc
         glaciogenic_perc, #only2 glaciogenic perc
         organic_perc, # very little organic coverage
         moraine_perc, # surficial geology
         lakeperc, 
         #water_perc, #waterbodies includes river water bodies, want to specify lakes
         scaledgpp, #scalednpp, not npp because of strong correlation to gpp, write that in methods
         meanelev_m, meanslope_deg, # geormorph
         barrenland_perc, forest_perc, grassland_perc, #landcover
         lichenmoss_perc, shrubland_perc, # landcover
         wetland_perc, #essentially no wetland cover
         percslump17act, percslump17all,
         slumpacccount, strahlerimpactacc,
         wmeanSOCC_100CM,
         gldistkm,
         meanrough,
         RainTot24,
         RainTot48,
         RainTot72,
         RainTot96)

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

## (1.7) Add Y/N for slump impact for RDA plotting ====================

d$slumpYN <- "N"
d$slumpYN[d$percslump17act>0] <- "Y"
d$slumpYN[d$site==29] <- "N"

##### ========== (2) Conduct linear analysis ==========================================================================

# (2.1) Step 1 : Outliers assessment ====================

dotchart(d$tssyield)
dotchart(d$tocyield) # site 0 might need to be removed as an outlier
dotchart(d$pocyield) # site 0 might need to be removed as an outlier
dotchart(d$docyield)

# continue for all variables .... # decided to log10(x+1) transform all variables

# (2.2) Step 2: Transformations ====================
# for outlier data and removal of variables lacking data 

x=2

d$logtocy <- log10(d$tocyield+x)
d$logpocy <- log10(d$pocyield+x)
d$logdocy <- log10(d$docyield+x)

d$logwy <- log10(d$wateryield+x)

d$logcond <- log10(d$mcond_uscm+x)
d$logcay <- log10(d$Cayield+x)
d$lognay <- log10(d$Nayield+x)
d$logmgy <- log10(d$Mgyield+x)
d$logso4y <- log10(d$SO4yield+x)
d$logsry <- log10(d$Sryield+x)
d$logfey <- log10(d$Feyield+x)

d$logspower <- log10(d$streampower+x)
d$logsslope <- log10(d$slope+x)

d$shale_log <- log10(d$percshale+x)
d$col_log <- log10(d$colluvial_perc+x)
d$pied_log <- log10(d$piedmont_perc+x)
d$moraine_log <- log10(d$moraine_perc+x)
d$barren_log <- log10(d$barrenland_perc+x)
d$grass_log <-  log10(d$grassland_perc+x)
d$lichen_log <- log10(d$lichenmoss_perc+x)
d$lake_log <- log10(d$lakeperc+x)
d$forest_log <- log10(d$forest_perc+x)
d$shrub_log <- log10(d$shrubland_perc+x)

d$pslumpact_log <- log10(d$percslump17act+x)
d$pslumpall_log <- log10(d$percslump17all+x)
d$logslumpcount <- log10(d$slumpacccount+x)
d$logslumpstrahler <- log10(d$strahlerimpactacc+x)

#d$log18O <- log10(d$delta18opermille+x)

d$loggpp <- log10(d$scaledgpp+x) 

d$logelev <- log10(d$meanelev_m+x)
d$logslope <- log10(d$meanslope_deg+x) 
d$logrough <- log10(d$meanrough+x) 

d$logsoc <- log10(d$wmeanSOCC_100CM+x)
d$loggldist <- log10(d$gldistkm +x)
d$lograin <- log10(d$RainTot96+x)


# (2.2) Select initial data ====================
darch <- d
d <- d %>% 
  select(site, slumpYN,
         logtocy, logpocy, logdocy, JDay,
         logwy, logcond, logcay, lognay, logmgy, logso4y, logsry, logfey,
         delta18opermille,
         logspower, logsslope,
         shale_log,
         col_log, pied_log, moraine_log,
         lake_log, loggpp, 
         logelev, logslope, 
         barren_log, forest_log, grass_log, lichen_log,
         shrub_log, 
         pslumpact_log, pslumpall_log, 
         logslumpcount, logslumpstrahler,
         logsoc, loggldist, logrough, lograin)

# (2.4) Examine and reduce variables ====================

pairs(d[, c("logtocy", "logpocy", "logdocy", 
            "logcay", "lognay", "logmgy", "logso4y", "logsry", "logfey",
            "delta18opermille")], 
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)
# choose so4 as representative of all ions and poc as rep of toc and tss


pairs(d[, c(13:ncol(d))], #12 or 2
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)

pairs(d[, c("logslope", 
            "barren_log",
            "logrough", 
            "logelev" )], #12 or 2
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)

pairs(d[, c("pslumpact_log", 
            "pslumpall_log",
            "logslumpcount", 
            "logslumpstrahler" )], #12 or 2
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)

pairs(d[, c("loggpp", 
            "shrub_log")], #12 or 2
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)

pairs(d[, c("logtocy", "logpocy", "logdocy",
            "logslumpcount", "pslumpact_log")], #12 or 2
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)

##### ========== (4) Plot POC:TOC, TOC, and TSS across watershed scales ==========================================================================
# need to bring in UP, IN, DN data
# need to have separate values for UP, IN, DN, sub-catchments, and transects

## (4.1) Read in chapter 1 values =====
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

## (4.2) Select from dall =====

dall$slumpYN <- "N"
dall$slumpYN[dall$percslump17act>0] <- "Y"
dall$slumpYN[dall$site==29] <- "N"

dall$poctoc <- (dall$pocflux/(dall$pocflux+dall$docflux))*100
dall$tocflux <- dall$pocflux + dall$docflux
dall <- dall %>% select(site, campaign, WatershedArea, strahlerstream, slumpYN,
                        poctoc, pocflux, docflux, tocflux, tssflux)

## (4.3) Merge chapter 1 with sub-catchments and transects for graph =====

poctoc <- merge(dall, c1, 
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

## (4.4) Graph =====
options(scipen = 100)

theme <-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.background = element_blank(),axis.line.x = element_line(colour="black"),
              axis.line.y = element_line(colour="black"),
              axis.text = element_text(colour="black",size=14),legend.background=element_blank(),
              text=element_text(size = 16),
              legend.title=element_blank(),
              legend.position = c(0,1),
              legend.direction="horizontal",
              legend.justification = c(0,1),
              aspect.ratio=1)

perpoc <- ggplot() + 
  geom_point(data=poctoc,
             aes(x=strahlerstream, y=poctoc, shape=campaign, fill=slumpYN),
             size=3) +
  scale_shape_manual(limits=c("2017synoptic", "2017transect"),
                     values=c(22, 21),
                     labels=c("subcatch.", "Stony main")) +
  scale_fill_manual(values=c("SkyBlue", "Orange"), 
                    labels=c("NO RTS", "RTS")) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6)) +
  labs(x=expression("Strahler Stream Order"), y="% of TOC as POC") +
  theme + theme(legend.position="none") 

poctoc$tocyield <- poctoc$tocflux/poctoc$WatershedArea

toc <- ggplot() + 
  geom_point(data=poctoc,
             aes(x=strahlerstream, y=tocyield, shape=campaign, fill=slumpYN),
             size=3) +
  scale_shape_manual(limits=c("2017synoptic", "2017transect"),
                     values=c(22, 21),
                     labels=c("subcatch.", "Stony main")) +
  scale_fill_manual(values=c("SkyBlue", "Orange"), 
                    labels=c("NO RTS", "RTS")) +
  
  guides(fill=guide_legend(override.aes = list(shape=22))) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6)) +
  scale_y_log10() +
  labs(x=expression("Strahler Stream Order"), y=expression("TOC Yield (mg km"^-2*" s"^-1*")"))+
  theme + 
  theme(legend.position="top", 
        legend.box="horizontal", 
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key = element_rect(fill = NA, color = NA),
        plot.tag = element_text(size = 18, face = "bold")) 

##### ============================== Section 5: Export plots ===================================================

# Ensure the Figures directory exists
dfg <- "Figures/"
if(!dir.exists(dfg)) dir.create(dfg)

## (5.1) Figure 1: TOC Yield and %POC across scales (Row 1) ======

# We use the patchwork library to combine the toc and perpoc plots row-wise
# guides = "collect" will pull the shared legend to the top.
fig1 <- toc + perpoc + 
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') 

if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Figures/Fig4/")) dir.create("Figures/Fig4/", recursive = TRUE)

# Save Figure 4
ggsave(filename = "Figures/Fig4/Figure4_TOCyields.pdf", 
       plot = fig1, width = 10, height = 5, units = "in", dpi = 300)

ggsave("Figures/Fig4/Figure4_TOCyields.png", plot=fig1, width = 10, height = 5, units="in", dpi = 300)



##### ========== (6) Linear model of TOC ==========================================================================
# maybe this needs to be moved to the RDA so that 
### (6.1) TOC ====
l <- lm(log10(tocyield) ~ percslump17act + scaledgpp + wateryield,
        data =darch)

l <- lm(log10(tocyield) ~ percslump17act + wateryield,
        data =darch)

darch$predtoc <- (darch$percslump17act*0.89315) + (darch$wateryield*0.01941) + 1.15933

ggplot() + geom_point(data=darch, 
                      aes(y=predtoc, x=log10(tocyield))) +
  geom_abline(slope=1) + theme

### (6.2) POC ====
l <- lm(log10(pocyield) ~ percslump17act + scaledgpp + wateryield,
        data =darch)

l <- lm(log10(pocyield) ~ percslump17act,
        data =darch)

darch$predtoc <- (darch$percslump17act*0.89315) + (darch$wateryield*0.01941) + 1.15933

ggplot() + geom_point(data=darch, 
                      aes(y=predtoc, x=log10(tocyield))) +
  geom_abline(slope=1) + theme


### (6.3) DOC ====
l <- lm(log10(docyield) ~ percslump17act + scaledgpp + wateryield,
        data =darch)

l <- lm(log10(docyield) ~ scaledgpp + wateryield,
        data =darch)

