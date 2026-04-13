# ===========================================================================================================#
# OCsynopticfluxcontrols.R
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

##### ========== (3) Run RDA analysis ==========================================================================

library(vegan)

d <- na.omit(d)

Y <- d %>% select(logtocy,
                  logpocy,
                  logdocy)

pairs(Y, 
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)

#mod0p <- rda(Y ~ 1, data = d, scale=TRUE)

mod0p <- rda(Y ~ Condition(JDay), data = d, scale=TRUE)

#d$slumpYNcode <- 0
#d$slumpYNcode[d$slumpYN=="Y"] <- 1

mod1p <- rda(Y ~ delta18opermille +
               logcay +
               logwy + 
               logspower+
               logsslope +
               shale_log + 
               #    shale_log*pslump_log +
               col_log +
               #   col_log*pslump_log +
               pied_log + 
               #   pied_log*pslump_log +
               moraine_log +
               #   moraine_log*pslump_log +
               lake_log +
               loggpp +
               #   loggpp*pslump_log +
               forest_log +
               grass_log +
               lichen_log +
               logsoc +
               # loggldist +
               # loggldist*pslump_log + # too many terms, seems redundant to elev and not weighted across the watershed
               logslope +
               lograin +
               #  lograin*pslump_log +
               pslumpact_log +
               # pslumpall_log +
               logslumpcount +
               # logslumpstrahler +
               #JDay, 
               # slumpYNcode +
               Condition(JDay),
             data = d, scale=TRUE)

rda <- mod1p

plot(rda, scaling=2)
(R2adj <- RsquareAdj(rda)$r.squared)
(R2adj <- RsquareAdj(rda)$adj.r.squared)

## Global test of the RDA result
anova(rda, permutations = how(nperm = 5000))
## Tests of all canonical axes
anova(rda, by = "axis", permutations = how(nperm = 5000))

# Apply Kaiser-Guttman criterion to residual axes
# this is not really necessary (pg 221 numerical ecology)
rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)]

vif.cca(rda)

step.p.forward <-
  ordiR2step(mod0p, 
             scope = formula(mod1p), 
             direction = "forward", 
             permutations = how(nperm = 5000),
             R2permutations = 5000)


rdasimp <- rda(Y ~ Condition(JDay) + pslumpact_log + logwy + loggpp,
               data = d, scale=TRUE)

## Global test of the RDA result
anova(rdasimp, permutations = how(nperm = 5000))
## Tests of all canonical axes
anova(rdasimp, by = "axis", permutations = how(nperm = 5000))
## Tests of all terms
anova(rdasimp, by = "terms", permutations = how(nperm = 5000))

plot(rdasimp)

out = varpart(Y, ~JDay, ~logwy + loggpp + pslumpact_log,
              data = d, scale=TRUE)

out = varpart(Y, ~JDay, ~logwy, ~ loggpp, ~pslumpact_log,
              data = d, scale=TRUE)
plot(out)
out
(R2adj <- RsquareAdj(rdasimp)$r.squared)
(R2adj <- RsquareAdj(rdasimp)$adj.r.squared)

adonis2(Y ~ JDay + pslumpact_log + logwy + loggpp, data=d, permutations = 5000)

##### ========== (4) Plot RDA analysis ==========================================================================

# scores and figure

scor = scores(rdasimp, display=c("sp", "cn", "bp", "lc"), scaling=2) 

sites <- data.frame(scor$constraints)
sites$site <- d$site
sites$slumpYN <- d$slumpYN

species_centroids <- data.frame(scor$species)
species_centroids
species_centroids$species_names <- c("TOC","POC","DOC")

arrows <- data.frame(scor$biplot)
arrows$pf_names <- c("log(%ActiveSlumps+2)","log(WaterYield+2)","log(GPP+2)")
arrows

mult <- attributes(scores(rdasimp))$const

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

rdagraph <- ggplot(species_centroids, aes(x = RDA1, y= RDA2)) +
  geom_vline(xintercept = 0,linetype="dashed", colour="grey")+ylab("RDA2 (14.8%)")+
  geom_hline(yintercept = 0,linetype="dashed", colour="grey")+xlab("RDA1 (65.9%)")+
  geom_text(data = sites, 
            aes(x= RDA1, y = RDA2-0.15, label=site), 
            size=2.5, colour="grey40") + 
  #geom_point(colour="blue",
  #          size = 2, shape=17)+
  geom_text(data = species_centroids, 
            aes(label = species_names),
            colour = "blue", size = 4)+
  coord_cartesian(x = c(-2, 1.5), y = c(-1.5, 1.5))+
  geom_point(data = sites, 
             aes(fill=slumpYN, shape=slumpYN), 
             size=4, colour="white") +
  #  geom_text(data = sites, 
  #            aes(x= RDA1-0.15, y = RDA2, label=site), 
  #           size=4) +
  scale_shape_manual(limits= c("Y", "N"),
                     breaks= c("Y", "N"),
                     values= c(21, 22)) +
  scale_fill_manual(limits= c("Y", "N"),
                    breaks= c("Y", "N"),
                    values= c("pink", "Sky Blue")) +
  geom_segment(data = arrows,
               aes(x = 0, xend = (RDA1),
                   y = 0, yend = (RDA2)),
               arrow = arrow(length = unit(0.4, "cm")), colour = "grey30")+
  geom_text(data = arrows,
            aes(x= 1.2*RDA1, y = 1.2*RDA2, #we add 10% to the text to push it slightly out from arrows
                label = pf_names), #otherwise you could use hjust and vjust. I prefer this option
            size = 4,
            hjust = 0.5)+
  theme

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
        legend.key = element_rect(fill = NA, color = NA)) 

##### ========== (5) Linear model of TOC ==========================================================================

### (5.1) TOC ====
l <- lm(log10(tocyield) ~ percslump17act + scaledgpp + wateryield,
        data =darch)

l <- lm(log10(tocyield) ~ percslump17act + wateryield,
        data =darch)

darch$predtoc <- (darch$percslump17act*0.89315) + (darch$wateryield*0.01941) + 1.15933

ggplot() + geom_point(data=darch, 
                      aes(y=predtoc, x=log10(tocyield))) +
  geom_abline(slope=1) + theme

### (5.2) POC ====
l <- lm(log10(pocyield) ~ percslump17act + scaledgpp + wateryield,
        data =darch)

l <- lm(log10(pocyield) ~ percslump17act,
        data =darch)

darch$predtoc <- (darch$percslump17act*0.89315) + (darch$wateryield*0.01941) + 1.15933

ggplot() + geom_point(data=darch, 
                      aes(y=predtoc, x=log10(tocyield))) +
  geom_abline(slope=1) + theme


### (5.3) DOC ====
l <- lm(log10(docyield) ~ percslump17act + scaledgpp + wateryield,
        data =darch)

l <- lm(log10(docyield) ~ scaledgpp + wateryield,
        data =darch)

##### ============================== Section 6: Export plots ===================================================

# Ensure the Figures directory exists
dfg <- "Figures/"
if(!dir.exists(dfg)) dir.create(dfg)

## (6.1) Figure 1: TOC Yield and %POC across scales (Row 1) ======

# We use the patchwork library to combine the toc and perpoc plots row-wise
# guides = "collect" will pull the shared legend to the top.
fig1 <- toc + perpoc + 
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') & 
  theme(legend.position = "top",
        plot.tag = element_text(size = 18, face = "bold"))

# Save Figure 1
ggsave(filename = paste0(dfg, "Figure1_TOC_POC.pdf"), 
       plot = fig1, width = 10, height = 5, units = "in", dpi = 300)


## (6.2) Figure 2: RDA Graphs (Rows 2 & 3) ======

# Note: The current script only generates one single RDA plot object (`rdagraph`).
# To recreate the multi-panel grid seen in your reference image (c, d, e), you will 
# need to assign your other RDA plots to objects (e.g., rdagraph2, rdagraph3).

# IF you have 3 RDA plots, you can uncomment and use this layout:
# fig2 <- (rdagraph + rdagraph2) / (rdagraph3 + plot_spacer()) +
#   plot_annotation(tag_levels = list(c('c', 'd', 'e')), tag_prefix = '(', tag_suffix = ')') &
#   theme(plot.tag = element_text(size = 18, face = "bold"))

# For now, exporting the single rdagraph generated in the script above as Figure 2:
fig2 <- rdagraph + 
  labs(tag = "(c)") +
  theme(plot.tag = element_text(size = 18, face = "bold"))

# Save Figure 2
ggsave(filename = paste0(dfg, "Figure2_RDAs.pdf"), 
       plot = fig2, width = 9, height = 9, units = "in", dpi = 300)