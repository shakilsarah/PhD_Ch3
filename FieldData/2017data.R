# ===========================================================================================================#
# OCfluxdataset.R
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
source(paste0(wd,"ThesisDrafts/Statistics/MixedEffectsModels/HighstatLibV10.R"))

##### ========== (1) DATA PREP ==========================================================================

# (1.1) OC and TSS fluxes ====================

tss <- read_excel(paste0(df, "2015-17 TSS Combine.xlsx"))
poc <- read_excel(paste0(df, "20162017POCPO14CTrans.xlsx"))
doc <- read_excel(paste0(df, "2017DOC.xlsx"))
dis <- read_excel(paste0(df, "Discharge Step 6_ Extrapolate Any partial Q.xlsx"))
disscott <- read_excel(paste0(df, "Zolkos Q area coords compiled.xlsx"), sheet="transectforR")
docscott <- read_excel(paste0(df, "c4_sc_doc.xlsx"), sheet="c4_sc_doc")

#match discott to poc code
disscott$loc <- "NA"
disscott$trans <- "SS-0000"
docscott$site <- paste0(docscott$Site, "-", docscott$Loc_ORIG)
docscott$loc <- "NA"
docscott$trans <- "SS-0000"
docscott <- docscott[!is.na(docscott$DOCuM),]
docscott$DOCmgL <- signif((docscott$DOCuM*12.0107)/1000, digits = 5)

# there are poc duplicates because of multiple calibration reps for poc, 
# average (will be same as individual numbers because identical)

# (1.1.1) POC ==========

poc <- poc %>% 
      filter(`Transect Location`!="PERI" &
             `Sample Type`=="Sample" & 
             (is.na(`MatCode`) | `MatCode`!="DOC")) %>%
      select(site = `Slump Site`,
             loc = `Stream Location`,
             date = `Sampling Date`,
             trans = `Transect Location`,
             POCmgL = `POC (mg/L)`,
             PO13C = d13C,
             POCTSSrat = `Per POC`)  %>%
      group_by(site, loc, date, trans) %>%
      summarize(POCmgL = mean(POCmgL, na.rm=TRUE),
                PO13C = mean(PO13C, na.rm=TRUE),
                POCTSSrat = mean(POCTSSrat, na.rm=TRUE))

# since we did filter by sample. sample length is 1 for all
# checked with lengthnona during summarize, removed after confirmation

# (1.1.2) DOC ==========
doc <- doc %>% 
  filter(`Sample Type` == "Sample") %>%
  select(site = `Slump Site`,
         loc = `Stream Location`,
         date = `Sampling Date`, 
         trans = `Transect Location`,
         DOCmgL = NPOC) 

docscott <- docscott %>% select(site,
                           loc,
                           date = sampdate, 
                           trans,
                           DOCmgL)

doc <- merge(doc, docscott, all=TRUE)

doc$DOCmgL <- as.numeric(doc$DOCmgL)

# (1.1.3) TSS ==========
tss <- tss %>% 
  select(site = `Slump or Confluence Site`,
         loc = `Stream Location`,
         date = `Sampling Date`,
         trans = `Transect Location`,
         tssmgL = `TSS (mg/L, Avg)`) %>%
  group_by(site, loc, date, trans) %>%
  summarize(tssmgL = mean(tssmgL, na.rm=TRUE))

# (1.1.3) Discharge ==========
dis <- dis %>% select(site = `Slump Site`,
                      loc = `Stream Location`,
                      date = `Sampling Date`, 
                      trans = `Transect Location`,
                      dism3s = `Full Est Q`)

disscott <- disscott %>% select(site = slumpsite, 
                               loc,
                               date = Date,
                               trans,
                               dism3s = `Full Q (m3/s)`)

dis$date2 <- as.Date(dis$date)
dis$dism3s[dis$site=="SC Outlet" & dis$date2=="2017-08-07"] <- NA
dis$date2 <- NULL

dis <- merge(dis, disscott, all=TRUE)

# (1.1.4) Merge ==========
oc <- merge(tss, poc, all=TRUE)
oc <- merge(oc, doc, all=TRUE)
oc <- merge(oc, dis)

## (1.2) gis data ====================
gis <- read.csv(paste0(df, "watmaster.csv"))

pairs(gis[,c(8, 13:ncol(gis))])

gis <- gis%>%
       select(site, loc, trans, date, #id variables
              lat, long, year, campaign, trannum, # useful indexes
              WatershedArea, 
              percshale, # bedrock,
              colluvial_perc, piedmont_perc, # surficial geology
              alluvial_perc, # few sites with alluvial and low percentages when present
              bedrock_perc, # no bedrock perc
              fluvial_perc, # only1 fluvial perc
              glaciogenic_perc, #only2 glaciogenic perc
              organic_perc, # very little organic coverage
              moraine_perc, # surficial geology
              lakeperc, water_perc, #waterbodies
              scaledgpp, #scalednpp, not npp because of strong correlation to gpp, write that in methods
              meanelev_m, meanslope_deg, # geormorph
              barrenland_perc, forest_perc, grassland_perc, #landcover
              lichenmoss_perc, shrubland_perc, # landcover
              wetland_perc, #essentially no wetland cover
              percslump17act, percslump17all, 
              slumpacccount, strahlerimpactacc, #slump
              wmeanSOCC_100CM, # SOC
              gldistkm,
              meanrough, strahlerstream
              )

## (1.3) Merge OC with gis ====================
oc$date <- as.Date(oc$date)
gis$date <- as.Date(gis$date)
t <- merge(oc, gis)
# only interested in 2017 year
t <- t[t$campaign=="2017synoptic" | t$campaign=="2017transect",]

## (1.4) Fill in 0s where necessary for gis spatial data ====================
# fill NAs with 0
t$colluvial_perc[is.na(t$colluvial_perc)] <- 0
t$piedmont_perc[is.na(t$piedmont_perc)] <- 0
t$alluvial_perc[is.na(t$alluvial_perc)] <- 0
t$bedrock_perc[is.na(t$bedrock_perc)] <- 0
t$fluvial_perc[is.na(t$fluvial_perc)] <- 0
t$organic_perc[is.na(t$organic_perc)] <- 0
t$glaciogenic_perc[is.na(t$glaciogenic_perc)] <- 0
t$moraine_perc[is.na(t$moraine_perc)] <- 0
t$lakeperc[is.na(t$lakeperc)] <- 0
t$water_perc[is.na(t$water_perc)] <- 0
t$barrenland_perc[is.na(t$barrenland_perc)] <- 0
t$forest_perc[is.na(t$forest_perc)] <- 0
t$grassland_perc[is.na(t$grassland_perc)] <- 0
t$lichenmoss_perc[is.na(t$lichenmoss_perc)] <- 0
t$shrubland_perc[is.na(t$shrubland_perc)] <- 0
t$wetland_perc[is.na(t$wetland_perc)] <- 0
t$percslump17act[is.na(t$percslump17act)] <- 0
t$percslump17all[is.na(t$percslump17all)] <- 0

## (1.5)  Rain data ====================

dfrain <- "ThesisDrafts/Chapter3/WeatherStationData/"

rain24 <- read.csv(paste0(df, "24hrsClimateWindow.csv"))
rain48 <- read.csv(paste0(df, "48hrsClimateWindow.csv"))
rain72 <- read.csv(paste0(df, "72hrsClimateWindow.csv"))
rain96 <- read.csv(paste0(df, "96hrsClimateWindow.csv"))

rain24 <- rain24 %>% select(date=Date, site, RainTot24)
rain48 <- rain48 %>% select(date=Date, site, RainTot48)
rain72 <- rain72 %>% select(date=Date, site, RainTot72)
rain96 <- rain96 %>% select(date=Date, site, RainTot96)

rain <- merge(rain24, rain48)
rain <- merge(rain, rain72)
rain <- merge(rain, rain96)

## (1.5.1)  Merge with 2017 synoptic data ====================

rain$date <- as.Date(rain$date)
t <- merge(t, rain, by=c("site", "date"), all.x=TRUE)

## (1.6)  Rain data for Stony mainstem ====================
rain24t <- read.csv(paste0(df, "24hrsClimateWindow_stonytransscott.csv"))
rain48t <- read.csv(paste0(df, "48hrsClimateWindow_stonytransscott.csv"))
rain72t <- read.csv(paste0(df, "72hrsClimateWindow_stonytransscott.csv"))
rain96t <- read.csv(paste0(df, "96hrsClimateWindow_stonytransscott.csv"))

rain24t <- rain24t %>% select(date=Date, site, RainTot24)
rain48t <- rain48t %>% select(date=Date, site, RainTot48)
rain72t <- rain72t %>% select(date=Date, site, RainTot72)
rain96t <- rain96t %>% select(date=Date, site, RainTot96)

raint <- merge(rain24t, rain48t)
raint <- merge(raint, rain72t)
raint <- merge(raint, rain96t)

## (1.6.1)  Merge with 2017 synoptic data ====================

raint$date <- as.Date(raint$date)

t <- merge(t, raint, by=c("site", "date"), all.x=TRUE)

t$RainTot24.x[is.na(t$RainTot24.x)] <- t$RainTot24.y[is.na(t$RainTot24.x)]
t$RainTot48.x[is.na(t$RainTot48.x)] <- t$RainTot48.y[is.na(t$RainTot48.x)]
t$RainTot72.x[is.na(t$RainTot72.x)] <- t$RainTot72.y[is.na(t$RainTot72.x)]
t$RainTot96.x[is.na(t$RainTot96.x)] <- t$RainTot96.y[is.na(t$RainTot96.x)]

t$RainTot24.y <- NULL; t$RainTot48.y <- NULL; t$RainTot72.y <- NULL; t$RainTot96.y <- NULL

t <- t %>% rename(RainTot24 = RainTot24.x, 
                  RainTot48 = RainTot48.x,
                  RainTot72 = RainTot72.x,
                  RainTot96 = RainTot96.x)

## (1.7)  Add in streambed particle size data ====================

ps <- read.csv(paste0(df, "streambedparticlesizedateadded.csv"))
ps <- ps %>% select(site, D50, psand)
ps$streambednotes <- NA

ps$streambednotes[ps$site=="1"] <- "no rocks/pebbles, moss layer over bedrock"
ps$streambednotes[ps$site=="12-1"] <- "silty, few small pebbles, twigs. Orange phytoplankton"
ps$streambednotes[ps$site=="13"] <- "algae bottom - not completed"
ps$streambednotes[ps$site=="62"] <- "no rocks/pebbles, moss layer over bedrock. <2 silt material around"

# merge with master

t <- merge(ps, t, all=TRUE)

## (1.8)  Add in cation data ===============================

cations <- read_excel(
  paste0(df, "Sarah Shakil 2017 reported on Feb 22, 18.xlsx"),
  sheet="forR")

catids <- read_excel(paste0(df, "SarahShakil_Peel2017_SampleList.xlsx"),
                     sheet="CationsforR")

# (1.8.1) Clean and merge with site ids =====
cations <- cations %>%
  select(site= Site, 
         NamgL = `Na (mg/L)`,
         KmgL = `K (mg/L)`,
         CamgL = `Ca (mg/L)`,
         MgmgL = `Mg (mg/L)`,
         FemgL = `Fe (mg/L)`,
         SrmgL = `Sr (?g/L)`) %>%
  filter(NamgL!=-1)

cations$NamgL[cations$NamgL=="<MDL"] <- 0.0160/2
cations$KmgL[cations$KmgL=="<MDL"] <- 0.0092/2
cations$CamgL[cations$CamgL=="<MDL"] <- 0.0053/2
cations$MgmgL[cations$MgmgL=="<MDL"] <- 0.0103/2
cations$FemgL[cations$FemgL=="<MDL"] <- 0.0161/2
cations$SrmgL[cations$SrmgL=="<MDL"] <- 0.0157/2

cat <-  merge(catids, cations, by.x="Bottle Number", by.y="site")

cat <- cat %>%
  filter(samptype=="Sample") %>%
  select(site, date=Date, NamgL, KmgL, CamgL, MgmgL, FemgL, SrmgL)

cat$NamgL <- as.numeric(cat$NamgL)
cat$KmgL <- as.numeric(cat$KmgL)
cat$CamgL <- as.numeric(cat$CamgL)
cat$MgmgL <- as.numeric(cat$MgmgL)
cat$FemgL <- as.numeric(cat$FemgL)
cat$SrmgL <- as.numeric(cat$SrmgL)

# (1.8.2)  Merge with master sheet =====

t <- merge(t, cat, by=c("site", "date"), all.x=TRUE)

## (1.9)  Add in sulphate data ===============================

so4cl <- read_excel(
  paste0(df, "Sarah Shakil 2017 reported on Feb 22, 18.xlsx"), 
  sheet="forR")

anids <- read_excel(paste0(df, "SarahShakil_Peel2017_SampleList.xlsx"),
                    sheet="AnionsforR")

## (1.9.1) Clean and merge with site ids =====
so4cl <- so4cl %>%
  select(site= Site, 
         ClmgL = `Cl (mg/L)`,
         SO4mgL = `SO4 (mg/L)`) %>%
  filter(ClmgL!=-1)

so4cl$ClmgL[so4cl$ClmgL=="<MDL"] <- 0.03/2
so4cl$SO4mgL[so4cl$SO4mgL=="<MDL"] <- 0.04/2

so4 <-  merge(anids, so4cl, by.x="Bottle Number", by.y="site")

so4 <- so4 %>%
  filter(samptype=="Sample") %>%
  select(site, date=Date, ClmgL, SO4mgL)

so4$ClmgL <- as.numeric(so4$ClmgL)
so4$SO4mgL <- as.numeric(so4$SO4mgL)

# (1.9.2)  Merge with master sheet =====

t <- merge(t, so4, by=c("site", "date"), all.x=TRUE)

## (1.10)  Add in sulphate data ===============================

o18 <- read_excel(
  paste0(df, "Sarah Shakil 2017 reported on Feb 22, 18.xlsx"),
  sheet="forR")

o18ids <- read_excel(paste0(df, "SarahShakil_Peel2017_SampleList.xlsx"),
                     sheet="18O_DHforR")

## (1.10.1) Clean and merge with site ids =====
o18 <- o18 %>%
  select(site= Site, 
         deltadpermille = `Delta D (0/00)`,
         delta18opermille = `Delta 18O (0/00)`) %>%
  filter(deltadpermille!=10000.0)

o18 <-  merge(o18ids, o18, by.x="Bottle Number", by.y="site")

o18 <- o18 %>%
  filter(samptype=="Sample") %>%
  select(site, date=Date, deltadpermille, delta18opermille)

# (1.10.2)  Merge with master sheet =====

t <- merge(t, o18, by=c("site", "date"), all.x=TRUE)

## (1.11)  Add in cleaned field data (slope, ysi, channel morphometry) ===============================

field <- read.csv(paste0(wd, df, "field2017datacleaned.csv"))

field <- field %>% 
    select(site, date, slope, time, lat, long, 
           melev_m, mtemp_degC, sdtemp, 
           mDO_perc, mDO_mgL, mcond_uscm,
           mpH, avgvelocity_ms, Q_m3s, maxdepth_m,
           avgdepth_m, bankfulwidth_m, wettedwidth_m)

t <- merge(t, field, by=c("site", "date"), all=TRUE)

## (1.11)  Add in Chla concentrations (ug/L) ===============================

chla <- read_excel(paste0(df, "ChlaCalc-SampleRun-Shakil-21SEP2017.xlsx"), sheet="forR")

# chla detection limit is 1ug/L for a 200 mL sample, will roughly equate that too 0.2ug
# can't do half DL because the volumes filtered are so variable it will inflate some values to unreal numbers
# so only use values above the chla on filter (ug) detection limit, will then compare this to 

chla$finalconc_ugL_dlcor <- NA
chla$finalconc_ugL_dlcor[chla$chlafilter_ug>0.2] <- chla$finalconc_ugL[chla$chlafilter_ug>0.2]
chla$DLnote <- NA
chla$DLnote[is.na(chla$finalconc_ugL_dlcor)] <- "Below Detection Limit of 1ugL per 200 mL interpreted as 0.2 ug Chla in extraction"

chla <- chla %>% filter(samptype=="Sample") %>% select(-samptype)

# merge with master 
t <- merge(t, chla, by=c("site", "date"), all=TRUE)

## (1.12) Flux and Yield calculations ===============================

t$tssflux <- t$tssmgL*t$dism3s*1000
t$pocflux <- t$POCmgL*t$dism3s*1000
t$docflux <- t$DOCmgL*t$dism3s*1000

t$tssyield <- t$tssflux/t$WatershedArea
t$pocyield <- t$pocflux/t$WatershedArea
t$docyield <- t$docflux/t$WatershedArea
t$tocyield <- t$pocyield + t$docyield
t$percPOCy <- (t$pocyield/t$tocyield)*100
t$percDOCy <- (t$docyield/t$tocyield)*100
t$Cayield <- (t$CamgL*t$dism3s*1000)/t$WatershedArea
t$Nayield <- (t$NamgL*t$dism3s*1000)/t$WatershedArea
t$Mgyield <- (t$MgmgL*t$dism3s*1000)/t$WatershedArea
t$SO4yield <- (t$SO4mgL*t$dism3s*1000)/t$WatershedArea
t$Clyield <- (t$ClmgL*t$dism3s*1000)/t$WatershedArea
t$Sryield <- (t$SrmgL*t$dism3s*1000)/t$WatershedArea
t$Feyield <- (t$FemgL*t$dism3s*1000)/t$WatershedArea
t$wateryield <- (t$dism3s*1000)/t$WatershedArea

t$streampower <- (1000)*(9.8)*(t$Q_m3s)*(t$slope)

## (1.11) Write CSV ====================

write.csv(t, paste0(df, "2017data.csv"))
