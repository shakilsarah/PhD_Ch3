# ===========================================================================================================#
# StonyCreekTransect.R
# By: Sarah Shakil
# Contact: shakil@ualberta.ca
# Background: Stony Creek Transect plots and D50 vs watershed area plots
# Last Updated: October 10 2021
# Modified: Combined into a single publication-ready figure using patchwork, 
#           added ggrepel for labels, fixed panel legends, changed flux to yield,
#           unified shape shading/positioning for panel j, and moved panel k legend to top.
#===========================================================================================================#

##### (i) Workspace PREP ==============================================================================

## Clear list 
list=rm(list=ls(all=TRUE))

## Set working directory
df <- "data/"

# load libraries
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(patchwork) # For combining ggplots into a single figure
library(ggrepel)   # Added to prevent overlapping text labels


# load functions
se <- function(x) {
  sqrt(var(x, na.rm=TRUE)/length(x[!is.na(x)])) 
}

## create standard error function for use later
lengthnona <- function(x) {
  length(x[!is.na(x)])
}

##### ========== (1) DATA PREP ==========================================================================

## (1.1) Read in OC yields, PO13C, and POCTSSrat ====================
d <- read.csv(paste0(df, "2017data.csv"))

## (1.11) Update the surficial geology =============================

# read in the updated watmaster file 
watmaster <- read.csv(paste0(df,"watmaster_wlakes_2026.csv"))

# replace NA in the surficial geology columns with 0

watmaster <- watmaster %>%
  mutate(
    across(
      c(
        colluvial_perc,
        morainal_perc,
        alluvial_perc,
        glaciofluvial_perc,
        organic_perc,
        glaciolacustrine_perc
      ),
      ~ replace_na(.x, 0)
    )
  )


# Columns to bring in from watmaster
watmaster_cols <- c(
  "site",
  "date",
  "trans",
  "percshale",
  "colluvial_perc",
  "morainal_perc",
  "alluvial_perc",
  "glaciofluvial_perc",
  "organic_perc",
  "glaciolacustrine_perc",
  "lakeperc"
)

# Columns to remove from d before adding updated columns from watmaster
cols_to_remove_from_d <- c(
  "percshale",
  "colluvial_perc",
  "piedmont_perc",
  "alluvial_perc",
  "bedrock_perc",
  "fluvial_perc",
  "glaciogenic_perc",
  "organic_perc",
  "moraine_perc",
  "lakeperc"
)

# check if there are any missing site-dates in watmaster
missing_site_dates <- d %>%
  distinct(site, date) %>%
  anti_join(
    watmaster %>% distinct(site, date),
    by = c("site", "date")
  )

# View missing combinations
missing_site_dates
# yes - there are two 

# Subset watmaster to only site-date combinations present in d
watmaster_subset <- watmaster %>%
  semi_join(
    d %>% select(site, date),
    by = c("site", "date")
  ) %>%
  select(all_of(watmaster_cols))

# Remove old columns from d, then add selected watmaster columns by site and date
d_updated <- d %>%
  select(-any_of(cols_to_remove_from_d)) %>%
  left_join(
    watmaster_subset,
    by = c("site", "date")
  )

# Check result
dim(d_updated)
head(d_updated)

# Optional: check whether watmaster has duplicate site-date combinations
# Doesn't matter because none of the SE samples are used here 
watmaster %>%
  count(site, date) %>%
  filter(n > 1)

## Rename d_updated to d and continue 

d <- d_updated

# Continue 

d <- d %>%
  select(campaign, site, date, D50, psand, tssmgL,
         POCmgL, DOCmgL, tssflux, pocflux, docflux, 
         tssyield, pocyield, docyield, tocyield, 
         PO13C, POCTSSrat,
         RainTot24, RainTot48, RainTot72, RainTot96,
         WatershedArea,
         delta18opermille, Cayield, wateryield, streampower,
         slope, percshale, 
         colluvial_perc, morainal_perc,  # updated 
         lakeperc, scaledgpp, forest_perc, 
         grassland_perc, lichenmoss_perc, wmeanSOCC_100CM,
         meanslope_deg, RainTot96, percslump17act, percslump17all,
         slumpacccount, strahlerimpactacc, strahlerstream, dism3s)

d$date <- as.character(as.Date(d$date))

## (1.2) Read in transect distances ====================
dista <- read_excel(paste0(df, "20162017POCPO14CTrans.xlsx"))
dista$`Sampling Date` <- as.character(as.Date(dista$`Sampling Date`))
dist <- dista %>%
  filter(`Sample Type`=="Sample" & !(`Stream Location`=="PERI")) %>%
  select(site = `Slump Site`,
         loc = `Stream Location`,
         trans = `Transect Location`,
         date = `Sampling Date`,
         distm) %>%
  group_by(site, loc, trans, date) %>%
  summarize(distm=mean(distm, na.rm=TRUE)) %>%
  ungroup

dist$distm[dist$site=="SC-2B"] <- 2858.267449

## (1.2) Reformat for 14C data ====================
po14c <- dista %>% 
  filter(MatCode=="POC" & !(is.na(F14C))) %>%
  select(site = `Slump Site`,
         loc = `Stream Location`,
         trans = `Transect Location`,
         date = `Sampling Date`, 
         F14C_poc=F14C, 
         F14Cerror_poc=`F14C error`) %>%
  group_by(site, date) %>%
  summarize(F14C_poc=mean(F14C_poc, na.rm=TRUE),
            F14Cerror_poc=mean(F14Cerror_poc, na.rm=TRUE)) %>% # just because there are duplicate rows
  ungroup

do14c <- dista %>% 
  filter(MatCode=="DOC" & !(is.na(F14C))) %>%
  select(site = `Slump Site`,
         loc = `Stream Location`,
         trans = `Transect Location`,
         date = `Sampling Date`,
         F14C_doc=F14C, 
         F14Cerror_doc=`F14C error`) %>%
  group_by(site, date) %>%
  summarize(F14C_doc=mean(F14C_doc, na.rm=TRUE),
            F14Cerror_doc=mean(F14Cerror_doc, na.rm=TRUE)) %>% # just because there are duplicate rows
  ungroup

## (1.3) Read in OC optics ====================
o <- read.csv(paste0(df, "masteroptics2.csv"))
o <- o %>%
  select(site, date, JDay, prcntC1_p, prcntC3_p, prcntC2_p, prcntC4_p, prcntC5_p,
         SR_d, SUVA254, slumpYN)
o$sum <- (o$prcntC1_p) + o$prcntC2_p + o$prcntC3_p + o$prcntC4_p
o$date <- as.character(as.Date(o$date))

## (1.4) Distance from nearest active slump ====================

sdist <- read_excel(paste0(df, "slumpstreamdist.xlsx"))

## (1.5) Merge together ====================

a <- merge(d, dist, by=c("site", "date"), all=TRUE)
a <- merge(a, o, by=c("site", "date"), all=TRUE)
a <- merge(a, sdist, by=c("site"), all.x=TRUE)
a <- merge(a, po14c, by=c("site", "date"), all.x=TRUE)
a <- merge(a, do14c, by=c("site", "date"), all.x=TRUE)
a$nearslumpdist_m<- as.numeric(a$nearslumpdist_m)

## (1.6) Select variables of interest ====================

aarch <- a %>%
  filter(!(site=="DC-3" | site=="DC-4" | site=="SC-1A")) %>%
  select(-D50, -psand, -note)

a <- a %>% filter(campaign=="2017transect" | site==3) %>%
  filter(!(site=="DC-3" | site=="DC-4" | site=="SC-1A")) %>%
  select(-D50, -psand, -note)

## (1.7) Read in slump UP and DN from Chapter 1 ====================
c1 <- read.csv(paste0(df, "PeelPlateau_RTSandstream_geochem.csv"))

cav <- c1 %>% 
  filter(sampletype=="Stream" & 
           (streamlocation=="UP"|streamlocation=="DN") ) %>%
  group_by(streamlocation) %>%
  summarize(P1av=mean(prcntC1, na.rm=TRUE), P1se=se(prcntC1),
            P3av=mean(prcntC3, na.rm=TRUE), P3se=se(prcntC3),
            percPOCav=mean(percPOC, na.rm=TRUE), percPOCse=se(percPOC),
            PO13Cav=mean(PO13C, na.rm=TRUE), PO13Cse=se(PO13C))

##### ========== (2) %fines vs watershed area ==========================================================================

## (2.1) Save theme ====================
theme <- theme(panel.grid.major = element_blank(),
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
               legend.key=element_blank(),
               panel.border = element_rect(colour = "black", fill=NA, size=1),
               legend.spacing = unit(0, "mm"),
               legend.key.size = unit(0.001, "cm"), # This hides keys in scatterplots
               plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

paneltheme <- theme(strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    panel.spacing.x = unit(0.4, "inches"))

scalex <- scale_x_continuous(limits=c(0, 85)) 

###  (2.2) Reformat constitutent data for graphing  ==================

along <- a%>%
  pivot_longer(cols=c("DOCmgL", "POCmgL", "tssmgL",
                      "docflux", "pocflux", "tssflux",
                      "docyield", "pocyield", "tssyield",
                      "POCTSSrat", "PO13C", "prcntC1_p", "prcntC3_p",
                      "SR_d", "SUVA254"),
               names_to="variables",
               values_to="values") %>%
  select(site, distm, date, F14C_poc, F14Cerror_poc, 
         F14C_doc, F14Cerror_doc,
         variables, values, slumpYN, WatershedArea,
         scaledgpp, percslump17act, wateryield, meanslope_deg,
         streampower, morainal_perc, colluvial_perc, Cayield, strahlerstream) # updated

along$panel <- "conc"
along$panel[along$variables=="docflux"|
              along$variables=="pocflux"|
              along$variables=="tssflux"] <- "flux"
along$panel[along$variables=="docyield"|
              along$variables=="pocyield"|
              along$variables=="tssyield"] <- "yield"
along$panel[along$variables=="POCTSSrat"|
              along$variables=="PO13C"|
              along$variables=="prcntC1_p"|
              along$variables=="prcntC3_p"|
              along$variables=="SR_d"|
              along$variables=="SUVA254"] <- "comp"

along$frac <- "doc"
along$frac[along$variables=="POCmgL"|
             along$variables=="pocflux"|
             along$variables=="pocyield"] <- "poc"
along$frac[along$variables=="tssmgL"|
             along$variables=="tssflux"|
             along$variables=="tssyield"] <- "tss"
along$frac[along$variables=="POCTSSrat"] <- "%poc"
along$frac[along$variables=="PO13C"] <- "PO13C"
along$frac[along$variables=="prcntC1_p"] <- "%P1"
along$frac[along$variables=="prcntC3_p"] <- "%P3"
along$frac[along$variables=="SR_d"] <- "SR"
along$frac[along$variables=="SUVA254"] <- "SUVA"

## (2.3) Remformat geospatial data for graphing  ====================

along2 <- a%>%
  pivot_longer(cols=c("scaledgpp", "percslump17act", "wateryield",
                      "meanslope_deg", "morainal_perc",  # updated
                      "colluvial_perc", "slumpacccount"),
               names_to="variables",
               values_to="values") %>%
  select(site, distm, date, variables, values, slumpYN, WatershedArea, strahlerstream)

## (2.4.1) Watershed Variables ======
along2a <- along2 %>% filter(variables=="scaledgpp")
along2b <- along2 %>% filter(variables=="percslump17act")
wdes1 <- ggplot() + 
  geom_point(data=along2a, 
             aes(x=distm/1000, 
                 y=values, shape=variables),
             colour="black", size=4) +
  geom_line(data=along2a, 
            aes(x=distm/1000, y=values), linetype=2) +
  geom_point(data=along2b, 
             aes(x=distm/1000, 
                 y=values/2, shape=variables),
             colour="black", size=4) +
  geom_line(data=along2b, 
            aes(x=distm/1000, y=values/2), linetype=2) +
  geom_text_repel(data=along2b, 
                  aes(x=distm/1000, y = values/2, label = rownames(along2a)), 
                  nudge_y = -0.04, min.segment.length = 0.2) + 
  scale_shape(labels=c("%RTS", "GPP")) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 2, 
                                         name=expression("% RTS"["active"]))) +
  labs(x="Distance (km)", y="Mean Watershed GPP") +
  scalex +
  theme + theme(legend.position="top")

var <- c("morainal_perc", "colluvial_perc")

along2c <- along2 %>% filter(variables %in% var)

var <- c("meanslope_deg")

along2d <- along2 %>% filter(variables %in% var)

wdes2 <- ggplot() + 
  geom_point(data=along2c, 
             aes(x=distm/1000, 
                 y=values/10, shape=variables),
             colour="black", size=4) +
  geom_line(data=along2c, 
            aes(x=distm/1000, y=values/10, 
                group=variables), linetype=2) +
  geom_point(data=along2d, 
             aes(x=distm/1000, 
                 y=values, shape=variables),
             colour="black", size=4) +
  geom_line(data=along2d, 
            aes(x=distm/1000, y=values, group=variables), linetype=2) +
  scale_shape(labels=c("% Colluvial", "Watershed Slope", "% Morainal")) +
  labs(x="Distance (km)", 
       y="% Colluvial/Morainal or \n Mean Watershed Slope (\u00B0)") +
  scalex +
  theme + theme(legend.position="top")

<<<<<<< HEAD:Analyses/StonyCreekTransect.R
# Close all active graphics devices
graphics.off()
wdes2

## (2.4.2) Conc orig ======
=======
>>>>>>> 88ba80aeb5fec992760acba6d25d3a434ec889fb:Analyses/Figure3-StonyTrans.R

## (2.4.3) Conc and Yield (yield_doc and yield_poc) ======

yield_doc <- ggplot() + 
  geom_point(data=along[along$panel=="conc" & along$frac=="doc",], 
             aes(x=distm/1000, y=values*10, fill=slumpYN, shape="Concentration"), colour="black", size=4) +
  geom_line(data=along[along$panel=="conc" & along$frac=="doc",], 
            aes(x=distm/1000, y=values*10), linetype=3) +
  geom_point(data=along[along$panel=="yield" & along$frac=="doc",], 
             aes(x=distm/1000, y=values, fill=slumpYN, shape="Yield"), colour="black", size=4) +
  geom_line(data=along[along$panel=="yield" & along$frac=="doc",], 
            aes(x=distm/1000, y=values), linetype=2) +
  scale_fill_manual(values=c("white", "grey"), labels=c("NO RTS", "RTS")) +
  scale_shape_manual(values=c("Concentration"=22, "Yield"=21)) +
  guides(
    shape = guide_legend(override.aes = list(fill = "grey", size = 4), order = 1, nrow = 1),
    fill = guide_legend(override.aes = list(shape = 22, size = 4), order = 2, nrow = 1)
  ) +
  annotate("text", x = Inf, y = Inf, label = "DOC", hjust = 1.5, vjust = 1.5, size = 6, fontface = "bold") +
  labs(x="", y=expression("Yield (mg km"^-2*" s"^-1*")")) +
  scalex+
  scale_y_continuous(sec.axis = sec_axis(~ . / 10, name=expression("Concentration mg L"^-1))) +
  theme +
  theme(
    legend.position = "top", 
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.size = unit(0.5, "cm")
  ) +
  paneltheme

yield_poc <- ggplot() + 
  geom_point(data=along[along$panel=="conc" & along$frac=="poc",], 
             aes(x=distm/1000, y=values*10, fill=slumpYN, shape="Concentration"), colour="black", size=4) +
  geom_line(data=along[along$panel=="conc" & along$frac=="poc",], 
            aes(x=distm/1000, y=values*10), linetype=3) +
  geom_point(data=along[along$panel=="yield" & along$frac=="poc",], 
             aes(x=distm/1000, y=values, fill=slumpYN, shape="Yield"), colour="black", size=4) +
  geom_line(data=along[along$panel=="yield" & along$frac=="poc",], 
            aes(x=distm/1000, y=values), linetype=2) +
  scale_fill_manual(values=c("white", "grey"), labels=c("NO RTS", "RTS")) +
  scale_shape_manual(values=c("Concentration"=22, "Yield"=21)) +
  annotate("text", x = Inf, y = Inf, label = "POC", hjust = 1.5, vjust = 1.5, size = 6, fontface = "bold") +
  labs(x="", y=expression("Yield (mg km"^-2*" s"^-1*")")) +
  scalex+
  scale_y_continuous(sec.axis = sec_axis(~ . / 10, name=expression("Concentration mg L"^-1))) +
  theme + 
  theme(legend.position="none") + 
  paneltheme

## (2.4.5) %POC (comp1a) ======

comp1a <- ggplot() +
  geom_rect(data = cav, aes(ymin = percPOCav-percPOCse, 
                            ymax = percPOCav+percPOCse,
                            xmin = -Inf, xmax = Inf, 
                            group=streamlocation),
            fill="grey", colour=NA, alpha=0.2)+
  geom_line(data=along[along$frac=="%poc",], 
            aes(x=distm/1000, y=values), linetype=1) +
  geom_point(data=along[along$frac=="%poc",], 
             aes(x=distm/1000, y=values, fill=slumpYN),
             colour="black", size=4, shape=22) +
  scale_fill_manual(values=c("white", "grey"),
                    labels=c("NO RTS", "RTS")) +
  scalex +
  labs(x="", y="%POC") +
  theme +
  theme(legend.position = c(1,1),
        legend.direction="horizontal",
        legend.justification = c(1,1))

## (2.4.6) PO13C (comp1b) ======

comp1b <- ggplot() +
  geom_rect(data = cav, aes(ymin = PO13Cav-PO13Cse, 
                            ymax = PO13Cav+PO13Cse,
                            xmin = -Inf, xmax = Inf, group=streamlocation),
            fill="grey", colour=NA, linetype=2, alpha=0.2)+
  geom_line(data=along[along$frac=="PO13C",], 
            aes(x=distm/1000, y=values), linetype=2) +
  geom_point(data=along[along$frac=="PO13C",], 
             aes(x=distm/1000, y=values, fill=slumpYN),
             colour="black", size=4, shape=22) +
  scale_fill_manual(values=c("white", "grey"),
                    labels=c("NO RTS", "RTS")) +
  scalex +
  labs(x="", y=expression(delta^13*"C"["POC"])) +
  theme + theme(legend.position="none")

## (2.4.8) P1 (comp2a) ======

comp2a <- ggplot() + 
  geom_rect(data = cav, aes(ymin = P1av-P1se, 
                            ymax = P1av+P1se,
                            xmin = -Inf, xmax = Inf, group=streamlocation),
            fill="grey", colour=NA, alpha=0.2)+
  geom_line(data=along[along$frac=="%P1",], 
            aes(x=distm/1000, y=values), linetype=2) +
  geom_point(data=along[along$frac=="%P1",], 
             aes(x=distm/1000, y=values, fill=slumpYN),
             colour="black", size=4, shape=22) +
  scale_fill_manual(values=c("white", "grey"),
                    labels=c("NO RTS", "RTS"),
                    guide = guide_legend(
                      override.aes = list(shape=22))) +
  labs(x="", y="% P1") +
  scalex +
  theme + theme(legend.position="none")

## (2.4.8) P3 (comp2b) ======

comp2b <- ggplot() + 
  geom_rect(data = cav, aes(ymin = P3av-P3se, 
                            ymax = P3av+P3se,
                            xmin = -Inf, xmax = Inf, group=streamlocation),
            fill="grey", colour=NA, alpha=0.2)+
  geom_line(data=along[along$frac=="%P3",], 
            aes(x=distm/1000, y=values), linetype=2) +
  geom_point(data=along[along$frac=="%P3",], 
             aes(x=distm/1000, y=values, fill=slumpYN),
             colour="black", size=4, shape=22) +
  scale_fill_manual(values=c("white", "grey"),
                    labels=c("NO RTS", "RTS"),
                    guide = guide_legend(
                      override.aes = list(shape=22))) +
  labs(x="", y="% P3") +
  scalex +
  theme + theme(legend.position="none")

## (2.4.9) C14 (c14) ======
docsed14c <- read_excel(paste0(df, "PO14C.xlsx"))

doc14c <- docsed14c %>%
  filter(MatCode=="DOC" & `Slump Site`=="SC Outlet")
doc14c$distm <- 70634.884
doc14c$distm[doc14c$`Sampling Date`=="2017-08-07"] <- 71781.598

sed14c <- docsed14c %>%
  filter(MatCode=="A" & `Slump Site`=="SC Outlet")
sed14c$distm <- 70634.884
sed14c$distm[sed14c$`Sampling Date`=="2017-08-07"] <- 71781.598

c14 <- ggplot() + 
  geom_line(data=a, 
            aes(x=distm/1000, y=F14C_poc), linetype=1) +
  geom_point(data=a, 
             aes(x=distm/1000, y=F14C_poc, fill=slumpYN, shape="POC"),
             colour="black", size=4) +
  geom_point(data=doc14c, 
             aes(x=distm/1000, y=F14C, shape="DOC"),
             colour="black", size=4, fill="grey") +
  geom_point(data=sed14c, 
             aes(x=distm/1000, y=F14C, shape="streambank"),
             colour="black", size=4, fill="grey") +
  geom_errorbar(data=a, 
                aes(x=distm/1000, ymin=F14C_poc-F14Cerror_poc, 
                    ymax=F14C_poc+F14Cerror_poc)) +
  geom_errorbar(data=doc14c, 
                aes(x=distm/1000, ymin=F14C-`F14C error`, 
                    ymax=F14C+`F14C error`)) +
  geom_errorbar(data=sed14c, 
                aes(x=distm/1000, ymin=F14C-`F14C error`, 
                    ymax=F14C+`F14C error`)) +
  scale_fill_manual(values=c("white", "grey"),
                    labels=c("NO RTS", "RTS")) +
  scale_shape_manual(values=c("DOC"=23, "POC"=22, "streambank"=24)) +
  guides(
    fill = "none",
    shape = guide_legend(override.aes = list(fill = "grey", size = 4), title=NULL, nrow=1)
  ) +
  labs(x="Distance (km)", y=expression("F"^14*"C")) +
  scalex +
  theme + 
  theme(
    legend.position = "top", 
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.key.size = unit(0.5, "cm")
  )

## (2.4.9) C14 and P1 P3 correlation (cor) ======
cor <- ggplot() + 
  geom_smooth(data=along[along$frac=="%P1"|along$frac=="%P3",],
              method='lm', 
              aes(y=values, x=F14C_poc, group=frac),
              formula= y~x, colour="black") +
  geom_point(data=along[along$frac=="%P1"|along$frac=="%P3",], 
             aes(x=F14C_poc, y=values, shape=frac, fill=slumpYN),
             colour="black",
             size=4) +
  scale_fill_manual(values=c("white", "grey"),
                    labels=c("NO RTS", "RTS")) + 
  scale_shape_manual(values = c(21, 22)) +
  # --- FIX 1: Set guide layouts to horizontal rows ---
  guides(
    fill = guide_legend(override.aes = list(shape = c(21,22)), nrow = 1, order = 1),
    shape = guide_legend(nrow = 1, order = 2, title=NULL)
  ) +
  labs(y="% fluorescence", x=expression("F"^14*"C")) +
  scale_x_log10() + 
  scale_y_log10() +
  theme +
  # --- FIX 2: Move legend to the top of the plot and layout horizontally ---
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.size = unit(0.5, "cm")
  )

## (2.4.10) Rainfall ======

a$site2 <- rownames(a)
arain <- a%>%
  pivot_longer(
    cols=c(RainTot24, RainTot48, RainTot72, RainTot96),
    names_to = "hours", values_to="rain") %>%
  select(site2, hours, rain)

rain <- ggplot() + 
  geom_bar(data=arain,
           aes(y=rain, x=site2, fill=hours), colour="black",
           stat="identity", position="dodge") +
  scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11)) +
  geom_point(data=a, aes(x=site2, y=dism3s, shape="Discharge"), 
             colour="black", fill="SkyBlue", size=4) +
  scale_fill_manual(limits=c("RainTot24", "RainTot48", "RainTot72", "RainTot96"),
                    labels=c("24", "48", "72", "96"),
                    values=c("black", "grey50", "grey90", "white")) +
  scale_shape_manual(values=c("Discharge"=23)) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size=5), order=1),
         shape = guide_legend(override.aes = list(fill="SkyBlue", colour="black", size=4), order=2)) +
  labs(x="Site", y=expression("Total Rainfall (mm) \n or Discharge (m"^3*"s"^-1*")")) +
  theme + 
  theme(legend.position="top",
        legend.key.size = unit(0.5, "cm")) 


## (2.5) PRINT SINGLE FIGURE ====================

<<<<<<< HEAD:Analyses/StonyCreekTransect.R
library(grid)
library(svglite)

# Function to save ggplot or grob as SVG
savePlotSVG <- function(myPlot, w, h, filename, out_dir = "Figures") {
  
  # Create output folder if it does not already exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Full output path
  out_file <- file.path(out_dir, paste0(filename, ".svg"))
  
  # Open SVG device
  svglite::svglite(file = out_file, width = w, height = h)
  
  # Draw plot/grob
  if (inherits(myPlot, "ggplot")) {
    print(myPlot)
  } else {
    grid::grid.draw(myPlot)
  }
  
  # Close device
  dev.off()
}

# Choose where you want the plots saved
fig_dir <- "D:/5_Projects/Shakil_Ch3/figures/transect"

# Save figures as SVGs
savePlotSVG(
  cbind(ggplotGrob(wdes1), ggplotGrob(wdes2), size = "first"),
  w = 8.75, h = 3,
  filename = "Stonytransa",
  out_dir = fig_dir
)

savePlotSVG(
  ggplotGrob(flux),
  w = 9, h = 2.45,
  filename = "Stonytransb",
  out_dir = fig_dir
)

savePlotSVG(
  cbind(ggplotGrob(comp1a), ggplotGrob(comp1b), size = "first"),
  w = 8.75 - 0.35, h = 2.45,
  filename = "Stonytransc",
  out_dir = fig_dir
)

savePlotSVG(
  cbind(ggplotGrob(comp2a), ggplotGrob(comp2b), size = "first"),
  w = 8.75 - 0.35, h = 2.45,
  filename = "Stonytransd",
  out_dir = fig_dir
)

savePlotSVG(
  cbind(ggplotGrob(c14), ggplotGrob(cor), size = "first"),
  w = 8.75 - 0.35, h = 2.5,
  filename = "Stonytranse",
  out_dir = fig_dir
)

savePlotSVG(
  ggplotGrob(rain),
  w = 9, h = 2.3,
  filename = "Stonytransrain",
  out_dir = fig_dir
)



=======
# Combine all plots using patchwork syntax:
# / denotes a new row, | denotes side-by-side
final_figure <- rain /
  (wdes1 | wdes2) /
  (yield_doc | yield_poc) /
  (comp1a | comp1b) /
  (comp2a | comp2b) /
  (c14 | cor) +
  plot_layout(heights = c(1, 1, 1, 1, 1, 1)) +
  plot_annotation(tag_levels = 'a', 
                  tag_prefix = '(', 
                  tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 14))

# Check for/Create the nested Fig3 directory before saving
if(!dir.exists("Figures")) {
  dir.create("Figures")
}
if(!dir.exists("Figures/Fig3")) {
  dir.create("Figures/Fig3")
}

ggsave("Figures/Fig3/StonyCreekTransect_FullFigure.pdf", final_figure, width = 10, height = 15)
ggsave("Figures/Fig3/StonyCreekTransect_FullFigure.png", final_figure, width = 10, height = 15, dpi = 300)

## (2.4) Stats ====================

library(car)
## (2.4.1) Test of linear model on POC yield =====

a$logpocy <- log10(a$pocyield)

# (Eq. 1)
a$predpoceq1 <- ((a$percslump17act)*0.893) + (0.019*(a$wateryield)) + 1.159
a$predratpoceq1 <- (a$predpoceq1)/a$logpocy

# (Eq. 2)
a$predpoc <- ((a$percslump17act)*1.34) + 0.336
a$predratpoceq2 <- (a$predpoceq2)/a$logpocy

View(a[,c("site", "percslump17act", "predpoc", "logpocy", "pocyield", "predratpoc", "RainTot96")])
# linear model seems really unstable

## (2.4.2) Test of linear model on TOC yield =====

a$logtocy <- log10(a$tocyield)
a$predtoc <- (a$percslump17act)*0.893 + 0.019*a$wateryield + 1.159
a$predrattoc <- (a$predtoc)/a$logtocy

View(a[,c("site", "percslump17act", "predtoc", "logtocy", "tocyield", "predrattoc", "RainTot96")])
# linear model seems really unstable

## (2.4.3) Test of linear model on F14C =====
#an <- a[!is.na(a$dism3s),] can't remember why site without discharge was removed, since no reason don't do it
l1 <- lm(log10(F14C_poc) ~ percslump17act + meanslope_deg + scaledgpp + JDay + RainTot96, data=a[-1,])
l2 <- lm(log10(F14C_poc) ~ meanslope_deg + scaledgpp + JDay + RainTot96, data=a[-1,])
anova(l1, l2)
l3 <- lm(log10(F14C_poc) ~ scaledgpp + JDay + RainTot96, data=a[-1,])
anova(l2, l3)
l4 <- lm(log10(F14C_poc) ~ scaledgpp + RainTot96, data=a[-1,])
anova(l3, l4)

plot(l4)
vif(l4)
summary(l4)

## (2.4.4) Stats for discussion =====

slumpallareastony <- (a$percslump17all[a$distm==max(a$distm)])*(a$WatershedArea[a$distm==max(a$distm)])
slumpactareastony <- (a$percslump17act[a$distm==max(a$distm)])*(a$WatershedArea[a$distm==max(a$distm)])

a$perctotareaslump <- (a$percslump17all*a$WatershedArea)/slumpallareastony
a$percactareaslump <- (a$percslump17act*a$WatershedArea)/slumpactareastony

megaslumps <- 536144.810233 + 83176.439906 + 431886.92817 + (160647.604697 + 44281.827707)
megaslumpskm2 <- megaslumps / (10^6)
megaslumpskm2/slumpallareastony
View(a[,c("site", "perctotareaslump", "percactareaslump", "percslump17all")])

