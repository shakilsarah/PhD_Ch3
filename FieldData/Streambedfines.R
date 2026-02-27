# ===========================================================================================================#
# Streambedfines.R
# By: Sarah Shakil
# Contact: shakil@ualberta.ca
# Background: Stony Creek Transect plots and D50 vs watershed area plots
# Last Updated: October 10 2021
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

d <- d %>%
  select(campaign, site, date, D50, psand, tssmgL,
         POCmgL, DOCmgL, tssflux, pocflux, docflux, 
         tssyield, pocyield, docyield, tocyield, 
         PO13C, POCTSSrat,
         RainTot24, RainTot48, RainTot72, RainTot96,
         WatershedArea,
         delta18opermille, Cayield, wateryield, streampower,
         slope, percshale, colluvial_perc, piedmont_perc,
         moraine_perc, lakeperc, scaledgpp, forest_perc, 
         grassland_perc, lichenmoss_perc, wmeanSOCC_100CM,
         meanslope_deg, RainTot96, percslump17act, percslump17all,
         slumpacccount, strahlerimpactacc, strahlerstream)

d$date <- as.character(as.Date(d$date))

## (1.2) Read in transect distances ====================
dist <- read_excel(paste0(df, "20162017POCPO14CTrans.xlsx"))
dist <- dist %>%
        filter(`Sample Type`=="Sample" & `Stream Location`=="PERI") %>%
        select(site = `Slump Site`,
               loc = `Stream Location`,
               trans = `Transect Location`,
               date = `Sampling Date`,
               distm, 
               F14C,
               F14error = `F14C error`)

dist$date <- as.character(as.Date(dist$date))

## (1.3) Read in OC optics ====================
o <- read.csv(paste0(df, "masteroptics2.csv"))
o <- o %>%
     select(site, date, JDay, prcntC1_p, prcntC3_p, SR_d, SUVA254, slumpYN)

o$date <- as.character(as.Date(o$date))

## (1.4) Distance from nearest active slump ====================

sdist <- read_excel(paste0(df, "slumpstreamdist.xlsx"))

## (1.5) Merge together ====================

a <- merge(d, dist, by=c("site", "date"), all=TRUE)
a <- merge(a, o, by=c("site", "date"), all=TRUE)
a <- merge(a, sdist, by=c("site"), all.x=TRUE)
a$nearslumpdist_m<- as.numeric(a$nearslumpdist_m)

##### ========== (2) %fines vs watershed area ==========================================================================

## (2.1) Save theme ====================
theme <-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.background = element_blank(),axis.line.x = element_line(colour="black"),
              axis.line.y = element_line(colour="black"),
              axis.text = element_text(colour="black",size=12),legend.background=element_blank(),
              text=element_text(size = 12),
              legend.title=element_blank(),
              legend.position = c(1,1),
              legend.direction="vertical",
              legend.justification = c(1,1),
              aspect.ratio=1)

## (2.2) Graph ====================
asyn <- a[a$campaign=="2017synoptic",]

finea <- ggplot() + 
  geom_point(data=asyn[!is.na(asyn$psand),], 
             aes(x=WatershedArea, y=(psand), fill=slumpYN),
             size=3.5,
             shape=21, colour="grey20") +
  scale_fill_manual(breaks=c("Y", "N"), 
                    values=c("pink", "Sky Blue")) + 
  labs(x=expression("Watershed Area (km"^2*")"),
       y="%fines") +
  theme

fineb <- ggplot() + 
  geom_point(data=asyn[!is.na(asyn$psand),], 
             aes(x=nearslumpdist_m, y=(psand)),
             size=3.5,
             shape=21, colour="grey20", fill="pink") +
  scale_x_log10() +
  labs(x=expression("Distance to nearest slump (m)"),
       y="%fines") +
  theme

finec <- ggplot() + 
  geom_point(data=asyn[!is.na(asyn$psand) & asyn$slumpYN=="Y",], 
             aes(x=percslump17act, y=(psand)),
             size=3.5,
             shape=21, colour="grey20", fill="pink") +
  scale_x_log10() +
  labs(x=expression("% area active slumps"),
       y="%fines") +
  theme
## (2.3) Print graphs ====================

library(grid)

savePlotpdf <- function(myPlot) {
  pdf(file = paste0(wd, "Figures/", filename, ".pdf"),
      width = 6, height = 3)
  print(myPlot)
  dev.off()
}

g1 <- ggplotGrob(finea)
g2 <- ggplotGrob(fineb)

ga <- cbind(g1, g2, size="first")
grid.draw(ga)
filename <- "streambedfines"
savePlotpdf(grid.draw(ga))

## (2.3) Stats analyses ====================
l <- lm((psand) ~ as.factor(slumpYN) + RainTot96 + JDay,
                  data=asyn)
l <- lm((psand) ~ as.factor(slumpYN) + JDay,
        data=asyn)
l <- lm((psand) ~ as.factor(slumpYN),
        data=asyn)
anova(l)

l <- lm(psand ~ nearslumpdist_m, data=asyn)
