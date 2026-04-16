# ==============================================================================#
# 2017RDA_Publication_StrictMatch.R
# Background: RDA stats and 3-panel publication figure
# Method: Strict adherence to original listwise deletion (Complete Case)
# ==============================================================================#

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

# continue for all variables .... # decided to log10(x+2) transformation for all variables

# (2.2) Step 2: Transformations ====================
# for outlier data and removal of variables lacking data 

library(dplyr)

x <- 2

d <- d %>%
  mutate(
    logtocy          = log10(tocyield + x),
    logpocy          = log10(pocyield + x),
    logdocy          = log10(docyield + x),
    logwy            = log10(wateryield + x),
    logcond          = log10(mcond_uscm + x),
    logcay           = log10(Cayield + x),
    lognay           = log10(Nayield + x),
    logmgy           = log10(Mgyield + x),
    logso4y          = log10(SO4yield + x),
    logsry           = log10(Sryield + x),
    logfey           = log10(Feyield + x),
    logspower        = log10(streampower + x),
    logsslope        = log10(slope + x), # careful: original code uses 'slope' here
    shale_log        = log10(percshale + x),
    col_log          = log10(colluvial_perc + x),
    pied_log         = log10(piedmont_perc + x),
    moraine_log      = log10(moraine_perc + x),
    barren_log       = log10(barrenland_perc + x),
    grass_log        = log10(grassland_perc + x),
    lichen_log       = log10(lichenmoss_perc + x),
    lake_log         = log10(lakeperc + x),
    forest_log       = log10(forest_perc + x),
    shrub_log        = log10(shrubland_perc + x),
    pslumpact_log    = log10(percslump17act + x),
    pslumpall_log    = log10(percslump17all + x),
    logslumpcount    = log10(slumpacccount + x),
    logslumpstrahler = log10(strahlerimpactacc + x),
    # log18O         = log10(delta18opermille + x),
    loggpp           = log10(scaledgpp + x), 
    logelev          = log10(meanelev_m + x),
    logslope         = log10(meanslope_deg + x), # careful: original code uses 'meanslope_deg' here
    logrough         = log10(meanrough + x), 
    logsoc           = log10(wmeanSOCC_100CM + x),
    loggldist        = log10(gldistkm + x),
    lograin          = log10(RainTot96 + x)
  )


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


rdasimp_ocyield <- rda(Y ~ Condition(JDay) + pslumpact_log + logwy + loggpp,
               data = d, scale=TRUE)

## Global test of the RDA result
anova(rdasimp_ocyield, permutations = how(nperm = 5000))
## Tests of all canonical axes
anova(rdasimp_ocyield, by = "axis", permutations = how(nperm = 5000))
## Tests of all terms
anova(rdasimp_ocyield, by = "terms", permutations = how(nperm = 5000))

plot(rdasimp_ocyield)

out = varpart(Y, ~JDay, ~logwy + loggpp + pslumpact_log,
              data = d, scale=TRUE)

out = varpart(Y, ~JDay, ~logwy, ~ loggpp, ~pslumpact_log,
              data = d, scale=TRUE)
plot(out)
out
(R2adj <- RsquareAdj(rdasimp_ocyield)$r.squared)
(R2adj <- RsquareAdj(rdasimp_ocyield)$adj.r.squared)

#adonis2(Y ~ JDay + pslumpact_log + logwy + loggpp, data=d, permutations = 5000)

##### ========== (4) Plot RDA analysis (PANEL A) ==========================================================================
library(ggrepel) # Load ggrepel for non-overlapping labels

# Extract scores and figure data
scor <- scores(rdasimp_ocyield, display = c("sp", "cn", "bp", "lc"), scaling = 2) 

sites <- data.frame(scor$constraints) %>%
  mutate(site = d$site, slumpYN = d$slumpYN)

species_centroids <- data.frame(scor$species) %>%
  mutate(species_names = c("TOC", "POC", "DOC"))

arrows <- data.frame(scor$biplot) %>%
  mutate(pf_names = c("log(%ActiveSlumps+2)", "log(WaterYield+2)", "log(GPP+2)"))

# Calculate and format the Adjusted R-Squared label
R2adj_val <- RsquareAdj(rdasimp_ocyield)$adj.r.squared
panel_a_label <- paste0("(a) RDA Adj R² = ", round(R2adj_val * 100), "%")

# ---------------------------------------------------------
# PUBLICATION-READY THEME
# ---------------------------------------------------------
pub_theme <- theme_classic(base_size = 12, base_family = "sans") +
  theme(
    # Text formatting
    text = element_text(color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    
    # Axis lines and ticks
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    
    # Legend formatting (inside plot, top left, fully transparent)
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    
    # Panel & Background
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Margins & Ratio
    plot.margin = margin(t = 15, r = 15, b = 10, l = 10),
    aspect.ratio = 1
  )

# Calculate R2
R2adj_val <- RsquareAdj(rdasimp_ocyield)$adj.r.squared

# 1. Extract the summary of the RDA model
summ_a <- summary(rdasimp_ocyield)

# 2. Extract the 'Proportion Explained' for RDA1 and RDA2 (Row 2 of the importance matrix)
# Multiply by 100 and round to 1 decimal place
rda1_var_a <- round(summ_a$cont$importance[2, "RDA1"] * 100, 1)
rda2_var_a <- round(summ_a$cont$importance[2, "RDA2"] * 100, 1)

# 3. Build the dynamic axis labels
xlab_dynamic_a <- paste0("RDA1 (", rda1_var_a, "%)")
ylab_dynamic_a <- paste0("RDA2 (", rda2_var_a, "%)")

# Use bquote to construct the math expression for the title outside the plot
panel_a_title <- bquote("(a) "*R[adj]^2*" = "*.(round(R2adj_val * 100))*"%")

set.seed(42) # for repeatability of ggrepl label placement

ocyield <- ggplot(species_centroids, aes(x = RDA1, y = RDA2)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.5) + 
  
  # automatic plugin of R2adj and proportion explained by each axis
  labs(title = panel_a_title, x = xlab_dynamic_a, y = ylab_dynamic_a) +

  # repelled labels
  geom_text_repel(data = sites, aes(x = RDA1, y = RDA2, label = site), size = 3, colour = "grey40", max.overlaps = Inf, box.padding = 0.3) + 
  geom_text_repel(data = species_centroids, aes(label = species_names), colour = "blue4", size = 4.5, fontface = "bold") +
 
   coord_cartesian(x = c(-2, 1.5), y = c(-1.5, 1.5)) +
  
  # Data points
  geom_point(data = sites, aes(fill = slumpYN, shape = slumpYN), size = 4, colour = "black", stroke = 0.5) +
  scale_shape_manual(values = c("Y" = 21, "N" = 22)) +
  scale_fill_manual(values = c("Y" = "#FFB6C1", "N" = "#87CEEB")) + 
  
  # Environmental Arrows
  geom_segment(data = arrows, aes(x = 0, xend = RDA1, y = 0, yend = RDA2), arrow = arrow(length = unit(0.3, "cm"), type = "closed"), colour = "black", linewidth = 0.6) +
  
  # Repelled arrow labels
  # Swapped text_repel for label_repel with a semi-transparent white background
  geom_label_repel(data = arrows, aes(x = 1.15 * RDA1, y = 1.15 * RDA2, label = pf_names), 
                   size = 4, fontface = "italic", colour = "black", 
                   fill = alpha("white", 0.7), label.size = NA, box.padding = 0.5) +
  
  pub_theme +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 14)) # Style the new title


# Print the plot
ocyield

##### ========== (1) DATA PREP FOR OPTICS RDAS ==========================================================================

## (1.1) Read in file ====================

d <- read.csv(paste0(df, "2017data.csv"))

d <- d %>% select(site, date, campaign, loc, trans, tssmgL, POCmgL, DOCmgL, tssyield,
                  pocyield, tocyield, docyield,
                  delta18opermille,
                  Cayield,
                  wateryield, 
                  streampower,
                  slope,
                  meanslope_deg,
                  percshale,
                  colluvial_perc,
                  piedmont_perc, 
                  moraine_perc,
                  lakeperc,
                  scaledgpp,
                  forest_perc,
                  grassland_perc,
                  lichenmoss_perc,
                  shrubland_perc,
                  wmeanSOCC_100CM,
                  gldistkm,
                  meanrough,
                  RainTot96,
                  percslump17act, percslump17all,
                  slumpacccount, strahlerimpactacc,
                  PO13C,
                  POCTSSrat,
                  FemgL)

o <- read.csv(paste0(df, "masteroptics.csv"))

o <- o %>% filter(o$samptype=="Sample")

o <- o %>% select(-X, -dilfac, -samptype, -repnum)

## (1.2) Merge optics with data needed from d (tss conc) ====================

d <- merge(d, o, 
           by.x=c("site", "date", "loc"),
           by.y=c("site", "date", "streamloc"), all.x=TRUE)

# (1.3) Standardize atot to TSS ====================

d$atot250450TSSstd <- d$atot250450_p/d$tssmgL

# (1.4) Calculate Fe corrected SUVA254 ====================

#calculate absorbance due to Fe
## (from Poulin B A, Ryan J N and Aiken G R 2014 Effects of Iron on Optical Properties of Dissolved Organic Matter Environmental Science & Technology 48 10098-106)
d$a254Fe <- (0.0653*d$FemgL)+0.002

d$a254decfecor <- d$a254dec_d-d$a254Fe

d$SUVA254 <- d$a254decfecor/d$DOCmgL

## (1.5) Calculate julian date ====================

d$JDay <- julian(as.Date(d$date), origin = as.Date("2016-12-31"))
# julian() sets origin=0 in jday counts, so start the day before Jan. 1st (i.e Dec. 31st of prev. year)

## (1.6) Add Y/N for slump impact for RDA plotting ====================

d$slumpYN <- "N"
d$slumpYN[d$percslump17act>0] <- "Y"
d$slumpYN[d$site==29] <- "N"

## (1.7) Write for other plotting ====================

write.csv(d, paste0(df, "masteroptics2.csv"))

## (1.8) Fix site codes for RDA plotting ====================

d$site <- as.character(d$site)

d$site[d$site=="40-alt1"] <- "40"
d$site[d$site=="12-1"] <- "12"
d$site[d$site=="9-alt"] <- "9"

## (1.9) Filter to sites of interest ====================

d <- d %>% filter(d$site!="SC Outlet")

d$campaign[is.na(d$campaign)] <- "2017synoptic"
d <- d[d$campaign=="2017synoptic",]

# (1.10) Log transform the variables you know need to be transformed ====================

x <- 2

d <- d %>%
  mutate(
    # Chemistry (Ca, Mg, and Na strongly correlated with SO4)
    logcay           = log10(Cayield + x),
    
    # Geology & Soils
    shale_log        = log10(percshale + x),
    col_log          = log10(colluvial_perc + x),
    pied_log         = log10(piedmont_perc + x),
    moraine_log      = log10(moraine_perc + x),
    logsoc           = log10(wmeanSOCC_100CM + x),
    
    # Landcover & Biology
    loggpp           = log10(scaledgpp + x), 
    forest_log       = log10(forest_perc + x),
    grass_log        = log10(grassland_perc + x),
    lichen_log       = log10(lichenmoss_perc + x),
    lake_log         = log10(lakeperc + x),
    
    # Slump Metrics
    pslumpact_log    = log10(percslump17act + x),
    pslumpall_log    = log10(percslump17all + x),
    logslumpcount    = log10(slumpacccount + x),
    logslumpstrahler = log10(strahlerimpactacc + x),
    
    # Hydrology & Geomorphology
    logspower        = log10(streampower + x),
    logsslope        = log10(slope + x),
    logslope         = log10(meanslope_deg + x), 
    logwy            = log10(wateryield + x),
    lograin          = log10(RainTot96 + x)
  )

# (1.11) Check the new variables ====================

attach(d)
dotchart(log10(atot250450TSSstd+x))
dotchart(log10(POCTSSrat+x))
dotchart(PO13C)
dotchart((prcntC1_p))
dotchart(prcntC3_p)
dotchart(SUVA254)
dotchart(SR_d)
detach(d)

# (1.12) Transform necessary new variables ====================

d$atot250450TSSstd_log <- log10(d$atot250450TSSstd+x)
d$POCTSSrat_log <- log10(d$POCTSSrat+x)
#d$HIXp_log <- log10(d$HIX_p+x)

# (1.13) Select variables of interest ====================
d <- d %>%
  select(site, slumpYN,
         JDay,
         delta18opermille,
         logcay,
         logwy, 
         logspower,
         logsslope,
         shale_log, 
         col_log,
         pied_log, 
         moraine_log,
         lake_log,
         loggpp,
         forest_log,
         grass_log,
         lichen_log,
         logsoc,
         logslope,
         lograin,
         pslumpact_log,
         pslumpall_log,
         logslumpcount,
         logslumpstrahler,
         atot250450TSSstd_log,
         POCTSSrat,
         POCTSSrat_log, 
         PO13C,
         prcntC1_p,
         prcntC3_p,
         SUVA254 , 
         SR_d)

pairs(d[, c(21:ncol(d))], 
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)


##### ========== (2) Run RDA analysis 1 ==========================================================================

library(vegan)

d <- na.omit(d)

Y <- d %>% select(
  #atot250450TSSstd_log, was really only in Ch1 PCA because no %POC
  POCTSSrat_log,
  PO13C,
  prcntC1_p,
  #prcntC2_p, 
  prcntC3_p,
  #  prcntC4_p,
  #  prcntC5_p,
  #BIX_d,
  #HIX_d,
  #prcntC1_d,
  #prcntC2_d,
  #prcntC3_d,
  SUVA254,
  SR_d)

pairs(Y, 
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)

mod0p <- rda(Y ~ Condition(JDay), data = d, scale=TRUE)

mod1p <- rda(Y ~ delta18opermille +
               logcay +
               logwy + 
               logspower+
               logsslope +
               shale_log + 
               col_log +
               pied_log + 
               moraine_log +
               lake_log +
               loggpp +
               forest_log +
               grass_log +
               lichen_log +
               logsoc +
               logslope +
               lograin +
               pslumpact_log +
               # pslumpall_log +
               logslumpcount +
               # logslumpstrahler +
               #JDay,
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

step.forward <-
  ordiR2step(mod0p, 
             scope = formula(mod1p), 
             direction = "forward", 
             permutations = how(nperm = 5000),
             R2permutations = 5000)

RsquareAdj(step.forward)

rdasimp <- rda(Y ~ Condition(JDay) + logcay + loggpp,
               data = d, scale=TRUE)

plot(rdasimp, scaling=2)

## Global test of the RDA result
anova(rdasimp, permutations = how(nperm = 5000))
## Tests of all canonical axes
anova(rdasimp, by = "axis", permutations = how(nperm = 5000))
## Tests of all terms
anova(rdasimp, by = "terms", permutations = how(nperm = 5000))

out = varpart(Y, ~JDay, ~logcay + loggpp,
              data = d, scale=TRUE)

out = varpart(Y, ~JDay, ~logcay, ~loggpp,
              data = d, scale=TRUE)
#plot(out)
out

##### ========== (3) Plot RDA analysis (PANEL B) ================================================================
library(ggrepel)

# Extract scores and data for Panel B
scor_b <- scores(rdasimp, display=c("sp", "cn", "bp", "lc"), scaling=2) 

sites_b <- data.frame(scor_b$constraints) %>%
  mutate(site = d$site, slumpYN = d$slumpYN)

species_centroids_b <- data.frame(scor_b$species) %>%
  mutate(
    species_names = c("log(%POC+2)", "PO13C", "P1", "P3", "DOM-SUVA", "DOM-SR"),
    OCfrac = c("POM", "POM", "POM", "POM", "DOM", "DOM")
  )

arrows_b <- data.frame(scor_b$biplot) %>%
  mutate(pf_names = c("log(CaYield+2)", "log(GPP+2)"))

# Calculate Adjusted R-Squared label for Panel B
R2adj_val_b <- RsquareAdj(rdasimp)$adj.r.squared
panel_b_label <- paste0("(b) RDA Adj R² = ", round(R2adj_val_b * 100), "%")

# 1. Extract the summary of the RDA model
summ_a <- summary(rdasimp)

# 2. Extract the 'Proportion Explained' for RDA1 and RDA2 (Row 2 of the importance matrix)
# Multiply by 100 and round to 1 decimal place
rda1_var_a <- round(summ_a$cont$importance[2, "RDA1"] * 100, 1)
rda2_var_a <- round(summ_a$cont$importance[2, "RDA2"] * 100, 1)

# 3. Build the dynamic axis labels
xlab_dynamic_a <- paste0("RDA1 (", rda1_var_a, "%)")
ylab_dynamic_a <- paste0("RDA2 (", rda2_var_a, "%)")

# Calculate R2 and construct math expression
R2adj_val_b <- RsquareAdj(rdasimp)$adj.r.squared
panel_b_title <- bquote("(b) "*R[adj]^2*" = "*.(round(R2adj_val_b * 100))*"%")

set.seed(42)

panel_b <- ggplot(species_centroids_b, aes(x = RDA1, y = RDA2)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.5) +
  
  # dynamic labels
  labs(title = panel_b_title, x = xlab_dynamic_a, y = ylab_dynamic_a) +
  
  geom_text_repel(data = sites_b, aes(x = RDA1, y = RDA2, label = site), size = 3, colour = "grey40", max.overlaps = Inf, box.padding = 0.3) + 
  geom_point(data = sites_b, aes(fill = slumpYN, shape = slumpYN), size = 4, colour = "black", stroke = 0.5) +
  geom_point(aes(colour = OCfrac), size = 2.5, shape = 17) +
  geom_text_repel(aes(label = species_names, colour = OCfrac), size = 4.5, fontface = "bold", show.legend = FALSE) +
  scale_shape_manual(values = c("Y" = 21, "N" = 22)) +
  scale_fill_manual(values = c("Y" = "#FFB6C1", "N" = "#87CEEB")) + 
  scale_colour_manual(values = c("POM" = "brown", "DOM" = "blue4")) +
  geom_segment(data = arrows_b, aes(x = 0, xend = RDA1, y = 0, yend = RDA2), arrow = arrow(length = unit(0.3, "cm"), type = "closed"), colour = "black", linewidth = 0.6) +
  
  # Use label_repel with background
  geom_label_repel(data = arrows_b, aes(x = 1.15 * RDA1, y = 1.15 * RDA2, label = pf_names), 
                   size = 4, fontface = "italic", colour = "black", 
                   fill = alpha("white", 0.7), label.size = NA, box.padding = 0.5) +
  
  pub_theme +
  theme(legend.position = "none", plot.title = element_text(hjust = 0, face = "bold", size = 14))

##### ========== (4) Run RDA analysis 2 ==========================================================================

d$slumpYNcode <- 0
d$slumpYNcode[d$slumpYN=="Y"] <- 1

mod0p2 <- rda(Y ~ Condition(JDay), data = d, scale=TRUE)

mod1p2 <- rda(Y ~ delta18opermille +
                logcay +
                logwy + 
                logspower+
                logsslope +
                shale_log + 
                col_log +
                pied_log + 
                moraine_log +
                lake_log +
                loggpp +
                forest_log +
                grass_log +
                lichen_log +
                logsoc +
                logslope +
                lograin +
                pslumpact_log +
                # pslumpall_log +
                logslumpcount +
                # logslumpstrahler +
                as.factor(slumpYN) + 
                #JDay,
                Condition(JDay),
              data = d, scale=TRUE)


## Global test of the RDA result
anova(mod1p2, permutations = how(nperm = 5000))
## Tests of all canonical axes
anova(mod1p2, by = "axis", permutations = how(nperm = 5000))

plot(mod1p2)

step.forward2 <-
  ordiR2step(mod0p2, 
             scope = formula(mod1p2), 
             direction = "forward", 
             permutations = how(nperm = 5000),
             R2permutations = 5000)


rdasimp2 <- rda(Y ~Condition(JDay) + as.factor(slumpYN) + 
                  logcay + moraine_log + logspower + loggpp ,
                data = d, scale=TRUE)

plot(rdasimp2)


## Global test of the RDA result
anova(rdasimp2, permutations = how(nperm = 5000))
## Tests of all canonical axes
anova(rdasimp2, by = "axis", permutations = how(nperm = 5000))

out2 = varpart(Y, ~JDay, ~as.factor(slumpYN) + logcay + moraine_log + logspower + loggpp,
               data = d, scale=TRUE)

slump = varpart(Y, ~JDay, ~as.factor(slumpYN), ~ logcay + moraine_log + logspower + loggpp,
                data = d, scale=TRUE)

cay = varpart(Y, ~JDay, ~logcay, ~as.factor(slumpYN) + moraine_log + logspower + loggpp,
              data = d, scale=TRUE)

spower = varpart(Y, ~JDay, ~logspower, ~as.factor(slumpYN) + logcay + moraine_log + loggpp,
                 data = d, scale=TRUE)

gpp= varpart(Y, ~JDay, ~ loggpp, ~as.factor(slumpYN) + logcay + moraine_log + logspower,
             data = d, scale=TRUE)

mor= varpart(Y, ~JDay, ~moraine_log, ~as.factor(slumpYN) + logcay + loggpp + logspower,
             data = d, scale=TRUE)


##### ========== (5) Plot RDA analysis 2 (PANEL C) ==============================================================

# Extract scores and data for Panel C
scor_c <- scores(rdasimp2, display=c("sp", "cn", "bp", "lc"), scaling=2) 

sites_c <- data.frame(scor_c$constraints) %>%
  mutate(site = d$site, slumpYN = d$slumpYN)

species_centroids_c <- data.frame(scor_c$species) %>%
  mutate(
    species_names = c("log(%POC+2)", "PO13C", "P1", "P3", "DOM-SUVA", "DOM-SR"),
    OCfrac = c("POM", "POM", "POM", "POM", "DOM", "DOM")
  )

arrows_c <- data.frame(scor_c$biplot) %>%
  mutate(pf_names = c("factor", "log(CaYield+2)", "log(%moraine+2)", "log(spower+2)", "log(GPP+2)")) %>%
  filter(pf_names != "factor") # Remove factor so it plots as centroid

factorcent_c <- data.frame(scor_c$centroids) %>%
  mutate(names = c("NO RTS", "RTS"))

# 1. Extract the summary of the RDA model
summ_a <- summary(rdasimp2)

# 2. Extract the 'Proportion Explained' for RDA1 and RDA2 (Row 2 of the importance matrix)
# Multiply by 100 and round to 1 decimal place
rda1_var_a <- round(summ_a$cont$importance[2, "RDA1"] * 100, 1)
rda2_var_a <- round(summ_a$cont$importance[2, "RDA2"] * 100, 1)

# 3. Build the dynamic axis labels
xlab_dynamic_a <- paste0("RDA1 (", rda1_var_a, "%)")
ylab_dynamic_a <- paste0("RDA2 (", rda2_var_a, "%)")

# Calculate R2 and construct math expression
R2adj_val_c <- RsquareAdj(rdasimp2)$adj.r.squared
panel_c_title <- bquote("(c) "*R[adj]^2*" = "*.(round(R2adj_val_c * 100))*"%")

set.seed(42)

panel_c <- ggplot(species_centroids_c, aes(x = RDA1, y = RDA2)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.5) +
  
  # Add the title
  labs(title = panel_c_title, x = "RDA1 (53.8%)", y = "RDA2 (11.0%)") +
  
  geom_text_repel(data = sites_c, aes(x = RDA1, y = RDA2, label = site), size = 3, colour = "grey40", max.overlaps = Inf, box.padding = 0.3) + 
  geom_point(data = sites_c, aes(fill = slumpYN, shape = slumpYN), size = 4, colour = "black", stroke = 0.5) +
  geom_point(aes(colour = OCfrac), size = 2.5, shape = 17) +
  geom_text_repel(aes(label = species_names, colour = OCfrac), size = 4.5, fontface = "bold", show.legend = FALSE) +
  geom_point(data = factorcent_c, aes(x = RDA1, y = RDA2), size = 3, shape = 17, colour = "black") +
  geom_text_repel(data = factorcent_c, aes(x = RDA1, y = RDA2, label = names), size = 4.5, fontface = "bold", colour = "black") +
  scale_shape_manual(values = c("Y" = 21, "N" = 22)) +
  scale_fill_manual(values = c("Y" = "#FFB6C1", "N" = "#87CEEB")) + 
  scale_colour_manual(values = c("POM" = "brown", "DOM" = "blue4")) +
  geom_segment(data = arrows_c, aes(x = 0, xend = RDA1, y = 0, yend = RDA2), arrow = arrow(length = unit(0.3, "cm"), type = "closed"), colour = "black", linewidth = 0.6) +
  
  # Use label_repel with background
  geom_label_repel(data = arrows_c, aes(x = 1.15 * RDA1, y = 1.15 * RDA2, label = pf_names), 
                   size = 4, fontface = "italic", colour = "black", 
                   fill = alpha("white", 0.7), label.size = NA, box.padding = 0.5) +
  
  pub_theme +
  theme(legend.position = "none", plot.title = element_text(hjust = 0, face = "bold", size = 14))

##### ============================== Section 6: Export Multi-Panel Plot =========================================
# Since you loaded the 'patchwork' library earlier, you can ditch the clunky 'grid' code entirely.
# This binds them side-by-side (A + B + C) and exports a single, perfectly aligned PDF.
# Note: 'ocyield' is assumed to be Panel A from your previous code block.

combined_plot <- ocyield + panel_b + panel_c + plot_layout(ncol = 1)

if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Figures/Fig5/")) dir.create("Figures/Fig5/", recursive = TRUE)

# Save the combined figure as a high-resolution PDF
ggsave("Figures/Fig5/Figure5_RDAs.pdf", plot = combined_plot, width = 6, height = 15, units = "in")
ggsave("Figures/Fig5/Figure5_RDAs.png", plot=combined_plot, width = 6, height = 15, units="in", dpi = 300)
