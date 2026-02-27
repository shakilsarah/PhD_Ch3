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

##### ========== (1) DATA PREP ==========================================================================

## (1.1) Read in file ====================
df <- "ThesisDrafts/Chapter3/Data/" 

d <- read.csv(paste0(wd, df, "synopticflux.csv"))

d$wateryield <- (d$dism3s*1000)/d$WatershedArea

## (1.2) Select variables ====================
d <- d%>%
  select(site, 
         date,
         tocyield,
         pocyield,
         docyield,
         percPOCy,
         percDOCy,
         wateryield,
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
         percslump, slumpacccount,
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

##### ========== (2) Conduct linear analysis ==========================================================================

# (2.1) Step 1 : Outliers

dotchart(d$tocyield) # site 0 might need to be removed as an outlier
dotchart(d$pocyield) # site 0 might need to be removed as an outlier
dotchart(d$docyield)

dotchart(d$percDOCy)
dotchart(d$percPOCy)

dotchart(d$WatershedArea) # site 10 potentially outlier

# geology
dotchart(d$percshale)
dotchart(d$colluvial_perc)
dotchart(d$piedmont_perc)
dotchart(d$alluvial_perc) # alot of zeros
dotchart(d$bedrock_perc) # all zeros
dotchart(d$fluvial_perc) # only 1 not zero
dotchart(d$glaciogenic_perc) # only 2 not zero
dotchart(d$organic_perc) # only 4 not zero, site 38 might be an outlier
dotchart(d$moraine_perc)


dotchart(d$lakeperc) # keep an eye on site 35 as an outlier, but honestly doesn't seem hugely concerning
dotchart(d$scaledgpp)
dotchart(d$meanelev_m)
dotchart(d$meanslope_deg)

#landcover
dotchart(d$barrenland_perc)
dotchart(d$forest_perc) # site 57 might be an outlier
dotchart(d$grassland_perc)
dotchart(d$lichenmoss_perc)
dotchart(d$shrubland_perc) # site 12-1 potentially an outlier
dotchart(d$wetland_perc) # I think remove wetland because essentially no wetland coverage (only 2 sites are non-zero)

# slump
dotchart(d$percslump)
dotchart(d$slumpacccount)

#SOC
dotchart(d$wmeanSOCC_100CM) # 45 slightly but log transforming doesn't make a difference and doesn't seem a big concern

#glacial distance
dotchart(d$gldistkm)

#terrain roughness
dotchart(d$meanrough)

#Rain
dotchart(log10(d$RainTot24+1))
dotchart(log10(d$RainTot48+1))
dotchart(log10(d$RainTot72+1))
dotchart(log10(d$RainTot96+1)) # very abnormal
plot(d$pocyield~d$RainTot96)

dotchart(d$wateryield)

# (2.1) Step 2: Transformations for outlier data and removal of variables lacking data

dotchart(log10(d$tocyield+1))
d$logtocy <- log10(d$tocyield)

dotchart(log10(d$pocyield))
d$logpocy <- log10(d$pocyield)

d$logdocy <- log10(d$docyield)

d$logwy <- log10(d$wateryield)

# arcsine percent transformations
d$shale_as <- asin(sqrt(((d$percshale)/100)))
d$col_as <- asin(sqrt(((d$colluvial_perc)/100)))
d$pied_as <- asin(sqrt(((d$piedmont_perc)/100)))
d$moraine_as <- asin(sqrt(((d$moraine_perc)/100)))
d$barren_as <- asin(sqrt(((d$barrenland_perc)/100)))
d$grass_as <-  asin(sqrt(((d$grassland_perc)/100)))
d$lichen_as <- asin(sqrt(((d$lichenmoss_perc)/100)))
d$pslump_as <- asin(sqrt(((d$percslump)/100)))
d$lake_as <- asin(sqrt(((d$lakeperc)/100)))
d$forest_as <- asin(sqrt(((d$forest_perc)/100)))
d$shrub_as <- asin(sqrt(((d$shrubland_perc)/100)))
d$pPOC_as <- asin(sqrt(((d$percPOCy)/100)))

d <- d %>% 
  select(site,
         logtocy, logpocy, logdocy, pPOC_as, JDay,
         logwy,
         shale_as,
         col_as, pied_as, moraine_as,
         lake_as, scaledgpp, 
         meanelev_m, meanslope_deg, 
         barren_as, forest_as, grass_as, lichen_as,
         shrub_as, pslump_as, slumpacccount,
         wmeanSOCC_100CM, gldistkm, meanrough, RainTot96)

pairs(d[, c(6:ncol(d))], 
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)

d <- d %>%
    select(-meanslope_deg, -barren_as, -meanelev_m, -scaledgpp)


corvif(d[, c(6:ncol(d))])

d <- d %>%
  select(-col_as, -RainTot96, -shale_as, -shrub_as) 

pairs(d[, c(6:ncol(d))], 
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)

corvif(d[, c(6:ncol(d))])

d <- d %>%
  select(-moraine_as) 

corvif(d[, c(6:ncol(d))])

d <- d %>%
  select(-forest_as) 

corvif(d[, c(6:ncol(d))])

pairs(d[, c(2:ncol(d))], 
      lower.panel=panel.smooth2, 
      upper.panel=panel.cor,
      diag.panel=panel.hist)

# run linear model

# model fit is singular for JDay as random effect, therefore can't use lmer

toc <- lm(logtocy ~ pied_as + lake_as + grass_as + lichen_as + pslump_as +
               slumpacccount + wmeanSOCC_100CM + gldistkm + meanrough + JDay + logwy,
          data=d)

step(toc, direct="backward")
drop1(toc)


poc <- lm(logpocy ~ pied_as + lake_as + grass_as + lichen_as + pslump_as +
            slumpacccount + wmeanSOCC_100CM + gldistkm + meanrough + JDay + logwy,
          data=d)

step(poc, direct="backward")
drop1(poc)
plot(poc)


doc <- lm(logdocy ~ pied_as + lake_as + grass_as + lichen_as + pslump_as +
            slumpacccount + wmeanSOCC_100CM + gldistkm + meanrough + JDay + logwy,
          data=d)
step(doc, direct="backward")
drop1(doc)
plot(doc)

# could potentially proceed after removing moraine and shrub since VIF is less than 10 then.....

#d <- d %>% 
#    select(-colluvial_perc, # 0.9 corr with elev and slope
#          -meanelev_m, # 0.9 cor with slope
#         -barrenland_perc, # 0.9 cor with slope
#        -logshrub, # 0.9 cor with gpp
#       -WatershedArea,
#      -meanslope_deg)

#d <- d %>% select(-meanslope_deg,
#                 -percshale,
#                -grassland_perc,
#               -lichenmoss_perc,
#              -logforest) # remove all land cover might be causing VIF to be too high

# too many variables are getting removed, try running a pca to assess gradients

##### ========== (2)  ==========================================================================

# run a PCA to find common gradients

exp <- d 

## DCA
DCA <- exp[,2:ncol(exp)]
DCA$gldistkm <- DCA$gldistkm+2
comp.dca <- decorana(DCA) 
summary(comp.dca)
## axis lengths are less than 1.5

# PCA based on a correlation matrix
# Argument scale=TRUE calls for a standardization of the variables

pca <- rda(exp[,2:ncol(exp)], scale=TRUE)
pca
summary(pca) # Default scaling is 2
summary(pca, scaling=1)
# to remove the sites and species scores from the summary() command, add axes=0 argument

# Examine and plot partial results from PCA output
?cca.object
# explains how an ordination object produced by vegan 
# produced by vegan is structured and how to extract its results

# Eigenvalues
ev <- pca$CA$eig # extract eigen values from pca

# Scree plot and broken stick model
screeplot(pca, bstick=TRUE, npcs=length(pca$CA$eig))
PCAsignificance(pca, axes=4)
# Plot PCA results - Ryan's version ###

# extract pca scores, sites gives scores per site, species gives scores per species
POMscors1 = scores(pca, display=c("sites", "species"), scaling=1, choices=c(1,2,3)) 
POMscors2 = scores(pca, display=c("sites", "species"), scaling=2, choices=c(1,2,3)) 


### PCA Arrows ======

# create arrows for PCA plot - scaling 1
POMarrows1 <- data.frame(POMscors1$species)
POMarrows1$comp_metr <- rownames(POMarrows1)
#POMarrows1$comp_metr <- colnames(compf[,3:ncol(compf)])
POMarrows1

# create arrows for PCA plot - scaling 2
POMarrows2 <- data.frame(POMscors2$species)
POMarrows2$comp_metr <- rownames(POMarrows2)
#POMarrows2$comp_metr <- colnames(compf[,3:ncol(compf)])
POMarrows2

### PCA Site Centroids ======

# create centroids for sites - scaling 1
POMsites1 <- data.frame(POMscors1$sites)
POMsites1$site <- exp$site
POMsites1

# create centroids for sites - scaling 2
POMsites2 <- data.frame(POMscors2$sites)
POMsites1$site <- exp$site
POMsites2


### Graph Themes =====
theme<-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             panel.background = element_blank(),axis.line.x = element_line(colour="black"),
             axis.line.y = element_line(colour="black"),
             axis.text = element_text(colour="black",size=14),legend.background=element_blank(),
             text=element_text(size = 18),
             legend.title=element_blank(),
             legend.position = c(1, 0), 
             legend.justification = c(1, 0),
             plot.title = element_text(margin = margin(b = -20))
)


### Plottin PCA results ====

library(ggplot2)
# to extend ggplot to 3D plotting
#library(gg3D)
# code below for 3D use of plotly: https://stackoverflow.com/questions/45052188/how-to-plot-3d-scatter-diagram-using-ggplot

## plot PCA results using scaling 1

### Prep Circle of equilibrium ====

## obtain stretching factor for vegan from cleanplot.pca (Numerical Ecology in R, 2018)
cleanplot.pca(pca, scaling=1, silent=FALSE) # silent = FALSE prints the computation results to provide radius

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

dat <- circleFun(c(0,0),((((2/19)^0.5)*4.714458)*2),npoints = 100)
#geom_path will do open circles, geom_polygon will do filled circles
# 3.855654 is the general scaling constant used internally by vegan (see summary(pca))

pcagraph1 <- ggplot(POMsites1, aes(x = PC1, y= PC2)) +
  geom_point(data = POMsites1) +
  geom_text(data = POMsites1, aes(label=site))+
  geom_segment(data = POMarrows1,
               aes(x = 0, xend = (PC1),
                   y = 0, yend = (PC2)),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black")+
  geom_text(data = POMarrows1,
            aes(x= 1.5*PC1, y = 1.5*PC2, #we add 10% to the text to push it slightly out from arrows
                label = comp_metr), #otherwise you could use hjust and vjust. I prefer this option
            size = 6,
            hjust = 0.5)+
  geom_path(data=dat, aes(x,y))+
  geom_vline(xintercept = 0,linetype="dashed")+ylab("PCA 2 (51%)")+
#  scale_y_continuous(limits=c(-2.5,2.5), breaks=c(-2, -1, 0, 1, 2)) + 
 # scale_x_continuous(limits=c(-2.5,2.5), breaks=c(-2, -1, 0, 1, 2)) +
  geom_hline(yintercept = 0,linetype="dashed")+xlab("PCA 1 (31%)") +
  ggtitle("B. distance-biplot") +
  theme + theme(legend.position="none")

