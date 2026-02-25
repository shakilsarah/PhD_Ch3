# ===========================================================================================================#
# 2017RDAoptics.R
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

d <- read.csv(paste0(wd, df, "2017data.csv"))

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

o <- read.csv(paste0(wd, df, "masteroptics.csv"))

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

write.csv(d, paste0(wd, df, "masteroptics2.csv"))

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

d$logcay <- log10(d$Cayield+x) # Ca, Mg, and Na were strongly correlated with SO4

d$shale_log <- log10(d$percshale+x)
d$col_log <- log10(d$colluvial_perc+x)
d$pied_log <- log10(d$piedmont_perc+x)
d$moraine_log <- log10(d$moraine_perc+x)
d$logsoc <- log10(d$wmeanSOCC_100CM+x)

d$loggpp <- log10(d$scaledgpp+x) 
d$forest_log <- log10(d$forest_perc+x)
d$grass_log <-  log10(d$grassland_perc+x)
d$lichen_log <- log10(d$lichenmoss_perc+x)
d$lake_log <- log10(d$lakeperc+x)

d$pslumpact_log <- log10(d$percslump17act+x)
d$pslumpall_log <- log10(d$percslump17all+x)
d$logslumpcount <- log10(d$slumpacccount+x)
d$logslumpstrahler <- log10(d$strahlerimpactacc+x)

d$logspower <- log10(d$streampower+x)
d$logsslope <- log10(d$slope+x)

d$logslope <- log10(d$meanslope_deg+x) 
d$logwy <- log10(d$wateryield+x)
d$lograin <- log10(d$RainTot96+x)

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


##### ========== (3) Plot RDA analysis ==========================================================================

# scores and figure

scor = scores(rdasimp, display=c("sp", "cn", "bp", "lc"), scaling=2) 

sites <- data.frame(scor$constraints)
sites$site <- d$site
sites$slumpYN <- d$slumpYN

species_centroids <- data.frame(scor$species)
species_centroids
species_centroids$species_names <- c(
                                      "log(%POC+2)",
                                      "PO13C","P1", "P3", "DOM-SUVA", "DOM-SR")
species_centroids$OCfrac <- c("POM", "POM", "POM", "POM", "DOM", "DOM")

arrows <- data.frame(scor$biplot)
arrows$pf_names <- c("log(CaYield+2)","log(GPP+2)")
arrows

mult <- attributes(scores(rdasimp))$const

theme<-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             panel.background = element_blank(),axis.line.x = element_line(colour="black"),
             axis.line.y = element_line(colour="black"),
             axis.text = element_text(colour="black",size=14),legend.background=element_blank(),
             text=element_text(size = 16),
             legend.position="none",
             #legend.title=element_blank(),
            # legend.position = c(0,1),
             #legend.direction="horizontal",
             #legend.justification = c(0,1),
             aspect.ratio=1)

rdagraph <- ggplot(data = species_centroids, 
                   aes(x = RDA1, y= RDA2)) +
  geom_vline(xintercept = 0,linetype="dashed", colour="grey")+ylab("RDA2 (4.4%)")+
  geom_hline(yintercept = 0,linetype="dashed", colour="grey")+xlab("RDA1 (39.9%)")+
  geom_point(data = sites, 
             aes(fill=slumpYN, shape=slumpYN), 
             size=4, colour="white") +
  geom_text(data = sites, 
            aes(x= RDA1, y = RDA2-0.15, label=site), 
            size=2.5, colour="grey40") +
  geom_point(aes(colour=OCfrac),
            size = 2, shape=17)+
  geom_text(data = species_centroids, 
            aes(x = RDA1*1.1, y= RDA2*1.1,
                label = species_names,  colour=OCfrac),
            size = 4)+
 # coord_cartesian(x = c(-2, 1.5), y = c(-1.5, 1.5))+
  scale_shape_manual(limits= c("Y", "N"),
                     breaks= c("Y", "N"),
                     values= c(21, 22)) +
  scale_fill_manual(limits= c("Y", "N"),
                      breaks= c("Y", "N"),
                      values= c("pink", "Sky Blue")) +
  scale_colour_manual(limits= c("POM", "DOM"),
                    breaks= c("POM", "DOM"),
                    values= c("brown", "blue")) +
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


##### ========== (5) Plot RDA analysis 2==========================================================================

# scores and figure

scor = scores(rdasimp2, display=c("sp", "cn", "bp", "lc"), scaling=2) 

sites <- data.frame(scor$constraints)
sites$site <- d$site
sites$slumpYN <- d$slumpYN

species_centroids <- data.frame(scor$species)
species_centroids
species_centroids$species_names <- c("log(%POC+2)","PO13C","P1", "P3", 
                                     "DOM-SUVA", "DOM-SR")
species_centroids$OCfrac <- c("POM", "POM", "POM", "POM", "DOM", "DOM")

arrows <- data.frame(scor$biplot)
arrows$pf_names <- c("factor", "log(CaYield+2)", 
                     "log(%moraine+2)","log(spower+2)",
                     "log(GPP+2)")
arrows <- arrows[arrows$pf_names!="factor",] # want to plot this as a centroid not an arrow

factorcent <- data.frame(scor$centroids)
factorcent$names <- c("NO RTS", "RTS")

mult <- attributes(scores(rdasimp2))$const

theme<-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             panel.background = element_blank(),axis.line.x = element_line(colour="black"),
             axis.line.y = element_line(colour="black"),
             axis.text = element_text(colour="black",size=14),legend.background=element_blank(),
             text=element_text(size = 16),
             legend.position="none",
             #legend.title=element_blank(),
             # legend.position = c(0,1),
             #legend.direction="horizontal",
             #legend.justification = c(0,1),
             aspect.ratio=1)

rdagraph2 <- ggplot(data = species_centroids, 
                   aes(x = RDA1, y= RDA2)) +
  geom_vline(xintercept = 0,linetype="dashed", colour="grey")+ylab("RDA2 (11.0%)")+
  geom_hline(yintercept = 0,linetype="dashed", colour="grey")+xlab("RDA1 (53.8%)")+
  geom_point(data = sites, 
             aes(fill=slumpYN, shape=slumpYN), 
             size=4, colour="white") +
  geom_text(data = sites, 
            aes(x= RDA1, y = RDA2-0.15, label=site), 
            size=2.5, colour="grey40") +
  geom_point(aes(colour=OCfrac),
             size = 2, shape=17)+
  geom_text(data = species_centroids, 
            aes(x = RDA1*1.1, y= RDA2*1.1,
                label = species_names,  colour=OCfrac),
            size = 4)+
  geom_point(data=factorcent, aes(x=RDA1, y=RDA2),
             size = 2.5, shape=17, colour="black")+
  geom_text(data = factorcent, 
            aes(x = RDA1*1.1, y= RDA2*1.1,
                label = names),
            size = 4, colour="black")+
  # coord_cartesian(x = c(-2, 1.5), y = c(-1.5, 1.5))+
  scale_shape_manual(limits= c("Y", "N"),
                     breaks= c("Y", "N"),
                     values= c(21, 22)) +
  scale_fill_manual(limits= c("Y", "N"),
                    breaks= c("Y", "N"),
                    values= c("pink", "Sky Blue")) +
  scale_colour_manual(limits= c("POM", "DOM"),
                      breaks= c("POM", "DOM"),
                      values= c("brown", "blue")) +
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



##### ============================== Section 4: Export plot ===================================================
library(grid)


savePlotpdf <- function(myPlot) {
  pdf(file = paste0(wd, "ThesisDrafts/Chapter3/Graphs/Fig3/", filename, ".pdf"),
      width = 4, height = 4)
  print(myPlot)
  dev.off()
}

g1 <- ggplotGrob(rdagraph)
grid.draw(g1)
filename <- "rdaoccomp"
savePlotpdf(grid.draw(g1))

g1 <- ggplotGrob(rdagraph2)
grid.draw(g1)
filename <- "rdaoccomp2"
savePlotpdf(grid.draw(g1))
