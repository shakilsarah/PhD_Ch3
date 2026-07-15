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

df <- "Data/"

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
source("functions/HighstatLibV10.R")

##### ========== (1) DATA PREP ==========================================================================

## (1.1) Read in file ====================

d <- read.csv(paste0(df, "2017data.csv"))

## (1.11) Update the surificial geology 


# read in the updated watmaster file 
watmaster <- read.csv("D:/5_Projects/git/repos/Phd_Ch3/data/watmaster_wlakes_2026.csv")

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
  distinct(site, date, trans) %>%
  anti_join(
    watmaster %>% distinct(site, date, trans),
    by = c("site", "date", "trans")
  )

# View missing combinations
missing_site_dates
# yes - there are two. 


# Subset watmaster to only site-date combinations present in d
watmaster_subset <- watmaster %>%
  semi_join(
    d %>% select(site, date, trans),
    by = c("site", "date", "trans")
  ) %>%
  select(all_of(watmaster_cols))

# Remove old columns from d, then add selected watmaster columns by site and date
d_updated <- d %>%
  select(-any_of(cols_to_remove_from_d)) %>%
  left_join(
    watmaster_subset,
    by = c("site", "date", "trans")
  )

# Check result
dim(d_updated)
head(d_updated)

# Optional: check whether watmaster has duplicate site-date combinations
# Doesn't matter because none of the SE samples are used here 
watmaster %>%
  count(site, date, trans) %>%
  filter(n > 1)

## Rename d_updated to d and continue 

d <- d_updated

# Continue 

d <- d %>% select(site, date, campaign, loc, trans, tssmgL, 
                  POCmgL, DOCmgL, tssyield,
                  pocyield, tocyield, docyield,
                  delta18opermille,
                  Cayield,
                  wateryield, 
                  streampower,
                  slope,
                  meanslope_deg,
                  percshale,
                  colluvial_perc,
                 # piedmont_perc, 
                  morainal_perc,
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

d$logcay <- log10(d$Cayield+x) # Ca, Mg, and Na were strongly correlated with SO4

d$shale_log <- log10(d$percshale+x)
d$colluvial_log <- log10(d$colluvial_perc+x)
#d$pied_log <- log10(d$piedmont_perc+x)
d$morainal_log <- log10(d$morainal_perc+x)
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
             colluvial_log,
           #  pied_log, 
             morainal_log,
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

# fit the  null model
mod0p <- rda(Y ~ Condition(JDay), data = d, scale=TRUE)

# fit the full model
mod1p <- rda(Y ~ delta18opermille +
               logcay +
               logwy + 
               logspower+
               logsslope +
               shale_log + 
               colluvial_log +
               #pied_log +    # bye bye 
               morainal_log +
               lake_log +
               loggpp +
               forest_log +
               grass_log +
               lichen_log +
               logsoc +
              # logslope +
               lograin +
               pslumpact_log +
               # pslumpall_log +
               logslumpcount +
               # logslumpstrahler +
               #JDay,
               Condition(JDay),
             data = d, scale=TRUE)


rda <- mod1p

x11()

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

rdasimp_morainal <- rda(Y ~ Condition(JDay) + logcay + morainal_log,
               data = d, scale=TRUE)


plot(rdasimp, scaling=2)
plot(rdasimp_morainal, scaling=2)

## Global test of the RDA result
anova(rdasimp, permutations = how(nperm = 5000))
anova(rdasimp_morainal, permutations = how(nperm = 5000))

## Tests of all canonical axes
anova(rdasimp, by = "axis", permutations = how(nperm = 5000))
anova(rdasimp_morainal, by = "axis", permutations = how(nperm = 5000))

## Tests of all terms
anova(rdasimp, by = "terms", permutations = how(nperm = 5000))
anova(rdasimp_morainal, by = "terms", permutations = how(nperm = 5000))

## Test marginal -- each predictor still considered sig when considered uniquely? 
anova(rdasimp, by = "margin", permutations = how(nperm = 5000))
anova(rdasimp_morainal, by = "margin", permutations = how(nperm = 5000))


# based on all the data, proceed with gpp. 

#out = varpart(Y, ~JDay, ~logcay + loggpp,
#              data = d, scale=TRUE)

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

#arrows$pf_names <- c("log(CaYield+2)", "log(%morainal+2)")

arrows

mult <- attributes(scores(rdasimp))$const

theme<-theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             axis.line.x = element_line(colour="black"),
             axis.line.y = element_line(colour="black"),
             axis.text = element_text(colour="black",size=14),
             legend.background=element_blank(),
             text=element_text(size = 16),
             legend.position="none",
             #legend.title=element_blank(),
            # legend.position = c(0,1),
             #legend.direction="horizontal",
             #legend.justification = c(0,1),
             aspect.ratio=1)

rdagraph <- ggplot(data = species_centroids, 
                   aes(x = RDA1, y= RDA2)) +
  geom_vline(xintercept = 0,linetype="dashed", 
             colour="grey")+ylab("RDA2 (4.4%)")+
  
  geom_hline(yintercept = 0,linetype="dashed", 
             colour="grey")+xlab("RDA1 (39.9%)")+
  
  geom_point(data = sites, 
             aes(fill=slumpYN, shape=slumpYN), 
             size=4, colour="black") +
  
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
               arrow = arrow(length = unit(0.2, "cm")), colour = "grey10")+
  geom_text(data = arrows,
            aes(x= 1.2*RDA1, y = 1.2*RDA2, #we add 10% to the text to push it slightly out from arrows
                label = pf_names), #otherwise you could use hjust and vjust. I prefer this option
            size = 4,
            hjust = 0.5)+
  theme


rdagraph

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
                colluvial_log +
                #pied_log + 
                morainal_log +
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

rda2 <- mod1p2 

(R2adj <- RsquareAdj(rda2)$r.squared)
(R2adj <- RsquareAdj(rda2)$adj.r.squared)

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

#original
rdasimp2 <- rda(Y ~Condition(JDay) + as.factor(slumpYN) + 
                  logcay + morainal_log + logspower + loggpp ,
                data = d, scale=TRUE)
# updated
rdasimp2_new <- rda(Y ~Condition(JDay) + as.factor(slumpYN) + 
                  logcay + logwy + loggpp,
                data = d, scale=TRUE)

RsquareAdj(rdasimp2_new)

plot(rdasimp2)
plot(rdasimp2_new)

## Global test of the RDA result
anova(rdasimp2, permutations = how(nperm = 5000))
anova(rdasimp2_new, permutations = how(nperm = 5000))

## Tests of all canonical axes
anova(rdasimp2, by = "axis", permutations = how(nperm = 5000))
anova(rdasimp2_new, by = "axis", permutations = how(nperm = 5000))

## Tests of all terms
anova(rdasimp2, by = "terms", permutations = how(nperm = 5000))
anova(rdasimp2_new, by = "terms", permutations = how(nperm = 5000))

## Test marginal -- each predictor still considered sig when considered uniquely? 
anova(rdasimp2, by = "margin", permutations = how(nperm = 5000))
anova(rdasimp2_new, by = "margin", permutations = how(nperm = 5000))

# conclusion is to proceed with the new model. 

rdasimp2 <- rdasimp2_new 

# get the uniuqe effect of each predictor


# Unique effect of calcium yield
rda_cay_unique <- rda(
  Y ~ logcay +
    Condition(JDay + as.factor(slumpYN) + logwy + loggpp),
  data = d,
  scale = TRUE
)

# Unique effect of water yield
rda_wy_unique <- rda(
  Y ~ logwy +
    Condition(JDay + as.factor(slumpYN) + logcay + loggpp),
  data = d,
  scale = TRUE
)

# Unique effect of GPP
rda_gpp_unique <- rda(
  Y ~ loggpp +
    Condition(JDay + as.factor(slumpYN) + logcay + logwy),
  data = d,
  scale = TRUE
)

rda_slump_unique # 0.0978
rda_cay_unique # 0.237
rda_wy_unique  # 0.0717
rda_gpp_unique # 0.0626

# update the following - predictors are now:
# as.factor(slumpYN) +  logcay + logwy + loggpp


# Multiple versions here since only 4 predictors can be supported and 
# we have 5 

out2 = varpart(Y, ~JDay, ~as.factor(slumpYN) + 
                 logcay + 
                 logwy + 
                 loggpp,
               data = d, scale=TRUE)

slump = varpart(Y, ~JDay, ~as.factor(slumpYN), ~ logcay + 
                  logwy + 
                  loggpp,
                data = d, scale=TRUE)

cay = varpart(Y, ~JDay, ~logcay, ~as.factor(slumpYN) + 
                logwy + 
                loggpp,
              data = d, scale=TRUE)

wy = varpart(Y, ~JDay, ~logwy, ~as.factor(slumpYN) +
                   logcay + 
                   loggpp,
                 data = d, scale=TRUE)

gpp= varpart(Y, ~JDay, ~ loggpp, ~as.factor(slumpYN) + 
               logcay +
               logwy,
             data = d, scale=TRUE)

out2
slump
cay
gpp
wy


# JDay always explains barely anything. Try without

no_day = varpart(Y, ~as.factor(slumpYN),
                 ~logcay, 
                  ~logwy, 
                  ~loggpp,
                        data = d, scale=TRUE)

no_day
# probably the most represnetative 


#mor= varpart(Y, ~JDay, ~morainal_log, ~as.factor(slumpYN) + 
#               logcay + 
#               loggpp +
#               logspower,
#             data = d, scale=TRUE)


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
                     "log(WaterYield+2)",
                     "log(GPP+2)")

arrows <- arrows[arrows$pf_names!="factor",] # want to plot this as a centroid not an arrow

factorcent <- data.frame(scor$centroids)
factorcent$names <- c("NO RTS", "RTS")

mult <- attributes(scores(rdasimp2))$const

theme<-theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             axis.line.x = element_line(colour="black"),
             axis.line.y = element_line(colour="black"),
             axis.text = element_text(colour="black",size=14),
             legend.background=element_blank(),
             text=element_text(size = 16),
             legend.position="none",
             #legend.title=element_blank(),
             # legend.position = c(0,1),
             #legend.direction="horizontal",
             #legend.justification = c(0,1),
             aspect.ratio=1)

rdagraph2 <- ggplot(data = species_centroids, 
                   aes(x = RDA1, y= RDA2)) +
  geom_vline(xintercept = 0,linetype="dashed", 
             colour="grey")+ylab("RDA2 (11.0%)")+
  
  geom_hline(yintercept = 0,linetype="dashed", 
             colour="grey")+xlab("RDA1 (53.8%)")+
  
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
               arrow = arrow(length = unit(0.4, "cm")), colour = "grey130")+
  geom_text(data = arrows,
            aes(x= 1.2*RDA1, y = 1.2*RDA2, #we add 10% to the text to push it slightly out from arrows
                label = pf_names), #otherwise you could use hjust and vjust. I prefer this option
            size = 4,
            hjust = 0.5)+
  theme



rdagraph2

##### ============================== Section 4: Export plot ===================================================
library(grid)


savePlotpdf <- function(myPlot) {
  pdf(file = paste0("Figures/", filename, ".pdf"),
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
