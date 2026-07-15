# ===========================================================================================================#
# Endmembertable.R
# By: Sarah Shakil
# Contact: shakil@ualberta.ca
# Background: Preliminary stats
# Last Updated: October 20 2020
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

## Call book code with necessary functions 
# (e.g. pairs function; Zuur 2009, Mixed Effects Models and Extensions)
source("functions/HighstatLibV10.R")

##### ========== (1) DATA PREP ==========================================================================

# (1.1) 2017 streambank and merge ====================

### 13C =====
sed13c <- read_excel(paste0(df, "Shakil_PeelHWSBPO13C.xlsx"))

# use this one (?) -MT
sed13c <- read_excel(paste0(df, "Sediment_PerPOC.xlsx"))

sed13c$cat <- "SB"
sed13c$cat[sed13c$`Transect Location`=="PLE"] <- "PLE"
sed13c$cat[sed13c$`Transect Location`=="HOL"] <- "HOL"
sed13c$cat[sed13c$`Transect Location`=="OHO"] <- "UAL"
sed13c$cat[sed13c$`Transect Location`=="AHO"|
           sed13c$`Transect Location`=="BHO"] <- "LAL"

sed13c <- sed13c %>%
  filter((`Stream Location`=="NA"|
            `Stream Location`=="erosion"|
            `Stream Location`=="streambank" |
            `Transect Location`=="PLE"|
            `Transect Location`=="HOL"|
            `Transect Location`=="OHO"|
            `Transect Location`=="AHO"|
            `Transect Location`=="BHO")&
           `Sampling Date` > 2017 & 
           `Sample Type`=="Sample") %>%
  select(cat, `Slump Site`, `Stream Location`, 
         `Transect Location`, `Sampling Date`,
         `Sample Type`, d13C)

### %POC ====

sedperpoc <- read_excel(paste0(df, "Sediment_PerPOC.xlsx"))

sedperpoc$cat <- "SB"
sedperpoc$cat[sedperpoc$`Transect Location`=="PLE"] <- "PLE"
sedperpoc$cat[sedperpoc$`Transect Location`=="HOL"] <- "HOL"
sedperpoc$cat[sedperpoc$`Transect Location`=="OHO"] <- "UAL"
sedperpoc$cat[sedperpoc$`Transect Location`=="AHO"|
              sedperpoc$`Transect Location`=="BHO"] <- "LAL"

sedperpoc <- sedperpoc %>%
  filter((`Stream Location`=="NA"|
          `Stream Location`=="erosion"|
          `Stream Location`=="streambank" |
          `Transect Location`=="PLE"|
          `Transect Location`=="HOL"|
          `Transect Location`=="OHO"|
          `Transect Location`=="AHO"|
          `Transect Location`=="BHO")&
          `Sampling Date` > 2017 & 
          `Sample Type`=="Sample") %>%
  select(cat, `Slump Site`, `Stream Location`, 
         `Transect Location`, `Sampling Date`,
         `Sample Type`, PercPOC)

# (1.2) F14C of sed ====================
sed14c <- read_excel(paste0(df, "20162017POCPO14CTrans.xlsx"))

sed14c <- sed14c %>%
          filter(MatCode=="A") %>%
          select(`Slump Site`, `Stream Location`,
                 `Transect Location`, `Sampling Date`, 
                 F14C_A = F14C, F14Cerror_A = `F14C error`)
sed14c$cat <- "SB"
          
# (1.3) Ch1 headwall and peri samples ====================
c1 <- read.csv(paste0(df, "PeelPlateau_RTSandstream_geochem.csv"))

## there is no data on the periphyton samples for the table


c1$PO13C[c1$endmembertype=="Periphyton"]
c1$cat <- as.character(c1$headwallcat)
c1$cat[c1$endmembertype=="Periphyton"] <- "PERI"
c1 <- c1 %>%
      filter(sampletype=="Endmember" &
             !is.na(F14C_A | F14C_DIC)) %>%
      select(cat, slumpsite, 
             streamlocation, samplingdate,
             F14C_A,
             F14Cerror_A, F14C_DIC)

c1$slumpsite <- as.character(c1$slumpsite)
c1$slumpsite[c1$slumpsite=="FM2"] <- "SB"
c1$slumpsite[c1$slumpsite=="FM3"] <- "SC"

c1$streamlocation <- "HW"

# (1.4) Merge files ====================

sed17 <- merge(sed13c, sedperpoc, all=TRUE)
sed <- merge(sed17, sed14c, all=TRUE)

sed$`Sampling Date` <- as.Date(sed$`Sampling Date`)
c1$samplingdate <- as.Date(c1$samplingdate)

em <- merge(c1, sed, 
            by.x=c("cat", "slumpsite", 
                   "streamlocation", "samplingdate",
                   "F14C_A", "F14Cerror_A"), 
            by.y=c("cat", "Slump Site", 
                   "Stream Location", "Sampling Date",
                   "F14C_A", "F14Cerror_A"),
            all=TRUE)

write.csv(em,paste0(df, "Figures/endmembersummaryall.csv") )

##### ========== (2) Calculate summaries ==========================================================================

# (2.1) Avg, Sd, Range ====================

emsum <- em %>%
  group_by(cat) %>%
  summarize(
    perpocav = mean(PercPOC, na.rm = TRUE),
    perpocsd = sd(PercPOC, na.rm = TRUE),
    perpocse = se(PercPOC),
    perpocmin = min(PercPOC, na.rm = TRUE),
    perpocmax = max(PercPOC, na.rm = TRUE),
    
    d13cav = mean(d13C, na.rm = TRUE),
    d13csd = sd(d13C, na.rm = TRUE),
    d13cse = se(d13C),
    d13cmin = min(d13C, na.rm = TRUE),
    d13cmax = max(d13C, na.rm = TRUE),
    
    F14Cav = mean(F14C_A, na.rm = TRUE),
    F14Csd = sd(F14C_A, na.rm = TRUE),
    F14Cse = se(F14C_A),
    F14Cmin = min(F14C_A, na.rm = TRUE),
    F14Cmax = max(F14C_A, na.rm = TRUE),
    .groups = "drop"
  )

# check structure 
# there is no periphyton category

str(emsum)

# Add the periphyton data manually, based on the table in the manuscript
## see the footnotes for Table 2 -MT 

emsum <- emsum %>%
  add_row(
    cat       = "PERI",
    perpocav  = 13.75,
    perpocsd  = NA_real_,
    perpocse  = NA_real_,
    perpocmin = NA_real_,
    perpocmax = NA_real_,
    d13cav    = -33.39,
    d13csd    = NA_real_,
    d13cse    = 2.10,
    d13cmin   = -38.83,
    d13cmax   = -30.02,
    F14Cav    = 0.80,
    F14Csd    = NA_real_,
    F14Cse    = NA_real_,
    F14Cmin   = 0.673,
    F14Cmax   = NA_real_
  )

# convert to long 
emsum_long <- emsum %>%
  select(
    cat,
    perpocav, perpocse, perpocmin, perpocmax,
    d13cav, d13cse, d13cmin, d13cmax,
    F14Cav, F14Cse, F14Cmin, F14Cmax
  ) %>%
  pivot_longer(
    cols = -cat,
    names_to = c("variable", ".value"),
    names_pattern = "(perpoc|d13c|F14C)(av|se|min|max)"
  ) %>%
  mutate(
    variable = recode(
      variable,
      perpoc = "SOC",
      d13c = "d13C",
      F14C = "F14C_A"
    ),
    variable = factor(
      variable,
      levels = c("SOC", "d13C", "F14C_A")
    )
  )

# check structure 
str(emsum_long)


# Plot: average and se  ---------------------------------------------------

# 13C x %POC plot

# reshape the summary data so each category has separate SOC and d13C columns
emsum_plot <- emsum_long %>%
  filter(variable %in% c("SOC", "d13C")) %>%
  select(cat, variable, av, se) %>%
  pivot_wider(
    names_from = variable,
    values_from = c(av, se),
    names_glue = "{variable}_{.value}"
  )

# Plot category avg with standard-error bars on both axes
c13_poc <- ggplot(
  emsum_plot,
  aes(
    x = d13C_av,
    y = SOC_av,
    color = cat
  )
) +
  geom_errorbar(
    aes(
      ymin = SOC_av - SOC_se,
      ymax = SOC_av + SOC_se
    ),
    width = 0
  ) +
  geom_errorbar(
    aes(
      xmin = d13C_av - d13C_se,
      xmax = d13C_av + d13C_se
    ),
    height = 0
  ) +
  geom_point(size = 3) +
  labs(
    x = expression(delta^13 * C),
    y = "% POC",
    color = "Category"
  ) +
  theme_classic()

# check 
c13_poc

# 13C x 14C plot 

# Reshape d13C and F14C values into separate columns

emsum_plot_14C <- emsum_long %>%
  filter(variable %in% c("d13C", "F14C_A")) %>%
  select(cat, variable, av, se) %>%
  pivot_wider(
    names_from = variable,
    values_from = c(av, se),
    names_glue = "{variable}_{.value}"
  )

# Plot mean d13C against mean F14C with SE on both axes
c13_c14 <-ggplot(
  emsum_plot_14C,
  aes(
    x = d13C_av,
    y = F14C_A_av,
    color = cat
  )
) +
  geom_errorbar(
    aes(
      ymin = F14C_A_av - F14C_A_se,
      ymax = F14C_A_av + F14C_A_se
    ),
    width = 0,
    na.rm = TRUE
  ) +
  geom_errorbarh(
    aes(
      xmin = d13C_av - d13C_se,
      xmax = d13C_av + d13C_se
    ),
    height = 0,
    na.rm = TRUE
  ) +
  geom_point(size = 3) +
  labs(
    x = expression(delta^13 * C~("\u2030")),
    y = expression(F^14 * C),
    color = "Category"
  ) +
  theme_classic()

c13_poc
c13_c14

library(patchwork)
# plot them together 

c13_poc + c13_c14

# remove legend
combined_plot <- c13_poc +
  theme(legend.position = "none") +
  c13_c14 +
  plot_layout(ncol = 2)

combined_plot



# Plot: all measurements --------------------------------------------------

## Try the plot with all measurements instead of average 


# d13C x %POC plot
c13_poc_all <- ggplot(
  em,
  aes(
    x = d13C,
    y = PercPOC,
    color = cat
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  labs(
    x = expression(delta^13 * C~("\u2030")),
    y = "% POC",
    color = "Category"
  ) +
  theme_classic()

# d13C x F14C plot
c13_c14_all <- ggplot(
  em,
  aes(
    x = d13C,
    y = F14C_A,
    color = cat
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8,
    position = position_jitter(width = 0.08, height = 0),
    na.rm = TRUE
  ) +
  labs(
    x = expression(delta^13 * C~("\u2030")),
    y = expression(F^14 * C),
    color = "Category"
  ) +
  theme_classic()

# Check individual plots
c13_poc_all
c13_c14_all

# there are several measurements of soc, 13c, 14c of each member
# but there are very few samples with paired measurements 

# there are only 8 samples with 14C measurements
# of which, only 2  of them have paired 13C measurements.. 
# streambank samples.. 

# Combine plots with one shared legend
combined_plot_all <- c13_poc_all +
  c13_c14_all +
  plot_layout(ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

combined_plot_all



