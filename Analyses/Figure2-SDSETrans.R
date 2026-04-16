# ===========================================================================================================#
# TransectGraphs.R
# By: Sarah Shakil
# Background: Publication-ready transect graphs
# ===========================================================================================================#

##### (i) Workspace PREP ==============================================================================

## Clear list 
rm(list=ls(all=TRUE))

## Set working directory
# setwd("YOUR_DIRECTORY_HERE")
df <- "Data/"

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

# colourblind palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2","#D55E00","#CC79A7")
cbbPalette <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

##### ========== (1) DATA PREP =========================================================================

# (1.1) Watershed Areas
wa <- read_excel(paste0(df, "watershedareas.xlsx"), sheet="watershedareasforR")

# (1.2) 2015 Transect
#(1.2.1) POC, Q
trans2015 <- read_excel(paste0(df, "2015Transects.xlsx"))
trans2015$POCflux <- trans2015$POCmgL * trans2015$Qm3s * 1000

trans2015 <- merge(trans2015, wa, by.y=c("site", "loc", "trans", "date"),
                   by.x=c("Slump Site", "Stream Location", "Transect Location", "Sampling Date"),
                   all.x=TRUE) # FIXED TYPO: changed a..x=TRUE to all.x=TRUE

trans2015$POCyieldmgLkm2 <- (trans2015$POCmgL * trans2015$Qm3s * 1000) / trans2015$WatershedArea

cols.num <- c("POCmgL", "Qm3s", "POCyieldmgLkm2", "Temp", "Cond")
trans2015[cols.num] <- sapply(trans2015[cols.num], as.numeric)

#(1.2.2) 13C
c13 <- read_excel(paste0(df, "2015_POC_Amount_13C.xlsx"))
c13 <- c13[c13$`Multiple Analysis Number`==1, 
           c("Slump Site", "Stream Location", "Transect Location",
             "Sampling Date", "Sample Type", "d13C")]

#(1.2.3) TSS for percPOC
pc <- read_excel(paste0(df, "2015_POC_TSS.xlsx"))
pc$TSSmgL <- (pc$`Filter Post-weight (mg)` - pc$`Filter Pre-weight (mg)`) / pc$`Filtered Volume (L)`
pc <- pc[,c("Slump Site", "Stream Location", "Transect Location",
            "Sampling Date", "Sample Type", "TSSmgL")]

#(1.2.4) merge
a <- merge(trans2015, c13, all.x=TRUE)
trans <- merge(a, pc, all.x=TRUE)
trans <- trans[trans$`Sample Type`=="Sample",]
trans$`Sample Type` <- NULL


#(1.2.5) Calculate %POC
trans$percPOC <- trans$POCmgL / trans$TSSmgL

# (1.2.6) Pivot longer
translong <- trans %>%
  pivot_longer(cols=c("percPOC", "d13C", "POCmgL", "POCflux", "POCyieldmgLkm2"),
               names_to="type", values_to = "value")

translong$`Sampling Date` <- as.character(translong$`Sampling Date`)
translong$`Stream Location` <- as.factor(translong$`Stream Location`)

##### ========== (2) 2015 Transect Graphing =================================================

options(scipen=999)
library(patchwork) 

# 2.1 Universal Theme
pub_theme <- theme_bw(base_size = 12, base_family = "sans") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# 2.2 Fix the Legend Dates
# This restores the missing Sampling_Period column!
translong <- translong %>%
  mutate(Sampling_Period = case_when(
    `Sampling Date` %in% c("2015-07-04", "2015-07-06") ~ "JUL-4 / JUL-6",
    `Sampling Date` %in% c("2015-07-21", "2015-07-22") ~ "JUL-21 / JUL-22",
    `Sampling Date` %in% c("2015-08-16", "2015-08-21") ~ "AUG-16 / AUG-21"
  )) %>%
  mutate(Sampling_Period = factor(Sampling_Period, 
                                  levels = c("JUL-4 / JUL-6", "JUL-21 / JUL-22", "AUG-16 / AUG-21")))

# Also update your factor levels right above your graphing section:
translong$type <- factor(translong$type, 
                         levels=c("percPOC", "d13C", "POCmgL", "POCflux", "POCyieldmgLkm2"))

# 2.3 Create Shared Plot Layers 
base_layers <- list(
  geom_point(size = 3.5, stroke = 0.8),
  scale_shape_manual(limits = c("UP", "IN", "DN"), values = c(24, 21, 25), name = "Location:"),
  scale_fill_manual(limits = c("JUL-4 / JUL-6", "JUL-21 / JUL-22", "AUG-16 / AUG-21"),
                    values = c("black", "grey50", "white"), name = "Date:",
                    guide = guide_legend(override.aes = list(shape = 21))),
  pub_theme,
  coord_cartesian(clip = "off"), # CRITICAL: Allows text to be drawn outside the plot boundaries
  theme(plot.margin = margin(t = 20, r = 5, b = 5, l = 5)) # Adds 20pt of space at the top for the labels
)

# Helper for hiding X-axis on the top 3 rows
no_x_axis <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

# --- ROW 1: % POC ---
p1_SD <- ggplot(translong %>% filter(type == "percPOC", `Slump Site` == "SD"),
                aes(x = distm, y = value, fill = Sampling_Period, shape = `Stream Location`)) +
  base_layers + no_x_axis + ggtitle("SD") + ylab("% POC") +
  annotate("text", x = -Inf, y = Inf, label = "(a)", hjust = 0, vjust = -0.8, size = 5)

p1_SE <- ggplot(translong %>% filter(type == "percPOC", `Slump Site` == "SE"),
                aes(x = distm, y = value, fill = Sampling_Period, shape = `Stream Location`)) +
  base_layers + no_x_axis + ggtitle("SE") + ylab(NULL) +
  annotate("text", x = -Inf, y = Inf, label = "(b)", hjust = 0, vjust = -0.8, size = 5)

# --- ROW 2: d13C ---
p2_SD <- ggplot(translong %>% filter(type == "d13C", `Slump Site` == "SD"),
                aes(x = distm, y = value, fill = Sampling_Period, shape = `Stream Location`)) +
  base_layers + no_x_axis + ylab(expression(delta^13*"C-POC")) +
  annotate("text", x = -Inf, y = Inf, label = "(c)", hjust = 0, vjust = -0.8, size = 5)

p2_SE <- ggplot(translong %>% filter(type == "d13C", `Slump Site` == "SE"),
                aes(x = distm, y = value, fill = Sampling_Period, shape = `Stream Location`)) +
  base_layers + no_x_axis + ylab(NULL) +
  annotate("text", x = -Inf, y = Inf, label = "(d)", hjust = 0, vjust = -0.8, size = 5)

# --- ROW 3: POC mg L-1 (WITH LOG SCALE) ---
p3_SD <- ggplot(translong %>% filter(type == "POCmgL", `Slump Site` == "SD"),
                aes(x = distm, y = value, fill = Sampling_Period, shape = `Stream Location`)) +
  base_layers + no_x_axis + scale_y_log10() + ylab(expression("POC (mg L"^-1*")")) +
  annotate("text", x = -Inf, y = Inf, label = "(e)", hjust = 0, vjust = -0.8, size = 5)

p3_SE <- ggplot(translong %>% filter(type == "POCmgL", `Slump Site` == "SE"),
                aes(x = distm, y = value, fill = Sampling_Period, shape = `Stream Location`)) +
  base_layers + no_x_axis + scale_y_log10() + ylab(NULL) +
  annotate("text", x = -Inf, y = Inf, label = "(f)", hjust = 0, vjust = -0.8, size = 5)

# --- ROW 4: POC flux ---
p4_SD <- ggplot(translong %>% filter(type == "POCflux", `Slump Site` == "SD"),
                aes(x = distm, y = value, fill = Sampling_Period, shape = `Stream Location`)) +
  base_layers + no_x_axis + ylab(expression("POC flux mg s"^-1)) +
  annotate("text", x = -Inf, y = Inf, label = "(g)", hjust = 0, vjust = -0.8, size = 5)

p4_SE <- ggplot(translong %>% filter(type == "POCflux", `Slump Site` == "SE"),
                aes(x = distm, y = value, fill = Sampling_Period, shape = `Stream Location`)) +
  base_layers + no_x_axis + ylab(NULL) +
  annotate("text", x = -Inf, y = Inf, label = "(h)", hjust = 0, vjust = -0.8, size = 5)

# --- ROW 5: POC yield ---
p5_SD <- ggplot(translong %>% filter(type == "POCyieldmgLkm2", `Slump Site` == "SD"),
                aes(x = distm, y = value, fill = Sampling_Period, shape = `Stream Location`)) +
  base_layers + xlab("Distance (m)") + ylab(expression("POC yield mg s"^-1*"km"^-2)) +
  annotate("text", x = -Inf, y = Inf, label = "(i)", hjust = 0, vjust = -0.8, size = 5)

p5_SE <- ggplot(translong %>% filter(type == "POCyieldmgLkm2", `Slump Site` == "SE"),
                aes(x = distm, y = value, fill = Sampling_Period, shape = `Stream Location`)) +
  base_layers + xlab("Distance (m)") + ylab(NULL) +
  annotate("text", x = -Inf, y = Inf, label = "(j)", hjust = 0, vjust = -0.8, size = 5)

##### ========== (3) COMPILE AND SAVE ====================================================================

if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Figures/Fig2")) dir.create("Figures/Fig2", recursive = TRUE)

# Stitch the 10 plots together
combined_plot <- (p1_SD | p1_SE) / 
  (p2_SD | p2_SE) / 
  (p3_SD | p3_SE) / 
  (p4_SD | p4_SE) / 
  (p5_SD | p5_SE) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top",
        legend.box = "horizontal",
        legend.margin = margin(b = 10)) 

# Save PDF (Increased height from 11 to 13.5 to accommodate the extra row)
ggsave(filename = "Figures/Fig2/2015_Transect_Combined.pdf", 
       plot = combined_plot, 
       width = 8.5, height = 13.5, units = "in")

# Save JPEG
ggsave(filename = "Figures/Fig2/2015_Transect_Combined.jpg", 
       plot = combined_plot, 
       width = 8.5, height = 13.5, units = "in",
       dpi = 300, bg = "white")