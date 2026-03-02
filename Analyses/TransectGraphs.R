# ===========================================================================================================#
# TransectGraphs.R
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

# load functions
se <- function(x) {
  sqrt(var(x, na.rm=TRUE)/length(x[!is.na(x)])) 
}

## create standard error function for use later
lengthnona <- function(x) {
  length(x[!is.na(x)])
}

# load functions
savePlotpdf <- function(myPlot) {
  pdf(file = paste0("Figures/", name),
      width = w, height = h)
  print(myPlot)
  dev.off()
}

# colourblind palette

# The palette with grey:
cbPalette<-c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2","#D55E00","#CC79A7")

# The palette with black:
cbbPalette<-c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

##### ========== (1) DATA PREP =========================================================================

# (1.2) 2015 Transect

#(1.2.1) POC, Q
trans2015 <- read_excel(paste0(df, "2015Transects.xlsx"))
trans2015$POCflux <- trans2015$POCmgL*trans2015$Qm3s*1000

cols.num <- c("POCmgL", "Qm3s", "POCflux", "Temp", "Cond")
trans2015[cols.num] <- sapply(trans2015[cols.num],as.numeric)
sapply(trans2015, class)

#(1.2.2) 13C
c13 <- read_excel(paste0(df, "2015_POC_Amount_13C.xlsx"))
c13 <- c13[c13$`Multiple Analysis Number`==1, 
           c("Slump Site", "Stream Location", "Transect Location",
             "Sampling Date", "Sample Type", "d13C")]

#(1.2.3) TSS for percPOC
pc <- read_excel(paste0(df, "2015_POC_TSS.xlsx"))
pc$TSSmgL <- (pc$`Filter Post-weight (mg)`-pc$`Filter Pre-weight (mg)`)/pc$`Filtered Volume (L)`
pc <- pc[,c("Slump Site", "Stream Location", "Transect Location",
            "Sampling Date", "Sample Type", "TSSmgL")]
#(1.2.4) merge
a <- merge(trans2015, c13, all.x=TRUE)
trans <- merge(a, pc, all.x=TRUE)
trans <- trans[trans$`Sample Type`=="Sample",]
trans$`Sample Type` <- NULL


#(1.2.5) New calcs (percPOC)
trans$percPOC <- trans$POCmgL/trans$TSSmgL
#trans$POCmgL <- trans$POCmgL/(trans$Qm3s*1000)

translong <- trans %>%
             pivot_longer(cols=c("percPOC", "d13C", "POCmgL", "POCflux"),
                          names_to="type", values_to = "value")

translong$`Sampling Date` <- as.character(translong$`Sampling Date`)
translong$`Stream Location` <- as.factor(translong$`Stream Location`)

#translongm <- translong %>%
 # group_by(`Slump Site`, `Stream Location`, 
  #         `Transect Location`, `distm`, `type`) %>%
  #summarize(n=lengthnona(value), 
   #         means=mean(value, na.rm=TRUE),
    #        sd=sd(value, na.rm=TRUE),
     #       sem=se(value))

#translongm$type <- as.factor(as.character(translongm$type))

# (1.3) 2016 Transect
data.all <- read_excel(paste0(df, "20162017POCPO14CTrans.xlsx"))
data.all$distkm <- data.all$distm/1000
data.samples <- data.all[data.all$`Sample Type`=="Sample",]

##### ========== (2) 2015 Transect Graph =========================================================================

## command to remove scientific notation from ggplot
options(scipen=999)

plain <- function(x) {
  signif(x, 1)
}

# save theme
theme <- theme(panel.background = element_rect(fill="white"), 
               panel.border = element_rect(colour="black", fill=NA, size=1),
               plot.margin = unit( c(0.05,0.05,0.05,0.05) , "in"),# fill must equal NA for pane.border or it will block out graph
               legend.key = element_rect(fill="transparent", colour="transparent"),
               axis.text=element_text(size=11, colour="black"),
               axis.title=element_text(size=11, colour="black"),
               axis.text.x = element_text(angle = 0),
               legend.position=c(0,1), 
               legend.justification=c(0, 0),
               legend.direction="horizontal",
               legend.background=element_blank(),
               legend.text=element_text(size=11),
               legend.spacing.y = unit(0.1, 'cm'),
               aspect.ratio=0.7,
               legend.margin = margin(0,0,0,0, unit="cm"),
               strip.background = element_blank(),
               strip.placement = "outside") 

# set facet order for graphing
translong$type <- factor(translong$type, 
                         levels=c("percPOC", "d13C", "POCmgL", "POCflux"))

# 2.1 SE 2015 Transect =====

SE <- ggplot() + 
  #geom_line(data=translong[translong$`Slump Site`=="SE",] ,
   #         aes(x=distm, y=value,
    #            colour=as.character(`Sampling Date`)),
     #       linetype=2) +
  geom_point(data = translong[translong$`Slump Site`=="SE",], 
             aes(x=distm, y=value, 
                 fill=as.character(`Sampling Date`), 
                 shape=`Stream Location`), 
             size=4) +
  scale_shape_manual(breaks=c("UP", "IN", "DN"),
                     limits=c("UP", "IN", "DN"),
                     values=c(24, 21, 25),
                     labels=c("UP", "IN", "DN"),
                     guide=guide_legend(title=NULL)) +
  scale_fill_manual(breaks=c("2015-07-06", "2015-07-21", "2015-08-21"),
                     limits=c("2015-07-06", "2015-07-21", "2015-08-21"),
                     values=c("black", "grey", "white"),
                     labels=c("JUL-6", "JUL-21", "AUG-21"),
                     guide=guide_legend(title=NULL)) +
  scale_colour_manual(breaks=c("2015-07-06", "2015-07-21", "2015-08-21"),
                    limits=c("2015-07-06", "2015-07-21", "2015-08-21"),
                    values=c("black", "grey30", "grey"),
                    labels=c("JUL-6", "JUL-21", "AUG-21"),
                    guide=guide_legend(title=NULL)) +
  ylab (NULL) + xlab("Distance (m)") + 
  guides(color = guide_legend(override.aes = list(shape = 21))) +
  facet_wrap(~type, scales="free_y", ncol=1,
             strip.position = "left", 
             labeller = as_labeller(c(percPOC = "%POC",
                                      d13C = expression(delta*"13C"),
                                      POCmgL = expression("POC conc (mg L"^-1*")"),
                                      POCflux = expression("POC flux (mg s"^-1*")"))
             ) ) + theme 

SE2 <- ggplot() + 
  geom_point(data = translong[translong$`Slump Site`=="SE"&
                              translong$type=="POCmgL",], 
             aes(x=distm, y=value, 
                 fill=as.character(`Sampling Date`), 
                 shape=`Stream Location`), 
                 size=4) +
  scale_shape_manual(breaks=c("UP", "IN", "DN"),
                     limits=c("UP", "IN", "DN"),
                     values=c(24, 21, 25),
                     labels=c("UP", "IN", "DN"),
                     guide=guide_legend(title=NULL)) +
  scale_fill_manual(breaks=c("2015-07-06", "2015-07-21", "2015-08-21"),
                    limits=c("2015-07-06", "2015-07-21", "2015-08-21"),
                    values=c("black", "grey", "white"),
                    labels=c("JUL-6", "JUL-21", "AUG-21"),
                    guide=guide_legend(title=NULL)) +
  scale_colour_manual(breaks=c("2015-07-06", "2015-07-21", "2015-08-21"),
                      limits=c("2015-07-06", "2015-07-21", "2015-08-21"),
                      values=c("black", "grey30", "grey"),
                      labels=c("JUL-6", "JUL-21", "AUG-21"),
                      guide=guide_legend(title=NULL)) +
  ylab (NULL) + xlab(NULL) + scale_y_log10() +
  theme + 
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# 2.1 SD 2015 Transect =====

# SD
transSDlong <- translong[translong$`Slump Site`=="SD",] 

SD <- ggplot() + 
  #geom_line(data=transSElong,
  #         aes(x=distm, y=value,
  #            colour=as.character(`Sampling Date`)),
  #       linetype=2) +
  geom_point(data = translong[translong$`Slump Site`=="SD",], 
             aes(x=distm, y=value, 
                 fill=as.character(`Sampling Date`), 
                 shape=`Stream Location`), 
             size=4) +
  scale_shape_manual(breaks=c("UP", "IN", "DN"),
                     limits=c("UP", "IN", "DN"),
                     values=c(24, 21, 25),
                     labels=c("UP", "IN", "DN"),
                     guide=guide_legend(title=NULL)) +
  scale_fill_manual(breaks=c("2015-07-04", "2015-07-22", "2015-08-16"),
                    limits=c("2015-07-04", "2015-07-22", "2015-08-16"),
                    values=c("black", "grey", "white"),
                    labels=c("JUL-4", "JUL-22", "AUG-16"),
                    guide=guide_legend(title=NULL)) +
  scale_colour_manual(breaks=c("2015-07-04", "2015-07-22", "2015-08-16"),
                      limits=c("2015-07-04", "2015-07-22", "2015-08-16"),
                      values=c("black", "grey30", "grey"),
                      labels=c("JUL-4", "JUL-22", "AUG-16"),
                      guide=guide_legend(title=NULL)) +
  ylab (NULL) + xlab("Distance (m)") + 
  guides(color = guide_legend(override.aes = list(shape = 21))) +
  facet_wrap(~type, scales="free_y", ncol=1,
             strip.position = "left", 
             labeller = as_labeller(c(percPOC = "%POC",
                                      d13C = expression(delta*"13C"),
                                      POCmgL = expression("POC conc (mg L"^-1*")"),
                                      POCflux = expression("POC flux (mg s"^-1*")"))
             ) )  + theme

SD2 <- ggplot() + 
  geom_point(data = translong[translong$`Slump Site`=="SD"&
                                translong$type=="POCmgL",], 
             aes(x=distm, y=value, 
                 fill=as.character(`Sampling Date`), 
                 shape=`Stream Location`), 
             size=4) +
  scale_shape_manual(breaks=c("UP", "IN", "DN"),
                     limits=c("UP", "IN", "DN"),
                     values=c(24, 21, 25),
                     labels=c("UP", "IN", "DN"),
                     guide=guide_legend(title=NULL)) +
  scale_fill_manual(breaks=c("2015-07-04", "2015-07-22", "2015-08-16"),
                    limits=c("2015-07-04", "2015-07-22", "2015-08-16"),
                    values=c("black", "grey", "white"),
                    labels=c("JUL-4", "JUL-22", "AUG-16"),
                    guide=guide_legend(title=NULL)) +
  scale_colour_manual(breaks=c("2015-07-04", "2015-07-22", "2015-08-16"),
                      limits=c("2015-07-04", "2015-07-22", "2015-08-16"),
                      values=c("black", "grey30", "grey"),
                      labels=c("JUL-4", "JUL-22", "AUG-16"),
                      guide=guide_legend(title=NULL)) +
  ylab (NULL) + xlab(NULL) + scale_y_log10() +
  theme + 
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

##### ========== (3) DATA PREP =========================================================================
### (3.1) load in graphing libraries =====
library(gtable)
library(grid)
library(gridExtra)


### (3.2) Graphs and dimensions =====
SE <- ggplotGrob(SE)
SD <- ggplotGrob(SD)

### (3.3) Print panel graphs =====
w = 4; h=11

name <- "2015SEtransect.pdf"
savePlotpdf(grid.draw(SE))

name <- "2015SDtransect.pdf"
savePlotpdf(grid.draw(SD))


### (3.4) Print concentration graphs =====
w = 3.7; h=3.8

name <- "2015SEtransectconc.pdf"
savePlotpdf(grid.draw(SE2))

w = 3.6; h=3.75
name <- "2015SDtransectconc.pdf"
savePlotpdf(grid.draw(SD2))