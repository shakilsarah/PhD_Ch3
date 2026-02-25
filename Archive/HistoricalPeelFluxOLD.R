# ===========================================================================================================#
# HistoricalPeelFlux.R
# By: Sarah Shakil
# Contact: shakil@ualberta.ca
# Background: Peel historical flux loadest outputs
# Last Updated: October 15 2021
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
library(lubridate)

# load functions
se <- function(x) {
  sqrt(var(x, na.rm=TRUE)/length(x[!is.na(x)])) 
}

## create standard error function for use later
lengthnona <- function(x) {
  length(x[!is.na(x)])
}


##### ========== (1) Check gaps in input Q data ==========================================================================

dfq <- "ThesisDrafts/Chapter3/peeldataset/Shakil Loadrunner/20210930_09.51.32.399_DOC/inputs/"

qin <- read.delim(paste0(wd, dfq, 
                             "qloadest.txt"), 
                      header = TRUE, sep = "\t", skip=5,
                      row.names=NULL)

qin$datetime <- as.Date(qin$datetime)
qgaps <- qin %>% 
  arrange(datetime) %>% 
  mutate(diff_row = datetime - lag(datetime)) %>%
  filter(diff_row>7)

qgaps$year <- year(qgaps$datetime)

#DateRange <- seq(min(qin$datetime), max(qin$datetime), by = 1)
#gaps <- as.data.frame(DateRange[!DateRange %in% qin$datetime])
##### ========== (1) Read in data ==========================================================================

df <- "ThesisDrafts/Chapter3/peeldataset/Shakil Loadrunner/"

### ===== (1.1) DOC =====

df2 <- "20210930_09.51.32.399_DOC/outputs/"
variable <- "DOC"

## (1.1.1) Monthly =====

docmnth <- read.delim(paste0(wd, df, df2, 
                             "allsite_", variable, "_monthflux.txt"), 
                      header = TRUE, sep = "\t", skip=15,
                      row.names=NULL)

data <- docmnth
names(data)[1:(ncol(data)-1)] <- names(data)[2:ncol(data)]
data[, ncol(data)] <- NULL
docmnth <- data

## (1.1.2) Annual =====

docann <- read.delim(paste0(wd, df, df2, 
                             "allsite_", variable, "_annualflux.txt"), 
                      header = TRUE, sep = "\t", skip=15,
                      row.names=NULL)

data <- docann
names(data)[1:(ncol(data)-1)] <- names(data)[2:ncol(data)]
data[, ncol(data)] <- NULL
docann <- data

## (1.1.3) Flux =====

docf <- read.delim(paste0(wd, df, df2, 
                            "allsite_", variable, "_flux.txt"), 
                     header = TRUE, sep = "\t", skip=20,
                     row.names=NULL)

### ===== (1.2) POC =====

df2 <- "20210930_09.53.28.453_POC/outputs/"
variable <- "POC"

## (1.2.1) Monthly =====

pocmnth <- read.delim(paste0(wd, df, df2, 
                             "allsite_", variable, "_monthflux.txt"), 
                      header = TRUE, sep = "\t", skip=15,
                      row.names=NULL)

data <- pocmnth
names(data)[1:(ncol(data)-1)] <- names(data)[2:ncol(data)]
data[, ncol(data)] <- NULL
pocmnth <- data

## (1.2.2) Annual =====

pocann <- read.delim(paste0(wd, df, df2, 
                            "allsite_", variable, "_annualflux.txt"), 
                     header = TRUE, sep = "\t", skip=15,
                     row.names=NULL)

data <- pocann
names(data)[1:(ncol(data)-1)] <- names(data)[2:ncol(data)]
data[, ncol(data)] <- NULL
pocann <- data

## (1.2.3) Flux =====

pocf <- read.delim(paste0(wd, df, df2, 
                          "allsite_", variable, "_flux.txt"), 
                   header = TRUE, sep = "\t", skip=20,
                   row.names=NULL)

### ===== (1.3) TSS =====

df2 <- "20210930_09.56.42.601_TSS/outputs/"
variable <- "TSS"

## (1.3.1) Monthly =====

tssmnth <- read.delim(paste0(wd, df, df2, 
                             "allsite_", variable, "_monthflux.txt"), 
                      header = TRUE, sep = "\t", skip=15,
                      row.names=NULL)

data <- tssmnth
names(data)[1:(ncol(data)-1)] <- names(data)[2:ncol(data)]
data[, ncol(data)] <- NULL
tssmnth <- data

## (1.3.2) Annual =====

tssann <- read.delim(paste0(wd, df, df2, 
                            "allsite_", variable, "_annualflux.txt"), 
                     header = TRUE, sep = "\t", skip=15,
                     row.names=NULL)

data <- tssann
names(data)[1:(ncol(data)-1)] <- names(data)[2:ncol(data)]
data[, ncol(data)] <- NULL
tssann <- data

## (1.3.3) Flux =====

tssf <- read.delim(paste0(wd, df, df2, 
                          "allsite_", variable, "_flux.txt"), 
                   header = TRUE, sep = "\t", skip=20,
                   row.names=NULL)

### ===== (1.4) TDS =====

df2 <- "20210930_09.54.53.799_TDS/outputs/"
variable <- "TDS"
## (1.4.1) Monthly =====

tdsmnth <- read.delim(paste0(wd, df, df2, 
                             "allsite_", variable, "_monthflux.txt"), 
                      header = TRUE, sep = "\t", skip=15,
                      row.names=NULL)

data <- tdsmnth
names(data)[1:(ncol(data)-1)] <- names(data)[2:ncol(data)]
data[, ncol(data)] <- NULL
tdsmnth <- data

## (1.4.2) Annual =====

tdsann <- read.delim(paste0(wd, df, df2, 
                            "allsite_", variable, "_annualflux.txt"), 
                     header = TRUE, sep = "\t", skip=15,
                     row.names=NULL)

data <- tdsann
names(data)[1:(ncol(data)-1)] <- names(data)[2:ncol(data)]
data[, ncol(data)] <- NULL
tdsann <- data

## (1.4.3) Flux =====

tdsf <- read.delim(paste0(wd, df, df2, 
                          "allsite_", variable, "_flux.txt"), 
                   header = TRUE, sep = "\t", skip=20,
                   row.names=NULL)

##### ========== (2) Merge data into one file ==========================================================================

### ===== (2.1) Add in identifying columns =====

pocmnth$var <- "poc"
docmnth$var <- "doc"
tssmnth$var <- "tss"
tdsmnth$var <- "tds"

pocann$var <- "poc"
docann$var <- "doc"
tssann$var <- "tss"
tdsann$var <- "tds"

pocf$var <- "poc"
docf$var <- "doc"
tssf$var <- "tss"
tdsf$var <- "tds"

### ===== (2.2) Merge the monthlies =====

mnth <- merge(docmnth, pocmnth, all=TRUE)
mnth <- merge(mnth, tssmnth, all=TRUE)
mnth <- merge(mnth, tdsmnth, all=TRUE)

mnth$timeunit <- "month"
mnth$month <- as.numeric(substr(as.character(mnth$Date), start = 1, stop = 2))
mnth$year <- substr(as.character(mnth$Date), start = 4, stop = 7)
mnth$Date2 <- as.Date(paste0(mnth$year, "-", mnth$month, "-01"))
mnth$year <- year(mnth$Date2)

### ===== (2.3) Merge the annuals and format anything necessary =====

ann <- merge(docann, pocann, all=TRUE)
ann <- merge(ann, tssann, all=TRUE)
ann <- merge(ann, tdsann, all=TRUE)

ann$timeunit <- "year"

### ===== (2.4) Merge fluxes =====

f <- merge(docf, pocf, all=TRUE)
f <- merge(f, tssf, all=TRUE)
f <- merge(f, tdsf, all=TRUE)

f$timeunit <- "day"

f$Date <- as.Date(f$Date, format="%m/%d/%Y")
f$year <- year(f$Date)
f$month <- month(f$Date)

### ===== (2.5) Remove commas from numeric columns of interest =====

mnth$Flow <- as.numeric(gsub(",","", mnth$Flow))
mnth$AMLE.Load <- as.numeric(gsub(",","", mnth$AMLE.Load))
mnth$AMLE.Conc <- as.numeric(gsub(",","", mnth$AMLE.Conc))

ann$Flow <- as.numeric(gsub(",","", ann$Flow))
ann$AMLE.Load <- as.numeric(gsub(",","", ann$AMLE.Load))
ann$AMLE.Conc <- as.numeric(gsub(",","", ann$AMLE.Conc))

f$Flow <- as.numeric(gsub(",","", f$Flow))
f$AMLE.Load <- as.numeric(gsub(",","", f$AMLE.Load))
f$AMLE.Conc <- as.numeric(gsub(",","", f$AMLE.Conc))

### ===== (2.4) Remove years with major month and annuals gaps =====

ann <- ann[!(ann$Date %in% qgaps$year),]
mnth <- mnth[!(mnth$Date %in% qgaps$year),]

### ===== (2.5) Annual: Add in days for scaling and account for leap years =====
ann$days <- 365
ann$days[ann$Date==1972 |
         ann$Date==1976 |
         ann$Date==1980 |
         ann$Date==1984 |
         ann$Date==1988 |
         ann$Date==1992 |
         ann$Date==1996 |
         ann$Date==2000 |
         ann$Date==2004 |
         ann$Date==2008 |
         ann$Date==2012 |
         ann$Date==2016] <- 366

# convert flow to km3/yr, original units were ft3/s
ann$diskm3yr <- ((ann$Flow/35.3147)/(10^9))*(60*60*24*ann$days)

# convert load to annual load in Gigamoles per year (original was kg/day)
ann$annualloadGMyr <- ((ann$AMLE.Load*1000)/12.0107)/(10^9)*ann$days
ann$milltonnesyr <- ((ann$AMLE.Load*ann$days)/1000)/(10^6)

### ===== (2.6) Monthly: Add in days for scaling and account for leap years =====

mdays <- as.data.frame(unique(mnth$month))
colnames(mdays) <- "month"
mdays$days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

mnth <- merge(mnth, mdays, by="month")

mnth$days[mnth$month==2 &
          (mnth$year==1972 |
           mnth$year==1976 |
           mnth$year==1980 |
           mnth$year==1984 |
           mnth$year==1988 |
           mnth$year==1992 |
           mnth$year==1996 |
           mnth$year==2000 |
           mnth$year==2004 |
           mnth$year==2008 |
           mnth$year==2012 |
           mnth$year==2016)] <- 29

##### ========== (3) Sum flux sheet and see if it's identical to annual ========

### ===== (3.1) Annual =====

fyr <- f %>%
  group_by(var, year) %>%
  summarize(annload=sum(AMLE.Load),
            annconc=mean(AMLE.Conc), 
            sdconc=mean(AMLE.Conc))

a <- merge(fyr, ann[,c("Date", "Flow", "AMLE.Load", "AMLE.Conc", "var", "days")],
              by.x=c("var", "year"), by.y=c("var", "Date"))

a$annloadscaled <- a$AMLE.Load*a$days
a$rat <- a$annloadscaled/a$annload

a <- a[a$rat < 1.0001,] 

### ===== (3.2) Monthly =====

fmn <- f %>%
  group_by(var, year, month) %>%
  summarize(mnthload=sum(AMLE.Load),
            mnthconc=mean(AMLE.Conc), 
            sdconc=mean(AMLE.Conc))

m <- merge(fmn, 
           mnth[,c("var", "year", "month", "Date", 
                   "Flow", "AMLE.Load", "AMLE.Conc", "var", "days")],
           by=c("var", "year", "month"))

m$mnthloadscaled <- m$AMLE.Load*m$days
m$rat <- m$mnthloadscaled/m$mnthload

m <- m[m$rat < 1.0001,] 

### ===== (3.3) Spring & Summer =====

m5to6 <- m[m$month==5|
             m$month==6,]

msum56 <- m5to6 %>%
  group_by(var, year) %>%
  summarize(springload=sum(mnthload)) 

m7to8 <- m[m$month==7|
           m$month==8,]

msum78 <- m7to8 %>%
        group_by(var, year) %>%
        summarize(summmerload=sum(mnthload)) 

a <- merge(a, msum56, by=c("var", "year"), all.x = TRUE)
a <- merge(a, msum78, by=c("var", "year"), all.x = TRUE)

a$summerperc <- (a$summmerload/a$annload)*100
a$springperc <- (a$springload/a$annload)*100

##### ========== (4) Do conversions necessary ==================================

### ===== (4.1) Annual =====

# convert flow to km3/yr, original units were ft3/s
a$diskm3yr <- ((a$Flow/35.3147)/(10^9))*(60*60*24*a$days)

# convert load to annual load in Gigamoles per year (original was kg/day)
a$annualloadGMyr <- ((a$AMLE.Load*1000)/12.0107)/(10^9)*a$days
a$milltonnesyr <- ((a$AMLE.Load*a$days)/1000)/(10^6)

### ===== (4.2) Month =====

# convert flow to km3/yr, original units were ft3/s
m$diskm3avm <- ((m$Flow/35.3147)/(10^9))*(60*60*24)

# convert load to annual load in Gigamoles per year (original was kg/day)
m$loadGMavm <- ((m$AMLE.Load*1000)/12.0107)/(10^9)
m$milltonnesavm <- ((m$AMLE.Load)/1000)/(10^6)

##### ========== (3) Graph data ================================================

## (3.1) Save theme ====================
theme <-theme(panel.grid.major = element_blank(),
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
              # aspect.ratio=1
              legend.key=element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              legend.spacing = unit(0, "mm"),
              legend.key.size = unit(0.001, "cm"),
              plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

paneltheme <- theme(strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    panel.spacing.x = unit(0.4, "inches"))

## (3.2) Discharge ====================
dis <- ggplot(data=a, 
       aes(x=year, y=diskm3yr)) +
       labs(y=expression("Annual Discharge (km"^3*" y"^-1*")"), x="Year") +
       geom_point(size=4, shape=21, fill="grey", colour="black") + theme

## (3.3) Other variables ====================
oc <- ggplot(data=a[a$var=="doc"|a$var=="poc",], 
       aes(x=year, y=annualloadGMyr)) +
       geom_point(aes(fill=var), size=4, shape=21,
                  #fill="grey", 
                  colour="black") +
       labs(y=expression("Annual OC Load (Gmol y"^-1*")"), x="Year") +
      # facet_wrap(.~var, scales="free_y") +
  theme #+ paneltheme


t <- ggplot(data=a[a$var=="tss"|a$var=="tds",], 
             aes(x=year, y=milltonnesyr)) +
  geom_point(size=4, shape=21, fill="grey", colour="black") +
  labs(y=expression("Annual Flux (Million Tonnes y"^-1*")"), x="Year") +
  facet_wrap(.~var, scales="free_y") + theme + paneltheme

## (3.4) Monthly changes ====================
mdis <- ggplot(data=m, 
              aes(x=year, y=diskm3avm)) +
  labs(y=expression("Monthly Discharge (km"^3*" d"^-1*")"), x="Year") +
  geom_point(size=4, shape=21, fill="grey", colour="black") + 
  facet_wrap(.~month, scales="free_y") + theme 

moc <- ggplot(data=m[m$var=="doc"|m$var=="poc",], 
             aes(x=year, y=loadGMavm)) +
  geom_point(aes(fill=var), size=4, shape=21, colour="black") +
  labs(y=expression("Average OC Load (Gmol d"^-1*")"), x="Year") +
  facet_wrap(.~month, scales="free_y") + theme 

mt <- ggplot(data=m[m$var=="tds"|m$var=="tss",], 
              aes(x=year, y=milltonnesavm)) +
  geom_point(aes(fill=var), size=4, shape=21, colour="black") +
  labs(y=expression("Annual Flux (Million Tonnes d"^-1*")"), x="Year") +
  facet_wrap(.~month, scales="free_y") + theme 

# increasing may and decreasing june is likely just showing a shift in hydrograph,
# potential disconnection between peak discharge thaw and thaw of sediments on the Plateau

## (3.5) Summer Proportion changes ====================

oc_summerperc <- ggplot(data=a[a$var=="doc"|a$var=="poc",], 
             aes(x=year, y=summerperc)) +
  geom_point(aes(fill=var), size=4, shape=21,
             #fill="grey", 
             colour="black") +
  labs(y="Proportion of Load in summer (June-Aug, %)", x="Year") +
  # facet_wrap(.~var, scales="free_y") +
  theme #+ paneltheme

oc_springperc <- ggplot(data=a[a$var=="doc"|a$var=="poc",], 
                        aes(x=year, y=springperc)) +
  geom_point(aes(fill=var), size=4, shape=21,
             #fill="grey", 
             colour="black") +
  labs(y="Proportion of Load in spring (May-June, %)", x="Year") +
  # facet_wrap(.~var, scales="free_y") +
  theme #+ paneltheme

t_summerperc <- ggplot(data=a[a$var=="tss"|a$var=="tds",], 
            aes(x=year, y=summerperc)) +
  geom_point(size=4, shape=21, fill="grey", colour="black") +
  labs(y="Proportion of Load in summer (June-Aug, %)", x="Year") +
  facet_wrap(.~var, scales="free_y") + theme + paneltheme

t_springperc <- ggplot(data=a[a$var=="tss"|a$var=="tds",], 
                       aes(x=year, y=springperc)) +
  geom_point(size=4, shape=21, fill="grey", colour="black") +
  labs(y="Proportion of Load in summer (June-Aug, %)", x="Year") +
  facet_wrap(.~var, scales="free_y") + theme + paneltheme

## (3.6) POC:DOC ratio ====================

# (3.6.1) Annual ====
awide <-  a %>%
  pivot_wider(id_cols=year,
              names_from=var,
              values_from=c(annualloadGMyr))

awide$pocdoc <- (awide$poc/awide$doc)

pocdocrat <- ggplot(data=awide, 
             aes(x=year, y=pocdoc)) +
  geom_point(size=4, shape=21,
             fill="grey", 
             colour="black") +
  labs(y="POC:DOC ratio", x="Year") +
  # facet_wrap(.~var, scales="free_y") +
  theme #+ paneltheme

# (3.6.2) Monthly ====

#mwide <-  m %>%
 # pivot_wider(id_cols=c(year, month),
  #            names_from=var,
   #           values_from=c(loadGMavm))

#mwide$pocdoc <- (mwide$poc/mwide$doc)

#pocdocrat <- ggplot(data=mwide, 
 #                   aes(x=year, y=pocdoc)) +
#  geom_point(size=4, shape=21,
 #            fill="grey", 
  #           colour="black") +
  #labs(y="POC:DOC ratio", x="Year") +
  # facet_wrap(.~month, scales="free_y") +
  #theme #+ paneltheme

##### ========== (4) Trend Analysis =================

library(zyp)

t <- lm(annualloadGMyr~year, trendpoc)

trendpoc <- a[a$var=="poc", c("year", "annualloadGMyr")]

zyp.trend.vector(trendpoc$annualloadGMyr, trendpoc$year,
                 method=c("zhang"), # need to figure out what method to use
                 conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)

trenddoc <- a[a$var=="doc", c("year", "annualloadGMyr")]

zyp.trend.vector(trenddoc$annualloadGMyr, trenddoc$year,
                 method=c("yuepilon", "zhang"),
                 conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)


trendtss <- a[a$var=="tss", c("year", "milltonnesyr")]

zyp.trend.vector(trendtss$milltonnesyr, trendtss$year,
                 method=c("yuepilon", "zhang"),
                 conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)

trendtss <- a[a$var=="tss" & a$year!=2013, c("year", "milltonnesyr")]

zyp.trend.vector(trendtss$milltonnesyr, trendtss$year,
                 method=c("yuepilon"),
                 conf.intervals=TRUE, preserve.range.for.sig.test=TRUE) # still not significant

trendtds <- a[a$var=="tds", c("year", "milltonnesyr")]

zyp.trend.vector(trendtds$milltonnesyr, trendtds$year,
                 method=c("yuepilon", "zhang"),
                 conf.intervals=TRUE, preserve.range.for.sig.test=TRUE)
