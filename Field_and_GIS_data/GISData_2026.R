# ===========================================================================================================#
# GISData.R
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
#library(ggplot2)
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

##### ========== (1) DATA PREP =========================================================================

# (1.1) Watershed area file ==========
wa <- read_excel(paste0(df, "watershedareas.xlsx"), sheet="watershedareasforR")

# (1.2) Watershed bedrock geology ==========

# Calculate the percentage of each watershed (sub-catchment) underlain by shale 
# bedrock and merges that value into the watershed master table. 

bgeo <- read.csv(paste0(df, "watershedbedgeol.csv"))
bgeo$shale <- NA

# Mark polygon as "Y" (yes) if any of the 4 litholgoy columns contain "Shale"
# Keep as NA if no "Shale"

bgeo$shale[bgeo$Lithology1=="Shale" | 
            bgeo$Lithology2=="Shale" |
            bgeo$Lithology3=="Shale" |
            bgeo$Lithology4=="Shale"] <- "Y"

# Keep only the required columns 
bgeo <- bgeo %>% select(Name, polygonareakm2, shale) 

# Calculate the wateshed area from the bedrock polygons, grouping polygons
# by name 

bgeowat <- bgeo %>%
            group_by(Name) %>%
            summarise(sumwatarea = sum(polygonareakm2))

# Calculate area by shale category
bgeoshale <- bgeo %>%
              group_by(Name, shale) %>%
              summarise(sumshalearea = sum(polygonareakm2))
# merge shale area with total watershed area
bgeo <- merge(bgeoshale, bgeowat, by=c("Name"))
# calculate the percentage of shale 
bgeo$percshale <- (bgeo$sumshalearea/bgeo$sumwatarea)*100
# keep only shale rows 
bgeo <- bgeo[!is.na(bgeo$shale),]
# keep only watershed name and percent shale
bgeo <- bgeo %>%
         select(Name, percshale)

# merge with watershed area
watmaster <- merge(wa, bgeo, by="Name")

# (1.3) Watershed surficial geology ==========

# ch3wssurfgeowlakes_v2.csv -- > with lakes 
# ch3wssurfgeowolakes_v2.csv --> lake area clipped out 
# just use the wo (without) lakes 

sgeo <- read.csv(paste0(df, "ch3wssurfgeowolakes_v2.csv"))

# update sgeo with a new column - Name_2 - which includes
# the new categories 

sgeo <- sgeo %>%
  mutate(
    Name_2 = case_when(
      Name %in% c("Af", "ApAt", "At") ~ "Alluvial",
      Name %in% c("Ct", "Cv", "Cx", "Cy") ~ "Colluvial",
      Name == "fOpOfPO" ~ "Organic",
      Name %in% c("Gp", "Gt") ~ "Glaciofluvial",
      Name == "Lb" ~ "Glaciolacustrine",
      Name %in% c("MhMrMm", "MpMbMpv") ~ "Morainal",
      Name == "Mv" ~ "Morainal/Bedrock",
      Name == "R/Cv" ~ "Bedrock/Colluvium",
      TRUE ~ NA_character_
    ) 
  ) %>%
  relocate(Name_2, .after = Name)
  

# rename columns and keep relevant ones: geoltype, site name, and area 

sgeo <- sgeo %>% select(geoltype=Name_2, Name=Name_1, polygonareakm2=Area_km2geodesic)

# Name_1 is the watershed name 
# Calculate the total area of  all surficial geology polygons within each watershed  

sgeowat <- sgeo %>%
            group_by(Name) %>%
            summarise (sumwatarea = sum(polygonareakm2))

# calculate the area by surficial geology type

sgeotype <- sgeo %>%
  group_by(Name, geoltype) %>%
  summarise (sumtypearea = sum(polygonareakm2))

#t <- merge(wsgeowat, wbgeowat, by="Name")
#t$t <- t$sumwatarea.x/t$sumwatarea.y
# odd that there are some minor discrepencies, but lakes weren't removed, should be fine

# merge type-specific area with total watershed area
sgeo <- merge(sgeotype, sgeowat, by="Name")
# now, each type of surfifial geology has the total area of that type of SG, 
# and the total area in the watershed (allows us to calculate the percentage)

# calculate the percentage 
sgeo$sgperc <- (sgeo$sumtypearea/sgeo$sumwatarea)*100

# clean up - keep only what we need 
sgeo <- sgeo %>% select(Name, geoltype, sgperc)

# convert to long format so each surficial geology type becomes its own column
# note that there are many NAs 

# also note - these are the names that are used later. The df must be modified 
# with new SG names from the beginning, then the rest of the code must be updated 
sgeo <- sgeo %>%
         pivot_wider(id_cols=Name, names_from = geoltype, values_from = sgperc)

# merge with watershed area
watmaster <- merge(watmaster, sgeo, by="Name")

# quick check to see that the surficial geology unit % sum to 100

surf_cols <- c("Bedrock/Colluvium","Colluvial","Alluvial","Organic",  
               "Glaciofluvial","Glaciolacustrine","Morainal","Morainal/Bedrock")

surf_check <- watmaster %>%
  mutate(
    surf_sum = rowSums(across(all_of(surf_cols)), na.rm = TRUE),
    surf_diff_from_100 = surf_sum - 100
  ) %>%
  select(Name, all_of(surf_cols), surf_sum, surf_diff_from_100)

# make a quick version of watmaster (unmodified) as a backup

watmaster_clean <- watmaster

# (1.4) Watershed nhnmodified lake area km2 ==========

# calculate the percentrage of each watershed covered by lakes, and 
# merge that value into watmaster

wlake <- read.csv(paste0(df, "watershedlakearea_updated20210726.csv"))

wlakewat <- wlake %>%
  select (Name, polygonareakm2) %>%
  group_by(Name) %>%
  summarise (sumwatarea = sum(polygonareakm2))

# compare the watershed areas with the surifical geology watershed areas 
# they are almost identical 

t <- merge(wlakewat, sgeowat, by="Name")
t$t <- t$sumwatarea.x/t$sumwatarea.y
# sarah original note: odd that there are some minor discrepencies, but lakes weren't removed, should be fine

# calculate the total lake area per the watershed 

lakearea <- wlake %>%
  filter (waterDefinitionText=="Lake") %>%
  group_by(Name, waterDefinitionText) %>%
  summarise (sumlakearea = sum(polygonareakm2))

# merge the total watershed area and lake area 
wlake <- merge(wlakewat, lakearea, by="Name")

# calculate the percent of lake area 
wlake$lakeperc <- (wlake$sumlakearea/wlake$sumwatarea)*100
# the highest is 9%, most <3 % 

length(unique(watmaster$Name))
# 83 - good 

# merge with watershed area
watmaster <- merge(watmaster, wlake, by="Name", all.x=TRUE)

# replace missing lake percentages with 0
watmaster$lakeperc[is.na(watmaster$lakeperc)] <- 0

# (1.5) Watershed GPP mean ==========

# this section add watershed-level GPP to watmaster 

#wgppmean1 <- read.csv(paste0(df, "ShedsGPPmean_20170712_2.csv"))
#wgppmean1 <- wgppmean1 %>% select(Name, mean20170712=mean)
#wgppmean2 <- read.csv(paste0(df, "ShedsGPPmean_20170720_2.csv"))
#wgppmean2 <- wgppmean2 %>% select(Name, mean20170720=mean)
#wgppmean3 <- read.csv(paste0(df, "ShedsGPPmean_20170805_2.csv"))
#wgppmean3 <- wgppmean3 %>% select(Name, mean20170805=mean)

#wgpp <- merge(wgppmean1, wgppmean2, by="Name")
#wgpp <- merge(wgpp, wgppmean3, by="Name")

#wgpp$mean20170712scaled <- wgpp$mean20170712*0.0001
#wgpp$mean20170720scaled <- wgpp$mean20170720*0.0001
#wgpp$mean20170805scaled <- wgpp$mean20170805*0.0001

#wgppmeancalc <- wgpp %>%
# select(Name, mean20170712scaled, 
#       mean20170720scaled, mean20170805scaled) %>%
#pivot_longer(cols=c(mean20170712scaled, 
#                   mean20170720scaled,
#                  mean20170805scaled),
#          names_to="gpp8d_date",
#         values_to="gpp8d") %>%
#  group_by(Name) %>%
# summarize(mean8dgpp = mean(gpp8d))

wgppmean <- read.csv(paste0(df, 
                            "ShedsGPPmean_stackedsum20170704to20170813.csv"))
# scale the GPP 
wgppmean$scaledgpp <- wgppmean$mean*0.0001

# simplify to only include the columns needed
wgppmean_clean <- wgppmean[, c("Name", "scaledgpp")]

# merge with master file
watmaster <- merge(watmaster, wgppmean_clean, by="Name")

# (1.6) Watershed NPP mean ==========

# add NPP mean 

wnppmean <- read.csv(paste0(df, "ShedsannualNPPmean_2017.csv"))
wnppmean$scalednpp <- wnppmean$mean*0.0001

# again - clean 
wnppmean_clean <- wnppmean[, c("Name", "scalednpp")]

# merge with master file
watmaster <- merge(watmaster, wnppmean_clean, by="Name")

# (1.7) Watershed elevation ==========

wele <- read.csv(paste0(df, "meanwatershedelevation.csv"))
wele <- wele %>% select(Name, meanelev_m=MEAN)

# merge with master file
watmaster <- merge(watmaster, wele, by="Name")

# (1.8) Watershed slope ==========

wslope <- read.csv(paste0(df, "meanwatershedslopes.csv"))
wslope <- wslope %>% select(Name, meanslope_deg=MEAN)

# merge with master file
watmaster <- merge(watmaster, wslope, by="Name")

# (1.9) Land cover ==========

# calculate the percentage of each watershed covered by each land-cover category

wslc <- read.csv(paste0(df, "watershedlandcovernolakes.csv"))

wslc$gridcode[wslc$FID_canlc2015_nolakes_nad83zn8==-1] <- NA
wslc <- wslc %>% select(Name, gridcode, polygonareakm2)

length(unique(wslc$gridcode))
# there are 14 gridcodes 
unique(wslc$gridcode)
# 10, 8, 11,  6,  5, 12,  1, 13, 18, 16, NA, 17,  2, 14

sum(is.na(wslc$gridcode))
# 153 
# There are 153 with no landcover distinction 

# read in the code for the landcover 

# the original LandCoverCircaCode.xlsx does not exist
# the code waqs found online and added to the data 

#wslc_legend <- read_excel(paste0(df, "LandCoverCircaCode.xlsx"), sheet="2015")
wslc_legend <- read_excel(paste0(df, "LandCover2015.xls"))

# subset to only include the landcover groups that are in wslc
wslc_legend <- wslc_legend %>%
  filter(Value %in% unique(wslc$gridcode))

# Keep only the Value and Class_EN
wslc_legend <- wslc_legend[, c("Value", "Class_EN")]

# reassign the Class_EN using a simplified 8 groups 
# rename Class_EN to gridcode aswell 

wslc_legend <- wslc_legend %>%
  mutate(
    Class_EN = case_when(
      Value %in% c(1, 2, 5, 6) ~ "forest",
      Value == 8 ~ "shrubland",
      Value == 10 ~ "grassland",
      Value %in% c(11, 12, 13) ~ "lichenmoss",
      Value == 14 ~ "wetland",
      Value == 16 ~ "barrenland",
      Value == 17 ~ "urbanbuilt",
      Value == 18 ~ "water",
      TRUE ~ Class_EN
    )
  ) %>%
  rename(
    gridcode = Value,
    landcover = Class_EN
  )


wslc <- merge(wslc, wslc_legend, by="gridcode", all.x=TRUE)

# these lines no longer apply 
#wslc$level1desc[wslc$FID_canlc2015_nolakes_nad83zn8==-1] <- NA

#wslc$lccat_sarah <- wslc$level1desc
#wslc$lccat_sarah[wslc$level1desc=="forest_needleleaf" |
#                 wslc$level1desc=="forest_broadleaf" |
#                 wslc$level1desc=="forest_mixed"] <- "forest"

# code for this section has been summarized given the new structure of wslc 

# summarise land cover area by watershed and landcover category 
wslc <- wslc %>% 
        select(Name, landcover, polygonareakm2) %>%
        group_by(Name, landcover) %>%
        summarize(
          polygonareakm2=sum(polygonareakm2))

# calculate total wateshed area based on landcover polygons 
wslcwat <- wslc %>%
           group_by(Name) %>%
           summarize(
             sumwatarea=sum(polygonareakm2))

# compare the landcover watershe area with bedrock         
t <- merge(wslcwat, bgeowat, by="Name")
t$t <- t$sumwatarea.x/t$sumwatarea.y

summary(t$t)
# looks good 

# merge landcover categeory with totla wateshed area 

wslc <- merge(wslc, wslcwat, by="Name")

# calculate the percentage cover 
wslc$landcover_perc <- (wslc$polygonareakm2/wslc$sumwatarea)*100
        
# reshape landcover categroies into separate columns 
wslc <- wslc %>% 
             select(Name, landcover, landcover_perc) %>%
               pivot_wider(
                          id_cols = Name,
                          names_from = landcover,
                          values_from = landcover_perc
  )

# merge with master 
watmaster <- merge(watmaster, wslc, by="Name")

# remove the "NA" column - we don't need it. 
# this is undistinguishable land cover type that does not impact the data 

watmaster$'NA' <- NULL

# (1.10) Slump percentages ==========

# (1.10.1) 2016 Slump percentages ==========
slumps2016 <- read.csv(paste0(df, 
                   "peelslumps2016_union_ch3archydrowatersheds20210616_clip.csv"))


slumpsum2016 <- slumps2016 %>%
  select(FID_Peel_Slumps_2016, Name, polygonareakm2)%>%
  filter(FID_Peel_Slumps_2016!=-1) %>%
  group_by(Name) %>%
  summarize(sumslump16=sum(polygonareakm2))

slumpsumwat2016 <- slumps2016 %>%
  select(Name, polygonareakm2) %>%
  group_by(Name) %>%
  summarize(sumwatarea=sum(polygonareakm2))

t <- merge(slumpsumwat2016, bgeowat, by="Name")
t$t <- t$sumwatarea.x/t$sumwatarea.y

slumps2016 <- merge(slumpsum2016, slumpsumwat2016, by="Name")

slumps2016$percslump2016 <- (slumps2016$sumslump16/slumps2016$sumwatarea)*100

# keep only necessary columns before merge 
slumps2016 <- slumps2016[, c("Name", "percslump2016")]

#merge with master
watmaster <- merge(watmaster, slumps2016, by="Name", all.x=TRUE)

# (1.10.1) 2017 Slump percentages SEE NOTE ==========

# SKIP -- NOT IN THE FINAL WATMASTER

#slumps17old <- read.csv(paste0(df, 
 #                             "rtsdelineationsstony_union_watersheds2021061.csv"))

## UH-OH!! ## 
slumps17act <- read.csv(paste0(df, 
                              "rtsdelineationsstonysarahaddCB_union_watersheds2021061.csv"))

# ^^ this file does not exist.
# it is not "rtsdelineationsstony_union_watersheds2021061"
# not sure - proceed without for now 

slumps17all <- read.csv(paste0(df, 
                               "allslumpsclippedtowatershed_union_watersheds2021061.csv"))

# slumps 2017 old =====


#slumpsum17old <- slumps17old %>%
#  select(FID_rts_delineations_stony, Name, Shape_Area)%>%
#  filter(FID_rts_delineations_stony!=-1) %>%
#  group_by(Name) %>%
#  summarize(sumslump17=sum(Shape_Area))

#slumpsumwat <- slumps17old %>%
#  select(Name, Shape_Area) %>%
#  group_by(Name) %>%
#  summarize(sumwatarea=sum(Shape_Area))
#t <- merge(slumpsumwat, bgeowat, by="Name")
#t$t <- t$sumwatarea.x/t$sumwatarea.y

#slumps17 <- merge(slumpsum17old, slumpsumwat, by="Name")
#slumps17$percslump17old <- (slumps17$sumslump17/slumps17$sumwatarea)*100

#slumps17old <- slumps17 %>% select(Name, percslump17old)

#merge with master
#watmaster <- merge(watmaster, slumps17old, by="Name", all.x=TRUE)

# slumps 2017active =====

## SKIPPED FOR NOWW ## 
slumpsum17act <- slumps17act %>%
  select(FID_rts_delineations_stony_sarahaddCB, Name, Shape_Area)%>%
  filter(FID_rts_delineations_stony_sarahaddCB!=-1) %>%
  group_by(Name) %>%
  summarize(sumslump17=sum(Shape_Area))

slumpsumwat <- slumps17act %>%
  select(Name, Shape_Area) %>%
  group_by(Name) %>%
  summarize(sumwatarea=sum(Shape_Area))
t <- merge(slumpsumwat, bgeowat, by="Name")
t$t <- t$sumwatarea.x/t$sumwatarea.y

slumps17 <- merge(slumpsum17act, slumpsumwat, by="Name")
slumps17$percslump17act <- (slumps17$sumslump17/slumps17$sumwatarea)*100

slumps17act <- slumps17 %>% select(Name, percslump17act)

#merge with master
watmaster <- merge(watmaster, slumps17act, by="Name", all.x=TRUE)

# slumps2017all =====

slumpsum17all <- slumps17all %>%
  select(FID_allslumpsclippedtowatershed, Name, Shape_Area)%>%
  filter(FID_allslumpsclippedtowatershed!=-1) %>%
  group_by(Name) %>%
  summarize(sumslump17=sum(Shape_Area))

slumpsumwat <- slumps17all %>%
  select(Name, Shape_Area) %>%
  group_by(Name) %>%
  summarize(sumwatarea=sum(Shape_Area))
t <- merge(slumpsumwat, bgeowat, by="Name")
t$t <- t$sumwatarea.x/t$sumwatarea.y

slumps17 <- merge(slumpsum17all, slumpsumwat, by="Name")
slumps17$percslump17all <- (slumps17$sumslump17/slumps17$sumwatarea)*100

slumps17all <- slumps17 %>% select(Name, percslump17all)

#merge with master
watmaster <- merge(watmaster, slumps17all, by="Name", all.x=TRUE)

# (1.11) Watercourse slump accumulations ==========

slumpacc <- read.csv(paste0(df, "watercourseacccountslumpjk_union_watersheds2021061.csv"))

# acc = strahler type stream order for slumping, acc_count=cumulative # of slumps
slumpacc <- slumpacc%>%
  select(Name, acc, acc_count, Strahler) %>%
  group_by(Name)%>%
  summarize(strahlerimpactacc=max(acc), 
            slumpacccount=max(acc_count),
            strahlerstream=max(Strahler))

#merge with master
watmaster <- merge(watmaster, slumpacc, by="Name", all.x=TRUE)



# make manual adjustments to the data 

watmaster$strahlerimpactacc[watmaster$Name=="site42"] <- 1
watmaster$strahlerimpactacc[watmaster$Name=="site48"] <- 1

watmaster$slumpacccount[watmaster$Name=="site42"] <- 1
watmaster$slumpacccount[watmaster$Name=="site48"] <- 1

watmaster$strahlerimpactacc[watmaster$year==2017 &
                            is.na(watmaster$strahlerimpactacc)] <- 0
watmaster$slumpacccount[watmaster$year==2017 &
                              is.na(watmaster$slumpacccount)] <- 0


# (1.12) 100 cm organic carbon stocks ==========

oc <- read.csv(paste0(df, "NCSCDv2Canada_int_watersheds2021061.csv"))

oc <- oc %>%
  select(Name, Shape_Area, SOCC_100CM)

ocwat <- oc %>%
  select(Name, Shape_Area) %>%
  group_by(Name) %>%
  summarize(sumwat=sum(Shape_Area))

oc <- merge(oc, ocwat)

oc$fracarea <- oc$Shape_Area/oc$sumwat

oc <- oc %>% 
  group_by(Name) %>%
  summarize(wmeanSOCC_100CM = weighted.mean(SOCC_100CM, fracarea)) %>%
  select(Name, wmeanSOCC_100CM)

# merge with watmaster

watmaster <- merge(watmaster, oc, all=TRUE)

# (1.13) Distance from western glacial limit ==========

gl <- read.csv(paste0(df, "distfromglaciallimit.csv"))

gl$NEAR_DIST[gl$site==1] <- gl$NEAR_DIST[gl$site==1]*-1
gl$NEAR_DIST[gl$site==45] <- gl$NEAR_DIST[gl$site==45]*-1
gl$NEAR_DIST[gl$site==30] <- gl$NEAR_DIST[gl$site==30]*-1
gl$NEAR_DIST[gl$site==31] <- gl$NEAR_DIST[gl$site==31]*-1
gl$NEAR_DIST[gl$site==32] <- gl$NEAR_DIST[gl$site==32]*-1

gl$gldistkm <- gl$NEAR_DIST/1000

gl <- gl %>%
     select(site, loc, trans, date, 
            lat, long, year, campaign,
            trannum, loc2, gldistkm)

gl$loc <- as.character(gl$loc)
gl$loc2 <- as.character(gl$loc2)

gl$loc[is.na(gl$loc)] <- "NA"
gl$loc2[is.na(gl$loc2)] <- "NA"

gl$date <- as.Date(gl$date)

watmaster$date <- as.Date(watmaster$date)

watmaster <- merge(watmaster, gl, all=TRUE)

watmaster$gldistkm[watmaster$site==27 & !is.na(watmaster$Name)] <- watmaster$gldistkm[watmaster$site==27 & is.na(watmaster$Name)]

watmaster <- watmaster[!(watmaster$site==27 & is.na(watmaster$Name)),]

# (1.14) Terrain roughness ==========

rough <-  read.csv(paste0(df, "watersheds2021061roughness_zonalstatistics.csv"))

rough <- rough %>% select(Name, meanrough=MEAN)

#merge with watmaster
watmaster <- merge(watmaster, rough, all=TRUE)

# (1.15) Revised site 10 details ==========

watmaster$WatershedArea[watmaster$site==10] <- 105.7384634815
watmaster$meanrough[watmaster$site==10] <- 2.84464457549354
watmaster$meanelev_m[watmaster$site==10] <- 449.828126537548
watmaster$meanslope_deg[watmaster$site==10] <- 5.64077872296588

#SOC site 10 revised =====
oc <- read.csv(paste0(df, "NCSCDv2Canada_int_revisedsite10wpoly.csv"))

oc <- oc %>%
  select(Shape_Area, SOCC_100CM)

ocwat <- oc %>%
  summarize(sumwat=sum(Shape_Area))

oc <- merge(oc, ocwat)

oc$fracarea <- oc$Shape_Area/oc$sumwat

oc <- oc %>% 
  summarize(wmeanSOCC_100CM = weighted.mean(SOCC_100CM, fracarea)) %>%
  select(wmeanSOCC_100CM)

# edit in watmaster 
watmaster$wmeanSOCC_100CM[watmaster$site==10] <- oc$wmeanSOCC_100CM

#slumpacccount site 10 revised ======

slumpacc <- read.csv(paste0(df, "watercourseacccountslumpjk_union_revisedsite10wpoly.csv"))

# acc = strahler type stream order for slumping, acc_count=cumulative # of slumps
slumpacc <- slumpacc%>%
  select(acc, acc_count, Strahler) %>%
  summarize(strahlerimpactacc=max(acc), 
            slumpacccount=max(acc_count),
            strahlerstream=max(Strahler))

# edit in watmaster 
watmaster$strahlerimpactacc[watmaster$site==10] <- slumpacc$strahlerimpactacc
watmaster$slumpacccount[watmaster$site==10] <- slumpacc$slumpacccount
watmaster$strahlerstream[watmaster$site==10] <- slumpacc$strahlerstream

# slumps 2017active site 10 revised ===== 

## SKIP THIS SECTION SINCE SOURCE FILE FOR slumps17act IS MISSING 

slumps17act <- read.csv(paste0(df, 
                               "rtsdelineationsstonysarahaddCB_union_revisedsite10wpoly_clip.csv"))

slumpsum17act <- slumps17act %>%
  select(FID_rts_delineations_stony_sarahaddCB, Shape_Area)%>%
  filter(FID_rts_delineations_stony_sarahaddCB!=-1) %>%
  summarize(sumslump17=sum(Shape_Area))

slumpsumwat <- slumps17act %>%
  select(Shape_Area) %>%
  summarize(sumwatarea=sum(Shape_Area))

slumps17 <- merge(slumpsum17act, slumpsumwat)
slumps17$percslump17act <- (slumps17$sumslump17/slumps17$sumwatarea)*100

# edit in watmaster 
watmaster$percslump17act[watmaster$site==10] <- slumps17$percslump17act

# slumps2017all site 10 revised =====

slumps17all <- read.csv(paste0(df, 
                               "allslumpsclippedtowatershed_union_revisedsite10wpoly_clip.csv"))
slumpsum17all <- slumps17all %>%
  select(FID_allslumpsclippedtowatershed, Shape_Area)%>%
  filter(FID_allslumpsclippedtowatershed!=-1) %>%
  summarize(sumslump17=sum(Shape_Area))

slumpsumwat <- slumps17all %>%
  select(Shape_Area) %>%
  summarize(sumwatarea=sum(Shape_Area))

slumps17 <- merge(slumpsum17all, slumpsumwat)
slumps17$percslump17all <- (slumps17$sumslump17/slumps17$sumwatarea)*100

# edit in watmaster 
watmaster$percslump17all[watmaster$site==10] <- slumps17$percslump17all

# landcover site 10 revised =====

# This section  modified since LandCoverCircaCode does not exist.
# New LandCover2025 file used to update the codes 

wslc <- read.csv(paste0(df, "canlc2015nolakesclip_union_revsites10.csv"))

# See sectionj 1.9 (Landcover). Read in the file and run the clean-up 
# to get the simplified categories 

# Do not run this line - file does not exist 
#wslc_legend <- read_excel(paste0(df, "LandCoverCircaCode.xlsx"), sheet="2015")


wslc$gridcode[wslc$FID_canlc2015_nolakes_nad83zn8==-1] <- NA

wslc <- wslc %>% select(gridcode, Shape_Area)

# Note that the first row has "0" value. Not sure why. Proceed anyway. 

wslc <- merge(wslc, wslc_legend, by="gridcode", all.x=TRUE)

# this replaces the -1 value with NA from the row that had a gridcode value of 0
# this does not apply anymore - it is already NA 
#wslc$level1desc[wslc$FID_canlc2015_nolakes_nad83zn8==-1] <- NA

#wslc$lccat_sarah <- wslc$level1desc

# not needed since the landcover categories are already simplfiied 

#wslc$lccat_sarah[wslc$level1desc=="forest_needleleaf" |
#                   wslc$level1desc=="forest_broadleaf" |
#                   wslc$level1desc=="forest_mixed"] <- "forest"

# assume that "lccat_sarah" is "land cover category" 
# this adjusted - lccat_sarah became landcover
wslc <- wslc %>% 
  select(landcover, Shape_Area) %>%
  group_by(landcover) %>%
  summarize(polygonareakm2=sum(Shape_Area))

# calculate the total watershed area
wslcwat <- wslc %>%
  summarize(sumwatarea=sum(polygonareakm2))

# merge together, then calculate the percentage 
# lccat_sarah_perc changed to landcover_perc
wslc <- merge(wslc, wslcwat)
wslc$landcover_perc <- (wslc$polygonareakm2/wslc$sumwatarea)*100


wslc$site <- "10"

# convert to long format - now it is ready to add back into the watmaster 
wslc <- wslc %>% 
  select(site, landcover, landcover_perc)%>%
  pivot_wider(names_from=landcover,
              values_from=landcover_perc)

#correct in watmaster file
# replace these values directly, skipping this block of code (it is not 
# necessary withe new setup) 

# watmaster has only one row where site == 10, so it is an easy replace. 

#watmaster$`barren land`[watmaster$site==10] <- wslc$`barren land`
#watmaster$forest[watmaster$site==10] <- wslc$forest
#watmaster$grassland[watmaster$site==10] <- wslc$grassland
#watmaster$`lichen/moss`[watmaster$site==10] <- wslc$`lichen/moss`
#watmaster$shrubland[watmaster$site==10] <- wslc$shrubland
#watmaster$`urban and built-up`[watmaster$site==10] <- wslc$`urban and built-up`
#watmaster$water[watmaster$site==10] <- wslc$water
#watmaster$wetland[watmaster$site==10] <- wslc$wetland


# Omit the "NA" category - -
cols_to_replace <- c(
  "barrenland", "forest", "grassland", "lichenmoss",
  "shrubland", "urbanbuilt", "water", "wetland"
)

watmaster[watmaster$site == 10, cols_to_replace] <- 
  wslc[wslc$site == 10, cols_to_replace]


# Watershed bedrock geology site 10 revised ==========
bgeo <- read.csv(paste0(df, "watershedbedgeol_nad83zn8_revsite10.csv"))

bgeo$shale <- NA
bgeo$shale[bgeo$Lithology1=="Shale" | 
             bgeo$Lithology2=="Shale" |
             bgeo$Lithology3=="Shale" |
             bgeo$Lithology4=="Shale"] <- "Y"

bgeo <- bgeo %>% select(Shape_Area, shale) 

bgeowat <- bgeo %>%
  summarise(sumwatarea = sum(Shape_Area))

bgeoshale <- bgeo %>%
  group_by(shale) %>%
  summarise(sumshalearea = sum(Shape_Area))

bgeo <- merge(bgeoshale, bgeowat)
bgeo$percshale <- (bgeo$sumshalearea/bgeo$sumwatarea)*100
bgeo <- bgeo[!is.na(bgeo$shale),]

bgeo <- bgeo %>%
  select(percshale)

#correct in watmaster file
watmaster$percshale[watmaster$site==10] <- bgeo$percshale

# Watershed surficial geology site 10 revised ==========

## SKIP THIS SECTION SINCE THE SURFICIAL GEOLOGY WAS ALREADY RE-DONE 

sgeo <- read.csv(paste0(df, "watershedgeol_nad83zn8_revsite10.csv"))

sgeo <- sgeo %>% select(geoltype=Name, Shape_Area)
sgeowat <- sgeo %>%
  summarise (sumwatarea = sum(Shape_Area))
sgeotype <- sgeo %>%
  group_by(geoltype) %>%
  summarise (sumtypearea = sum(Shape_Area))

sgeo <- merge(sgeotype, sgeowat)

sgeo$sgperc <- (sgeo$sumtypearea/sgeo$sumwatarea)*100
sgeo <- sgeo %>% select(geoltype, sgperc)
sgeo$site <- "10"
sgeo <- sgeo %>%
  pivot_wider(id_cols=site, names_from = geoltype, values_from = sgperc)

#correct in watmaster file
watmaster$Alluvial[watmaster$site==10] <- sgeo$Alluvial
watmaster$Colluvial[watmaster$site==10] <- sgeo$Colluvial
watmaster$Moraine[watmaster$site==10] <- sgeo$Moraine
watmaster$Organic[watmaster$site==10] <- sgeo$Organic
watmaster$Piedmont[watmaster$site==10] <- sgeo$Piedmont

# Watershed NHN lake site 10 revised ==========

wlake <- read.csv(paste0(df, "NHNLakes_union_revsite10.csv"))

wlakewat <- wlake %>%
  select (Shape_Area) %>%
  summarise (sumwatarea = sum(Shape_Area))

lakearea <- wlake %>%
  filter (waterDefinitionText=="Lake") %>%
  group_by(waterDefinitionText) %>%
  summarise (sumlakearea = sum(Shape_Area))

wlake <- merge(wlakewat, lakearea)
wlake$lakeperc <- (wlake$sumlakearea/wlake$sumwatarea)*100

#correct in watmaster file
watmaster$lakeperc[watmaster$site==10] <- wlake$lakeperc

# GPP and NPP site 10 revised ==========

watmaster$scaledgpp[watmaster$site==10] <- 2505.40446373062*0.0001 # see google drive file
watmaster$scalednpp[watmaster$site==10] <- 2955.97681159213*0.0001

# (2.2) Create master file ==========


# The code has been adjusted since most clean-up is already done 
# Use this opportunity to rename the columns to match with the original "watmaster"
# file. 

#watmaster <- watmaster %>%
#             select(Name,
#                    site, loc, trans,
#                    date, lat, long,
#                    year, campaign, trannum, loc2,
#                    WatershedArea, percshale,
#                    colluvial_perc=Colluvial, 
#                    piedmont_perc=Piedmont,
#                    alluvial_perc=Alluvial,
#                    bedrock_perc=Bedrock, 
#                    fluvial_perc=Fluvial, 
#                    glaciogenic_perc=Glaciogenic,
#                    moraine_perc=Moraine, 
#                    organic_perc=Organic, 
#                    lakeperc,
#                    scaledgpp, scalednpp,
#                    meanelev_m, meanslope_deg,
#                    barrenland_perc=`barren land`,
#                    #forestneedle_perc=forest_needleleaf,
#                   #forestbroad_perc=forest_broadleaf,
#                    #forestmixed_perc=forest_mixed,
#                    forest_perc=forest,
#                    grassland_perc=grassland,
#                    lichenmoss_perc=`lichen/moss`,
#                    shrubland_perc=shrubland,
#                    water_perc=water,
#                    urbanbuilt_perc=`urban and built-up`,
#                    wetland_perc=wetland,
#                    percslump2016,
#                    percslump17act,
#                    percslump17all,
#                    slumpacccount,
#                    strahlerimpactacc, 
#                    wmeanSOCC_100CM,
#                    gldistkm, meanrough,
#                    strahlerstream)

# check the column name of watmaster
colnames(watmaster)

# remove columns that are unnecessary 
cols_to_remove <- c("notes", "Shape_Length", "Shape_Area", "sumwatarea",
                    "waterDefinitionText", "sumlakearea")

watmaster <- watmaster %>% 
  select(-all_of(cols_to_remove))

# rename some of the columns 
watmaster <- watmaster %>%
  rename(
    bedrockcolluvium_perc = `Bedrock/Colluvium`,
    colluvial_perc = Colluvial,
    morainal_perc = Morainal,
    morainalbedrock_perc = `Morainal/Bedrock`,
    alluvial_perc = Alluvial,
    glaciofluvial_perc = Glaciofluvial,
    organic_perc = Organic,
    glaciolacustrine_perc = Glaciolacustrine,
    barrenland_perc = barrenland,
    forest_perc = forest,
    grassland_perc = grassland,
    lichenmoss_perc = lichenmoss,
    shrubland_perc = shrubland,
    water_perc = water,
    urbanbuilt_perc = urbanbuilt,
    wetland_perc = wetland
  )

# check how variables compare to each other
#pairs(watmaster[,12:ncol(watmaster)])
# consider just removing land cover...

# print watmaster
write.csv(watmaster, paste0(df, "watmaster_2026.csv"))
