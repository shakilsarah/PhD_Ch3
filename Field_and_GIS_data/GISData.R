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
bgeo <- read.csv(paste0(df, "watershedbedgeol.csv"))
bgeo$shale <- NA
bgeo$shale[bgeo$Lithology1=="Shale" | 
            bgeo$Lithology2=="Shale" |
            bgeo$Lithology3=="Shale" |
            bgeo$Lithology4=="Shale"] <- "Y"
bgeo <- bgeo %>% select(Name, polygonareakm2, shale) 
bgeowat <- bgeo %>%
            group_by(Name) %>%
            summarise(sumwatarea = sum(polygonareakm2))
bgeoshale <- bgeo %>%
              group_by(Name, shale) %>%
              summarise(sumshalearea = sum(polygonareakm2))
bgeo <- merge(bgeoshale, bgeowat, by=c("Name"))
bgeo$percshale <- (bgeo$sumshalearea/bgeo$sumwatarea)*100
bgeo <- bgeo[!is.na(bgeo$shale),]
bgeo <- bgeo %>%
         select(Name, percshale)

# merge with watershed area
watmaster <- merge(wa, bgeo, by="Name")

# (1.3) Watershed surficial geology ==========
sgeo <- read.csv(paste0(df, "surficialgeolareacote.csv"))
sgeo <- sgeo %>% select(geoltype=Name, Name=Name_1, polygonareakm2)
sgeowat <- sgeo %>%
            group_by(Name) %>%
            summarise (sumwatarea = sum(polygonareakm2))
sgeotype <- sgeo %>%
  group_by(Name, geoltype) %>%
  summarise (sumtypearea = sum(polygonareakm2))

#t <- merge(wsgeowat, wbgeowat, by="Name")
#t$t <- t$sumwatarea.x/t$sumwatarea.y
# odd that there are some minor discrepencies, but lakes weren't removed, should be fine
sgeo <- merge(sgeotype, sgeowat, by="Name")

sgeo$sgperc <- (sgeo$sumtypearea/sgeo$sumwatarea)*100
sgeo <- sgeo %>% select(Name, geoltype, sgperc)
sgeo <- sgeo %>%
         pivot_wider(id_cols=Name, names_from = geoltype, values_from = sgperc)

# merge with watershed area
watmaster <- merge(watmaster, sgeo, by="Name")

# (1.4) Watershed nhnmodified lake area km2 ==========
wlake <- read.csv(paste0(df, "watershedlakearea_updated20210726.csv"))

wlakewat <- wlake %>%
  select (Name, polygonareakm2) %>%
  group_by(Name) %>%
  summarise (sumwatarea = sum(polygonareakm2))

t <- merge(wlakewat, sgeowat, by="Name")
t$t <- t$sumwatarea.x/t$sumwatarea.y
# odd that there are some minor discrepencies, but lakes weren't removed, should be fine

lakearea <- wlake %>%
  filter (waterDefinitionText=="Lake") %>%
  group_by(Name, waterDefinitionText) %>%
  summarise (sumlakearea = sum(polygonareakm2))

wlake <- merge(wlakewat, lakearea, by="Name")
wlake$lakeperc <- (wlake$sumlakearea/wlake$sumwatarea)*100

# merge with watershed area
watmaster <- merge(watmaster, wlake, by="Name", all.x=TRUE)
watmaster$lakeperc[is.na(watmaster$lakeperc)] <- 0

# (1.5) Watershed GPP mean ==========

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
wgppmean$scaledgpp <- wgppmean$mean*0.0001

# merge with master file
watmaster <- merge(watmaster, wgppmean, by="Name")

# (1.6) Watershed NPP mean ==========
wnppmean <- read.csv(paste0(df, "ShedsannualNPPmean_2017.csv"))
wnppmean$scalednpp <- wnppmean$mean*0.0001

# merge with master file
watmaster <- merge(watmaster, wnppmean, by="Name")

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

wslc <- read.csv(paste0(df, "watershedlandcovernolakes.csv"))
wslc_legend <- read_excel(paste0(df, "LandCoverCircaCode.xlsx"), sheet="2015")
wslc$gridcode[wslc$FID_canlc2015_nolakes_nad83zn8==-1] <- NA
wslc <- wslc %>% select(Name, gridcode, polygonareakm2)
wslc <- merge(wslc, wslc_legend, by="gridcode", all.x=TRUE)
wslc$level1desc[wslc$FID_canlc2015_nolakes_nad83zn8==-1] <- NA
wslc$lccat_sarah <- wslc$level1desc
wslc$lccat_sarah[wslc$level1desc=="forest_needleleaf" |
                 wslc$level1desc=="forest_broadleaf" |
                 wslc$level1desc=="forest_mixed"] <- "forest"

wslc <- wslc %>% 
        select(Name, lccat_sarah, polygonareakm2) %>%
        group_by(Name, lccat_sarah) %>%
        summarize(polygonareakm2=sum(polygonareakm2))

wslcwat <- wslc %>%
           group_by(Name) %>%
           summarize(sumwatarea=sum(polygonareakm2))
        t <- merge(wslcwat, bgeowat, by="Name")
        t$t <- t$sumwatarea.x/t$sumwatarea.y

wslc <- merge(wslc, wslcwat, by="Name")
wslc$lccat_sarah_perc <- (wslc$polygonareakm2/wslc$sumwatarea)*100
        
wslc <- wslc %>% 
        pivot_wider(id_cols=Name,
                    names_from=lccat_sarah,
                    values_from=lccat_sarah_perc)

# merge with master 
watmaster <- merge(watmaster, wslc, by="Name")

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


#merge with master
watmaster <- merge(watmaster, slumps2016, by="Name", all.x=TRUE)

# (1.10.1) 2017 Slump percentages ==========
slumps17old <- read.csv(paste0(df, 
                              "rtsdelineationsstony_union_watersheds2021061.csv"))


slumps17act <- read.csv(paste0(df, 
                              "rtsdelineationsstonysarahaddCB_union_watersheds2021061.csv"))

slumps17all <- read.csv(paste0(df, 
                               "allslumpsclippedtowatershed_union_watersheds2021061.csv"))

# slumps 2017 old =====
slumpsum17old <- slumps17old %>%
  select(FID_rts_delineations_stony, Name, Shape_Area)%>%
  filter(FID_rts_delineations_stony!=-1) %>%
  group_by(Name) %>%
  summarize(sumslump17=sum(Shape_Area))

slumpsumwat <- slumps17old %>%
  select(Name, Shape_Area) %>%
  group_by(Name) %>%
  summarize(sumwatarea=sum(Shape_Area))
t <- merge(slumpsumwat, bgeowat, by="Name")
t$t <- t$sumwatarea.x/t$sumwatarea.y

slumps17 <- merge(slumpsum17old, slumpsumwat, by="Name")
slumps17$percslump17old <- (slumps17$sumslump17/slumps17$sumwatarea)*100

slumps17old <- slumps17 %>% select(Name, percslump17old)

#merge with master
watmaster <- merge(watmaster, slumps17old, by="Name", all.x=TRUE)

# slumps 2017active =====

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

# merge with master
watmaster <- merge(watmaster, gl, all=TRUE)

watmaster$gldistkm[watmaster$site==27 & !is.na(watmaster$Name)] <- watmaster$gldistkm[watmaster$site==27 & is.na(watmaster$Name)]

watmaster <- watmaster[!(watmaster$site==27 & is.na(watmaster$Name)),]

# (1.14) Terrain roughness ==========

rough <-  read.csv(paste0(df, "watersheds2021061roughness_zonalstatistics.csv"))

rough <- rough %>% select(Name, meanrough=MEAN)

#merge with watmaster
watmaster <- merge(watmaster, rough, all=TRUE)

# (1.13) Revised site 10 details ==========

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

wslc <- read.csv(paste0(df, "canlc2015nolakesclip_union_revsites10.csv"))
wslc_legend <- read_excel(paste0(df, "LandCoverCircaCode.xlsx"), sheet="2015")
wslc$gridcode[wslc$FID_canlc2015_nolakes_nad83zn8==-1] <- NA
wslc <- wslc %>% select(gridcode, Shape_Area)
wslc <- merge(wslc, wslc_legend, by="gridcode", all.x=TRUE)
wslc$level1desc[wslc$FID_canlc2015_nolakes_nad83zn8==-1] <- NA
wslc$lccat_sarah <- wslc$level1desc
wslc$lccat_sarah[wslc$level1desc=="forest_needleleaf" |
                   wslc$level1desc=="forest_broadleaf" |
                   wslc$level1desc=="forest_mixed"] <- "forest"

wslc <- wslc %>% 
  select(lccat_sarah, Shape_Area) %>%
  group_by(lccat_sarah) %>%
  summarize(polygonareakm2=sum(Shape_Area))

wslcwat <- wslc %>%
  summarize(sumwatarea=sum(polygonareakm2))

wslc <- merge(wslc, wslcwat)
wslc$lccat_sarah_perc <- (wslc$polygonareakm2/wslc$sumwatarea)*100

wslc$site <- "10"
wslc <- wslc %>% 
  select(site, lccat_sarah, lccat_sarah_perc)%>%
  pivot_wider(names_from=lccat_sarah,
              values_from=lccat_sarah_perc)

#correct in watmaster file
watmaster$`barren land`[watmaster$site==10] <- wslc$`barren land`
watmaster$forest[watmaster$site==10] <- wslc$forest
watmaster$grassland[watmaster$site==10] <- wslc$grassland
watmaster$`lichen/moss`[watmaster$site==10] <- wslc$`lichen/moss`
watmaster$shrubland[watmaster$site==10] <- wslc$shrubland
watmaster$`urban and built-up`[watmaster$site==10] <- wslc$`urban and built-up`
watmaster$water[watmaster$site==10] <- wslc$water
watmaster$wetland[watmaster$site==10] <- wslc$wetland

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
watmaster <- watmaster %>%
             select(Name,
                    site, loc, trans,
                    date, lat, long,
                    year, campaign, trannum, loc2,
                    WatershedArea, percshale,
                    colluvial_perc=Colluvial, 
                    piedmont_perc=Piedmont,
                    alluvial_perc=Alluvial,
                    bedrock_perc=Bedrock, 
                    fluvial_perc=Fluvial, 
                    glaciogenic_perc=Glaciogenic,
                    moraine_perc=Moraine, 
                    organic_perc=Organic, 
                    lakeperc,
                    scaledgpp, scalednpp,
                    meanelev_m, meanslope_deg,
                    barrenland_perc=`barren land`,
                    #forestneedle_perc=forest_needleleaf,
                    #forestbroad_perc=forest_broadleaf,
                    #forestmixed_perc=forest_mixed,
                    forest_perc=forest,
                    grassland_perc=grassland,
                    lichenmoss_perc=`lichen/moss`,
                    shrubland_perc=shrubland,
                    water_perc=water,
                    urbanbuilt_perc=`urban and built-up`,
                    wetland_perc=wetland,
                    percslump2016,
                    percslump17act,
                    percslump17all,
                    slumpacccount,
                    strahlerimpactacc, 
                    wmeanSOCC_100CM,
                    gldistkm, meanrough,
                    strahlerstream)

# check how variables compare to each other

pairs(watmaster[,12:ncol(watmaster)])
# consider just removing land cover...

# print watmaster
write.csv(watmaster, paste0(df, "watmaster.csv"))
