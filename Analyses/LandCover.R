## Calculate % surficial geology and vegetation in the watershed 

library(dplyr)
library(tidyr)


## Clear list 
list=rm(list=ls(all=TRUE))

# read in the updated watmaster file 
watmaster <- read.csv("D:/5_Projects/git/repos/Phd_Ch3/data/watmaster_wlakes_2026.csv")

colnames(watmaster)

# Surficial geology percentage columns
surf_cols <- c(
  "percshale",
  "colluvial_perc",
  "morainal_perc",
  "alluvial_perc",
  "glaciofluvial_perc",
  "organic_perc",
  "glaciolacustrine_perc"
)

# Calculate total area of all subcatchments combined
total_watershed_area <- sum(watmaster$WatershedArea, na.rm = TRUE)

# Calculate total area and percent cover for each surficial geology type
surficial_summary <- watmaster %>%
  mutate(
    across(all_of(surf_cols), ~ replace_na(.x, 0))
  ) %>%
  summarise(
    across(
      all_of(surf_cols),
      ~ sum(WatershedArea * (.x / 100), na.rm = TRUE),
      .names = "{.col}_total_area"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "surficial_geology",
    values_to = "total_area"
  ) %>%
  mutate(
    surficial_geology = gsub("_total_area$", "", surficial_geology),
    total_watershed_area = total_watershed_area,
    percent_of_total_watershed_area = (total_area / total_watershed_area) * 100
  )

surficial_summary


# each (check) 
#lichen/moss
# shrubland 
#forests  
# barren-land 

# Vegetation percentage columns
veg_cols <- c(
  "barrenland_perc",
  "forest_perc",
  "grassland_perc",
  "lichenmoss_perc",
  "shrubland_perc"
)

# Calculate total area of all subcatchments combined
total_watershed_area <- sum(watmaster$WatershedArea, na.rm = TRUE)

# Calculate total area and percent cover for each vegetation type
vegetation_summary <- watmaster %>%
  mutate(
    across(all_of(veg_cols), ~ replace_na(.x, 0))
  ) %>%
  summarise(
    across(
      all_of(veg_cols),
      ~ sum(WatershedArea * (.x / 100), na.rm = TRUE),
      .names = "{.col}_total_area"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "vegetation_type",
    values_to = "total_area"
  ) %>%
  mutate(
    vegetation_type = gsub("_total_area$", "", vegetation_type),
    total_watershed_area = total_watershed_area,
    percent_of_total_watershed_area = (total_area / total_watershed_area) * 100
  )

vegetation_summary
