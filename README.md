# PhD_Ch3
Code for Chapter 3 of PhD paper

# Order for Analyses

# (1) Compiling Raw Data

# (1.1) Data
# Contains all data for running analyses for the paper. Some data was pulled in from previously published work.

# Note: masteroptics_p.csv is the file for optical data based on BEPOM data from the published study 
# Shakil S, Tank S E, Kokelj S V, Vonk J E and Zolkos S 2020 Particulate dominance of organic carbon mobilization from thaw slumps on the Peel Plateau, NT: Quantification and implications for stream systems and permafrost carbon release Environ. Res. Lett. 15 114019

# The Ocean Optics the within the data folder is the data obtained using the Ocean Optics equipment
# The Aqualog folder within the is the data obtained using the Aqualog

# (1.2) Optics 
# Contains code for processing all raw measurement data and/or sampleID data and/or data pulled in from previous publications in the folder data
# Data processing order with R files
# SpectralSlope_Zolkos.R -> DOMAbsorbance.R -> OpticsMaster.R

# (1.3) Field Site, Chemistry, and Watershed (GIS) Data -> Contained within Field_and_GIS_data Folder
# Contains all raw data from chemistry reports and field data collection and the R code to stitch it together for analyses
# Data processing order with R files
# SiteLocations.R -> ParticleSize.R -> Streambedfines.R -> Streamslope.R -> GISData.R -> 2017data.R

# (2) Statistical Analyses and Graphs

# Contains Analyses reported within the papers
# There shouldn't be any dependencies between files, should be able to do each analysis independently based on the data generated previously

# (3) Final Database 
# Database.R - Final data uploaded for the paper
# TBD

# Additional Folders
# Functions -> Any R functions that are needed to be called in
# Archive -> Code created by Sarah near the beginning of the project but doesn't seem to be needed any longer

