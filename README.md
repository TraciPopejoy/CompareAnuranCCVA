# Comparison of CCVAs of Southeast US anurans
Traci P. DuBose Updated 01/25/2023

-------------------------------

# Purpose

This repository contains R code for comparing climate change vulnerability assessments of anurans found within the Southeast United States.  

-------------------------------

# Intended Uses

Code is currently in draft form and preliminary.

-------------------------------

# Files


### Folder: data
This folder contains data that is used to parameterize or build scripts within this folder. Most notably is a phylogeny and image for trait imputation, a trait table to designate different activities for anurans (from ATraiU, Moore et al. 2020), and micro climate data for the mechanistic niche models. 

### Folder: results
This folder contains figures and tables produced from the code within this folder

#### 0_SettingTaxSpaceExtents.Rmd
This Rmarkdown identifies what we consider the Southeast United States. It also determines which species have > 30% of their area of occupancy in the region and looks at how the remaining species vary in their ecological traits. There is a PCA of trait space in this document. 

#### 1_TraitsBasedCCVA.Rmd
This Rmarkdown has plots and other analyses that could be used to evaluate climate change vulnerability based on species traits. It is very preliminary and not complete at all. 

#### 2_RCS_CCVA.Rmd
This Rmarkdown contains code to calculate and visualize the results of a Rarity and Climate Sensitivity index calculated at our regional extent. 2_RCS_CCVA.R contains paired down code from this markdown that can be run on the VT high powered computer.

#### 3_SDM_CCVA.Rmd
This Rmarkdown contains a for loop that could be used to evaluate habitat fragmentation, area, and potential climate niche breadth based on species distribution models. The species distribution model files should be within a folder in data (./data/SDMs/) and follow the following naming convention: "[species]_topmodel.rds is the species raw model, [species_pred.asc is the predicted model across the species range, and [species]_bi_pred.asc is the predicted model transformed to either present/absent using a threshold value. 

#### 4_MNM_CCVA.Rmd
This Rmarkdown contains code that brings in data created using code from [AnuranMNMs/4a_MNMAcrossSpaceCode.Rmd](https://github.com/TraciPopejoy/AnuranMNMs/blob/main/4a_MNMAcrossSpaceCode.Rmd). The data reports minimum, 25th quartile, and median daily minimum warming tolerance for species in out species list. Is then plotted to display thermal sensitivity.

#### 5_Compare_CCVAs.Rmd
This Rmarkdown brings in results from Rmds 2 through 4 above and plots their (1) scaled result, (2) difference between the regional RCS and MNMs, and (3) taxonomic differences in CCVA convergence. Figures in the outline document are generally created here.  


