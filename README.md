# LimitingSimilarity_vs_Preadaptation
R sripts for paper "Evidence for environmental filtering and limiting similarity depend on spatial scale and dissimilarity metrics", authored by: Perez-Navarro MA, Shepherd HER, Brian JI, Clark AT and Catford JA. The author of  these scripts is Perez-Navarro MA. The scripts allows to estimate species functional dissimilarity and phylogenetic distances and perform the paper statistical analyses and reproduce tables and figures.

## Getting Started

Before running the scripts, ensure you have R and RStudio installed on your computer. Additionally, some scripts require specific data files and dependencies. Follow the setup instructions below to prepare your environment.

### Prerequisites

- R and RStudio
- Required R packages: (ade4, ape, circlize, dendextend, devtools, DHARMa, dplyr, emmeans, geiger, ggplot, HH, interactions, lme4, lmerTest, MEtBrewer, NSR,  performance, phangorn, phylogram, phytools, radiant.data, RcolorBrewer, readr, request, stringr, Taxonstand, taxize, tibble, tidyr, tidyverse, V.PhyloMaker.)

### Installation

1. Clone this repository to your local machine:
```
git clone https://github.com/mangeles582/LimitingSimilarity_vs_Preadaptation.git
```
2. Navigate into the project directory:
```
cd LimitingSimilarity_vs_Preadaptation
```
3. Install dependencies:
```
install.packages(c("ade4", "ape", "circlize", "dendextend", "devtools", "DHARMa", "dplyr", "emmeans", "geiger", "ggplot2", "HH", "interactions", "lme4", "lmerTest", "MEtBrewer", "NSR", "performance", "phangorn", "phylogram", "phytools", "radiant.data", "RColorBrewer", "readr", "request", "stringr", "Taxonstand", "taxize", "tibble", "tidyr", "tidyverse", "V.PhyloMaker"))

4. Download data

- Plant cover data Experiment E014 https://doi.org/10.6073/pasta/f0a2a36d882d9c937f817b15553e69ee. 
- Environmental data Experiment E014 https://cedarcreek.umn.edu/research/experiments/e014
- Trait data upon request to m.angeles582@gmail.com and from experiment E133 https://doi.org/10.6073/pasta/5814d764417a40dd80e56434b13a10b9
- Species origin status
check the origin of introduced/or native species manually 
https://plants.usda.gov/home/plantProfile?symbol=POPR 
or using NSR -Native status resolver- within BIEN database project
https://bien.nceas.ucsb.edu/bien/tools/nsr/batch-mode/
https://github.com/EnquistLab/RNSR/blob/master/NSR/vignettes/NSR.Rmd

## Script Descriptions

### Script 01_a_multivariate_dissim

- **Purpose**: Estimate functional multivariate dissimilarity at neighbourhood scale (i.e, plot level 0.5m2)
- **Required data**: plant cover and trait data

### Script 01_b_multivariate_dissim_transect

- **Purpose**: Estimate functional multivariate dissimilarity at site scale (i.e, transect level ~40m2)
- **Required data**: plant cover and trait data

### Script 01_c_univariate_dissim

- **Purpose**: Estimate functional dissimilarity for single functional traits at neighbourhood and site scale
- **Required data**: plant cover and trait data

### Script 02_a_phylogenetic_dist

- **Purpose**: Estimate phylogenetic distance for each species within the community at neighbourhood scale (i.e, plot level 0.5m2)
- **Required files**: Plant cover

### Script 02_b_phylogenetic_dist_transect

- **Purpose**: Estimate phylogenetic distance for each species within the community at site scale (i.e, transect level ~40m2) 
- **Required files**: Plant cover

### Script 03_final_dataset_preparation

- **Purpose**: Prepare dataset for statistical analyses. Includes estimation of plant colonization time, interpolation of soil nutrient content when data were not available at plot or transect level and accurate renaming of transect codes among other transformations.
- **Required files**: Plant cover and functional and phylogenetic datasets obtained from code 01 and 02.

### Script 04_analyses_and_figures

- **Purpose**: Run paper final statistical analyses and prepare main text figures.
- **Required files**: Dataset with all plant cover, functional dissimilarities, phylogenetic distances and environmental data obtained from code 03. Final dataset contains data that are not property of the main author so these are not freely available in this repository. Data are available upon request from m.angeles582@gmail.com

### Script 05_supplmentary_figures

- **Purpose**: Run paper final statistical analyses and prepare Appendix S1 figures.
- **Required files**: Dataset with all plant cover, functional dissimilarities, phylogenetic distances and environmental data obtained from code 03. Final dataset contains data that are not property of the main author so these are not freely available in this repository. Data are available upon request from m.angeles582@gmail.com

## Running the Scripts

Each script can be run independently, provided the required files are in place. It is advised to follow the script order as listed for a smooth workflow. 

## License

This work is licensed under a Creative Commons Attribution 1.0 International License (CC BY 4.0). This license allows others to share, copy and redistribute the material in any medium or format, and adapt, remix, transform, and build upon the material for any purpose. 

You can view a copy of this license at [http://creativecommons.org/licenses/by/4.0/](http://creativecommons.org/licenses/by/4.0/).

### How to Cite This Work

If you use this work, or any part of it, in your research or project, please provide appropriate credit to the authors and a link to the original source and the license. Here is a suggested format for citation:

- Evidence for environmental filtering and limiting similarity depend on spatial scale and dissimilarity metrics. 2025. Perez-Navarro MA, Shepherd HER, Brian JI, Clark AT and Catford JA

### Previous version of the paper

https://www.biorxiv.org/content/10.1101/2023.08.18.553820v1.abstract




