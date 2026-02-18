# Hole et al. (2026) 
Hole, G.M., Büntgen, U., Buchwal, A., Rees, W.G., Wheeler, H.C. Thermal thresholds constrain lignification and frost injury in Low Arctic shrubs of northwestern Canada. [In prep].

Code and data for analysis of taxon-specific thermal limitation of cell wall lignification and frost injury in Low Arctic shrubs.

This repository contains the analysis code and derived datasets required to reproduce the results of Hole et al. (2026). Primary datasets are archived separately (see Data Availability section). 


# Repository Structure
<pre>
shrub_climate_anomaly_response/
│   
└── R/   # All R scripts (analysis workflow)  
│   00_paths_setup.R  
│   01_chronology_formation.R  
│   02_climate_growth_analysis.R  
│   03_BRFR_plots.R  
│   04_dailyclimateplots.R  
│   05_Alnus_GLMM_LWBR.R  
│   06_Salix_GLMM_LWBR.R  
│   07_Alnus_GLMM_LWFR.R  
│   08_Salix_GLMM_LWFR.R  
│   09_Alnus_GLMM_EWFR.R  
│   10_Salix_GLMM_EWFR.R  
│   11_GLMM_output_heatmaps.R  
│   12_LWBR_Tthreshold_Alnus.R  
│   13_LWBR_Tthreshold_Salix.R  
│   utils_packages.R  
│  
├── data/  
    │
    ├── anomalies/               #anatomical anomaly occurrence datasets
    │   alnusBRFRdata.csv  
    │   salixBRFRdata.csv  
    │
    ├── rwi/                     #  RWI datasets
    │   ISR_alnus_RWI.csv
    │   ISR_salix_RWI.csv
    │
    ├── sample_metadata/         # sample metadata including latitude, longitude
    │   ISR_alnus_subset_samples.csv
    │   ISR_salix_subset_samples.csv
    |
    ├── id_crosswalk/             # mapping to UK PDC identifiers
    │   ISR_salix_subset_samples_pdc_key.csv
    │   ISR_alnus_subset_samples_pdc_key.csv
    │
    ├── climate_raw/              #external climate data inputs
    │   CRU4_tmp.dat  
    │   ERA5_prcp_daily.dat  
    │   ERA5_t2m_daily.dat  
    │   ERA5_tmin_daily.dat  
    │
    ├── derived/                  # site-level chronologies derived from sample-level RWI
    │   alspl10SV.csv 
    │   saspl10SV.csv
│  
├── figures/    #output figures (created when scripts run)  
├── output/      #output tables (created when scripts run)  
│  
├── shrub_climate_anomaly_response.Rproj  
└── README.md  
</pre>

# How to Run the Project

1.Open the project

Double-click:   
~~~
shrub_climate_anomaly_response.Rproj
~~~

2.Install and load required packages  
~~~
source(here::here("R/utils_packages.R"))
load_project_packages()
source(here::here("R/00_paths_setup.R"))
~~~

3. Each script is standalone and can be run independently after opening the project. Scripts may be executed in numerical order to reproduce the full workflow:
~~~  
00_paths_setup.R  - project setup  
01_chronology_formation.R  - form growth ring width index (RWI) chronologies
02_climate_growth_analysis.R  - climate/RWI relationships
03_BRFR_plots.R  - create blue ring and frost ring frequency plots
04_dailyclimateplots.R  - create climate plots
05_Alnus_GLMM_LWBR.R  - Alnus latewood blue ring - climate model  
06_Salix_GLMM_LWBR.R  - Salix latewood blue ring - climate model  
07_Alnus_GLMM_LWFR.R  - Alnus latewood frost ring - climate model  
08_Salix_GLMM_LWFR.R  - Salix latewood frost ring - climate model  
09_Alnus_GLMM_EWFR.R  - Alnus earlywood frost ring - climate model  
10_Salix_GLMM_EWFR.R  - frost earlywood frost ring - climate model  
11_GLMM_output_heatmaps.R  - heatmap plot of model outputs
12_LWBR_Tthreshold_Alnus.R - calculate Temp threshold for  Alnus LWBR
13_LWBR_Tthreshold_Salix.R - calculate Temp threshold for Salix LWBR  
~~~  

# Data availability  
The ```data/``` directory contains analysis-ready datasets.
Anatomical anomaly and RWI datasets archived at Zenodo (review link; DOI pending publication)
sample inventory data archived at UK Polar Data Centre https://doi.org/10.5285/b0c6fdb0-2bb5-435c-93e2-e309481ceaf1.

The UK PDC dataset uses unique machine-readable sample identifiers. For analysis, we use shorter alias identifiers (analysis_id) as R packages and plotting workflows are sensitive to long IDs. A crosswalk between analysis_id and the UK PDC pdc_sample_id is provided in ISR_salix_subset_samples_pdc_key.csv and ISR_alnus_subset_samples_pdc_key.csv (archived with the analysis dataset DOI). All derived tables can be joined back to the UK PDC dataset using this crosswalk.
