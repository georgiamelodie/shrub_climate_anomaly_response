# Hole et al. (2026) 
Hole, G.M., Büntgen, U., Buchwal, A., Rees, W.G., Wheeler, H.C. Thermal thresholds constrain lignification and frost injury in Low Arctic shrubs of northwestern Canada. [In prep].

The code in this repository reproduces the results of Hole et al. (2026) using dataset [DOI]. 
This is a subset of data that is available on the UK Polar Data Centre (UK PDC) as part of the dataset “Shrub ring width measurements of Alnus alnobetula and Salix spp. collected from the Inuvialuit Settlement region, Northwest Territories, Canada, 2022-2024”: https://doi.org/10.5285/bd62a79f-473b-4c00-a111-3bbe6bd446fd. Sample IDs can be cross-referenced using files 'ISR_salix_subset_samples_pdc_key.csv','ISR_alnus_subset_samples_pdc_key.csv'.

Code and data for analysis of taxon-specific thermal limitation of cell wall lignification in Low Arctic shrubs.

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
│   alnusBRFRdata.csv  
│   alspl10SV.csv  
│   CRU4_tmp_.dat  
│   ERA5_prcp_daily.dat  
│   ERA5_t2m_daily.dat  
│   ERA5_tmin_daily.dat  
│   ISR_alnus_RWI.csv  
│   ISR_alnus_subset_samples.csv
│   ISR_alnus_subset_sampleS_pdc_key.csv
│   ISR_salix_RWI.csv  
│   ISR_salix_subset_samples.csv
│   ISR_salix_subset_samples_pdc_key.csv
│   salixBRFRdata.csv  
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
The ```Data/``` directory contains analysis-ready datasets.
