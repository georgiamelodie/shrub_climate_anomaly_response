# Hole et al. (2026) 
Hole, G.M., Büntgen, U., Buchwal, A., Rees, W.G., Wheeler, H.C. Blue Rings reveal thermal thresholds of wood formation in Low Arctic tundra shrubs. Scientific Reports [In press].
<a href="https://doi.org/10.5281/zenodo.22285900"><img src="https://zenodo.org/badge/1158478095.svg" alt="DOI"></a>


Code and data for analysis of taxon-specific thermal limitation of cell wall lignification and frost injury in Low Arctic shrubs.

This repository contains the analysis code and derived datasets required to reproduce the results of Hole et al. (2026). Primary datasets are archived separately (see Data Availability section). 


## Repository Structure
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
│   14_LWBR_Tthreshold_taxon_comparison.R
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

## How to Run the Project

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
10_Salix_GLMM_EWFR.R  - Salix earlywood frost ring - climate model  
11_GLMM_output_heatmaps.R  - heatmap plot of model outputs
12_LWBR_Tthreshold_Alnus.R - calculate Temp threshold for  Alnus LWBR
13_LWBR_Tthreshold_Salix.R - calculate Temp threshold for Salix LWBR
14_LWBR_Tthreshold_taxon_comparison.R - calculate uncertainty overlap in taxon Temp threshold
~~~  


## Data availability
The ```data/``` directory contains analysis-ready datasets. 
A complete archive (code and data required to run the analysis) is available via a private Zenodo deposit for peer review. The public DOI will be added upon publication.

## Data sources
### Anatomical anomaly datasets
Anatomical anomaly datasets (BRFR) contain newly generated manually identified blue ring and frost ring occurrences.  

### Ring width datasets
Ring Width Index (RWI) data include both previously archived material and newly processed measurements.  
Partially overlapping RWI data archived at UK Polar Data Centre: https://doi.org/10.5285/b0c6fdb0-2bb5-435c-93e2-e309481ceaf1.

### Derived
RWI site-level chronologies derived from sample RWI data.  

### Sample identifiers  
To ensure compatibility with R modelling workflows, shortened analysis identifiers are used in this repository. Mapping between analysis_ID and pdc_sample_ID is provided in:  
~~~
data/id_crosswalk/
~~~
These crosswalk tables allow full traceability to UK PDC archived identifiers.  
Sample inventory data archived at UK Polar Data Centre:  https://doi.org/10.5285/b0c6fdb0-2bb5-435c-93e2-e309481ceaf1.  


### Climate Data  
These files are provided in data/climate_raw/ to ensure reproducibility of model outputs and temperature threshold analyses.  
Original data sources:  
ERA5: Copernicus Climate Data Store.  
CRU: Climatic Research Unit, University of East Anglia.  

## License  
MIT
