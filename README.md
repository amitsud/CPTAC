# CPTAC
Script for processing data from the Clinical Proteomic Tumor Analysis Consortium (CPTAC)

Description

This repository contains an R-based bioinformatics pipeline designed for the integration and analysis of multi-omic data from the Clinical Proteomic Tumor Analysis Consortium (CPTAC). The pipeline processes three primary data modalities—Copy Number Variation (CNV), Transcriptomics (RNA-seq), and Proteomics—across multiple cancer cohorts (including BRCA, ccRCC, COAD, GBM, and others). It performs cross-platform integration by implementing  quality control, covariate regression, and kernel density normalization.


Data access

The code is designed to work with the CPTAC Pancancer harmonized datasets which are available at the [PDC Data Commons](https://proteomic.datacommons.cancer.gov/pdc/cptac-pancancer). 

Download the data and update the base_path variable in the script to point to your local data directory.

This project is available for academic and research use.
