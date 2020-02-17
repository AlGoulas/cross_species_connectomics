# cross-species-connectomics
Perform cross-species and cross-scales analysis of mammalian brain connectomes

# description
The repository hosts the functions for the analysis (network metrics, statistics) (aux), the data (brain connectomes, cytological data for braina areas, spatial distances between areas and certain meta-data, e.g., species name) (data) and the main function to perform the analysis (main).

The main function makes use of the Brain Connectivity Toolbox to compute node-wise network metrics.
https://sites.google.com/site/bctnet/

Add the folder with all the subfolders to your MATLAB path. Load the data in you workspace and run the main function (CrossSpecies_NetworkMetrics.m).

For examining all 7 species/datasets, run the following command:

species_index=[1:7];


