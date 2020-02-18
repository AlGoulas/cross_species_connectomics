# Cross-species, cross-scales connectomics
Perform cross-species and cross-scales analysis of mammalian brain connectomes

![Cross-species, cross-scales connectomics](species_connectomes_fig.png)

# Description
The repository hosts the functions for the analysis (network metrics, statistics) (aux), the data (brain connectomes, cytological data for braina areas, spatial distances between areas and certain meta-data, e.g., species name) (data) and the main function to perform the analysis (main).

The main function performs the analysis as described in:

Goulas A, Majka P, Rosa MGP, Hilgetag CC (2019) A blueprint of mammalian cortical connectomes. PLoS Biol 17(3): e2005346. https://doi.org/10.1371/journal.pbio.2005346

The main function makes use of the Brain Connectivity Toolbox to compute node-wise network metrics.
https://sites.google.com/site/bctnet/

Add the folder with all the subfolders to your MATLAB path. Load the data in your workspace and run the main function (CrossSpecies_NetworkMetrics.m).

# Example

For examining all 7 species/datasets, run the following command:

```
species_index=[1:7];
[Stats, Plotting]=CrossSpecies_NetworkMetrics(Species, species_index, 1, 100, 'core-periphery', 'ks');
```
 
If only dataset e.g., 3 and 5 needs to be analyzed run:

```
[Stats, Plotting]=CrossSpecies_NetworkMetrics(Species, [3 5], 1, 100, 'core-periphery', 'ks');
```
See description of CrossSpecies_NetworkMetrics.m for a complete description of the inputs and the outputs.

# Citations

The following papers should be cited if any functions or data are used:

Goulas A, Majka P, Rosa MGP, Hilgetag CC (2019) A blueprint of mammalian cortical connectomes. PLoS Biol 17(3): e2005346. https://doi.org/10.1371/journal.pbio.2005346

Scannell JW, Young MP. The connectional organization of neural systems in the cat cerebral cortex. Current Biology. 1993;3(4):191–200. https://doi.org/10.1016/0960-9822(93)90331-H

Markov NT, Ercsey-Ravasz MM, Ribeiro Gomes AR, Lamy C, Magrou L, Vezoli J, et al. A Weighted and Directed Interareal Connectivity Matrix for Macaque Cerebral Cortex. Cerebral Cortex. 2014;24(1):17–36. https://doi.org/10.1093/cercor/bhs270

Zingg B, Hintiryan H, Gou L, Song M, Bay M, Bienkowski M, et al. Neural Networks of the Mouse Neocortex. Cell. 2014;156(5):1096–1111. https://doi.org/10.1016/j.cell.2014.02.023

Oh SW, Harris JA, Ng L, Winslow B, Cain N, Mihalas S, et al. A mesoscale connectome of the mouse brain. Nature. 2014;508:207–214. http://dx.doi.org/10.1038/nature13186

Atapour N, Majka P, Wolkowicz IH, Malamanova D, Worthy KH, Rosa MGP. Neuronal Distribution Across the Cerebral Cortex of the Marmoset Monkey (Callithrix jacchus). Cerebral Cortex. 2018; 29(9):3836-3863. https://dx.doi.org/10.1093/cercor/bhy263
   
Goulas A, Uylings HBM, Hilgetag CC. Principles of ipsilateral and contralateral cortico-cortical connectivity in the mouse. Brain Structure and Function. 2017;222(3):1281–1295. https://doi.org/10.1007/s00429-016-1277-y

Majka P, Chaplin TA, Yu HH, Tolpygo A, Mitra PP, Wójcik DK, et al. Towards a comprehensive atlas of cortical connections in a primate brain: Mapping tracer injection studies of the common marmoset into a reference digital template. Journal of Comparative Neurology. 2016;524(11):2161–2181. https://doi.org/10.1002/cne.24023

Beul SF, Barbas H, Hilgetag CC. A Predictive Structural Model of the Primate Connectome. Scientific Reports. 2017;7:43176. https://doi.org/10.1038/srep43176 


