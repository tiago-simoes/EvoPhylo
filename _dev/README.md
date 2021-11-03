# MorphoEvol_Tetrapods

R scripts  to replicate all analyses prior to and subsequent to evolutionary tree inference from the paper: 

[Simões, T. R. & Pierce, S. E. 2021. Sustained High Rates of Morphological Evolution During the Rise of Tetrapods. Nature Ecology & Evolution.](https://doi.org/10.1038/s41559-021-01532-x)

All data sets used for this paper (including the ones not relevant to the protocols provided in the R scripts) are available at [Harvard Dataverse](https://doi.org/10.7910/DVN/NNVTTD)

Each R script includes detailed comments to guide any reader interested on replicating the analysis herein, or utilizing our protocols for their own datasets (e.g., guides on performing character partitioning and assessing the strength of selection). If you find any of the protocols below useful and use them for your own studies, please cite the paper referenced above.

## Codes

### 1. Automated character partitioning: CharPartitioning_Tetrapods.R

Customized script for the new protocol presented in this study to automatically detect character partitions based on the Gower distance metric and graphic cluster visualization using t-SNEs. The phylogenetic data matrix can be imported by saving a copy of the dataset as a .csv file (as in the example provided), or using functions available in the package Claddis to import the phylogenetic data matrix into R and subsequently converting it into matrix file.

*Example files: available in Examples > CharPartitioning folder*

### 2. NewFUNCTIONS.R

New functions developed by Joëlle Barido-Sottani (see references therein) or by us to import tree files from BEAST2 and Mr. Bayes with metadata embedded (e.g., evolutionary rates), and to remove "dummy taxa" introduced by the "treeWoffset" operator from BEAST2.

### 3. WTY_Tetrapods.R

Script for the RWTY back to extract diagnostic parameters from all analyses.

*Example files: available in Examples > RWTY folder*

### 4. StatTestsy&Plots_Rates_Tetrapods.R

Script to perform statistical tests, output summary statistic values, and plot graphs related to morphological evolutionary rates obtained from the final analyses.

*Example files: available in Examples > Rates & FBD parameter stats folder*

### 5. StatTests&Plots_FBD_Tetrapods.R

Script to perform statistical tests, output summary statistic values, and plot graphs related to FBD parameters outputs (net diversification, relative extinction, and relative fossilization) obtained from the final analyses.

*Example files: available in Examples > Rates&FBD_Stats folder*

### 6. Detecting selection strength: StatsTests&Plots_SelectionStrength_Tetrapods.R

Customized script for the new protocol presented in this study to assess the strength of natural selection inferred from clock-based evolutionary rates using morphological data (applicable to living and/or extinct species).

*Example files: available in Examples > SelectionStrength folder*

The script for selection strength runs with files already available at the *Rates&FBD_Stats folder* and other files that are produced by the StatTests&Plots_FBD_Tetrapods.R script. Those files are repeated or already provided in the SelectionStrength folder so this script can be run independently of the others and to exemplify the kind of data input necessary for this protocol.
