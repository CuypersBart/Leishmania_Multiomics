# Leishmania Multiomics

This repository contains the two jupyter notebooks that were used to generate the results from the manuscript:
"Four layer multi-omics reveals molecular responses to aneuploidy in Leishmania" by Cuypers Bart et al.
Additionally, the repository contains the R scripts that were used to generate the main plots. 

## Leish_Aneuploidy_Multiomics.ipynb
- Imports and processes annotation, aneuploidy, CNV, SNP, transcriptomics and proteomics data of 6 aneuploid Leishmania strains. 
- Calculates contrasts between all possible pairs of strains
- Outputs transcript and protein fold changes (CNV-containing genes excluded), annotated with genome annotation metadata.

## Multiomics_Correlator.ipynb
- Imports and processes annotation, aneuploidy, CNV, SNP, transcriptomics and proteomics data of 6 aneuploid Leishmania strains
- Normalises the expression of each transcript and protein by its expression in a strain that is disomic for its encoding chromosome
- Averages these values per chromosome

## Data Folder
Contains the preprocessed 'omics data required by both jupyter notebooks. 
