# Gene Family Evolution and Copy Number Variation in Brassicaceae

This repository contains the analysis scripts, selected data, and output files associated with the manuscript:

**"Copy number variation drives gene family evolution and ecological adaptation in Brassicaceae."**

---

## Overview

Copy number variation (CNV) is a major source of genomic diversity, shaping gene family evolution and contributing to adaptive differences. In this study, we investigated CNV across four Brassicaceae species (*Arabidopsis thaliana*, *A. lyrata*, *A. halleri*, and *Arabis alpina*) to identify gene families undergoing rapid expansion or contraction.

Key findings include:

- Identification of **231 rapidly evolving gene families** enriched in defense, stress response, metabolism, and genetic information processing.
- Population-level CNV analyses in *A. thaliana* and *A. lyrata* using **160 long-read assemblies**.
- Evidence for **pseudogenization and functional divergence** among paralogs despite strong purifying selection.
- **Geographic structuring** of CNV patterns, mirroring allele-frequency–based population structure.
- Ecological modeling showing **CNV–climate associations**, particularly with temperature gradients, highlighting CNV as a driver of local adaptation.

---

## Repository Contents

This repository is structured around the main analyses described in the paper:

- **`family_discovery.R`**  
  Identification of rapidly evolving gene families using birth–death models of gene family evolution.

- **`DefStress.R`**  
  Focused analyses on defense- and stress-related gene families, including pseudogenization, functional divergence, and TE-association.

- **`pop_structure.R`**  
  Analyses of genetic divergence across populations of *A. thaliana* and *A. lyrata*, including geographic clustering and comparison with allele-frequency–based structure.

- **`ecology.R`**  
  Ecological modeling linking CNV patterns to climate variables. Highlights associations between CNV and bioclimatic gradients, particularly temperature.

- **`uniprot_query.py`**  
  Script to query UniProt for gene family functional annotations.

- **Data and Outputs**  
  Selected intermediate and final results are provided to ensure reproducibility of key figures and analyses presented in the manuscript.

---

## Citation

If you use the scripts or data in this repository, please cite the associated publication:

> Karunarathne P., et al. (2025). *Copy number variation drives gene family evolution and ecological adaptation in Brassicaceae.*  
> [Journal, Volume, Pages, DOI — to be added upon publication]

---