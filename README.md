# Applying the Dempster-Shafer Fusion Theory to combine independent land cover products: a case study on the mapping of oil palm plantations in Sumatra, Indonesia

![R](https://img.shields.io/badge/Language-R-blue)
![R Packages](https://img.shields.io/badge/Packages-sf%2C%20terra%2C%20pracma-blue?logo=r)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.envsoft.2023.123456-orange)](Link to the paper)
![Open Source](https://img.shields.io/badge/Open--Source-Yes-brightgreen)

## Table of Contents
- [Context](#Context)
- [Usage](#Usage)
- [Features](#Features)
- [Examples](#Examples)
- [References](#References)

## Context

This project provides a tool to process spatial data and perform classification using Dempster-Shafer Theory (DST). 
Dempster-Shafer Theory (DST) was develop by **A. Dempster in 1967** and later extended by **G. Shafer in 1976**. It generalizes the Bayesian theory of probability to assess the likelihood of an available evidence. This method is well known and widely used with remote sensing data for its management of uncertainty through the functions of belief and plausibility. It is based on a concept of ignorance to provide sound analysis by not interpreting the lack of information as evidence against an hypothesis (**Lein, 2003**). Its application is based on 4 steps :

- **1. The discernment framework (DF)** sets the hypothesis of the fusion process, i.e. it defines all the possible classes potentially assigned to a pixel:
- **2. The Mass Functions Assignment (MFA)** determines the belief level in each input source, i.e. it assigns numeric mass functions to each pixel depending on its original class.
- **3. The Dempster-Shafer fusion rule** combines the mass functions to estimate the conflict (i.e. uncertainty) and the belief (i.e. confidence) for each pixel to belong to each hypothesis.
- **4. The decision rule (DR)** relies on these metrics to make a final decision, i.e. assign a hypothesis to each pixel.

For more information, see the published article of this work :

## Usage

**If you use these codes, please cite our work:**

To use the codes, you can clone the directory

Clone the repository:
   ```bash
   git clone https://github.com/username/project-name.git
```

Attention R packages are required for the use of this code. 

```bash 
install.packages(c("sf", "terra", "tidyverse"))
```

## Features
This work is based on data from the scientific literature: 

Descals 
IIASA
Xu
Mapbiomas
Gaveau 

## Examples

An example of an application is available in the following article: 

## References
- Dempster, A. P. 1967. “Upper and Lower Probabilities Induced by a Multivalued Mapping.” The Annals of Math. Stat. 38 (2): 325–39. https://doi.org/10.1214/aoms/1177698950.
- Lein, J. K. 2003. “Applying Evidential Reasoning Methods to Agricultural Land Cover Classification.” International Journal of Remote Sensing, January. https://doi.org/10.1080/0143116031000095916.
- Shafer, Glenn. 1976. A Mathematical Theory of Evidence. Princeton Univ. Press.

