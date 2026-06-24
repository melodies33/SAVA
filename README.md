# SAVA: Safe, Always-Valid Alpha-Investing Rules For Doubly Sequential Online Inference
This repository contains the code and reproducibility materials for the paper "Safe, Always-Valid Alpha-Investing Rules For Doubly Sequential Online Inference". The repository includes:

- The R package SAVA, which implements the proposed SAVA methodology.
- Scripts for reproducing the simulation studies reported in Appendix B.
- Scripts for reproducing the real-data analysis reported in Section 4.
- The data files and plotting scripts required to reproduce the numerical results and figures in the paper.

## Repository Structure
- SAVA: source files for the SAVA R package.
- main_sava.R: reproduces the simulation results in Appendix B.
- plot.R: generates the corresponding simulation figures.
- real-amazon.R: reproduces the real-data analysis in Section 4.
- realdata_plot.R: generates the corresponding real-data figures.
- AMAZON_All_Beauty.zip, AMAZON_FASHION.zip and AMAZON_Luxury_Beauty.zip: contains the data for real data analysis in Section 4.

## Installation
The package source files are contained in the SAVA directory.

To install the package locally, open a terminal in the repository directory and run

`R CMD INSTALL SAVA`

Alternatively, install directly from GitHub:

```r
library(devtools)
devtools::install_github(
  "melodies33/SAVA",
  subdir = "SAVA"
)
```
## Reproducing the Results
Run the following scripts in order:

1. real-amazon.R
2. main_sava.R

These scripts generate CSV files containing the numerical results.

Then run

1. realdata_plot.R
2. plot.R

to reproduce the figures reported in the paper.

## Amazon Datasets
The Amazon review datasets used in the numerical experiments are available at

https://nijianmo.github.io/amazon/index.html

Under K-cores and ratings-only data, download the ratings-only datasets for:
- Amazon Fashion
- All Beauty
- Luxury Beauty

Rename the downloaded files as

- AMAZON_FASHION
- AMAZON_All_Beauty
- AMAZON_Luxury_Beauty

These datasets are required to reproduce the real-data analysis in Section 4.

For convenience and reproducibility, the exact dataset files used in our experiments are also included in this repository.

