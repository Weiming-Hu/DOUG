# DOUG

### Overview

**Please make sure you have installed [Git Large File Storage](https://git-lfs.github.com/) when you clone this repository to be able to down the data files.**

This is the repository for the Dynamically Optimized Unstructured Grid algorithm. The associated paper is 

Weiming Hu, Guido Cervone, Dynamically Optimized Unstructured Grid (DOUG) for Analog Ensemble of Numerical Weather Predictions Using Evolutionary Algorithms.

The repository also provides the source code for the R package `RAnEn` which is used to generate part of the test data. The original data are too big to be included as a test. Therefore, precomputed data are provided.

### Requirement

Please make sure R is installed with the following packages available:

- sp
- GA
- maps
- gstat
- ncdf4
- rgdal
- fields
- raster
- deldir
- spatstat
- maptools
- stringr
- RColorBrewer
- velox
- RAnEn

To install `RAnEn`, please run the following command in R:

```
install.packages("https://github.com/Weiming-Hu/AnalogsEnsemble/raw/master/RAnalogs/releases/RAnEn_latest.tar.gz", repos = NULL)
```

More information about `RAnEn` can be found [here](https://weiming-hu.github.io/AnalogsEnsemble/).

**Please make sure you have installed [Git Large File Storage](https://git-lfs.github.com/) when you clone this repository to be able to down the data files.**

### Usage

Please source the R script `main.R` to run a test case. Parameters can be set in several files separately:

- different-DOUG.R
- GA_core.R
- GA_real-data.R
- main.R

This `main.R` runs three different algorithms (Pure random, Random + elitism, and DOUG) to compare the performance. The final step of the script create a figure for the optimized unstructured grid.
