Welcome to the GitHub repository for the following publication: [MoCHI: neural networks to fit interpretable models and quantify energies, energetic couplings, epistasis and allostery from deep mutational scanning data (Faure AJ & Lehner B, 2024)](https://www.biorxiv.org/content/10.1101/2024.01.21.575681)

Here you'll find an R package with all scripts to reproduce the figures and results from the computational analyses described in the paper.

# Table Of Contents

* **1. [Required Software](#required-software)**
* **2. [Required Data](#required-data)**
* **3. [Installation Instructions](#installation-instructions)**
* **4. [Usage](#usage)**

# Required Software

To run the mochims pipeline you will need the following software and associated packages:

* **[_R_](https://www.r-project.org/)** (Cairo, bio3d, data.table, ggplot2, ggrepel, hexbin, plot3D, reshape2)

# Required Data

Fitness scores, inferred free energy changes and required miscellaneous files should be downloaded from **[here]()** and unzipped in your project directory (see 'base_dir' option) i.e. where output files should be written.

# Installation Instructions

Make sure you have git and conda installed and then run (expected install time <5min):

```
# Install dependencies manually (preferably in a fresh conda environment)
conda install -c conda-forge cairo r-base>4.0.0 r-bio3d r-cairo r-data.table r-devtools r-ggplot2 r-ggrepel r-hexbin r-plot3d r-reshape2 r-roxygen2

# Open an R session and install the mochims R package
devtools::install_github("lehner-lab/mochims")
```

# Usage

The top-level function **mochims()** is the recommended entry point to the pipeline and by default reproduces the figures and results from the computational analyses described in the following publication: [MoCHI: neural networks to fit interpretable models and quantify energies, energetic couplings, epistasis and allostery from deep mutational scanning data (Faure AJ & Lehner B, 2024)](https://www.biorxiv.org/content/10.1101/2024.01.21.575681). See [Required Data](#required-data) for instructions on how to obtain all required data and miscellaneous files before running the pipeline. Expected run time <20min.

```
library(mochims)
mochims(base_dir = "MY_PROJECT_DIRECTORY")
```
