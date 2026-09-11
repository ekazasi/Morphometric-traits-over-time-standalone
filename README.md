# Simulating Morphometric Traits Across Evolutionary Time
This repository contains a computational pipeline and R/Shiny web tool designed to simulate the evolution of 3D bone landmarks across a phylogenetic tree.

## Features
Evolutionary Models: supports **Brownian Motion** (Neutral Drift) and **Ornstein-Uhlenbeck** (Natural Selection).

3D Visualization: Interactive Plotly wireframes for ancestral bone structures, tip configurations and transitional timeinterval shapes.

Statistical Diagnostics: Built-in Mantel tests and PCA morphospace expansion plots.

## Requirements
To run this pipeline, you need R installed along with the following packages:
```geomorph```, ```phytools```, ```mvMORPH```, ```Rphylopars```, ```vegan```, ```viridis```, ```plotly```, ```shiny```, ```ggplot2```, ```dplyr```, ```ggtree```.

In order to be able to reproduce the same results, make sure you acquire the same R version, as well as, tool versions:
- **R version:** 4.5.2 (2025-10-31)
- ```geomorph```: geomorph_4.1.0
- ```mvMORPH```: mvMORPH_1.2.1
- ```phytools```: phytools_2.5-2
- ```Rphylopars```: Rphylopars_0.3.10
- ```vegan```: vegan_2.7-3
- ```viridis```: viridis_0.6.5
- ```plotly```: plotly_4.12.0
- ```shiny```: shiny_1.13.0
- ```ggplot2```: ggplot2_4.0.2
- ```dplyr```: dplyr_1.2.0
- ```ggtree```: ggtree_4.0.4

Other important packages and their versions:
- ```BiocManager```: BiocManager_1.30.27
- ```RRPP```: RRPP_2.1.2

<br>

To ensure reproducibility and match the specific R version (4.5.2) and package dependencies used in this project (above), I recommend using the provided conda environment.

<br>

**-> Prerequisites:**
- Miniconda or Anaconda installed on your system.

<br>

*See below how to create the project's environment*

<br>

-------------------------------------------------------------------------

<br>

**Pipeline:** Represent bone 3D shapes by simulating landmark configurations. 

Stages:
1. Data Acquisition and Pre-processing:
    - Input: A Morphologika (.txt) file containing $k$ landmarks in $p$ directions for $n$ specimens
    - Generalized Procrustes Analysis (GPA) on the data
    - Create the consensus bone structure for the ancestral and target taxa

2. Dimensionality Reduction via PCA
3. Phylogenetic Tree Construction (based on user's evolutionary scenario)
4. Evolutionary Simulation (BM vs. OU)
5. Slicing of the tree in fixed time-intervals and simulating the landmark configurations in the specific timestamps
6. Diagnostic Evaluation and Visualization:
    - Statistical testing (Mandel Tests)
    - Morphospace plots
    - Interactive Tree Explorer
    - Shape vs Time distance correlation plots


## How to Replicate This Work
Follow these steps to run the simulation using the provided example dataset:
1. Clone the repository:
```bash
git clone https://github.com/ekazasi/Morphometric-traits-over-time
```
2. Create the environment
```bash
conda env create -f environment.yml
```
3. Activate the environment
```bash
conda activate morphometrics_pipeline
```

4. Open RStudio and set the working directory to the project folder
    - In the top menu bar, go to Session
    - Select Set Working Directory
    - Click Choose Directory...
    - Navigate to your folder (where you cloned the git repository) and click Open 

5. Run the .R script and access the Shiny web app.

6. Load Data: Within the app, upload the *data/example_data.txt* file to see the primate cranial configuration simulation.

<br>
<br>

Directory's tree structure:
```
project/
├── app.R 
├── data/
│   └── example_data.txt
└── README.md
```
<br>
<br>


*This work was made for the purposes of a project assignment in the BINP29 course, "DNA Sequencing Informatics II" from the Bioinformatics master of Lund University.*
