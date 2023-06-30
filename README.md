# Introduction to Systems Biology
**Hanzehogeschool Groningen: Bioinformatics Project Year 2, Period 8**

Learning about the relations between reality, experiments and models, how to describe biological systems as mathematical models that can be simulated, and how to interpret the results.


## About the project
To first get to know (more) about modelling biological systems a series of four weekly assignments had to be made, with each week expanding on the subjects and becoming more difficult. 
After these four weeks knowledge about the following topics should have been acquired;
- Modelling and simulation.
- Differential equations.
- Dynamics, convergence and equilibrium.
- Setting up and running simulations using R.

After gathering the knowledge from the previous assignments, it was time for the final assignment that consists of writing a scientific article about the process of replicating a chosen scientific research article while also expanding the model used in the chosen research.
The chosen research had to preferably use the `DeSolve` package for the R programming language because this is what had been used in the weekly assignments.


## Endreport
The chosen research article ["Energetic mismatch induced by warming decreases leaf litter decomposition by aquatic detritivores"](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13710),
studies the impact of rising temperatures and changes in Gammarus body size induced by warming on population dynamics and benthic organic matter dynamics in freshwater systems.  
A biological model was created in R to simulate scenarios with different temperatures that helps to get to the desired conclusions, and with that improving the ability to predict the impact of climate change on carbon stocks and ecosystem functions.

For this project the supplied R code from the chosen research article was adapted greatly, improving reproducibility by making the code modular and dynamic instead of terribly hard-coded. 
A scientific article was written in Rmarkdown and LaTeX, about the subject and the adaptations in the code. 


## Repository File Structure
### Project Tree
```bash
BIN2-P8
├── README.md
├── Chosen Research Files
│   ├── jane13710-sup-0001-data s1.docx
│   ├── Journal of Animal Ecology - ... by aquatic detrivores.pdf
│   ├── Réveillon_et_al._(2022)_R_Code.R
│   └── doi_10.5061_dryad.jh9w0vtdj__v2
│       ├── Data_Mismatch.txt
│       └── README_Data.pdf
├── endreport
│   ├── Data_Mismatch.txt
│   ├── end_report.pdf
│   ├── end_report.Rmd
│   ├── report_subfiles
│   │   ├── abbreviations.tex
│   │   ├── abstract.tex
│   │   ├── after_body.Rmd
│   │   ├── before_body.tex
│   │   ├── equations.tex
│   │   ├── ieee.csl
│   │   ├── import.tex
│   │   └── references.bib
│   ├── figures
│   │   ├── gammarus_fossarum.jpeg
│   │   ├── Model Diagram Annotated.png
│   │   ├── Model Diagram.png
│   │   ├── original_plots
│   │   │   └── *
│   │   └── reproduced_plots
│   │       └── *
│   └── src
│       ├── functions.R
│       ├── mainAnalysis.R
│       ├── model.R
│       └── simulateScenario.R
├── week1 mRNA assignment
│   └── *
├── week2 Corticoide assignment
│   └── *
├── week3 Corticoide assignment
│   └── *
└── week4 Feedback assignment
    └── *
```

### / Chosen Research Files
The publicly available files from the chosen research have been downloaded and put into this directory for easy access.
Which are the following:
* A subdirectory with the research dataset and it's supporting readme pdf.
* The research article as a pdf
* R code from the research article
* Supplementary materials document

### / endreport
This is the final assignment where a report is written and the biological model recreated (and expanded) should almost be considered as a separate project in the repository.
To knit the endreport or run the R code the working directory has to be set to the endreport directory.  
Everything is written and subdivided into multiple files, the `end_report.Rmd` is the main file and can be knitted to pdf. 
The endreport article is split into separate sections as files, which are subsequently located in the `report_subfiles` directory. This is so the the main file only contains the actual article text and not other things such as the title page. 
The R code in the `src` directory is separate and modular but has dependencies built in referring to each other. 
`model.R` is the basis of the project, `functions.R` has utility functions uses for the simulations, `simulateScenario.R` has a single function that is used as entry point for simulating scenarios, this function is used by the `mainAnalysis.R` which is the only one that actually executes code and simulates the scenarios

### / weekX assignments
The four weekly assignment directories contain Rmarkdown files, their knitted pdf's and some supplementary figures. 
There is not much to do with the weekly assignments except to read through the pdf's, they do not have anything substantial to reproduce. Unlike the endreport, which is properly reproducible because it has actual code files.


## Installation
*This project was written on MacOS in RStudio (version 2022.02.0) with R version 4.1.3 for Apple silicon arm64.*  

First, a working R environment is needed, which can be installed from [the CRAN website](https://cran.r-project.org/) by carefully following the instructions there.  
Second, either [RStudio](https://www.rstudio.com/products/rstudio/download/), or another editor of choice should be installed.
It should be noted that to be able to knit the documents LaTeX is needed.   
Finally, the required packages listed below should be installed, after which the project should be reproducible.

When running the Rmarkdown files from the weekly assignments the working directory should be the repository directory but when knitting or running files from the endreport the working directory should be the endreport directory itself.


## Required Packages
The following R packages are required for the endreport and should be installed through an R console using the `install.packages()` function.
It is recommended to install these exact versions of the libraries to not get issues because of differences between versions.
- data.table (1.14.2)
- deSolve (1.34)
- ggpubr (0.4.0)
- lme4 (1.1-29)
- quantmod (0.4.20)
- reshape2 (1.4.4)
- tidyverse (1.3.1)

To easily install any missing packages the code below can be used instead, which should be pasted and run an R console:
```r
install.packages("versions")
required_packages <- c("data.table", "deSolve", "ggpubr", "lme4",   "quantmod", "reshape2", "tidyverse")
package_versions  <- c("1.14.2",     "1.34",    "0.4.0",  "1.1-29", "0.4.20",   "1.4.4",    "1.3.1")

missing <- !(required_packages %in% installed.packages()[,"Package"])
if(any(missing)) install.versions(required_packages[missing], package_versions[missing])
```


## Useful links
* [Project Manual](https://bioinf.nl/~fennaf/thema08/)
* [Chosen Research Article (Webpage)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13710)
* [Chosen Research Article (Online PDF)](https://besjournals.onlinelibrary.wiley.com/doi/pdfdirect/10.1111/1365-2656.13710)
* [Chosen Research Dataset](https://datadryad.org/stash/dataset/doi:10.5061/dryad.jh9w0vtdj)
* [Chosen Research R Code](https://zenodo.org/record/6408937)
