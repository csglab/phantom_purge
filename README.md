# Phantom Purge

Summary
-------

This repository hosts code and R Markdown notebooks implementing the statistical modeling approach to estimating the sample index hopping rate and purging phantom molecules in multiplexed droplet-based single-cell RNA-seq data. 

Paper
-------

The bioRxiv preprint describing the approach can be accessed [here](https://www.biorxiv.org/content/10.1101/617225v4).

The optimal purging of phantom molecules in multiplexed droplet-based single-cell RNA-seq data. ([R. Farouni](http://rfarouni.github.io/) and [H. S. Najafabadi](http://csg.lab.mcgill.ca/), 2019)


Website
---------

The paper's [website](https://csglab.github.io/phantom_purge/index.html) hosts code and R Markdown notebooks implementing the statistical modeling approach.


Validation Data
---------

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3267922.svg)](https://doi.org/10.5281/zenodo.3267922)


The preprocessed read count datatable that was used for experimentally validating the method can be found on Zenodo. 


Installation Requirements
----------
The code was tested using R 3.6.0 running on Ubuntu 18.04.2 LTS and R 3.5.2 running on CentOS Linux 7.

To run the code, the following R packages are required.

- [rhdf5](https://www.bioconductor.org/packages/release/bioc/html/rhdf5.html)
- [DropletUtils](https://www.bioconductor.org/packages/release/bioc/html/DropletUtils.html)
- tidyverse
- matrixStats
- broom
- furrr
- tictoc
- data.table
- cowplot
- scales


Command Line Example
---------

To run the workflow from the command line, just use the following command to generate a reproducible report

```
Rscript -e 'library(rmarkdown); rmarkdown::render("./code/run_workflow.Rmd", "all", output_file=sprintf("./output/notebooks/workflow_%s.nb.html",commandArgs(trailingOnly=T)[1]))' ${dataset_name}
```
where `${dataset_name}` is the name of the folder containing all the *molecule_info.h5* files for all the samples that were multiplexed on the same lane. The files should be renamed by the sample name as such **sample_name.h5** 

Contact
---------

We use GitHub [issues](https://github.com/csglab/phantom_purge/issues) for tracking requests and bugs. Please submit an new issue if you any comment or you would like to report a software bug.
 
