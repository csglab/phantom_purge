# Phantom Purge


Paper Title
-------

The optimal purging of phantom molecules in multiplexed droplet-based single-cell RNA-seq data. ([R. Farouni](http://rfarouni.github.io/) and [H. S. Najafabadi](http://csg.lab.mcgill.ca/), 2019)


Summary
-------

This repository hosts code and R Markdown notebooks implementing the statistical modeling approach to estimating the sample index hopping rate and purging phantom molecules in multiplexed droplet-based single-cell RNA-seq data. The preprint detailing the method can be accessed here.

Code
_____

Code to run the phantom purge workflow.

To run the workflow from the command line, just use the following command to generate a reproducible report

```
Rscript -e 'library(rmarkdown); rmarkdown::render("./code/run_workflow.Rmd", "all", output_file=sprintf("./output/notebooks/workflow_%s.nb.html",commandArgs(trailingOnly=T)[1]))' ${dataset_name}
```
where `${dataset_name}` is the name of the folder containing all the *molecule_info.h5* files for all the samples that were multiplexed on the same lane. The files should be renamed by the sample name as such **sample_name.h5** 


 
