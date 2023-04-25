# Generalized Correlates Analysis Reporting

## Summary

This repository houses modular workflows for statistical analyses of correlates
of risk / protection and automated reporting of analytic results. It serves as
a generalized suite of tools, based on the analyses originally designed for the
USG Biostatistics Response Team's analysis of COVID-19 vaccine efficacy trials
(archived
[here](https://github.com/CoVPN/correlates_reporting_usgcove_archive/)). See
below for brief descriptions of each of the analysis modules. This repository is
designed as the second part of an analytic pipeline, with the [correlates
processing](https://github.com/CoVPN/correlates_processing) module serving as an
upstream component.


## Design notes

The correlates_reporting3 repo is more flexible compared to correlates_reporting2 (design doc). Reporting2 was designed such that each time the “button” was pressed, a correlates analysis was run for one endpoint/analysis population/dataset combination, e.g., correlates for the severe endpoint starting seven days post the D29 visit in the US part of the ENSEMBLE trial. A new analysis can be specific by adding a new TRIAL section and a new COR section in config.yml with minimal coding. Reporting3 is designed to facilitate implementation of more complex objectives, e.g., to assess whether a correlate is modified by naïve/non-naïve status. This is primarily achieved by pushing the logic concerning COR from _common.R into a new top layer within an analysis module, see cor_coxph_cove_boost.R for an example. 

correlates_reporting3 makes some improvements over correlates_reporting2:

* Assay metadata such as LLOQs, ULOQs, llox_label etc are stored in a csv file and used in data mapping, correlates_processing, and correlates_reportings. 

The basic structure hasn’t changed from the correlates_reporting2 repo. 

* config.yml Each section in the file describes one data file. 

* _common.R It acts as like a package and contains functions that implement shared logic such as making bootstrap samples, trichotomizing assay markers, define followup days for risk etc. Keeping these functions in _common.R is more nimble than spinning them off to a package and faster development cycle. The one exception is the cove.boost.collapse.strata function, which is used by both reporting3 and processing repos and is contained in the kyotil package.


## Quick start

All analysis code are written in R. renv is used to manage package versions.

* Guides on forking repo for development can be found in the next section.

* After cloning the forked repo, start R in the root directory. Enter 
```r
renv::restore()
```
to install packages, which takes a few hours, but only has to be performed once.


* To generate tables and figures to include in reports, run the following command in the module directory. For this command to work, there needs to be a Makefile in the module that specifies the action. (See cor_coxph as an example.)
```bash
make cor_coxph
```

* To generate pdf reports for a module, run the following command in the root directory of the repository. For this command to work, there needs to be a report.Rmd in the module directory.
```bash
bash _build_chapter.sh cor_coxph
```

* Troubleshooting: Error in file(filename, "r", encoding = encoding) :  cannot open the connection
  * Check if there is a .Rprofile in the home directory. Remove it if it exists
  


## Collaboration Guide

* Getting started: see our [contribution
   guidelines](https://github.com/CoVPN/correlates_reporting2/blob/master/CONTRIBUTING.md).
   
* [Code style guide](https://style.tidyverse.org/), with some modifications;
  this will largely be enforcd with [`styler`](https://styler.r-lib.org/).

* Project organization: _mostly_ independent subdirectories, each incorporating
  [`here`](https://here.r-lib.org/) for path resolution.

* Package version control and virtual environments using
  [`renv`](https://rstudio.github.io/renv/).


## List of Analysis Modules

* Correlates of Risk (CoR) Analyses
  * `cor_tabular`: Tabular descriptions of correlates of risk.
  * `cor_graphical`: Graphical descriptions of correlates of risk.
  * `cor_coxph`: Cox proportional hazards modeling of risk.
  * `cor_threshold`: Risk modeling based on correlate thresholds.
  * `cor_nonlinear`: Nonlinear modeling and evaluation.
  * `cor_surrogates`: Optimal surrogates analyses.
* Correlates of Protection (CoP) Analyses
  * `cop_prinstrat`: Principal stratification analyses.
  * `cop_stochastic`: Stochastic risk and vaccine efficacy evaluation.
  * `cop_mediation`: Correlate-mediated vaccine efficacy and risk.


---

## License

The contents of this repository are distributed under the GPL-3 license. See
file [`LICENSE.md`](https://github.com/CoVPN/correlates_reporting2/blob/master/LICENSE.md)
for details.
