#if (exists(".DEF.COMMON")) stop ("_common.R has already been loaded") else .DEF.COMMON=TRUE
library(methods)
library(dplyr)
library(marginalizedRisk)
library(survival)
library(parallel)
library(kyotil)
library(copcor)
library(glue)

if(Sys.getenv("TRIAL")=="") stop("Environmental variable TRIAL not defined!!!!!!!!!!!!!!")
TRIAL=Sys.getenv("TRIAL")

set.seed(98109)

if(!exists("verbose")) verbose=0
if (Sys.getenv("VERBOSE") %in% c("T","TRUE")) verbose=1
if (Sys.getenv("VERBOSE") %in% c("1", "2", "3")) verbose=as.integer(Sys.getenv("VERBOSE"))



###################################################################################################
# disable lower level parallelization in favor of higher level of parallelization

library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1L)
#stopifnot(blas_get_num_procs() == 1L) # Commented this out as it does not work as expected any more!
omp_set_num_threads(1L)
stopifnot(omp_get_max_threads() == 1L)
    
# for mclapply etc
numCores <- unname(ifelse(Sys.info()["sysname"] == "Windows", 1, as.integer(future::availableCores()/2)))

###################################################################################################
# load config, assay metadata, define common variables and labels

config <- config::get(config = Sys.getenv("TRIAL"))
study_name=config$study_name

# created named lists for assay metadata to easier access, e.g. assay_labels_short["bindSpike"]
assay_metadata = read.csv(paste0(dirname(attr(config,"file")),"/",config$assay_metadata))
assays=assay_metadata$assay
assay_labels=assay_metadata$assay_label; names(assay_labels)=assays
assay_labels_short=assay_metadata$assay_label_short; names(assay_labels_short)=assays
llox_labels=assay_metadata$llox_label; names(llox_labels)=assays
lloqs=assay_metadata$lloq; names(lloqs)=assays
lods=assay_metadata$lod; names(lods)=assays
lloxs=ifelse(llox_labels=="lloq", lloqs, lods)


form.s = Surv(EventTimePrimary, EventIndPrimary) ~ 1
form.0 = update (form.s, as.formula(config$covariates))
print(form.0)

# read analysis ready data
dat = read.csv(config$data_cleaned)


###############################################################################
# compute tfinal.tpeak

# returns smaller of the two: 1) time of the last case, 2) last time to have 15 at risk
# use case: 1) dat.ph2 is the ph2 dataset from a case control study
get.tfinal.tpeak.1 = function(dat.ph2, event.ind.col="EventIndPrimary", event.time.col="EventTimePrimary") {
  min(
    max (dat.ph2[[event.time.col]] [dat.ph2[[event.ind.col]]==1]),
    sort(dat.ph2[[event.time.col]], decreasing=T)[15]-1
  )
}

if (TRIAL=="moderna_boost") {
  # compute tfinal.tpeak as the minimum of the four quadrants and no larger than 105 days
  tfinal.tpeaks=c(
    get.tfinal.tpeak.1(subset(dat, Trt==1 &  naive & ph2.BD29), event.ind.col="EventIndOmicronBD29", event.time.col="EventTimeOmicronBD29"),
    get.tfinal.tpeak.1(subset(dat, Trt==1 &  naive & ph2.BD29), event.ind.col="EventIndOmicronBD29", event.time.col="EventTimeOmicronBD29"),
    get.tfinal.tpeak.1(subset(dat, Trt==1 & !naive & ph2.BD29), event.ind.col="EventIndOmicronBD29", event.time.col="EventTimeOmicronBD29"),
    get.tfinal.tpeak.1(subset(dat, Trt==1 & !naive & ph2.BD29), event.ind.col="EventIndOmicronBD29", event.time.col="EventTimeOmicronBD29"),
    105)
  myprint(tfinal.tpeaks)
  tfinal.tpeak=min(tfinal.tpeaks)
  myprint(tfinal.tpeak)
}







###############################################################################
# theme options

# fixed knitr chunk options
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  out.width = "80%",
  out.extra = "",
  fig.pos = "H",
  fig.show = "hold",
  fig.align = "center",
  fig.width = 6,
  fig.asp = 0.618,
  fig.retina = 0.8,
  dpi = 600,
  echo = FALSE,
  message = FALSE,
  warning = FALSE
)

# global options
options(
  digits = 6,
  #scipen = 999,
  dplyr.print_min = 6,
  dplyr.print_max = 6,
  crayon.enabled = FALSE,
  bookdown.clean_book = TRUE,
  knitr.kable.NA = "NA",
  repos = structure(c(CRAN = "https://cran.rstudio.com/"))
)

# no complaints from installation warnings
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

# overwrite options by output type
if (knitr:::is_html_output()) {
  #options(width = 80)
  
  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}
if (knitr:::is_latex_output()) {
  #knitr::opts_chunk$set(width = 67)
  #options(width = 67)
  options(cli.unicode = TRUE)
  
  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}

# create and set global ggplot theme
# borrowed from https://github.com/tidymodels/TMwR/blob/master/_common.R
theme_transparent <- function(...) {
  # use black-white theme as base
  ret <- ggplot2::theme_bw(...)
  
  # modify with transparencies
  trans_rect <- ggplot2::element_rect(fill = "transparent", colour = NA)
  ret$panel.background  <- trans_rect
  ret$plot.background   <- trans_rect
  ret$legend.background <- trans_rect
  ret$legend.key        <- trans_rect
  
  # always have legend below
  ret$legend.position <- "bottom"
  return(ret)
}

library(ggplot2)
theme_set(theme_transparent())
theme_update(
  text = element_text(size = 25),
  axis.text.x = element_text(colour = "black", size = 30),
  axis.text.y = element_text(colour = "black", size = 30)
)

# custom ggsave function with updated defaults
ggsave_custom <- function(filename = default_name(plot),
                          height= 15, width = 21, ...) {
  ggsave(filename = filename, height = height, width = width, ...)
}



###############################################################################

# race labeling
labels.race <- c(
  "White", 
  "Black or African American",
  "Asian", 
  if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & startsWith(attr(config, "config"),"janssen_la")) "Indigenous South American" else "American Indian or Alaska Native",
  "Native Hawaiian or Other Pacific Islander", 
  "Multiracial",
  if ((study_name=="COVE" | study_name=="MockCOVE")) "Other", 
  "Not reported and unknown"
)

# ethnicity labeling
labels.ethnicity <- c(
  "Hispanic or Latino", "Not Hispanic or Latino",
  "Not reported and unknown"
)

# baseline stratum labeling
if (study_name=="COVEBoost") {
  demo.stratum.labels <- c(
    "Age >= 65, URM",
    "Age < 65, At risk, URM",
    "Age < 65, Not at risk, URM",
    "Age >= 65, White non-Hisp",
    "Age < 65, At risk, White non-Hisp",
    "Age < 65, Not at risk, White non-Hisp"
  )
} else stop("unknown study_name 1")

