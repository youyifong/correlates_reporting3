# Sys.setenv(TRIAL = "moderna_boost")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

# load required libraries
suppressPackageStartupMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(tidyr, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(purrr, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(stringr, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(quadprog, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(here, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(methods, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(SuperLearner, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(e1071, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(glmnet, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(kyotil, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(argparse, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(vimp, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(nloptr, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(RhpcBLASctl, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(reticulate, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(FSDAM, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(ranger, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(xgboost, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(conflicted, warn.conflicts = FALSE))
suppressMessages(conflicted::conflict_prefer("filter", "dplyr"))
suppressMessages(conflict_prefer("summarise", "dplyr"))

source(here::here("..", "_common.R"))
#-----------------------------------------------

# This script does all of the setup for variable importance analyses.
# It is shared across run_cv_sl_varsets.R (which runs the Super Learners)
# and get_vimp.R (which compiles the estimates into variable importance)
# This file contains all of the study-specific code:
#   setting up variable sets of interest;
#   setting up the weight vector (if needed);
# After running this file, all objects have standard names and can be passed to
# the appropriate downstream file.

# Define code version to run
# the demo version is simpler and runs faster!
# the production version runs SL with a diverse set of learners
run_prod <- !grepl("Mock", study_name)
#run_prod <- FALSE

# load required files and functions
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs

# SL optimal surrogate analysis: Superlearner code requires computing environment with more than 10 cores!
num_cores <- parallel::detectCores()
# if(num_cores < 10) stop("Number of cores on this computing environment are less than 10! Superlearner code needs atleast 11 cores to run smoothly.")

# Study-specific setup ---------------------------------------------------------
# Each dataset requires the following:
#   1. any preprocessing (e.g., adding BAMA data)
#   2. define baseline risk factors into object briskfactors, and define a correction
#      in object briskfactors_correction
#   3. define marker variables of interest into object markerVars
#   4. identify the endpoint of interest (call it "endpoint")
#   5. identify the weights vector (call it "wt")
#   6. define dataset dat.ph1 with all phase 1 data (all columns measured on all study participants)
#   7. define dataset dat.ph2 with phase 2 data (all study participants with phase 2 variables measured)

# Read in data from ENSEMBLE: TRIAL = "janssen_pooled_partA"
if (study_name %in% c("ENSEMBLE")) {
  
  # Create combined new dataset which has imputed values of demographics (for phase 1 data) from dat.covar.imp AND
  # imputed values for markers (for phase 2 data) from dat.wide.v
  dat.ph1 <- dat.mock %>%
    filter(Perprotocol == 1 & Trt == 1 & Bserostatus == 0) %>%
    mutate(Day29pseudoneutid50minusbindRBD = Day29pseudoneutid50 - Day29bindRBD,
           Day29ADCPminusbindRBD = Day29ADCP - Day29bindRBD,
           Delta29overBbindSpike_2fold = ifelse(Day29bindSpike > (BbindSpike + log10(2)), 1, 0),
           Delta29overBbindSpike_4fold = ifelse(Day29bindSpike > (BbindSpike + log10(4)), 1, 0),
           Delta29overBbindRBD_2fold = ifelse(Day29bindRBD > (BbindRBD  + log10(2)), 1, 0),
           Delta29overBbindRBD_4fold = ifelse(Day29bindRBD > (BbindRBD  + log10(4)), 1, 0),
           Delta29overBpseudoneutid50_2fold = ifelse(Day29pseudoneutid50 > (Bpseudoneutid50 + log10(2)), 1, 0),
           Delta29overBpseudoneutid50_4fold = ifelse(Day29pseudoneutid50 > (Bpseudoneutid50 + log10(4)), 1, 0),
           Delta29overBADCP_2fold = ifelse(Day29ADCP > (BADCP + log10(2)), 1, 0),
           Delta29overBADCP_4fold = ifelse(Day29ADCP > (BADCP + log10(4)), 1, 0))
  
  if(Sys.getenv("TRIAL") == "janssen_pooled_partA"){
    # baseline risk factors
    briskfactors <- c("risk_score", "HighRiskInd", "LatinAmerica", "SouthAfrica", "Columbia", 
                      "EnrollPeriod2", "EnrollPeriod3", "EnrollPeriod4", 
                      "EnrollPeriod5", "EnrollPeriod6", "EnrollPeriod7")
    # briskfactors_correction <- "Y ~ x + X$risk_score + X$HighRiskInd + X$LatinAmerica + X$SouthAfrica + X$Columbia + X$EnrollPeriod2 + X$EnrollPeriod3 + X$EnrollPeriod4 + X$EnrollPeriod5 + X$EnrollPeriod6 + X$EnrollPeriod7"
    
    dat.ph1 <- dat.ph1 %>% 
      mutate(
        # Code Region and Columbia as binary indicator variables
        LatinAmerica = ifelse(Region == 1 & Country != 4, 1, 0),
        SouthAfrica = ifelse(Region == 2, 1, 0),
        Columbia = ifelse(Country == 4, 1, 0),
        # Code Enroll Period as binary indicator variables
        EnrollPeriod2 = ifelse(BinnedEnrollmentBiweekly == 1, 1, 0),
        EnrollPeriod3 = ifelse(BinnedEnrollmentBiweekly == 2, 1, 0),
        EnrollPeriod4 = ifelse(BinnedEnrollmentBiweekly == 3, 1, 0),
        EnrollPeriod5 = ifelse(BinnedEnrollmentBiweekly == 4, 1, 0),
        EnrollPeriod6 = ifelse(BinnedEnrollmentBiweekly == 5, 1, 0),
        EnrollPeriod7 = ifelse(BinnedEnrollmentBiweekly == 6, 1, 0))
    }else if(Sys.getenv("TRIAL") == "janssen_la_partA"){
      # baseline risk factors
      briskfactors <- c("risk_score", "HighRiskInd", "Columbia", 
                        "EnrollPeriod2", "EnrollPeriod3", "EnrollPeriod4", 
                        "EnrollPeriod5", "EnrollPeriod6", "EnrollPeriod7")

      dat.ph1 <- dat.ph1 %>% 
        mutate(
          # Code Columbia as binary indicator variables
          Columbia = ifelse(Country == 4, 1, 0),
          # Code Enroll Period as binary indicator variables
          EnrollPeriod2 = ifelse(BinnedEnrollmentBiweekly == 1, 1, 0),
          EnrollPeriod3 = ifelse(BinnedEnrollmentBiweekly == 2, 1, 0),
          EnrollPeriod4 = ifelse(BinnedEnrollmentBiweekly == 3, 1, 0),
          EnrollPeriod5 = ifelse(BinnedEnrollmentBiweekly == 4, 1, 0),
          EnrollPeriod6 = ifelse(BinnedEnrollmentBiweekly == 5, 1, 0),
          EnrollPeriod7 = ifelse(BinnedEnrollmentBiweekly == 6, 1, 0))
    }
  
  individualMarkers <- c("Day29bindSpike",
                         "Day29bindRBD",
                         "Day29pseudoneutid50",
                         "Day29ADCP",
                         "Day29pseudoneutid50minusbindRBD",
                         "Day29ADCPminusbindRBD")
  
  # markers of interest
  markerVars <- c("Day29bindSpike", #"Delta29overBbindSpike",
                  "Delta29overBbindSpike_2fold", "Delta29overBbindSpike_4fold",
                  "Day29bindRBD", #"Delta29overBbindRBD",
                  "Delta29overBbindRBD_2fold", "Delta29overBbindRBD_4fold",
                  "Day29pseudoneutid50", #"Delta29overBpseudoneutid50",
                  "Delta29overBpseudoneutid50_2fold", "Delta29overBpseudoneutid50_4fold",
                  "Day29ADCP", #"Delta29overBADCP",
                  "Delta29overBADCP_2fold", "Delta29overBADCP_4fold",
                  "Day29pseudoneutid50minusbindRBD",
                  "Day29ADCPminusbindRBD")
  
  # Identify the endpoint variable
  endpoint <- "EventIndPrimaryIncludeNotMolecConfirmedD29"
  wt <- "wt.D29"
  ptidvar <- "Ptid"
  
  # Drop any observation with NA values in Ptid, Trt, briskfactors, endpoint and wt.D57
  dat.ph1 <- dat.ph1 %>%
    drop_na(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(wt)) %>%
    arrange(desc(get(endpoint)))
  
  # phase two data (all variables measured on all participants in phase 2 data)
  dat.ph2_init <- dat.ph1 %>%
    filter(TwophasesampIndD29 == TRUE) %>%
    select(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(wt), any_of(markerVars)) %>%
    arrange(desc(get(endpoint)))
}

# Read in data from either COVE or MockCOVE studies
if (study_name %in% c("COVE", "MockCOVE")) {
  # baseline risk factors
  briskfactors <- c("risk_score", "HighRiskInd", "MinorityInd")
  briskfactors_correction <- "Y ~ x + X$risk_score + X$HighRiskInd + X$MinorityInd"

  individualMarkers <- c("Day57bindSpike",
                         "Day57bindRBD",
                         "Day57pseudoneutid50",
                         "Day57pseudoneutid80",
                         "Day57liveneutmn50",
                         "Day29bindSpike",
                         "Day29bindRBD",
                         "Day29pseudoneutid50",
                         "Day29pseudoneutid80",
                         "Day29liveneutmn50")

  # markers of interest
  markerVars <- c("Day57bindSpike", #"Delta57overBbindSpike",
                  "Delta57overBbindSpike_2fold", "Delta57overBbindSpike_4fold",
                  "Day57bindRBD", #"Delta57overBbindRBD",
                  "Delta57overBbindRBD_2fold", "Delta57overBbindRBD_4fold",
                  "Day57pseudoneutid50", #"Delta57overBpseudoneutid50",
                  "Delta57overBpseudoneutid50_2fold", "Delta57overBpseudoneutid50_4fold",
                  "Day57pseudoneutid80", #"Delta57overBpseudoneutid80",
                  "Delta57overBpseudoneutid80_2fold", "Delta57overBpseudoneutid80_4fold",
                  "Day57liveneutmn50", #"Delta57overBliveneutmn50",
                  "Delta57overBliveneutmn50_2fold", "Delta57overBliveneutmn50_4fold",

                  "Day29bindSpike", #"Delta29overBbindSpike",
                  "Delta29overBbindSpike_2fold", "Delta29overBbindSpike_4fold",
                  "Day29bindRBD", #"Delta29overBbindRBD",
                  "Delta29overBbindRBD_2fold", "Delta29overBbindRBD_4fold",
                  "Day29pseudoneutid50", #"Delta29overBpseudoneutid50",
                  "Delta29overBpseudoneutid50_2fold", "Delta29overBpseudoneutid50_4fold",
                  "Day29pseudoneutid80", #"Delta29overBpseudoneutid80",
                  "Delta29overBpseudoneutid80_2fold", "Delta29overBpseudoneutid80_4fold",
                  "Day29liveneutmn50", #"Delta29overBliveneutmn50",
                  "Delta29overBliveneutmn50_2fold", "Delta29overBliveneutmn50_4fold")

  # Identify the endpoint variable
  endpoint <- "EventIndPrimaryD57"
  wt <- "wt.D57"
  ptidvar <- "Ptid"

  # Create combined new dataset which has imputed values of demographics (for phase 1 data) from dat.covar.imp AND
  # imputed values for markers (for phase 2 data) from dat.wide.v
  dat.ph1 <- dat.mock %>%
    filter(Perprotocol == 1 & Trt == 1 & Bserostatus == 0) %>%
    mutate(Delta57overBbindSpike_2fold = ifelse(Day57bindSpike > (BbindSpike + log10(2)), 1, 0),
           Delta57overBbindSpike_4fold = ifelse(Day57bindSpike > (BbindSpike + log10(4)), 1, 0),
           Delta57overBbindRBD_2fold = ifelse(Day57bindRBD > (BbindRBD  + log10(2)), 1, 0),
           Delta57overBbindRBD_4fold = ifelse(Day57bindRBD > (BbindRBD  + log10(4)), 1, 0),
           Delta57overBpseudoneutid50_2fold = ifelse(Day57pseudoneutid50 > (Bpseudoneutid50  + log10(2)), 1, 0),
           Delta57overBpseudoneutid50_4fold = ifelse(Day57pseudoneutid50 > (Bpseudoneutid50  + log10(4)), 1, 0),
           Delta57overBpseudoneutid80_2fold = ifelse(Day57pseudoneutid80 > (Bpseudoneutid80  + log10(2)), 1, 0),
           Delta57overBpseudoneutid80_4fold = ifelse(Day57pseudoneutid80 > (Bpseudoneutid80  + log10(4)), 1, 0),
           Delta57overBliveneutmn50_2fold = ifelse(Day57liveneutmn50 > (Bliveneutmn50  + log10(2)), 1, 0),
           Delta57overBliveneutmn50_4fold = ifelse(Day57liveneutmn50 > (Bliveneutmn50  + log10(4)), 1, 0),

           Delta29overBbindSpike_2fold = ifelse(Day29bindSpike > (BbindSpike + log10(2)), 1, 0),
           Delta29overBbindSpike_4fold = ifelse(Day29bindSpike > (BbindSpike + log10(4)), 1, 0),
           Delta29overBbindRBD_2fold = ifelse(Day29bindRBD > (BbindRBD  + log10(2)), 1, 0),
           Delta29overBbindRBD_4fold = ifelse(Day29bindRBD > (BbindRBD  + log10(4)), 1, 0),
           Delta29overBpseudoneutid50_2fold = ifelse(Day29pseudoneutid50 > (Bpseudoneutid50  + log10(2)), 1, 0),
           Delta29overBpseudoneutid50_4fold = ifelse(Day29pseudoneutid50 > (Bpseudoneutid50  + log10(4)), 1, 0),
           Delta29overBpseudoneutid80_2fold = ifelse(Day29pseudoneutid80 > (Bpseudoneutid80  + log10(2)), 1, 0),
           Delta29overBpseudoneutid80_4fold = ifelse(Day29pseudoneutid80 > (Bpseudoneutid80  + log10(4)), 1, 0),
           Delta29overBliveneutmn50_2fold = ifelse(Day29liveneutmn50 > (Bliveneutmn50  + log10(2)), 1, 0),
           Delta29overBliveneutmn50_4fold = ifelse(Day29liveneutmn50 > (Bliveneutmn50  + log10(4)), 1, 0)) %>%
    # Drop any observation with NA values in Ptid, Trt, briskfactors, endpoint and wt.D57
    drop_na(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(wt)) %>%
    arrange(desc(get(endpoint)))

  # phase two data (all variables measured on all participants in phase 2 data)
  dat.ph2_init <- dat.ph1 %>%
    filter(TwophasesampIndD57 == TRUE) %>%
    select(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(wt), any_of(markerVars)) %>%
    #drop_na(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80) %>%
    arrange(desc(get(endpoint)))
}


# Read in data from COVEBoost study
if (study_name == "COVEBoost") {
  
  # baseline risk factors
  briskfactors <- c("risk_score", "HighRiskInd", "MinorityInd")
  briskfactors_correction <- "Y ~ x + X$risk_score + X$HighRiskInd + X$MinorityInd"
  SLrunBaseline = "briskfactors"
  SLcohort = "naive"
  
  # # baseline risk factors
  # briskfactors <- c("risk_score", "HighRiskInd", "MinorityInd", "BeforeMedianUnblind",
  #                   "CalendarBD1Interval2", "CalendarBD1Interval3", "CalendarBD1Interval4")
  # briskfactors_correction <- "Y ~ x + X$risk_score + X$HighRiskInd + X$MinorityInd + X$BeforeMedianUnblind +
  #                                     X$CalendarBD1Interval2 + X$CalendarBD1Interval3 + X$CalendarBD1Interval4"
  # SLrunBaseline = "briskfactors"
  # SLcohort = "naive"
  
  # # baseline risk factors + Trt
  # briskfactors <- c("risk_score", "HighRiskInd", "MinorityInd", "BeforeMedianUnblind",
  #                       "CalendarBD1Interval2", "CalendarBD1Interval3", "CalendarBD1Interval4", "Trt")
  # briskfactors_correction <- "Y ~ x + X$risk_score + X$HighRiskInd + X$MinorityInd + X$BeforeMedianUnblind +
  #                       X$CalendarBD1Interval2 + X$CalendarBD1Interval3 + X$CalendarBD1Interval4 + X$Trt"
  # SLrunBaseline = "briskfactors_Trt"
  # SLcohort = "naive"
  # 
  # # baseline risk factors + Trt + Days between Dose 2 and Dose 3
  # briskfactors <- c("risk_score", "HighRiskInd", "MinorityInd", "BeforeMedianUnblind",
  #                               "CalendarBD1Interval2", "CalendarBD1Interval3", "CalendarBD1Interval4", "Trt", "NumberdaysDS2toDS3")
  # briskfactors_correction <- "Y ~ x + X$risk_score + X$HighRiskInd + X$MinorityInd + X$BeforeMedianUnblind +
  #                               X$CalendarBD1Interval2 + X$CalendarBD1Interval3 + X$CalendarBD1Interval4 + X$Trt + X$NumberdaysDS2toDS3"
  # SLrunBaseline = "briskfactors_Trt_DaysDS2toDS3"
  # SLcohort = "naive"
  
  individualMarkers <- c("BD1bindSpike_BA.1", "BD29bindSpike_BA.1", "DeltaBD29overBD1bindSpike_BA.1",
                         "BD1pseudoneutid50_BA.1", "BD29pseudoneutid50_BA.1", "DeltaBD29overBD1pseudoneutid50_BA.1",
                         
                         "BD1bindSpike_D614G", "BD29bindSpike_D614G", "DeltaBD29overBD1bindSpike_D614G",
                         "BD1bindRBD_D614G", "BD29bindRBD_D614G", "DeltaBD29overBD1bindRBD_D614G",
                         "BD1pseudoneutid50_D614G", "BD29pseudoneutid50_D614G", "DeltaBD29overBD1pseudoneutid50_D614G")
  
  # Markers present in data
  # pseudoneutid50_D614G, pseudoneutid50_BA.1, bindRBD_D614G, bindSpike_D614G, 
  # bindSpike_Delta, bindSpike_Alpha, bindSpike_Beta, bindSpike_BA.1, bindSpike_Gamma
  
  # markers of interest
  markerVars <- c("BD1bindSpike_BA.1",
                  "BD29bindSpike_BA.1",
                  "DeltaBD29overBD1bindSpike_BA.1",
                  "DeltaBD29overBD1bindSpike_BA.1_2fold",
                  "DeltaBD29overBD1bindSpike_BA.1_4fold",
                  
                  "BD1pseudoneutid50_BA.1",             
                  "BD29pseudoneutid50_BA.1",
                  "DeltaBD29overBD1pseudoneutid50_BA.1",
                  "DeltaBD29overBD1pseudoneutid50_BA.1_2fold",
                  "DeltaBD29overBD1pseudoneutid50_BA.1_4fold",
                  
                  "BD1bindSpike_D614G",
                  "BD29bindSpike_D614G",
                  "DeltaBD29overBD1bindSpike_D614G",
                  "DeltaBD29overBD1bindSpike_D614G_2fold",
                  "DeltaBD29overBD1bindSpike_D614G_4fold",
                  
                  "BD1bindRBD_D614G",
                  "BD29bindRBD_D614G",
                  "DeltaBD29overBD1bindRBD_D614G",
                  "DeltaBD29overBD1bindRBD_D614G_2fold",
                  "DeltaBD29overBD1bindRBD_D614G_4fold",
                  
                  "BD1pseudoneutid50_D614G",             
                  "BD29pseudoneutid50_D614G",
                  "DeltaBD29overBD1pseudoneutid50_D614G",
                  "DeltaBD29overBD1pseudoneutid50_D614G_2fold",
                  "DeltaBD29overBD1pseudoneutid50_D614G_4fold")
  
  # Identify the endpoint variable
  endpoint <- "EventIndPrimaryOmicronBD29"
  wt <- "wt.BD29"
  ptidvar <- "Ptid"
  
  # Create combined new dataset which has imputed values of demographics (for phase 1 data) from dat.covar.imp AND
  # imputed values for markers (for phase 2 data) from dat.wide.v
  dat.ph1 <- dat %>%
    mutate(pooled = ifelse(naive == 1 | nnaive == 1, 1, 0)) %>%
    rename(BD1bindSpike_D614G = BD1bindSpike,
           BD29bindSpike_D614G = BD29bindSpike,
           DeltaBD29overBD1bindSpike_D614G = DeltaBD29overBD1bindSpike,
           BD1bindRBD_D614G = BD1bindRBD, 
           BD29bindRBD_D614G = BD29bindRBD,
           DeltaBD29overBD1bindRBD_D614G = DeltaBD29overBD1bindRBD,
           BD1pseudoneutid50_D614G = BD1pseudoneutid50,
           BD29pseudoneutid50_D614G = BD29pseudoneutid50,
           DeltaBD29overBD1pseudoneutid50_D614G = DeltaBD29overBD1pseudoneutid50) %>%
    filter(ph1.BD29 == TRUE, get(SLcohort) == 1) %>%  
    mutate(#CalendarBD1Interval2 = ifelse(CalendarBD1Interval == 2, 1, 0),
           #CalendarBD1Interval3 = ifelse(CalendarBD1Interval == 3, 1, 0),
           #CalendarBD1Interval4 = ifelse(CalendarBD1Interval == 4, 1, 0), 
           DeltaBD29overBD1bindSpike_BA.1_2fold = ifelse(BD29bindSpike_BA.1 > (BD1bindSpike_BA.1 + log10(2)), 1, 0),
           DeltaBD29overBD1bindSpike_BA.1_4fold = ifelse(BD29bindSpike_BA.1 > (BD1bindSpike_BA.1 + log10(4)), 1, 0),
           DeltaBD29overBD1pseudoneutid50_BA.1_2fold = ifelse(BD29pseudoneutid50_BA.1 > (BD1pseudoneutid50_BA.1  + log10(2)), 1, 0),
           DeltaBD29overBD1pseudoneutid50_BA.1_4fold = ifelse(BD29pseudoneutid50_BA.1 > (BD1pseudoneutid50_BA.1  + log10(4)), 1, 0),
           
           DeltaBD29overBD1bindSpike_D614G_2fold = ifelse(BD29bindSpike_D614G > (BD1bindSpike_D614G + log10(2)), 1, 0),
           DeltaBD29overBD1bindSpike_D614G_4fold = ifelse(BD29bindSpike_D614G > (BD1bindSpike_D614G + log10(4)), 1, 0),
           DeltaBD29overBD1bindRBD_D614G_2fold = ifelse(BD29bindRBD_D614G > (BD1bindRBD_D614G + log10(2)), 1, 0),
           DeltaBD29overBD1bindRBD_D614G_4fold = ifelse(BD29bindRBD_D614G > (BD1bindRBD_D614G + log10(4)), 1, 0),
           DeltaBD29overBD1pseudoneutid50_D614G_2fold = ifelse(BD29pseudoneutid50_D614G > (BD1pseudoneutid50_D614G  + log10(2)), 1, 0),
           DeltaBD29overBD1pseudoneutid50_D614G_4fold = ifelse(BD29pseudoneutid50_D614G > (BD1pseudoneutid50_D614G  + log10(4)), 1, 0)
           ) %>%
    # Drop any observation with NA values in Ptid, Trt, briskfactors, endpoint and wt.D57
    drop_na(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(wt)) %>%
    arrange(desc(get(endpoint)))
  
  # phase two data (all variables measured on all participants in phase 2 data)
  dat.ph2_init <- dat.ph1 %>%
    filter(ph2.BD29 == TRUE) %>%
    select(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(wt), any_of(markerVars)) %>%
    #drop_na(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80) %>%
    arrange(desc(get(endpoint)))
}


# Check if output directory exists; if not, create it!
if(study_name == "COVEBoost"){
  if(!dir.exists(paste0("output/", Sys.getenv("TRIAL"), "/", SLcohort, "/", SLrunBaseline))){
    if(!dir.exists(paste0("output/", Sys.getenv("TRIAL"), "/", SLcohort))){
      dir.create(paste0("output/", Sys.getenv("TRIAL"), "/", SLcohort))
    }
    dir.create(paste0("output/", Sys.getenv("TRIAL"), "/", SLcohort, "/", SLrunBaseline))
  }
} else {
  if(!dir.exists(paste0("output/", Sys.getenv("TRIAL")))){
    dir.create(paste0("output/", Sys.getenv("TRIAL")))
  }
}


# Read in data from HVTN705
if (study_name == "HVTN705") {
  # Add BAMA antigen data to the original dataset
  dat.mock <- dat.mock %>%
    select(-matches("Day1")) 
  # baseline risk factors
  briskfactors <- c("RSA", "Age", "BMI", "Riskscore")
  briskfactors_correction <- "Y ~ x + X$Riskscore + X$Age + X$BMI + X$RSA"
  # markers of interest
  markerVars <- c("Day210ELCZ", "Day210ELMo", "Day210ELISpotPTEEnv",                         
                  "Day210ADCPgp140C97ZAfib","Day210ADCPgp140Mos1fib","Day210IgG3gp140C97ZAfibritin40delta","Day210IgG3gp140Mos1fibritin40delta",          
                  "Day210IgG340mdw_gp120","Day210IgG340mdw_gp140","Day210IgG340mdw_V1V2","Day210IgG3gp4140delta",                       
                  "Day210IgG340mdw_multi","Day210IgG340mdw_gp120_gp140_vm","Day210IgG50mdw_V1V2","Day210mdw_xassay",                            
                  "Day210ADCCCAP8_pk","Day210ADCCCH58_pk","Day210ADCCWITO_pk","Day210ADCCCAP8_pAUC",                         
                  "Day210ADCCCH58_pAUC","Day210ADCCWITO_pAUC","Day210ICS4AnyEnvIFNg_OR_IL2","Day210ICS8AnyEnvIFNg_OR_IL2",                 
                  "Day210mdw_xassay_overall","Day210IgG3AE.A244.V1V2.Tags_293F40delta","Day210IgG3C.1086C.V1.V2.Tags40delta","Day210IgG3gp70.001428.2.42.V1V240delta",      
                  "Day210IgG3gp70.1012.11.TC21.3257.V1V240delta", "Day210IgG3gp70.1394C9G1.V1V240delta","Day210IgG3gp70.BF1266.431a.V1V240delta","Day210IgG3gp70.Ce1086.B2.V1V240delta",        
                  "Day210IgG3gp70.B.CaseA.V1.V240delta","Day210IgGAE.A244.V1V2.Tags_293F50delta","Day210IgGC.1086C.V1.V2.Tags50delta","Day210IgGgp70_001428.2.42.V1V250delta",       
                  "Day210IgGgp70_1012.11.TC21.3257.V1V250delta","Day210IgGgp70_1394C9G1.V1V250delta","Day210IgGgp70_9004SS.A3.4.V1V250delta","Day210IgGgp70_BF1266.431a.V1V250delta",       
                  "Day210IgGgp70_Ce1086.B2.V1V250delta","Day210IgGgp70.B.CaseA.V1.V250delta" 
  )
  
  individualMarkers <- markerVars
  # Identify the endpoint variable
  endpoint <- "Delta.D210"
  wt <- "wt.D210"
  ptidvar <- "Ptid"

  # Create combined new dataset which has imputed values of demographics (for phase 1 data) from dat.covar.imp AND
  # imputed values for markers (for phase 2 data) from dat.wide.v
  dat.ph1 <- dat.mock %>%
    rename(Ptid = Subjectid) %>%
    filter(Ph1ptids.D210 == 1 & Trt == 1) %>%
    # Drop any observation with NA values in Ptid, Trt, briskfactors, endpoint and wt.D210
    drop_na(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(wt)) %>%
    arrange(desc(get(endpoint)))

  dat.ph2_init <- dat.ph1 %>%
    filter(Ph2ptids.D210 == 1) %>%
    select(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(wt), any_of(markerVars)) %>%
    arrange(desc(get(endpoint)))
}
# Study-agnostic processing ----------------------------------------------------

# Limit total variables that will be included in models
nv <- sum(dat.ph2_init %>% select(matches(endpoint)))

# Remove any predictor variables that are indicator variables and have fewer than 10  0's or 1's (if study_name == "COVE")
# Remove a variable if the number of cases in the variable = 1 subgroup is <= 3 or the number of cases in the variable = 0 subgroup is <= 3 (if study_name != "COVE")
dat.ph2_drop_rare <- drop_predVars_with_fewer_0s_or_1s(dat.ph2_init, c(briskfactors, markerVars))

# Update predictor variables
pred_vars <- dat.ph2_drop_rare %>%
  select(-all_of(ptidvar), -Trt, -all_of(endpoint), -all_of(wt)) %>%
  colnames()

# Remove any baseline risk factors with more than 5% missing values. 
dat.ph2 <- drop_predVars_with_high_total_missing_values(dat.ph2_drop_rare, pred_vars)

# Update pred_vars
pred_vars <- dat.ph2 %>%
  select(-all_of(ptidvar), -Trt, -all_of(endpoint), -all_of(wt)) %>%
  colnames()

# Save ptids to merge with predictions later
ph2_vacc_ptids <- dat.ph2 %>%
  select(all_of(ptidvar), all_of(endpoint), all_of(wt))

# Update briskfactors
briskfactors <- briskfactors[briskfactors %in% pred_vars]
briskfactors_correction <- paste0("Y ~ x + X$", paste0(briskfactors, collapse = " + X$"))

# create "Z" matrix to use for (A)IPW efficient influence function computation
Z_plus_weights <- dat.ph1 %>%
  select(all_of(ptidvar), all_of(endpoint), all_of(wt), Trt, all_of(briskfactors)) %>%
  # Drop any observation with NA values in Ptid, Trt, briskfactors, endpoint or wt.D57
  drop_na(all_of(ptidvar), Trt, all_of(briskfactors), all_of(endpoint), all_of(wt))

# Study-specific combination score creation ------------------------------------
# This results in a dataset called "markers" with the final marker variables

# For COVE study
if (study_name %in% c("COVE", "MockCOVE")) {
  # Create combination scores across the 5 markers
  dat.ph2 <- dat.ph2 %>%
    # generate combination scores for d57
    left_join(get.pca.scores(dat.ph2 %>%
                               select(all_of(ptidvar), Day57bindSpike, Day57bindRBD,
                                 Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50)) %>%
                rename(comb_PC1_d57 = PC1,
                       comb_PC2_d57 = PC2),
              by = ptidvar) %>%
    mutate(
      comb_maxsig.div.score_d57 = get.maxSignalDivScore(
        dat.ph2 %>%
          select(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50),
          "Day57"
      )
    ) %>%
    # generate combination scores for d29
    left_join(get.pca.scores(dat.ph2 %>%
                               select(all_of(ptidvar), Day29bindSpike, Day29bindRBD,
                                 Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50)) %>%
                rename(comb_PC1_d29 = PC1,
                       comb_PC2_d29 = PC2),
              by = ptidvar) %>%
    mutate(
      comb_maxsig.div.score_d29 = get.maxSignalDivScore(
        dat.ph2 %>%
          select(Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50),
          "Day29"
      )
    ) %>%
    # generate combination scores for both d57 and d29
    left_join(get.pca.scores(dat.ph2 %>%
                               select(all_of(ptidvar), Day57bindSpike, Day57bindRBD,
                                 Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50,
                                 Day29bindSpike, Day29bindRBD,
                                 Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50)) %>%
                rename(comb_PC1_d57_d29 = PC1,
                       comb_PC2_d57_d29 = PC2),
              by = ptidvar) %>%
    mutate(
      comb_maxsig.div.score_d57_d29 = get.maxSignalDivScore(
        dat.ph2 %>% select(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50,
                           Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50),
          "Day57_29"
      )
    )

  if (run_prod) {
    dat.ph2 <- dat.ph2 %>%
      # mutate(
      #   comb_nlPCA1_d57 = comb_PC1_d57, comb_nlPCA2_d57 = comb_PC2_d57,
      #   comb_nlPCA1_d29 = comb_PC1_d29, comb_nlPCA2_d29 = comb_PC2_d29,
      #   comb_nlPCA1_d57_d29 = comb_PC1_d57_d29, comb_nlPCA2_d57_d29 = comb_PC1_d57_d29
      # )
      # generate non-linear combination scores for d57
      left_join(get.nonlinearPCA.scores(dat.ph2 %>%
                                          select(all_of(ptidvar), Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50)) %>%
                  rename(comb_nlPCA1_d57 = nlPCA1,
                         comb_nlPCA2_d57 = nlPCA2),
                by = ptidvar) %>%
      # generate non-linear combination scores for d29
      left_join(get.nonlinearPCA.scores(dat.ph2 %>%
                                          select(all_of(ptidvar), Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50)) %>%
                  rename(comb_nlPCA1_d29 = nlPCA1,
                         comb_nlPCA2_d29 = nlPCA2),
                by = ptidvar) %>%
      # generate non-linear combination scores for both d57 and d29
      left_join(get.nonlinearPCA.scores(dat.ph2 %>%
                                          select(all_of(ptidvar), Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50,
                                                 Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50)) %>%
                  rename(comb_nlPCA1_d57_d29 = nlPCA1,
                         comb_nlPCA2_d57_d29 = nlPCA2),
                by = ptidvar)
  }

  # finalize marker data
  markers <- dat.ph2 %>%
    select(-all_of(ptidvar), -Trt, -risk_score, -HighRiskInd, -MinorityInd, -EventIndPrimaryD57, -wt.D57) %>%
    colnames()
}



# For COVEBoost study
if (study_name == "COVEBoost") {
  # # Create combination scores across the 5 markers
  # dat.ph2 <- dat.ph2 %>%
  #   # generate combination scores for d57
  #   left_join(get.pca.scores(dat.ph2 %>%
  #                              select(all_of(ptidvar), 
  #                                     Day57bindSpike, Day57bindRBD,
  #                                     Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50)) %>%
  #               rename(comb_PC1_d57 = PC1,
  #                      comb_PC2_d57 = PC2),
  #             by = ptidvar) %>%
  #   mutate(
  #     comb_maxsig.div.score_d57 = get.maxSignalDivScore(
  #       dat.ph2 %>%
  #         select(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50),
  #       "Day57"
  #     )
  #   ) %>%
  #   # generate combination scores for d29
  #   left_join(get.pca.scores(dat.ph2 %>%
  #                              select(all_of(ptidvar), Day29bindSpike, Day29bindRBD,
  #                                     Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50)) %>%
  #               rename(comb_PC1_d29 = PC1,
  #                      comb_PC2_d29 = PC2),
  #             by = ptidvar) %>%
  #   mutate(
  #     comb_maxsig.div.score_d29 = get.maxSignalDivScore(
  #       dat.ph2 %>%
  #         select(Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50),
  #       "Day29"
  #     )
  #   ) %>%
  #   # generate combination scores for both d57 and d29
  #   left_join(get.pca.scores(dat.ph2 %>%
  #                              select(all_of(ptidvar), Day57bindSpike, Day57bindRBD,
  #                                     Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50,
  #                                     Day29bindSpike, Day29bindRBD,
  #                                     Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50)) %>%
  #               rename(comb_PC1_d57_d29 = PC1,
  #                      comb_PC2_d57_d29 = PC2),
  #             by = ptidvar) %>%
  #   mutate(
  #     comb_maxsig.div.score_d57_d29 = get.maxSignalDivScore(
  #       dat.ph2 %>% select(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50,
  #                          Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50),
  #       "Day57_29"
  #     )
  #   )
  # 
  # if (run_prod) {
  #   dat.ph2 <- dat.ph2 %>%
  #     # mutate(
  #     #   comb_nlPCA1_d57 = comb_PC1_d57, comb_nlPCA2_d57 = comb_PC2_d57,
  #     #   comb_nlPCA1_d29 = comb_PC1_d29, comb_nlPCA2_d29 = comb_PC2_d29,
  #     #   comb_nlPCA1_d57_d29 = comb_PC1_d57_d29, comb_nlPCA2_d57_d29 = comb_PC1_d57_d29
  #     # )
  #     # generate non-linear combination scores for d57
  #     left_join(get.nonlinearPCA.scores(dat.ph2 %>%
  #                                         select(all_of(ptidvar), Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50)) %>%
  #                 rename(comb_nlPCA1_d57 = nlPCA1,
  #                        comb_nlPCA2_d57 = nlPCA2),
  #               by = ptidvar) %>%
  #     # generate non-linear combination scores for d29
  #     left_join(get.nonlinearPCA.scores(dat.ph2 %>%
  #                                         select(all_of(ptidvar), Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50)) %>%
  #                 rename(comb_nlPCA1_d29 = nlPCA1,
  #                        comb_nlPCA2_d29 = nlPCA2),
  #               by = ptidvar) %>%
  #     # generate non-linear combination scores for both d57 and d29
  #     left_join(get.nonlinearPCA.scores(dat.ph2 %>%
  #                                         select(all_of(ptidvar), Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutmn50,
  #                                                Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80, Day29liveneutmn50)) %>%
  #                 rename(comb_nlPCA1_d57_d29 = nlPCA1,
  #                        comb_nlPCA2_d57_d29 = nlPCA2),
  #               by = ptidvar)
  # }
  
  # finalize marker data
  markers <- dat.ph2 %>%
    select(-all_of(ptidvar), -Trt, -risk_score, -HighRiskInd, -MinorityInd, #-CalendarBD1Interval2, -CalendarBD1Interval3, -CalendarBD1Interval4, -NumberdaysDS2toDS3, 
           -EventIndPrimaryOmicronBD29, -wt.BD29) %>%
    colnames()
}




# For HVTN705
if (study_name == "HVTN705") {
  # Create combination scores 
  dat.ph2 <- dat.ph2 %>%
    # generate combination scores for 6 primary markers
    left_join(get.pca.scores(dat.ph2 %>%
                               select(all_of(ptidvar), Day210ELCZ, Day210ELISpotPTEEnv, 
                                      Day210ADCPgp140C97ZAfib, Day210IgG340mdw_V1V2, 
                                      Day210IgG340mdw_gp120_gp140_vm, Day210mdw_xassay_overall)) %>%
                rename(comb_PC1_prim = PC1,
                       comb_PC2_prim = PC2),
              by = ptidvar) %>%
    mutate(
      comb_maxsig.div.score_prim = get.maxSignalDivScore(
        dat.ph2 %>% select(Day210ELCZ, Day210ELISpotPTEEnv, 
                           Day210ADCPgp140C97ZAfib, Day210IgG340mdw_V1V2, 
                           Day210IgG340mdw_gp120_gp140_vm, Day210mdw_xassay_overall),
        "705-Day210-primaryMarkers"
      )
    ) %>%
    # generate combination scores for all markers
    left_join(get.pca.scores(dat.ph2 %>%
                               select(all_of(ptidvar), all_of(markerVars))) %>%
                rename(comb_PC1_all = PC1,
                       comb_PC2_all = PC2),
              by = ptidvar) %>%
    mutate(
      comb_maxsig.div.score_all = get.maxSignalDivScore(
        dat.ph2 %>% select(all_of(markerVars)),
        "705-Day210-allMarkers"
      )
    ) 
  
  
  if (run_prod) {
    dat.ph2 <- dat.ph2 %>%
      # mutate(
      #   comb_nlPCA1_prim = comb_PC1_prim,
      #   comb_nlPCA2_prim = comb_PC2_prim,
      #   comb_nlPCA1_all = comb_PC1_all,
      #   comb_nlPCA2_all = comb_PC2_all
      # )
    # generate non-linear combination scores for 6 primary markers
    left_join(get.nonlinearPCA.scores(dat.ph2 %>%
                                        select(all_of(ptidvar), Day210ELCZ, Day210ELISpotPTEEnv,
                                               Day210ADCPgp140C97ZAfib, Day210IgG340mdw_V1V2,
                                               Day210IgG340mdw_gp120_gp140_vm, Day210mdw_xassay_overall)) %>%
                rename(comb_nlPCA1_prim = nlPCA1,
                       comb_nlPCA2_prim = nlPCA2),
              by = ptidvar) %>%
    # generate non-linear combination scores for all 41 markers
    left_join(get.nonlinearPCA.scores(dat.ph2 %>%
                                        select(all_of(ptidvar), all_of(markerVars))) %>%
                rename(comb_nlPCA1_all = nlPCA1,
                       comb_nlPCA2_all = nlPCA2),
              by = ptidvar)
  }
  
  # finalize marker data
  markers <- c(markerVars, "comb_PC1_prim", "comb_PC2_prim", "comb_maxsig.div.score_prim",
               "comb_PC1_all", "comb_PC2_all", "comb_maxsig.div.score_all", "comb_nlPCA1_prim",
               "comb_nlPCA2_prim", "comb_nlPCA1_all", "comb_nlPCA2_all")
}


if (study_name %in% c("ENSEMBLE")) {
  # Create combination scores across the 5 markers
  dat.ph2 <- dat.ph2 %>%
    # generate combination scores for d29
    left_join(get.pca.scores(dat.ph2 %>%
                               select(all_of(ptidvar), Day29bindSpike, Day29bindRBD,
                                      Day29pseudoneutid50, Day29ADCP)) %>%
                rename(comb_PC1_d29 = PC1,
                       comb_PC2_d29 = PC2),
              by = ptidvar) %>%
    mutate(
      comb_maxsig.div.score_d29 = get.maxSignalDivScore(
        dat.ph2 %>%
          select(Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29ADCP),
        "Day29"
      )
    )  
  
  if (run_prod) {
    dat.ph2 <- dat.ph2 %>%
      # generate non-linear combination scores for d29
      left_join(get.nonlinearPCA.scores(dat.ph2 %>%
                                          select(all_of(ptidvar), Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29ADCP)) %>%
                  rename(comb_nlPCA1_d29 = nlPCA1,
                         comb_nlPCA2_d29 = nlPCA2),
                by = ptidvar) 
  }
  
  # finalize marker data
  markers <- dat.ph2 %>%
    select(-all_of(ptidvar), -Trt, -all_of(briskfactors), -all_of(endpoint), -all_of(wt)) %>%
    colnames()
}

# Define study-specific variable sets of interest ------------------------------
# Results in two objects:
#   1. varset_names: the names of the variable sets
#   2. varset_matrix: a binary matrix. Each row is a variable set,
#     and a 1 in column j indicates that marker j should be included
# set up COVE-specific variable sets
if (study_name %in% c("COVE", "MockCOVE")) {
  # 1. None (No markers; only baseline risk variables), phase 2 data
  varset_baselineRiskFactors <- rep(FALSE, length(markers))

  # 2-12 (Day57)
  varset_bAbSpike_D57 <- create_varsets(markers, grep("(?=.*57)(?=.*bindSpike)", markers, value=TRUE, perl=TRUE))
  varset_bAbRBD_D57 <- create_varsets(markers, grep("(?=.*57)(?=.*bindRBD)", markers, value=TRUE, perl=TRUE))
  varset_pnabID50_D57 <- create_varsets(markers, grep("(?=.*57)(?=.*id50)", markers, value=TRUE, perl=TRUE))
  varset_pnabID80_D57 <- create_varsets(markers, grep("(?=.*57)(?=.*id80)", markers, value=TRUE, perl=TRUE))
  varset_lnabMN50_D57 <- create_varsets(markers, grep("(?=.*57)(?=.*mn50)", markers, value=TRUE, perl=TRUE))
  varset_bAb_pnabID50_D57 <- create_varsets(markers, c(grep("(?=.*57)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*57)(?=.*id50)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_pnabID80_D57 <- create_varsets(markers, c(grep("(?=.*57)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*57)(?=.*id80)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_lnabMN50_D57 <- create_varsets(markers, c(grep("(?=.*57)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*57)(?=.*mn50)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_combScores_D57 <- create_varsets(markers, c(grep("(?=.*57)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                         grep("(?=.*57)(?=.*comb)(^((?!29).)*$)", markers, value=TRUE, perl=TRUE)))
  varset_allMarkers_D57 <- create_varsets(markers, grep("(?=.*57)(^((?!comb).)*$)", markers, value=TRUE, perl=TRUE))
  varset_allMarkers_combScores_D57 <- create_varsets(markers, grep("(?=.*57)(^((?!29).)*$)", markers, value=TRUE, perl=TRUE))

  # 13-23 (Day29)
  varset_bAbSpike_D29 <- create_varsets(markers, grep("(?=.*29)(?=.*bindSpike)", markers, value=TRUE, perl=TRUE))
  varset_bAbRBD_D29 <- create_varsets(markers, grep("(?=.*29)(?=.*bindRBD)", markers, value=TRUE, perl=TRUE))
  varset_pnabID50_D29 <- create_varsets(markers, grep("(?=.*29)(?=.*id50)", markers, value=TRUE, perl=TRUE))
  varset_pnabID80_D29 <- create_varsets(markers, grep("(?=.*29)(?=.*id80)", markers, value=TRUE, perl=TRUE))
  varset_lnabMN50_D29 <- create_varsets(markers, grep("(?=.*29)(?=.*mn50)", markers, value=TRUE, perl=TRUE))
  varset_bAb_pnabID50_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*29)(?=.*id50)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_pnabID80_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*29)(?=.*id80)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_lnabMN50_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*29)(?=.*mn50)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_combScores_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                         grep("(?=.*29)(?=.*comb)(^((?!57).)*$)", markers, value=TRUE, perl=TRUE)))
  varset_allMarkers_D29 <- create_varsets(markers, grep("(?=.*29)(^((?!comb).)*$)", markers, value=TRUE, perl=TRUE))
  varset_allMarkers_combScores_D29 <- create_varsets(markers, grep("(?=.*29)(^((?!57).)*$)", markers, value=TRUE, perl=TRUE))

  # 24-34 (Day57 and Day29)
  varset_bAbSpike_D29_D57 <- create_varsets(markers, grep("(?=.*bindSpike)", markers, value=TRUE, perl=TRUE))
  varset_bAbRBD_D29_D57 <- create_varsets(markers, grep("(?=.*bindRBD)", markers, value=TRUE, perl=TRUE))
  varset_pnabID50_D29_D57 <- create_varsets(markers, grep("(?=.*id50)", markers, value=TRUE, perl=TRUE))
  varset_pnabID80_D29_D57 <- create_varsets(markers, grep("(?=.*id80)", markers, value=TRUE, perl=TRUE))
  varset_lnabMN50_D29_D57 <- create_varsets(markers, grep("(?=.*mn50)", markers, value=TRUE, perl=TRUE))
  varset_bAb_pnabID50_D29_D57 <- create_varsets(markers, c(grep("(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*id50)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_pnabID80_D29_D57 <- create_varsets(markers, c(grep("(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                           grep("(?=.*id80)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_lnabMN50_D29_D57 <- create_varsets(markers, c(grep("(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                            grep("(?=.*mn50)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_combScores_D29_D57 <- create_varsets(markers, c(grep("(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                         grep("(?=.*comb)(?=.*d57_d29)", markers, value=TRUE, perl=TRUE)))
  varset_allMarkers_D29_D57 <- create_varsets(markers, grep("(^((?!comb).)*$)", markers, value=TRUE, perl=TRUE))
  varset_allMarkers_combScores_D29_D57 <- create_varsets(markers, c(grep("(^((?!comb).)*$)", markers, value=TRUE, perl=TRUE),
                                                                    grep("(?=.*comb)(?=.*d57_d29)", markers, value=TRUE, perl=TRUE)))

  varset_names <- c("1_baselineRiskFactors",
                    "2_bAbSpike_D57", "3_bAbRBD_D57", "4_pnabID50_D57",
                    "5_pnabID80_D57", "6_lnabMN50_D57", "7_bAb_pnabID50_D57", "8_bAb_pnabID80_D57", "9_bAb_lnabMN50_D57",
                    "10_bAb_combScores_D57", "11_allMarkers_D57", "12_allMarkers_combScores_D57",

                    "13_bAbSpike_D29", "14_bAbRBD_D29", "15_pnabID50_D29",
                    "16_pnabID80_D29", "17_lnabMN50_D29", "18_bAb_pnabID50_D29", "19_bAb_pnabID80_D29", "20_bAb_lnabMN50_D29",
                    "21_bAb_combScores_D29", "22_allMarkers_D29", "23_allMarkers_combScores_D29",

                    "24_bAbSpike_D29_D57", "25_bAbRBD_D29_D57", "26_pnabID50_D29_D57",
                    "27_pnabID80_D29_D57", "28_lnabMN50_D29_D57", "29_bAb_pnabID50_D29_D57", "30_bAb_pnabID80_D29_D57", "31_bAb_lnabMN50_D29_D57",
                    "32_bAb_combScores_D29_D57", "33_allMarkers_D29_D57", "34_allMarkers_combScores_D29_D57")

  # set up a matrix of all
  varset_matrix <- rbind(varset_baselineRiskFactors,
                         varset_bAbSpike_D57, varset_bAbRBD_D57, varset_pnabID50_D57,
                         varset_pnabID80_D57, varset_lnabMN50_D57, varset_bAb_pnabID50_D57, varset_bAb_pnabID80_D57, varset_bAb_lnabMN50_D57,
                         varset_bAb_combScores_D57, varset_allMarkers_D57, varset_allMarkers_combScores_D57,

                         varset_bAbSpike_D29, varset_bAbRBD_D29, varset_pnabID50_D29,
                         varset_pnabID80_D29, varset_lnabMN50_D29, varset_bAb_pnabID50_D29, varset_bAb_pnabID80_D29, varset_bAb_lnabMN50_D29,
                         varset_bAb_combScores_D29, varset_allMarkers_D29, varset_allMarkers_combScores_D29,

                         varset_bAbSpike_D29_D57, varset_bAbRBD_D29_D57, varset_pnabID50_D29_D57,
                         varset_pnabID80_D29_D57, varset_lnabMN50_D29_D57, varset_bAb_pnabID50_D29_D57, varset_bAb_pnabID80_D29_D57, varset_bAb_lnabMN50_D29_D57,
                         varset_bAb_combScores_D29_D57, varset_allMarkers_D29_D57, varset_allMarkers_combScores_D29_D57)
}



# Define study-specific variable sets of interest ------------------------------
# Results in two objects:
#   1. varset_names: the names of the variable sets
#   2. varset_matrix: a binary matrix. Each row is a variable set,
#     and a 1 in column j indicates that marker j should be included
# set up COVE-specific variable sets
if (study_name == "COVEBoost") {
  # 1. None (No markers; only baseline risk variables), phase 2 data
  varset_briskfactors <- rep(FALSE, length(markers))
  # varset_briskfactors_trt <- rep(FALSE, length(markers))
  # varset_briskfactors_trt_DS2toS3 <- rep(FALSE, length(markers))
  
  ################# BA.1 #############################
  varset_BD1_bSpike_BA.1 <- create_varsets(markers, grep("(?=.*BD1)(?=.*bindSpike)(?=.*BA.1)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bSpike_BA.1 <- create_varsets(markers, grep("(?=.*BD29)(?=.*bindSpike)(?=.*BA.1)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bSpike_BA.1 <- create_varsets(markers, grep("(?=.*Delta)(?=.*BD29)(?=.*BD1)(?=.*bindSpike)(?=.*BA.1)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_nAb_BA.1 <- create_varsets(markers, grep("(?=.*BD1)(?=.*pseudoneut)(?=.*BA.1)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_nAb_BA.1 <- create_varsets(markers, grep("(?=.*BD29)(?=.*pseudoneutid50)(?=.*BA.1)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_nAb_BA.1 <- create_varsets(markers, grep("(?=.*Delta)(?=.*BD29)(?=.*BD1)(?=.*pseudoneutid50)(?=.*BA.1)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_BD29_bSpike_BA.1 <- create_varsets(markers, grep("(?=.*bindSpike)(?=.*BA.1)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bSpike_BA.1 <- create_varsets(markers, grep("(?=.*BD1)(?=.*bindSpike)(?=.*BA.1)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bSpike_BA.1 <- create_varsets(markers, grep("(?=.*BD29)(?=.*bindSpike)(?=.*BA.1)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_BD29_nAb_BA.1 <- create_varsets(markers, grep("(?=.*pseudoneutid50)(?=.*BA.1)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_nAb_BA.1 <- create_varsets(markers, grep("(?=.*BD1)(?=.*pseudoneutid50)(?=.*BA.1)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_nAb_BA.1 <- create_varsets(markers, grep("(?=.*BD29)(?=.*pseudoneutid50)(?=.*BA.1)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_bSpikenAb_BA.1 <- create_varsets(markers, grep("(?=.*BD1)(?=.*BA.1)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bSpikenAb_BA.1 <- create_varsets(markers, grep("(?=.*BD29)(?=.*BA.1)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bSpikenAb_BA.1 <- create_varsets(markers, grep("(?=.*Delta)(?=.*BD29)(?=.*BD1)(?=.*BA.1)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_BD29_bSpikenAb_BA.1 <- create_varsets(markers, grep("(?=.*BA.1)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bSpikenAb_BA.1 <- create_varsets(markers, grep("(?=.*BD1)(?=.*BA.1)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bSpikenAb_BA.1 <- create_varsets(markers, grep("(?=.*BD29)(?=.*BA.1)", markers, value=TRUE, perl=TRUE))
  
  ################# D614G #############################
  varset_BD1_bSpike_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*bindSpike)(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bSpike_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*bindSpike)(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bSpike_D614G <- create_varsets(markers, grep("(?=.*Delta)(?=.*BD29)(?=.*BD1)(?=.*bindSpike)(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_bRBD_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*bindRBD)(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bRBD_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*bindRBD)(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bRBD_D614G <- create_varsets(markers, grep("(?=.*Delta)(?=.*BD29)(?=.*BD1)(?=.*bindRBD)(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_nAb_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*pseudoneut)(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_nAb_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*pseudoneutid50)(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_nAb_D614G <- create_varsets(markers, grep("(?=.*Delta)(?=.*BD29)(?=.*BD1)(?=.*pseudoneutid50)(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  ##
  varset_BD1_BD29_bSpike_D614G <- create_varsets(markers, grep("(?=.*bindSpike)(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bSpike_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*bindSpike)(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bSpike_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*bindSpike)(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_BD29_bRBD_D614G <- create_varsets(markers, grep("(?=.*bindRBD)(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bRBD_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*bindRBD)(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bRBD_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*bindRBD)(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_BD29_nAb_D614G <- create_varsets(markers, grep("(?=.*pseudoneutid50)(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_nAb_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*pseudoneutid50)(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_nAb_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*pseudoneutid50)(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  ##
  varset_BD1_bSpikebRBD_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*D614G)(^((?!over).(?!pseudoneutid50))*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bSpikebRBD_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*D614G)(^((?!over).(?!pseudoneutid50))*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bSpikebRBD_D614G <- create_varsets(markers, grep("(?=.*Delta)(?=.*BD29)(?=.*BD1)(?=.*D614G)(^((?!pseudoneutid50).)*$)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_bSpikenAb_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*D614G)(^((?!RBD).(?!over))*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bSpikenAb_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*D614G)(^((?!RBD).(?!over))*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bSpikenAb_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*D614G)(?=.*over)(^((?!RBD).)*$)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_bRBDnAb_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*D614G)(^((?!Spike).(?!over))*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bRBDnAb_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*D614G)(^((?!Spike).(?!over))*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bRBDnAb_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*D614G)(?=.*over)(^((?!Spike).)*$)", markers, value=TRUE, perl=TRUE))
  ##
  varset_BD1_BD29_bSpikebRBD_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*RBD))(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bSpikebRBD_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*RBD))((?=.*BD1)|(?=.*Delta))(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bSpikebRBD_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*RBD))((?=.*BD29)|(?=.*Delta))(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_BD29_bSpikenAb_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*pseudoneutid50))(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bSpikenAb_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*pseudoneutid50))((?=.*BD1)|(?=.*Delta))(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bSpikenAb_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*pseudoneutid50))((?=.*BD29)|(?=.*Delta))(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_BD29_bRBDnAb_D614G <- create_varsets(markers, grep("((?=.*RBD)|(?=.*pseudoneutid50))(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bRBDnAb_D614G <- create_varsets(markers, grep("((?=.*RBD)|(?=.*pseudoneutid50))((?=.*BD1)|(?=.*Delta))(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bRBDnAb_D614G <- create_varsets(markers, grep("((?=.*RBD)|(?=.*pseudoneutid50))((?=.*BD29)|(?=.*Delta))(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  ##
  varset_BD1_BD29_bSpikebRBDnAb_D614G <- create_varsets(markers, grep("(?=.*D614G)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bSpikebRBDnAb_D614G <- create_varsets(markers, grep("((?=.*BD1)|(?=.*Delta))(?=.*D614G)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bSpikebRBDnAb_D614G <- create_varsets(markers, grep("((?=.*BD29)|(?=.*Delta))(?=.*D614G)", markers, value=TRUE, perl=TRUE))

  ################# BA.1 and D614G #############################
  varset_BD1_bSpike_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*bindSpike)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bSpike_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*bindSpike)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bSpike_BA.1_D614G <- create_varsets(markers, grep("(?=.*Delta)(?=.*BD29)(?=.*BD1)(?=.*bindSpike)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_nAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*pseudoneut)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_nAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*pseudoneutid50)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_nAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*Delta)(?=.*BD29)(?=.*BD1)(?=.*pseudoneutid50)", markers, value=TRUE, perl=TRUE))
  ##
  varset_BD1_BD29_bSpike_BA.1_D614G <- create_varsets(markers, grep("(?=.*bindSpike)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bSpike_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*bindSpike)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bSpike_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*bindSpike)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_BD29_nAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*pseudoneutid50)(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_nAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD1)(?=.*pseudoneutid50)", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_nAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*pseudoneutid50)", markers, value=TRUE, perl=TRUE))
  ##
  varset_BD1_bSpikebRBD_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD1)(^((?!over).(?!pseudoneutid50))*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bSpikebRBD_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD29)(^((?!over).(?!pseudoneutid50))*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bSpikebRBD_BA.1_D614G <- create_varsets(markers, grep("(?=.*Delta)(?=.*BD29)(?=.*BD1)(^((?!pseudoneutid50).)*$)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_bSpikenAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD1)(^((?!RBD).(?!over))*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bSpikenAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD29)(^((?!RBD).(?!over))*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bSpikenAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*over)(^((?!RBD).)*$)", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_bRBDnAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD1)(^((?!Spike).(?!over))*$)", markers, value=TRUE, perl=TRUE))
  varset_BD29_bRBDnAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD29)(^((?!Spike).(?!over))*$)", markers, value=TRUE, perl=TRUE))
  varset_DeltaBD29overBD1_bRBDnAb_BA.1_D614G <- create_varsets(markers, grep("(?=.*BD29)(?=.*over)(^((?!Spike).)*$)", markers, value=TRUE, perl=TRUE))
  ##
  varset_BD1_BD29_bSpikebRBD_BA.1_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*RBD))(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bSpikebRBD_BA.1_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*RBD))((?=.*BD1)|(?=.*Delta))", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bSpikebRBD_BA.1_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*RBD))((?=.*BD29)|(?=.*Delta))", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_BD29_bSpikenAb_BA.1_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*pseudoneutid50))(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bSpikenAb_BA.1_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*pseudoneutid50))((?=.*BD1)|(?=.*Delta))", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bSpikenAb_BA.1_D614G <- create_varsets(markers, grep("((?=.*Spike)|(?=.*pseudoneutid50))((?=.*BD29)|(?=.*Delta))", markers, value=TRUE, perl=TRUE))
  
  varset_BD1_BD29_bRBDnAb_BA.1_D614G <- create_varsets(markers, grep("((?=.*RBD)|(?=.*pseudoneutid50))(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bRBDnAb_BA.1_D614G <- create_varsets(markers, grep("((?=.*RBD)|(?=.*pseudoneutid50))((?=.*BD1)|(?=.*Delta))", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bRBDnAb_BA.1_D614G <- create_varsets(markers, grep("((?=.*RBD)|(?=.*pseudoneutid50))((?=.*BD29)|(?=.*Delta))", markers, value=TRUE, perl=TRUE))
  ##
  varset_BD1_BD29_bSpikebRBDnAb_BA.1_D614G <- create_varsets(markers, grep("(^((?!over).)*$)", markers, value=TRUE, perl=TRUE))
  varset_BD1_Delta_bSpikebRBDnAb_BA.1_D614G <- create_varsets(markers, grep("((?=.*BD1)|(?=.*Delta))", markers, value=TRUE, perl=TRUE))
  varset_BD29_Delta_bSpikebRBDnAb_BA.1_D614G <- create_varsets(markers, grep("((?=.*BD29)|(?=.*Delta))", markers, value=TRUE, perl=TRUE))
  ## complete set
  varset_BD1_BD29_Delta_bSpikebRBDnAb_BA.1_D614G <- create_varsets(markers, grep("", markers, value=TRUE, perl=TRUE))
  
  ###############################################################
  varset_names <- data.frame(varsets = c("varset_briskfactors", #"varset_briskfactors_trt", "varset_briskfactors_trt_DS2toS3", 
                                             
                                             "varset_BD1_bSpike_BA.1", "varset_BD29_bSpike_BA.1", "varset_DeltaBD29overBD1_bSpike_BA.1",
                                             "varset_BD1_nAb_BA.1", "varset_BD29_nAb_BA.1", "varset_DeltaBD29overBD1_nAb_BA.1",
                                             "varset_BD1_BD29_bSpike_BA.1", "varset_BD1_Delta_bSpike_BA.1", "varset_BD29_Delta_bSpike_BA.1",
                                             "varset_BD1_BD29_nAb_BA.1", "varset_BD1_Delta_nAb_BA.1", "varset_BD29_Delta_nAb_BA.1", 
                                             "varset_BD1_bSpikenAb_BA.1", "varset_BD29_bSpikenAb_BA.1", "varset_DeltaBD29overBD1_bSpikenAb_BA.1",
                                             "varset_BD1_BD29_bSpikenAb_BA.1", "varset_BD1_Delta_bSpikenAb_BA.1", "varset_BD29_Delta_bSpikenAb_BA.1",
                                             
                                             "varset_BD1_bSpike_D614G", "varset_BD29_bSpike_D614G", "varset_DeltaBD29overBD1_bSpike_D614G", 
                                             "varset_BD1_bRBD_D614G", "varset_BD29_bRBD_D614G", "varset_DeltaBD29overBD1_bRBD_D614G", 
                                             "varset_BD1_nAb_D614G", "varset_BD29_nAb_D614G", "varset_DeltaBD29overBD1_nAb_D614G", 
                                             "varset_BD1_BD29_bSpike_D614G", "varset_BD1_Delta_bSpike_D614G", "varset_BD29_Delta_bSpike_D614G", 
                                             "varset_BD1_BD29_bRBD_D614G", "varset_BD1_Delta_bRBD_D614G", "varset_BD29_Delta_bRBD_D614G", 
                                             "varset_BD1_BD29_nAb_D614G", "varset_BD1_Delta_nAb_D614G", "varset_BD29_Delta_nAb_D614G", 
                                             "varset_BD1_bSpikebRBD_D614G", "varset_BD29_bSpikebRBD_D614G", "varset_DeltaBD29overBD1_bSpikebRBD_D614G", 
                                             "varset_BD1_bSpikenAb_D614G", "varset_BD29_bSpikenAb_D614G", "varset_DeltaBD29overBD1_bSpikenAb_D614G", 
                                             "varset_BD1_bRBDnAb_D614G", "varset_BD29_bRBDnAb_D614G", "varset_DeltaBD29overBD1_bRBDnAb_D614G", 
                                             "varset_BD1_BD29_bSpikebRBD_D614G", "varset_BD1_Delta_bSpikebRBD_D614G", "varset_BD29_Delta_bSpikebRBD_D614G", 
                                             "varset_BD1_BD29_bSpikenAb_D614G", "varset_BD1_Delta_bSpikenAb_D614G", "varset_BD29_Delta_bSpikenAb_D614G", 
                                             "varset_BD1_BD29_bRBDnAb_D614G", "varset_BD1_Delta_bRBDnAb_D614G", "varset_BD29_Delta_bRBDnAb_D614G", 
                                             "varset_BD1_BD29_bSpikebRBDnAb_D614G", "varset_BD1_Delta_bSpikebRBDnAb_D614G", "varset_BD29_Delta_bSpikebRBDnAb_D614G",
                                             
                                             "varset_BD1_bSpike_BA.1_D614G", "varset_BD29_bSpike_BA.1_D614G", "varset_DeltaBD29overBD1_bSpike_BA.1_D614G", 
                                             "varset_BD1_nAb_BA.1_D614G", "varset_BD29_nAb_BA.1_D614G", "varset_DeltaBD29overBD1_nAb_BA.1_D614G",
                                             "varset_BD1_BD29_bSpike_BA.1_D614G", "varset_BD1_Delta_bSpike_BA.1_D614G", "varset_BD29_Delta_bSpike_BA.1_D614G",
                                             "varset_BD1_BD29_nAb_BA.1_D614G", "varset_BD1_Delta_nAb_BA.1_D614G", "varset_BD29_Delta_nAb_BA.1_D614G", 
                                             "varset_BD1_bSpikebRBD_BA.1_D614G", "varset_BD29_bSpikebRBD_BA.1_D614G", "varset_DeltaBD29overBD1_bSpikebRBD_BA.1_D614G",
                                             "varset_BD1_bSpikenAb_BA.1_D614G", "varset_BD29_bSpikenAb_BA.1_D614G", "varset_DeltaBD29overBD1_bSpikenAb_BA.1_D614G",
                                             "varset_BD1_bRBDnAb_BA.1_D614G", "varset_BD29_bRBDnAb_BA.1_D614G", "varset_DeltaBD29overBD1_bRBDnAb_BA.1_D614G",
                                             "varset_BD1_BD29_bSpikebRBD_BA.1_D614G", "varset_BD1_Delta_bSpikebRBD_BA.1_D614G", "varset_BD29_Delta_bSpikebRBD_BA.1_D614G",
                                             "varset_BD1_BD29_bSpikenAb_BA.1_D614G", "varset_BD1_Delta_bSpikenAb_BA.1_D614G", "varset_BD29_Delta_bSpikenAb_BA.1_D614G",
                                             "varset_BD1_BD29_bRBDnAb_BA.1_D614G", "varset_BD1_Delta_bRBDnAb_BA.1_D614G", "varset_BD29_Delta_bRBDnAb_BA.1_D614G",
                                             "varset_BD1_BD29_bSpikebRBDnAb_BA.1_D614G", "varset_BD1_Delta_bSpikebRBDnAb_BA.1_D614G", "varset_BD29_Delta_bSpikebRBDnAb_BA.1_D614G",
                                             "varset_BD1_BD29_Delta_bSpikebRBDnAb_BA.1_D614G")) %>%
    tibble::rownames_to_column() %>%
    rename(varsetNo = rowname) %>%
    mutate(varsets2 = str_replace(varsets, "varset", varsetNo)) %>%
    .$varsets2
  
  # set up a matrix of all
  varset_matrix <- rbind(varset_briskfactors, #varset_briskfactors_trt, varset_briskfactors_trt_DS2toS3, 
                         
                         varset_BD1_bSpike_BA.1, varset_BD29_bSpike_BA.1, varset_DeltaBD29overBD1_bSpike_BA.1,
                         varset_BD1_nAb_BA.1, varset_BD29_nAb_BA.1, varset_DeltaBD29overBD1_nAb_BA.1,
                         varset_BD1_BD29_bSpike_BA.1, varset_BD1_Delta_bSpike_BA.1, varset_BD29_Delta_bSpike_BA.1,
                         varset_BD1_BD29_nAb_BA.1, varset_BD1_Delta_nAb_BA.1, varset_BD29_Delta_nAb_BA.1, 
                         varset_BD1_bSpikenAb_BA.1, varset_BD29_bSpikenAb_BA.1, varset_DeltaBD29overBD1_bSpikenAb_BA.1,
                         varset_BD1_BD29_bSpikenAb_BA.1, varset_BD1_Delta_bSpikenAb_BA.1, varset_BD29_Delta_bSpikenAb_BA.1,
                         
                         varset_BD1_bSpike_D614G, varset_BD29_bSpike_D614G, varset_DeltaBD29overBD1_bSpike_D614G, 
                         varset_BD1_bRBD_D614G, varset_BD29_bRBD_D614G, varset_DeltaBD29overBD1_bRBD_D614G, 
                         varset_BD1_nAb_D614G, varset_BD29_nAb_D614G, varset_DeltaBD29overBD1_nAb_D614G, 
                         varset_BD1_BD29_bSpike_D614G, varset_BD1_Delta_bSpike_D614G, varset_BD29_Delta_bSpike_D614G, 
                         varset_BD1_BD29_bRBD_D614G, varset_BD1_Delta_bRBD_D614G, varset_BD29_Delta_bRBD_D614G, 
                         varset_BD1_BD29_nAb_D614G, varset_BD1_Delta_nAb_D614G, varset_BD29_Delta_nAb_D614G, 
                         varset_BD1_bSpikebRBD_D614G, varset_BD29_bSpikebRBD_D614G, varset_DeltaBD29overBD1_bSpikebRBD_D614G, 
                         varset_BD1_bSpikenAb_D614G, varset_BD29_bSpikenAb_D614G, varset_DeltaBD29overBD1_bSpikenAb_D614G, 
                         varset_BD1_bRBDnAb_D614G, varset_BD29_bRBDnAb_D614G, varset_DeltaBD29overBD1_bRBDnAb_D614G, 
                         varset_BD1_BD29_bSpikebRBD_D614G, varset_BD1_Delta_bSpikebRBD_D614G, varset_BD29_Delta_bSpikebRBD_D614G, 
                         varset_BD1_BD29_bSpikenAb_D614G, varset_BD1_Delta_bSpikenAb_D614G, varset_BD29_Delta_bSpikenAb_D614G, 
                         varset_BD1_BD29_bRBDnAb_D614G, varset_BD1_Delta_bRBDnAb_D614G, varset_BD29_Delta_bRBDnAb_D614G, 
                         varset_BD1_BD29_bSpikebRBDnAb_D614G, varset_BD1_Delta_bSpikebRBDnAb_D614G, varset_BD29_Delta_bSpikebRBDnAb_D614G,
                         
                         varset_BD1_bSpike_BA.1_D614G, varset_BD29_bSpike_BA.1_D614G, varset_DeltaBD29overBD1_bSpike_BA.1_D614G, 
                         varset_BD1_nAb_BA.1_D614G, varset_BD29_nAb_BA.1_D614G, varset_DeltaBD29overBD1_nAb_BA.1_D614G,
                         varset_BD1_BD29_bSpike_BA.1_D614G, varset_BD1_Delta_bSpike_BA.1_D614G, varset_BD29_Delta_bSpike_BA.1_D614G,
                         varset_BD1_BD29_nAb_BA.1_D614G, varset_BD1_Delta_nAb_BA.1_D614G, varset_BD29_Delta_nAb_BA.1_D614G, 
                         varset_BD1_bSpikebRBD_BA.1_D614G, varset_BD29_bSpikebRBD_BA.1_D614G, varset_DeltaBD29overBD1_bSpikebRBD_BA.1_D614G,
                         varset_BD1_bSpikenAb_BA.1_D614G, varset_BD29_bSpikenAb_BA.1_D614G, varset_DeltaBD29overBD1_bSpikenAb_BA.1_D614G,
                         varset_BD1_bRBDnAb_BA.1_D614G, varset_BD29_bRBDnAb_BA.1_D614G, varset_DeltaBD29overBD1_bRBDnAb_BA.1_D614G,
                         varset_BD1_BD29_bSpikebRBD_BA.1_D614G, varset_BD1_Delta_bSpikebRBD_BA.1_D614G, varset_BD29_Delta_bSpikebRBD_BA.1_D614G,
                         varset_BD1_BD29_bSpikenAb_BA.1_D614G, varset_BD1_Delta_bSpikenAb_BA.1_D614G, varset_BD29_Delta_bSpikenAb_BA.1_D614G,
                         varset_BD1_BD29_bRBDnAb_BA.1_D614G, varset_BD1_Delta_bRBDnAb_BA.1_D614G, varset_BD29_Delta_bRBDnAb_BA.1_D614G,
                         varset_BD1_BD29_bSpikebRBDnAb_BA.1_D614G, varset_BD1_Delta_bSpikebRBDnAb_BA.1_D614G, varset_BD29_Delta_bSpikebRBDnAb_BA.1_D614G,
                         varset_BD1_BD29_Delta_bSpikebRBDnAb_BA.1_D614G)
}



# set up HVTN705-specific variable sets
if (study_name == "HVTN705") {
  # 1. None (No markers; only baseline risk variables), phase 2 data
  varset_baselineFactors <- rep(FALSE, length(markers))  # markers[c(varset_baselineFactors)]
  
  # 2-16
  varset_M7ELISA <- create_varsets(markers, grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE))
  varset_M7Tcells <- create_varsets(markers, grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE))
  varset_M7ADCP <- create_varsets(markers, grep("ADCP", markers, value=TRUE, perl=TRUE))
  varset_M7IgG3 <- create_varsets(markers, grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE))
  varset_M7IgG3_gp140 <- create_varsets(markers, grep("(?=.*IgG3)(?=.*gp140)", markers, value=TRUE, perl=TRUE))
  varset_M7IgG3_gp120 <- create_varsets(markers, grep("(?=.*IgG3)(?=.*gp120)", markers, value=TRUE, perl=TRUE))
  varset_M7IgG3_V1V2 <- create_varsets(markers, grep("(?=.*IgG3)(?=.*V1)(?=.*V2)", markers, value=TRUE, perl=TRUE))
  varset_M7IgG3_gp41 <- create_varsets(markers, grep("(?=.*IgG3)(?=.*gp41)", markers, value=TRUE, perl=TRUE))
  varset_M7IgG3_breadthScore <- create_varsets(markers, grep("(?=.*IgG3)(?=.*mdw)(^((?!multi).)*$)", markers, value=TRUE, perl=TRUE))
  varset_M7IgG3_multi <- create_varsets(markers, grep("(?=.*IgG3)(?=.*mdw)(?=.*multi)", markers, value=TRUE, perl=TRUE))
  varset_M7IgGt_V1V2 <- create_varsets(markers, grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE))
  varset_M7ADCC <- create_varsets(markers, grep("ADCC", markers, value=TRUE, perl=TRUE))
  varset_M7_combPrimaryMarkers <- create_varsets(markers, grep("(?=.*prim)(?=.*comb)", markers, value=TRUE, perl=TRUE))
  varset_M7_combAllMarkers <- create_varsets(markers, grep("(?=.*all)(?=.*comb)", markers, value=TRUE, perl=TRUE))
  varset_M7_overallScore <- create_varsets(markers, grep("(?=.*mdw)(?=.*xassay)", markers, value=TRUE, perl=TRUE))
  
  varset_M7_2_3 <- create_varsets(markers, c(grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE),
                                             grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE)))
  varset_M7_2_4 <- create_varsets(markers, c(grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE),
                                             grep("ADCP", markers, value=TRUE, perl=TRUE)))
  varset_M7_2_5_12 <- create_varsets(markers, c(grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE),
                                                grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE),
                                                grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE)))
  varset_M7_2_13 <- create_varsets(markers, c(grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE),
                                              grep("ADCC", markers, value=TRUE, perl=TRUE)))
  varset_M7_3_4 <- create_varsets(markers, c(grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                             grep("ADCP", markers, value=TRUE, perl=TRUE)))
  varset_M7_3_5_12 <- create_varsets(markers, c(grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                                grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE),
                                                grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE)))
  varset_M7_3_13 <- create_varsets(markers, c(grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                              grep("ADCC", markers, value=TRUE, perl=TRUE)))
  varset_M7_4_5_12 <- create_varsets(markers, c(grep("ADCP", markers, value=TRUE, perl=TRUE),
                                                grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE),
                                                grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE)))
  varset_M7_4_13 <- create_varsets(markers, c(grep("ADCP", markers, value=TRUE, perl=TRUE),
                                              grep("ADCC", markers, value=TRUE, perl=TRUE)))
  varset_M7_5_12_13 <- create_varsets(markers, c(grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE),
                                                 grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE),
                                                 grep("ADCC", markers, value=TRUE, perl=TRUE)))
  varset_M7_2_3_4 <- create_varsets(markers, c(grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE),
                                               grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                               grep("ADCP", markers, value=TRUE, perl=TRUE)))
  varset_M7_2_3_5_12 <- create_varsets(markers, c(grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE),
                                                  grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                                  grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE),
                                                  grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE)))
  varset_M7_2_3_13 <- create_varsets(markers, c(grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE),
                                                grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                                grep("ADCC", markers, value=TRUE, perl=TRUE)))
  varset_M7_3_4_5_12 <- create_varsets(markers, c(grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                                  grep("ADCP", markers, value=TRUE, perl=TRUE),
                                                  grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE),
                                                  grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE)))
  varset_M7_3_4_13 <- create_varsets(markers, c(grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                                grep("ADCP", markers, value=TRUE, perl=TRUE),
                                                grep("ADCC", markers, value=TRUE, perl=TRUE)))
  varset_M7_4_5_12_13 <- create_varsets(markers, c(grep("ADCP", markers, value=TRUE, perl=TRUE),
                                                   grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE),
                                                   grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE),
                                                   grep("ADCC", markers, value=TRUE, perl=TRUE)))
  varset_M7_2_3_4_5_12 <- create_varsets(markers, c(grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE),
                                                    grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                                    grep("ADCP", markers, value=TRUE, perl=TRUE),
                                                    grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE),
                                                    grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE)))
  varset_M7_2_3_4_13 <- create_varsets(markers, c(grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE),
                                                  grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                                  grep("ADCP", markers, value=TRUE, perl=TRUE),
                                                  grep("ADCC", markers, value=TRUE, perl=TRUE)))
  varset_M7_3_4_5_12_13 <- create_varsets(markers, c(grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                                     grep("ADCP", markers, value=TRUE, perl=TRUE),
                                                     grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE),
                                                     grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE),
                                                     grep("ADCC", markers, value=TRUE, perl=TRUE)))
  varset_M7_2_3_4_5_12_13 <- create_varsets(markers, c(grep("ELCZ|ELMo", markers, value=TRUE, perl=TRUE),
                                                       grep("ELISpot|ICS", markers, value=TRUE, perl=TRUE),
                                                       grep("ADCP", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*IgG3)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*IgG)(?=.*V1)(^((?!IgG3).)*$)", markers, value=TRUE, perl=TRUE),
                                                       grep("ADCC", markers, value=TRUE, perl=TRUE)))
  
  varset_names <- c("1_baselineFactors", "2_M7ELISA", "3_M7Tcells","4_M7ADCP","5_M7IgG3",
                    "6_M7IgG3_gp140","7_M7IgG3_gp120","8_M7IgG3_V1V2","9_M7IgG3_gp41",
                    "10_M7IgG3_breadthScore","11_M7IgG3_multi","12_M7IgGt_V1V2","13_M7ADCC",
                    "14_M7_combPrimaryMarkers", "15_M7_combAllMarkers","16_M7_overallScore", 
                    "17_M7_2_3", "18_M7_2_4", "19_M7_2_5_12", "20_M7_2_13", "21_M7_3_4", 
                    "22_M7_3_5_12", "23_M7_3_13", "24_M7_4_5_12", "25_M7_4_13", 
                    "26_M7_5_12_13","27_M7_2_3_4", "28_M7_2_3_5_12","29_M7_2_3_13","30_M7_3_4_5_12",
                    "31_M7_3_4_13","32_M7_4_5_12_13","33_M7_2_3_4_5_12","34_M7_2_3_4_13",
                    "35_M7_3_4_5_12_13","36_M7_2_3_4_5_12_13")
  
  # set up a matrix of all
  varset_matrix <- rbind(varset_baselineFactors, varset_M7ELISA, varset_M7Tcells,varset_M7ADCP,varset_M7IgG3,
                         varset_M7IgG3_gp140,varset_M7IgG3_gp120,varset_M7IgG3_V1V2,varset_M7IgG3_gp41,
                         varset_M7IgG3_breadthScore,varset_M7IgG3_multi,varset_M7IgGt_V1V2,varset_M7ADCC,
                         varset_M7_combPrimaryMarkers, varset_M7_combAllMarkers,varset_M7_overallScore, 
                         varset_M7_2_3, varset_M7_2_4, varset_M7_2_5_12, varset_M7_2_13, varset_M7_3_4, 
                         varset_M7_3_5_12, varset_M7_3_13, varset_M7_4_5_12, varset_M7_4_13, 
                         varset_M7_5_12_13, varset_M7_2_3_4, varset_M7_2_3_5_12,varset_M7_2_3_13,varset_M7_3_4_5_12,
                         varset_M7_3_4_13,varset_M7_4_5_12_13,varset_M7_2_3_4_5_12,varset_M7_2_3_4_13,
                         varset_M7_3_4_5_12_13,varset_M7_2_3_4_5_12_13)
}




if (study_name %in% c("ENSEMBLE")) {
  # 1. None (No markers; only baseline risk variables), phase 2 data
  varset_baselineRiskFactors <- rep(FALSE, length(markers))
  
  # 2-16 (Day29)
  varset_bAbSpike_D29 <- create_varsets(markers, grep("(?=.*29)(?=.*bindSpike)", markers, value=TRUE, perl=TRUE))
  varset_bAbRBD_D29 <- create_varsets(markers, grep("(?=.*29)(?=.*bindRBD)", markers, value=TRUE, perl=TRUE))
  varset_pnabID50_D29 <- create_varsets(markers, grep("(?=.*29)(?=.*id50)", markers, value=TRUE, perl=TRUE))
  varset_ADCP_D29 <- create_varsets(markers, grep("(?=.*29)(?=.*ADCP)", markers, value=TRUE, perl=TRUE))
  varset_FunctionalminusbAb_D29 <- create_varsets(markers, grep("(?=.*29)(?=.*minus)", markers, value=TRUE, perl=TRUE))
  varset_bAb_pnabID50_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*29)(?=.*id50)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_ADCP_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*29)(?=.*ADCP)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_pnabID50_diff_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)(^((?!ADCP).)*$)", markers, value=TRUE, perl=TRUE),
                                                       grep("(?=.*29)(?=.*id50)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_ADCP_diff_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)(^((?!id50).)*$)", markers, value=TRUE, perl=TRUE),
                                                        grep("(?=.*29)(?=.*ADCP)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_pnabID50_ADCP_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)(^((?!minus).)*$)", markers, value=TRUE, perl=TRUE),
                                                            grep("(?=.*29)(?=.*id50)(^((?!minus).)*$)", markers, value=TRUE, perl=TRUE),
                                                            grep("(?=.*29)(?=.*ADCP)(^((?!minus).)*$)", markers, value=TRUE, perl=TRUE)))
  varset_allMarkers_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)", markers, value=TRUE, perl=TRUE),
                                                     grep("(?=.*29)(?=.*id50)", markers, value=TRUE, perl=TRUE),
                                                     grep("(?=.*29)(?=.*ADCP)", markers, value=TRUE, perl=TRUE)))
  varset_bAb_combScores_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*bind)(^((?!minus).)*$)", markers, value=TRUE, perl=TRUE),
                                                         grep("(?=.*29)(?=.*comb)", markers, value=TRUE, perl=TRUE)))
  varset_pnabID50_combScores_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*id50)(^((?!minus).)*$)", markers, value=TRUE, perl=TRUE),
                                                         grep("(?=.*29)(?=.*comb)", markers, value=TRUE, perl=TRUE)))
  varset_ADCP_combScores_D29 <- create_varsets(markers, c(grep("(?=.*29)(?=.*ADCP)(^((?!minus).)*$)", markers, value=TRUE, perl=TRUE),
                                                         grep("(?=.*29)(?=.*comb)", markers, value=TRUE, perl=TRUE)))
  varset_allMarkers_combScores_D29 <- create_varsets(markers, grep("(?=.*29)", markers, value=TRUE, perl=TRUE))
  
  varset_names <- c("1_baselineRiskFactors",
                    "2_bAbSpike_D29", "3_bAbRBD_D29", "4_pnabID50_D29",
                    "5_ADCP_D29", "6_FunctionalminusbAb_D29", "7_bAb_pnabID50_D29", "8_bAb_ADCP_D29", 
                    "9_bAb_pnabID50_diff_D29", "10_bAb_ADCP_diff_D29", 
                    "11_bAb_pnabID50_ADCP_D29", "12_allMarkers_D29",
                    "13_bAb_combScores_D29", "14_pnabID50_combScores_D29", "15_ADCP_combScores_D29",
                    "16_allMarkers_combScores_D29")
  
  # set up a matrix of all
  varset_matrix <- rbind(varset_baselineRiskFactors,
                         varset_bAbSpike_D29, varset_bAbRBD_D29, varset_pnabID50_D29,
                         varset_ADCP_D29, varset_FunctionalminusbAb_D29, varset_bAb_pnabID50_D29, varset_bAb_ADCP_D29, 
                         varset_bAb_pnabID50_diff_D29, varset_bAb_ADCP_diff_D29, 
                         varset_bAb_pnabID50_ADCP_D29, varset_allMarkers_D29, 
                         varset_bAb_combScores_D29, varset_pnabID50_combScores_D29, varset_ADCP_combScores_D29, 
                         varset_allMarkers_combScores_D29)
}

# add on all of the individual marker variables
for (i in seq_len(length(markers))) {
  this_varset <- grepl(paste0("\\b", markers[i], "\\b"), markers, perl = TRUE)
  varset_matrix <- rbind(varset_matrix, this_varset)
  varset_names <- c(varset_names, markers[i])
}

# Save varset_names for running batch job on cluster!
if(job_id == 1){
  if(study_name == "COVEBoost"){
    varset_names %>% as.data.frame() %>% rename(varset_names = ".") %>%
      mutate(varset_no = seq.int(nrow(.))) %>%
      select(varset_no, varset_names) %>%
      write.csv(here("output", Sys.getenv("TRIAL"), SLcohort, SLrunBaseline, "varset_names.csv"))
  } else {
    varset_names %>% as.data.frame() %>% rename(varset_names = ".") %>%
      mutate(varset_no = seq.int(nrow(.))) %>%
      select(varset_no, varset_names) %>%
      write.csv(here("output", Sys.getenv("TRIAL"), "varset_names.csv"))
  }
  
}

# Study-agnostic set up of final data to pass to Super Learner -----------------

# all possible covariates for the Super Learner (baseline risk factors + all markers)
X_covars2adjust_ph2_init <- dat.ph2 %>% select(all_of(c(briskfactors, markers)))
X_covars2adjust_ph2 <- X_covars2adjust_ph2_init
# scale all variables to have mean 0, sd 1
for (a in colnames(X_covars2adjust_ph2)) {
  X_covars2adjust_ph2[[a]] <- scale(X_covars2adjust_ph2_init[[a]],
                                    center = mean(X_covars2adjust_ph2_init[[a]], na.rm = TRUE),
                                    scale = sd(X_covars2adjust_ph2_init[[a]], na.rm = TRUE))
}

# the Super Learner library
sl_lib <- SL_library

# the outcome
Y <- dat.ph2 %>% pull(endpoint)

# Study-specific final weights, treatment dataset
if (study_name %in% c("COVE", "MockCOVE")) {
  weights <- dat.ph2$wt.D57
  treatmentDAT <- dat.ph2 %>%
    select(all_of(ptidvar), Trt, wt.D57, EventIndPrimaryD57, all_of(c(briskfactors, markers))) %>%
    filter(Trt == 1) %>%
    select(-Trt)
} else if (study_name == "COVEBoost") {
  weights <- dat.ph2$wt.BD29
  treatmentDAT <- dat.ph2 %>%
    select(all_of(ptidvar), Trt, wt.BD29, EventIndPrimaryOmicronBD29, all_of(c(briskfactors, markers))) 
} else if (study_name == "HVTN705") {
  weights <- dat.ph2$wt.D210
  treatmentDAT <- dat.ph2 %>%
    select(all_of(ptidvar), Trt, wt.D210, Delta.D210, all_of(c(briskfactors, markers))) %>%
    filter(Trt == 1) %>%
    select(-Trt)
} else if (study_name == "ENSEMBLE") {
  weights <- dat.ph2$wt.D29
  treatmentDAT <- dat.ph2 %>%
    select(all_of(ptidvar), Trt, wt.D29, all_of(c(endpoint, briskfactors, markers))) %>%
    filter(Trt == 1) %>%
    select(-Trt)
}

# match the rows in treatmentDAT to get Z, C
all_cc_treatment <- Z_plus_weights %>%
  filter(Ptid %in% treatmentDAT$Ptid)
# pull out the participants who are NOT in the cc cohort and received the vaccine
all_non_cc_treatment <- Z_plus_weights %>%
  filter(!(Ptid %in% treatmentDAT$Ptid))
# put them back together
phase_1_data_treatmentDAT <- dplyr::bind_rows(all_cc_treatment, all_non_cc_treatment) %>%
  select(-Trt)
# indicator of observed in phase 2
C <- (phase_1_data_treatmentDAT$Ptid %in% treatmentDAT$Ptid)
# all outcomes; first come outcomes in phase 2 sample, then the remainder
full_y <- phase_1_data_treatmentDAT %>%
  pull(!!endpoint)

if (study_name %in% c("COVE", "MockCOVE")) {
  Z_treatmentDAT <- phase_1_data_treatmentDAT %>%
    select(-Ptid, -wt.D57)
  all_ipw_weights_treatment <- phase_1_data_treatmentDAT %>%
    pull(wt.D57)
} else if (study_name == "COVEBoost") {
  Z_treatmentDAT <- phase_1_data_treatmentDAT %>%
    select(-Ptid, -wt.BD29)
  all_ipw_weights_treatment <- phase_1_data_treatmentDAT %>%
    pull(wt.BD29)
} else if (study_name == "HVTN705") {
  Z_treatmentDAT <- phase_1_data_treatmentDAT %>%
    select(-Ptid, -wt.D210)
  all_ipw_weights_treatment <- phase_1_data_treatmentDAT %>%
    pull(wt.D210)
} else if (study_name == "ENSEMBLE") {
  Z_treatmentDAT <- phase_1_data_treatmentDAT %>%
    select(-Ptid, -all_of(wt))
  all_ipw_weights_treatment <- phase_1_data_treatmentDAT %>%
    pull(all_of(wt))
}

# set up outer folds for cv variable importance; do stratified sampling
if(Sys.getenv("TRIAL") %in%  c("janssen_pooled_partA", "janssen_la_partA")){
  V_outer <- 10
}else{
  V_outer <- 5
}

if (sum(dat.ph2 %>% pull(endpoint)) <= 25) {
  V_inner <- length(Y) - 1
  maxVar <- 5
} else if(sum(dat.ph2 %>% pull(endpoint)) > 25){
  V_inner <- 5
  maxVar <- floor(nv/6)
}

if (study_name %in% c("COVE", "COVEBoost")){
  V_inner <- length(Y) - 1
  maxVar <- 5
}

# Use inner leave-one-out CV even though there are enough number of HIV cases in phase 2 (51 HIV cases)
# Based off Peter's email, Thu 11/11/2021 3:08 PM
if (study_name == "HVTN705"){
  V_inner <- length(Y) - 1
  maxVar <- floor(nv/6)
}

