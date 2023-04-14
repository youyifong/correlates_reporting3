#Sys.setenv(TRIAL = "vat08m_naive"); COR="D43"; Sys.setenv(VERBOSE = 1)
#Sys.setenv(TRIAL = "moderna_mock"); COR="D29"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_pooled_partA"); COR="D29IncludeNotMolecConfirmedstart1"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "azd1222"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "moderna_real"); COR="D57over29"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "moderna_real"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "azd1222_bAb"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_pooled_EUA"); COR="D29IncludeNotMolecConfirmedstart1"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "profiscov"); COR="D91"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "profiscov_lvmn"); COR="D43start48"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "azd1222"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_la_partAsenior"); COR="D29IncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "prevent19"); COR="D35"; Sys.setenv(VERBOSE = 1)
#Sys.setenv(TRIAL = "janssen_la_partA"); COR="D29IncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "hvtn705second"); COR="D210"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_pooled_partA"); COR="D29SevereIncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 

renv::activate(project = here::here(".."))     
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))    
source(here::here("..", "_common.R"))

library(kyotil) # p.adj.perm, getFormattedSummary
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(plotrix) # weighted.hist
library(parallel)
library(forestplot)
library(Hmisc) # wtd.quantile, cut2
library(xtable) # this is a dependency of kyotil
time.start=Sys.time()
myprint(study_name)
myprint(verbose)


dat = read.csv(config$data_cleaned)
assay_metadata = read.csv(config$assay_metadata)
assays=assay_metadata$assay

source(here::here("code", "params.R"))




# path for figures and tables etc
save.results.to = here::here("output");                                if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", attr(config,"config")); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");               if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B)
myprint(numPerm)

# uloq censoring, done here b/c should not be done for immunogenicity reports
# note that if delta are used, delta needs to be recomputed
for (a in assays) {
  uloq=assay_metadata$uloq[assay_metadata$assay==a]
  for (t in c("BD1", "BD29")  ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloq), log10(uloq), dat.mock[[t %.% a]])
  }
}    


myprint(tfinal.tpeak)
write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study_name))

    
dat.vac.naive=subset(dat,  Trt==1 & naive & ph1)
dat.pla.naive=subset(dat,  Trt==0 & naive & ph1)
dat.vac.nnaive=subset(dat, Trt==1 & !naive & ph1)
dat.pla.nnaive=subset(dat, Trt==0 & !naive & ph1)

# for use in competing risk estimation
dat.vac.naive.ph2=subset(dat.vac.naive, ph2)


# define trichotomized markers
dat.vac.naive = add.trichotomized.markers (dat.vac.naive, all.markers, wt.col.name="wt")
marker.cutpoints=attr(dat.vac.naive, "marker.cutpoints")
for (a in all.markers) {        
    q.a=marker.cutpoints[[a]]
    if (startsWith(a, "Day")) {
        # not fold change
        write(paste0(labels.axis[1,get.assay.from.name(a)], " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
    } else {
        # fold change
        write(paste0(a, " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
    }
}


#create twophase design object
design.vacc.naive<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.vac.naive)
with(dat.vac.naive, table(Wstratum, ph2))
    

# table of ph1 and ph2 cases
tab=with(dat.vac.naive, table(ph2, EventIndPrimary))
names(dimnames(tab))[2]="Event Indicator"
print(tab)
mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)


begin=Sys.time()
print(date())

# some checks

#with(dat.vac.naive.ph2, weighted.mean(Day35bindRBD<log10(100), wt))

#with(dat.vac.naive.ph2, table(EventIndPrimary, is.na(seq1.spike.weighted.hamming)))
#with(dat.vac.naive.ph2, table(EventIndPrimary, is.na(seq1.log10vl)))

#with(dat.vac.naive, table(EventIndPrimary, sieve.status, EventTimePrimary>0))
#with(dat.vac.naive, table(EventIndPrimaryD29, sieve.status, EventTimePrimary>0))
#with(dat.vac.naive, table(EventIndPrimaryIncludeNotMolecConfirmedD29, sieve.status, EventTimePrimary>0))
#with(dat.vac.naive, plot(EventTimePrimary, sieve.time, cex=.2)); abline(0,1)
#subset(dat.vac.naive, EventIndPrimary==0 & sieve.status==1)

#with(subset(dat.vac.naive, EventIndPrimaryIncludeNotMolecConfirmedD29==1), table(is.na(seq1.variant), EventIndPrimaryHasVLD29))


###################################################################################################
# estimate overall VE in the placebo and vaccine arms
###################################################################################################

source(here::here("code", "cor_coxph_marginalized_risk_no_marker.R"))

if(Sys.getenv("COR_COXPH_NO_MARKER_ONLY")==1) q("no")



###################################################################################################
# run PH models
###################################################################################################
    
source(here::here("code", "cor_coxph_ph.R"))


# unit testing of coxph results
if (Sys.getenv("TRIAL") == "janssen_pooled_EUA" & COR=="D29IncludeNotMolecConfirmedstart1") {
    tmp.1=c(rv$tab.1[,4], rv$tab.2[,"overall.p.0"])
    tmp.2=c("0.162","0.079","0.006",      "0.498","   ","   ","0.162","   ","   ","0.003","   ","   ")
    assertthat::assert_that(all(tmp.1==tmp.2), msg = "failed cor_coxph unit testing")    
    
} else if (attr(config, "config")=="moderna_real" & COR=="D57") {
    assertthat::assert_that(all(abs(p.unadj-c(0.004803168, 0.002172787, 0.000129743, 0.000202068, 0.064569846, 0.005631520, 0.009016447, 0.051800145, 0.011506959, 0.579164657))<1e-6), msg = "failed cor_coxph unit testing")    
    
} else if (attr(config, "config")=="prevent19" & COR=="D35") {
    assertthat::assert_that(all(abs(p.unadj-c(0.000453604, 0.0023274, 0.013258206))<1e-6), msg = "failed cor_coxph unit testing")    
    
}
print("Passed cor_coxph unit testing")    

###################################################################################################
# marginalized risk and controlled VE
###################################################################################################
    
# load ylims.cor[[1]] from D29 analyses, which is a list of two: 1 with placebo lines, 2 without placebo lines.
tmp=paste0(here::here(), paste0("/output/", attr(config,"config"), "/", COR, "/ylims.cor.", study_name, ".Rdata"))
if (file.exists(tmp)) load(tmp) # if it does not exist, the code will find alternative ylim

source(here::here("code", "cor_coxph_marginalized_risk_bootstrap.R"))

source(here::here("code", "cor_coxph_marginalized_risk_plotting.R"))

if (attr(config, "config") %in% c("moderna_real", "janssen_pooled_EUA")) source(here::here("code", "cor_coxph_samplesizeratio.R"))



print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
