#Sys.setenv(TRIAL = "moderna_boost"); Sys.setenv(VERBOSE = 1)

renv::activate(project = here::here(".."))     
source(here::here("..", "_common.R"))
source("code/params.R")

library(survey)
library(plotrix) # weighted.hist
library(forestplot)
library(Hmisc) # wtd.quantile, cut2

time.start=Sys.time()
myprint(TRIAL)
myprint(verbose)
begin=Sys.time()
print(date())

# need this function b/c svycoxh may error due to singularity if, e.g. all cases have the same marker value
run.svycoxph=function(f, design) {
  fit=try(svycoxph(f, design=design), silent=T)
  if (class(fit)[1]=="try-error") NA else fit
}


# path for figures and tables etc
save.results.to.0 = here::here("output");                                if (!dir.exists(save.results.to.0))  dir.create(save.results.to.0)
save.results.to.0 = paste0(save.results.to.0, "/", attr(config,"config")); if (!dir.exists(save.results.to.0))  dir.create(save.results.to.0)

# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B, numPerm)


# uloq censoring, done here b/c should not be done for immunogenicity reports
# note that if delta are used, delta needs to be recomputed
for (a in assays) {
  uloq=assay_metadata$uloq[assay_metadata$assay==a]
  for (t in c("BD1", "BD29")  ) {
    dat[[t %.% a]] <- ifelse(dat[[t %.% a]] > log10(uloq), log10(uloq), dat[[t %.% a]])
  }
}    


dat$ph1=dat$ph1.BD29
dat$ph2=dat$ph2.BD29
dat$wt=dat$wt.BD29
dat$EventIndPrimary =dat$EventIndOmicronBD29
dat$EventTimePrimary=dat$EventTimeOmicronBD29
dat$yy=dat$EventIndPrimary


# define subsets of data    
dat.vac.naive=subset(dat,  Trt==1 & naive & ph1.BD29)
dat.pla.naive=subset(dat,  Trt==0 & naive & ph1.BD29)
dat.vac.nnaive=subset(dat, Trt==1 & !naive & ph1.BD29)
dat.pla.nnaive=subset(dat, Trt==0 & !naive & ph1.BD29)


# compute overall risk in each quadrant
prev.vac.naive = get.marginalized.risk.no.marker(form.0, dat.vac.naive, tfinal.tpeak)
prev.pla.naive = get.marginalized.risk.no.marker(form.0, dat.pla.naive, tfinal.tpeak)
prev.vac.nnaive= get.marginalized.risk.no.marker(form.0, dat.vac.nnaive, tfinal.tpeak)
prev.pla.nnaive= get.marginalized.risk.no.marker(form.0, dat.pla.nnaive, tfinal.tpeak)
overall.risks=list(prev.vac.naive, prev.pla.naive, prev.vac.nnaive, prev.pla.nnaive)
myprint(prev.vac.naive, prev.pla.naive, prev.vac.nnaive, prev.pla.nnaive)

# a crude boxplot to get some idea of the distributions
# par(mfrow=c(2,2))
# for (idat in 1:4) {
#   myprint(idat)
#   if (idat==1) {dat.ph1 = dat.vac.naive;  ilabel="vac_naive"}
#   if (idat==2) {dat.ph1 = dat.pla.naive;  ilabel="pla_naive"}
#   if (idat==3) {dat.ph1 = dat.vac.nnaive; ilabel="vac_nnaive"}
#   if (idat==4) {dat.ph1 = dat.pla.nnaive; ilabel="pla_nnaive"}
#   myboxplot(BD29bindSpike~EventIndPrimary, dat.ph1, main=ilabel)
# }  


###################################################################################################
# loop through each quadrant
# 4 mock data not working yet
for (idat in 1:4) {
  # idat=1
  myprint(idat)
  if (idat==1) {dat.ph1 = dat.vac.naive;  ilabel="vac_naive"}
  if (idat==2) {dat.ph1 = dat.pla.naive;  ilabel="pla_naive"}
  if (idat==3) {dat.ph1 = dat.vac.nnaive; ilabel="vac_nnaive"}
  if (idat==4) {dat.ph1 = dat.pla.nnaive; ilabel="pla_nnaive"}
  
  save.results.to = paste0(save.results.to.0, "/", ilabel, "/") 
  if (!dir.exists(save.results.to))  dir.create(save.results.to)
  print(paste0("save results to ", save.results.to))
  
  #??? move out of this loop
  # add trichotomized markers. use the same cutpoints for all 4 quadrants
  dat.ph1 = add.trichotomized.markers (dat.ph1, "BD29"%.%assays)
  marker.cutpoints=attr(dat.ph1, "marker.cutpoints")
  for (a in "BD29"%.%assays) {        
    q.a=marker.cutpoints[[a]]
    to.write = paste0(gsub("_", "\\\\_",a), " [", concatList(round(q.a, 2), ", "), ")%")
    write(to.write, file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
  }
  
  dat.ph2 = subset(dat.ph1, ph2)
  
  write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study_name))
  
  # table of ph1 and ph2 cases
  tab=with(dat.ph1, table(ph2, EventIndPrimary))
  names(dimnames(tab))[2]="Event Indicator"
  print(tab)
  mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)
  
  
  
  ###################################################################################################
  # estimate overall marginalized risk (no markers) and VE
  
  source(here::here("code", "cor_coxph_marginalized_risk_no_marker.R"))
  

  ###################################################################################################
  # prepare for Cox models runs and margialized risks

  tpeak=29
  # the origin of followup days, may be different from tpeak, e.g., D43start48
  tpeak1 = 29
  
  all.markers = paste0("BD29", assays); names(all.markers)=all.markers
  
  
  ###################################################################################################
  # run PH models
  
  design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)
  with(dat.ph1, table(Wstratum, ph2, useNA="ifany"))
  
  source(here::here("code", "cor_coxph_ph.R"))
  
  # unit testing
  if (TRIAL == "") {
    tmp.1=c(rv$tab.1[,4], rv$tab.2[,"overall.p.0"])
    tmp.2=c("0.162","0.079","0.006",      "0.498","   ","   ","0.162","   ","   ","0.003","   ","   ")
    assertthat::assert_that(all(tmp.1==tmp.2), msg = "failed cor_coxph unit testing")    
    print("Passed cor_coxph unit testing")    
  } 
  
    
  ###################################################################################################
  # marginalized risk and controlled VE

  comp.risk=F
  COR=idat # only used in table figure labels
  
  source(here::here("code", "cor_coxph_marginalized_risk_bootstrap.R"))
  
  source(here::here("code", "cor_coxph_marginalized_risk_plotting.R"))

  # source(here::here("code", "cor_coxph_samplesizeratio.R"))
  
}


print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
