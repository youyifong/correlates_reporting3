#Sys.setenv(TRIAL = "moderna_boost"); Sys.setenv(VERBOSE = 1)

renv::activate(project = here::here(".."))     
source(here::here("..", "_common.R"))

library(survey)
library(plotrix) # weighted.hist
library(forestplot)
library(Hmisc) # wtd.quantile, cut2
library(kyotil)

time.start=Sys.time()
TRIAL=Sys.getenv("TRIAL")
myprint(TRIAL)
myprint(verbose)
begin=Sys.time()
print(date())


# read analysis ready data
dat = read.csv(config$data_cleaned)

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

# define subsets of data    
dat.vac.naive=subset(dat,  Trt==1 & naive & ph1.BD29)
dat.pla.naive=subset(dat,  Trt==0 & naive & ph1.BD29)
dat.vac.nnaive=subset(dat, Trt==1 & !naive & ph1.BD29)
dat.pla.nnaive=subset(dat, Trt==0 & !naive & ph1.BD29)


# loop through each quadrant
for (idat in 1:4) {
  if (idat==1) {dat.ph1 = dat.vac.naive;  ilabel="vac_naive"}
  if (idat==2) {dat.ph1 = dat.pla.naive;  ilabel="pla_naive"}
  if (idat==3) {dat.ph1 = dat.vac.nnaive; ilabel="vac_nnaive"}
  if (idat==4) {dat.ph1 = dat.pla.nnaive; ilabel="pla_nnaive"}
  
  save.results.to = paste0(save.results.to.0, "/", ilabel, "/") 
  if (!dir.exists(save.results.to))  dir.create(save.results.to)
  print(paste0("save results to ", save.results.to))
  
  dat.ph1$ph1=dat.ph1$ph1.BD29
  dat.ph1$ph2=dat.ph1$ph2.BD29
  dat.ph1$wt=dat.ph1$wt.BD29
  dat.ph1$EventIndPrimary =dat.ph1$EventIndOmicronBD29
  dat.ph1$EventTimePrimary=dat.ph1$EventTimeOmicronBD29

  dat.ph2 = subset(dat.ph1, ph2)
  
  tfinal.tpeak=get.tfinal.tpeak.case.control.rule1 (dat.ph2) 
  myprint(tfinal.tpeak)
  write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study_name))

  # define trichotomized markers
  dat.ph1 = add.trichotomized.markers (dat.ph1, "BD29"%.%assays)
  marker.cutpoints=attr(dat.ph1, "marker.cutpoints")
  for (a in "BD29"%.%assays) {        
      q.a=marker.cutpoints[[a]]
      to.write = paste0(a, " [", concatList(round(q.a, 2), ", "), ")%")
      write(to.write, file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
  }
  
  #create twophase design object
  design.vacc.naive<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)
  with(dat.ph1, table(Wstratum, ph2))
      
  # table of ph1 and ph2 cases
  tab=with(dat.ph1, table(ph2, EventIndPrimary))
  names(dimnames(tab))[2]="Event Indicator"
  print(tab)
  mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)
  
  # estimate overall marginalized risk (no markers) and VE
  source(here::here("code", "cor_coxph_marginalized_risk_no_marker.R"))
  

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

}


print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
