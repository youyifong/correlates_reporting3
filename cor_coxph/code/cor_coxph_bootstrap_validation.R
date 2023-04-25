# checking if cove.boost.collapse.strata is sensible

Sys.setenv(TRIAL = "moderna_boost"); Sys.setenv(VERBOSE = 1)

renv::activate(project = here::here(".."))     
source(here::here("..", "_common.R"))
source("code/params.R")

library(survey)
library(plotrix) # weighted.hist
library(forestplot)
library(Hmisc) # wtd.quantile, cut2

time.start=Sys.time()
TRIAL=Sys.getenv("TRIAL")
myprint(TRIAL)
myprint(verbose)
begin=Sys.time()
print(date())

# need this function b/c svycoxh may error due to singularity if, e.g. all cases have the same marker value
run.svycoxph=function(f, design) {
  fit=try(svycoxph(f, design=design), silent=T)
  if (class(fit)[1]=="try-error") NA else fit
}

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



idat=1
myprint(idat)
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

dat.ph1$yy=dat.ph1$EventIndPrimary

dat.ph2 = subset(dat.ph1, ph2)

tfinal.tpeak=get.tfinal.tpeak.case.control.rule1 (dat.ph2) 
myprint(tfinal.tpeak)
write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study_name))

# define trichotomized markers
dat.ph1 = add.trichotomized.markers (dat.ph1, "BD29"%.%assays)
marker.cutpoints=attr(dat.ph1, "marker.cutpoints")
for (a in "BD29"%.%assays) {        
  q.a=marker.cutpoints[[a]]
  to.write = paste0(gsub("_", "\\\\_",a), " [", concatList(round(q.a, 2), ", "), ")%")
  write(to.write, file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
}

# table of ph1 and ph2 cases
tab=with(dat.ph1, table(ph2, EventIndPrimary))
names(dimnames(tab))[2]="Event Indicator"
print(tab)
mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)

design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)

f=Surv(EventTimePrimary, EventIndPrimary) ~ MinorityInd + HighRiskInd + risk_score + BD29bindSpike

fit=svycoxph(f, design=design.1); fit

dat.ph2=subset(dat.ph1, ph2)
coxph(f, dat.ph2, weights=dat.ph2$wt) # same est as fit

# when mc.cores is >1, would get some strang results. some fits will be error when it will fit when run individually
out=mclapply(1:1000, mc.cores = 10, FUN=function(seed) {  
  if (verbose>=2) myprint(seed) 
  dat.b = try(bootstrap.cove.boost(dat.ph1, seed))
  if (inherits (dat.b, "try-error")) return (NULL)
  
  dat.ph2.b=subset(dat.b, ph2)
  fit=try(do.call("coxph", list(formula=f, data=dat.ph2.b, weights=dat.ph2.b$wt)))
  if(!inherits(fit,"try-error")) fit$coefficients else NULL
})

boot=do.call(cbind, out)
sd(boot["BD29bindSpike",])

# comparing analytical stderr and bootstrap ones

  

print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
